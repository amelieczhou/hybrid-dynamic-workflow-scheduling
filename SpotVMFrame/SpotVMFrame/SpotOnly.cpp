#include "stdafx.h"
#include "Algorithms.h"
#include "PricingModel.h"
#include "ReadTrace.h"
#include <time.h>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <string>
#include <omp.h>

using namespace boost;

extern float lambda;

float recur_cal(int num_multiplier, int maxn,float* ffp){
	if(num_multiplier==1){
//		printf("mutiply %d\n",maxn);
		return ffp[maxn];
	}
	else{
		float result = 0.0;
		for(int i=1; i<=maxn-num_multiplier+1; i++){
			float tmp=recur_cal(num_multiplier-1,maxn-i,ffp);
			result += ffp[i]*tmp;
//			printf("%d \n",i);
		}
		return result;
	}	
}
float calcu_pt(float* firstfail,float n,int rowindex,int spottype){
	float result = 0.0;//returned probability
	int columns=600;
	float* ffp = (float*)malloc(columns*sizeof(float));
	for(int i=0; i<columns; i++)
		if(rowindex>=2520)
			ffp[i] = 0.0;
		else
		ffp[i] = firstfail[rowindex*columns*4+spottype*columns+i];
	int maxn = 3;
	if(n<3) maxn = n;
	for(int k=0; k<maxn; k++){//1 to n multipliers
		int starter = n-k;
		int multiplier = k+1;
//		printf("calculating %d multipliers\n",multiplier);
		/*for(int count=1; count<=multiplier; count++){
			int tmp = count;
			int tmp2 = k;
			while(tmp>0){
				recur_cal(tmp,tmp2,ffp);
				tmp--;
				tmp2--;
			}		
		}*/
		result += recur_cal(multiplier,n,ffp);
	}
	free(ffp);
	return result;
}

bool check_hitrate(Graph* graph, float deadline, float meetdl, float spotdistrgap){
	//do the convolcation
	float finalprob = 0.0;
	int spotdistrsize = 400;
	std::pair<vertex_iter,vertex_iter> vvp = vertices(*graph);
	int* sizes = (int*)malloc((*vvp.second-*vvp.first)*sizeof(int));//save the distribution size for each task
	for(; vvp.first != vvp.second; vvp.first++){
		(*graph)[*vvp.first].tag = false;
		sizes[*vvp.first] = spotdistrsize;		
		for(int i=0; i<4000; i++)
			(*graph)[*vvp.first].cumulativeTime[i]=0;//(*graph)[*vvp.first].randspot[i];
	}

	vvp = vertices(*graph);
	for(; vvp.first != vvp.second; vvp.first++){
		//find all its parents
		in_edge_iterator in_i, in_end;
		edge_descriptor e;
		boost::tie(in_i, in_end) = in_edges(*vvp.first, (*graph));
		if(in_i == in_end) {
			(*graph)[*vvp.first].tag = true;
			float sumtask = 0.0;
			for(int i=0; i<spotdistrsize; i++){
				(*graph)[*vvp.first].cumulativeTime[i] = (*graph)[*vvp.first].randspot[i];
				sumtask += (*graph)[*vvp.first].cumulativeTime[i];
			}
			if(sumtask>1.0)
				printf("");
		}
		else{//max operation, parents may have different cumulative sizes
			if(*vvp.first == 16)
				printf("");
			int maxparentsize = 0;
			int maxiter = 0;
			for(; in_i!=in_end; ++in_i){
				e = *in_i;
				int src = source(e,(*graph));
				if(maxparentsize<sizes[src]){
					maxparentsize=sizes[src];
					maxiter = src;
				}
			}
			float* maxdistr = (float*)malloc(maxparentsize*sizeof(float));
			float* tmpdistr = (float*)malloc(maxparentsize*sizeof(float));
			for(int i=0; i<maxparentsize; i++){
				maxdistr[i] = (*graph)[maxiter].cumulativeTime[i];
				tmpdistr[i] = 0.0;
			}
			float sumtask = 0.0;
			boost::tie(in_i, in_end) = in_edges(*vvp.first, (*graph));
			for (; in_i != in_end; ++in_i) {
				e = *in_i;
				Vertex src = source(e,(*graph));
				if(src != maxiter) {
					calmaxdistr(maxdistr,(*graph)[src].cumulativeTime,tmpdistr,maxparentsize,sizes[src]);//cumu<0
					for(int kk=0; kk<maxparentsize; kk++){
						//maxdistr[kk]=maxdistr[kk]>(*graph)[src].cumulativeTime[kk]?maxdistr[kk]:(*graph)[src].cumulativeTime[kk];
						maxdistr[kk]=tmpdistr[kk];
				//		printf("%f\t",tmpdistr[kk]);
						//sumtask += maxdistr[kk];
					}								
				}
			}
			for(int kk=0; kk<maxparentsize; kk++)
				sumtask += maxdistr[kk];
			if(sumtask>1.0)
				printf("");
			free(tmpdistr);
			int colsize = 400+maxparentsize-1;//sizes[*vvp.first]+maxparentsize-1;
			float* newdistr = (float*)malloc(colsize*sizeof(float));
			//conv(maxdistr,(*graph)[*vvp.first].cumulativeTime,newdistr,maxparentsize,sizes[*vvp.first]);
			conv(maxdistr,(*graph)[*vvp.first].randspot,newdistr,maxparentsize,400);
			sumtask = 0.0;
			for(int j=0; j<colsize; j++){
				(*graph)[*vvp.first].cumulativeTime[j]=newdistr[j];
				sumtask += newdistr[j];
			}
			if(sumtask>1.0) 
				printf("");
			sizes[*vvp.first] = colsize;
			(*graph)[*vvp.first].tag = true;
			free(newdistr);
			free(maxdistr);
		}
	}//for vvp.first
	for(int i=0; i<4000; i++){
		//printf("%f\t",(*graph)[*vvp.second-1].cumulativeTime[i]);
		if(i*spotdistrgap<deadline)
			finalprob += (*graph)[*vvp.second-1].cumulativeTime[i];
	}

	//finalprob /= randomsize;
	if(finalprob>=meetdl) return true;
	else return false;
}
void update_randspot(DAG* dag, float* firstfail, int spotdistrsize, float spotdistrgap){
	int columns = 2400;
	std::pair<vertex_iter,vertex_iter> vp = vertices(*dag->g);
	for(; vp.first!=vp.second; vp.first++){
		for(int j=0; j<spotdistrsize; j++){
			(*dag->g)[*vp.first].randspot[j]=0.0;
		}
		int assignedtype = (*dag->g)[*vp.first].assigned_type;
		int rowindex = (int)((*dag->g)[*vp.first].prices[0]*1000) -1;
		for(int rndindex=0; rndindex<randomsize; rndindex++){
			std::vector<std::pair<float,float> > distri;
			float shortesttime = (*dag->g)[*vp.first].probestTime[assignedtype*randomsize+rndindex];
			float pinitial = 0.0;
			for(int i=0; i<=shortesttime; i++){
				if(rowindex>=2520)
					pinitial += 0;
				else
				pinitial += firstfail[rowindex*columns+assignedtype*columns/4+i];
			}
			pinitial = 1-pinitial;
			distri.push_back(std::pair<float,float>(shortesttime,pinitial));
			float t = shortesttime+1;
			while(t<2*shortesttime){
				float pt = 0.0;//the probability of t
				pt = pinitial*calcu_pt(firstfail,t-shortesttime,rowindex,assignedtype);
				distri.push_back(std::pair<float,float>(t,pt));
				t ++;
			}
			//assign to spotrand
			for(int i=0; i<distri.size(); i++){
				int id = (int)std::ceil(distri[i].first/spotdistrgap);
				if(id >= spotdistrsize) 
					(*dag->g)[*vp.first].randspot[spotdistrsize-1] += distri[i].second;
				else (*dag->g)[*vp.first].randspot[id] += distri[i].second;
			}
		}//rndindex
		for(int i=0; i<spotdistrsize; i++)
			(*dag->g)[*vp.first].randspot[i] /= randomsize;
	}
}
void SpotOnly::bidpDeter(){
	//read in the ffp
	int columns = 2400;
	int spotdistrsize = 400;
	float spotdistrgap = 1;//50;//mongtage deadline 17000//ceil(dag.deadline / spotdistrsize);
	if(dag->type == ligo) spotdistrgap = 2;//125.0;//deadline 42000
	else if(dag->type == epigenome) spotdistrgap = 1;//40.0;//deadline 15000
	float *firstfail = new float[2520*2400];//[760*2400];
	ReadTrace* tr = new ReadTrace();
	tr->readtomem("newData/ffp_2.dat",firstfail);
	tr->closefile();
	delete tr;

	std::pair<vertex_iter,vertex_iter> vp = vertices(*dag->g);
	for(;vp.first!=vp.second;vp.first++) {		
		int assignedtype = (*dag->g)[*vp.first].assigned_type;
		(*dag->g)[*vp.first].prices = new float[1];
		(*dag->g)[*vp.first].prices[0] = 0.5*priceOnDemand[assignedtype];	
	}
	update_randspot(dag,firstfail,spotdistrsize,spotdistrgap);
	for(int i=0; i<400; i++)
		printf("%f\t",(*dag->g)[*vp.second-1].randspot[i]);
	//check if the deadline hit rate is satisfied
	bool satisfied = check_hitrate(dag->g,dag->deadline,dag->meet_dl,spotdistrgap);
	if(satisfied)//satisfied
		return;
	else {
		while(!satisfied){
			//for each task, search for bidding price
			vp = vertices(*dag->g);
			for(;vp.first!=vp.second;vp.first++){
				//increase bidding price
				(*dag->g)[*vp.first].prices[0]+= 0.02;
				float test=(*dag->g)[*vp.first].prices[0];
			}
			//update randspot and re-check
			update_randspot(dag,firstfail,spotdistrsize,spotdistrgap);
			satisfied = check_hitrate(dag->g,dag->deadline,dag->meet_dl,spotdistrgap);
		}
	}
}

void SpotOnly::Simulate(){

	float arrival_time = 0;
	dag->arrival_time = 0;	
	float violation = 0;
	float ave_cost = 0;
	float ioseq[types],iorand[types],net_up[types],net_down[types];
	omp_set_num_threads(20);
	float viol_private[20];
	float cost_private[20];
	float fail_private[20];
	for(int i=0; i<20; i++) fail_private[i]=viol_private[i]=cost_private[i]=0;

	std::pair<vertex_iter, vertex_iter> vp;
	time_t start,end;
	std::vector<DAG*> workflows; //continuous workflow
	workflows.push_back(dag);

	std::ifstream infile;
	std::string a = "arrivaltime_integer_";
	std::string b;
	std::ostringstream strlamda;
	strlamda << lambda;
	b = strlamda.str();
	std::string c = ".txt";
	std::string fname = a + b + c;
	char timee[256];
	infile.open(fname.c_str());
	if(infile==NULL){
		printf("cannot find input file!\n");
		return;
	}
	infile.getline(timee,256); //jump the lamda line
	infile.getline(timee,256); //jump the 0 line
	//incomming jobs
	//while(arrival_time < max_t){
	while(workflows.size()<(int)num_jobs){
		infile.getline(timee,256);
		arrival_time = atof(timee);

		DAG* job = new DAG(dag->deadline+arrival_time,dag->meet_dl);		
		job->g = dag->g; job->type = dag->type;
		job->arrival_time = arrival_time;
		vp = vertices(*job->g);
		for(int i=0; i<(*vp.second - *vp.first); i++)
			(*job->g)[i].sub_deadline += arrival_time;
		workflows.push_back(job);
	}
	infile.close();
	float* datatrace = new float[35971*4];
	ReadTrace* tr = new ReadTrace();
	tr->readtomem("data2.csv",datatrace);
	tr->closefile();
	delete tr;

	//start simulation
	time(&start);
	#pragma omp parallel
	{
		#pragma omp for
		for(int monte=0; monte < randomsize; monte++)
		{
			std::vector<DAG*> jobs;
			for(int i=0; i<workflows.size(); i++){
				DAG* newdag = new DAG(*workflows[i]);
				vp = vertices(*newdag->g);
				for(int j=0; j<(*vp.second - *vp.first); j++)
					(*newdag->g)[j].sub_deadline = (*workflows[i]->g)[j].sub_deadline;
				jobs.push_back(newdag);
			}
			std::vector<SpotVM*> sVMTP[types];

			//EDF scheduling
			double t = 0;
			bool condition = false;
			double moneycost = 0.0;
			float r[types];
			int totalspotfail =0;
			int totalspot = 0;
			do{	
				//accept workflows
				for(int i=0; i<jobs.size(); i++){
					if((int)t == (int)jobs[i]->arrival_time){
						if(dag->type == montage){
							for(int j=0; j<4; j++) {
								(*jobs[i]->g)[j].status = ready;
								(*jobs[i]->g)[j].readyCountdown = -1;
								(*jobs[i]->g)[j].restTime = 0;				
							}
							for(int j=4; j<20; j++){
								(*jobs[i]->g)[j].status = not_ready;
								(*jobs[i]->g)[j].readyCountdown = -1;
								(*jobs[i]->g)[j].restTime = 0;				
							}
						}else if(dag->type == ligo){
							for(int j=0; j<9; j++) {
								(*jobs[i]->g)[j].status = ready;
								(*jobs[i]->g)[j].readyCountdown = -1;
								(*jobs[i]->g)[j].restTime = 0;				
							}
							for(int j=9; j<40; j++){
								(*jobs[i]->g)[j].status = not_ready;
								(*jobs[i]->g)[j].readyCountdown = -1;
								(*jobs[i]->g)[j].restTime = 0;	
							}
						}else if(dag->type == epigenome){				
							(*jobs[i]->g)[0].status = ready;
							(*jobs[i]->g)[0].readyCountdown = -1;
							(*jobs[i]->g)[0].restTime = 0;	
							for(int j=1; j<20; j++){
								(*jobs[i]->g)[j].status = not_ready;
								(*jobs[i]->g)[j].readyCountdown = -1;
								(*jobs[i]->g)[j].restTime = 0;	
							}
						}else{
							printf("what is the dag type?");
							exit(1);
						}
					}
				}
				int spotfail = 0;
				//step 0, check spotVM
				if(fmod(t, 3) == 0){ //check every 3 minutes
					int tracelag = rand()%35000;
					for(int i=0; i<types; i++)
						r[i] = datatrace[tracelag*types+i];
					if(r != NULL) {			
					//		std::cout<<r[0]<<"\t"<<r[1]<<"\t"<<r[2]<<"\t"<<r[3]<<"\n";
					}
					else{
						std::cout<< "cannot read the price\n";	
					}
					for(int i=0; i<types; i++){
						int size = sVMTP[i].size();
						for(int j=0; j<size; j++){
							bool check;
							check = function(sVMTP[i][j]->price, r[i]);//if not fail		

							if(!check){		//failed	
								spotfail += 1;
								totalspotfail += 1;
								if(sVMTP[i][j]->tk != NULL)	{
									////this run fails, stop the run
									//condition = true;
									////add it to fail
									//int id= omp_get_thread_num();
									//fail_private[id] += 1;
									//goto finish;
									sVMTP[i][j]->canAlloc = false;
									sVMTP[i][j]->tk->status = ready;
									sVMTP[i][j]->tk->readyCountdown = -1;
									

									
									sVMTP[i][j]->tk = NULL;
								}
							}
						}
					}
				}
				//step 1
				std::vector<taskVertex*> ready_task;
				for(int ji=0; ji<jobs.size(); ji++){
					vp = vertices(*jobs[ji]->g);
					for(int i=0; i < (*vp.second - *vp.first ); i++)
					{
						bool tag = true;
						//get parent vertices
						in_edge_iterator in_i, in_end;
						edge_descriptor e;
						boost::tie(in_i, in_end) = in_edges(i, *jobs[ji]->g);
						if(in_i == in_end) tag = false;
						else{
							for (; in_i != in_end; ++in_i) 	{
								e = *in_i;
								Vertex src = source(e, *jobs[ji]->g);					
								if((*jobs[ji]->g)[src].status != finished)	{
									tag = false;
									//break;
								}
							}
						}
						if((*jobs[ji]->g)[i].status == ready || (tag && (*jobs[ji]->g)[i].status != scheduled && (*jobs[ji]->g)[i].status != finished)){
							ready_task.push_back(&(*jobs[ji]->g)[i]);							
						}
					}
				}
				
				//std::sort(ready_task.begin(),ready_task.end(), myfunction);
				for(int i=0; i<ready_task.size(); i++)//earliest deadline first
				{
					taskVertex* curr_task=ready_task[i];
					int _config = curr_task->assigned_type;
					if(curr_task->readyCountdown == -1)//
					{
						bool find = false;
						int foundindex = -1;
						//check VM/SpotVM list for available machine
						int size = sVMTP[_config].size();
						for(int j=0; j<size; j++)
						{
							if(sVMTP[_config][j]->tk == NULL&&sVMTP[_config][j]->life_time>=curr_task->probestTime[_config*randomsize+monte]&&sVMTP[_config][j]->canAlloc)
							{
								find = true;
								sVMTP[_config][j]->tk = curr_task;
								foundindex = j;
								break;
							}
						}
						if(find) {
							curr_task->status = scheduled;
							curr_task->tasktime = t;
							curr_task->restTime =  curr_task->probestTime[_config*randomsize+monte] ;
							curr_task->readyCountdown = 0;
							sVMTP[_config][foundindex]->tk = curr_task;
						}
						else 			
						{
							curr_task->readyCountdown = SpotLag;
							curr_task->tasktime = t;
						}
					}
					else if(curr_task->readyCountdown == 0)
					{
						curr_task->status = scheduled;
						curr_task->restTime = curr_task->probestTime[_config*randomsize+monte] ;

						SpotVM* vm = new SpotVM(curr_task->prices[0]); 
						vm->life_time = 60;
						vm->tk = curr_task;
						vm->type = _config;
						sVMTP[_config].push_back(vm);						
					}			
				}
				//delete VMs without task
				for(int i=0; i<types; i++)//////////
				{
					int size1 = sVMTP[i].size();
					
					for(int j=0; j<size1; j++)
					{
						if(sVMTP[i][j]->tk == NULL || !sVMTP[i][j]->canAlloc)
						{
							if(sVMTP[i][j]->tk == NULL && sVMTP[i][j]->canAlloc){
								moneycost += r[i];
							}

							SpotVM* vm = sVMTP[i][j];
							delete vm;
							sVMTP[i].erase(sVMTP[i].begin()+j);
							j--;
							size1--;
						}
					}
				}
				//step 2
				std::vector<taskVertex*> scheduled_task;
				for(int ji=0; ji<jobs.size(); ji++){
					vp = vertices((*jobs[ji]->g));
					for(int i=0; i<(*vp.second - *vp.first ); i++)
						if((*jobs[ji]->g)[i].status == scheduled)
							scheduled_task.push_back(&(*jobs[ji]->g)[i]);
				}
				
				for(int i=0; i<scheduled_task.size(); i++)
				{
					scheduled_task[i]->restTime -= 1;////////////////////////////
					if(scheduled_task[i]->restTime <= 0) 
					{
						scheduled_task[i]->status = finished;
						scheduled_task[i]->end_time = t;
						scheduled_task[i]->tasktime = t - scheduled_task[i]->tasktime;
						scheduled_task[i]->cost = scheduled_task[i]->tasktime * priceOnDemand[scheduled_task[i]->assigned_type] /60.0;
						//make the vm.task = NULL
						for(int j=0; j<sVMTP[scheduled_task[i]->assigned_type].size(); j++)
							if(sVMTP[scheduled_task[i]->assigned_type][j]->tk == scheduled_task[i])
							{
								sVMTP[scheduled_task[i]->assigned_type][j]->tk = NULL;
								break;
							}
					}
				}				
				//step 3
				for(int i=0; i<types; i++)
				{
					int size1 = sVMTP[i].size();			
					
					for(int j=0; j<size1; j++)
					{
						sVMTP[i][j]->life_time -= 1;//
						if(sVMTP[i][j]->life_time == 0){
							sVMTP[i][j]->life_time = 60;
							moneycost += r[i];
						}			
					}
				}
				for(int i=0; i<ready_task.size(); i++)//////////////////////////////////if >0
					if(ready_task[i]->readyCountdown > 0)
						ready_task[i]->readyCountdown -= 1;
				t += 1;

				condition = false;
				int unfinishednum = 0;
				for(int ji=0; ji<jobs.size(); ji++){
					vp = vertices((*jobs[ji]->g));
					for(int i=0; i < (*vp.second - *vp.first ); i++){
						if((*jobs[ji]->g)[i].status!= finished)
						{
							condition = true;
							unfinishednum += 1;
						}					
					}
				}
			}while(condition);//there is a task not finished	
			//step 4 finalize
			for(int i=0; i<types; i++)
			{
				int size = sVMTP[i].size();
				for(int j=0; j<size; j++)
				{				
					moneycost += r[i];				
				}
			}
			printf("Money Cost: %.4f, Time: %.2f\n", moneycost, t);
			int id = omp_get_thread_num();
			printf("thread id is %d\n",id);
			float ave_time = 0;
			for(int i=0; i<jobs.size(); i++){
				vp = vertices((*jobs[i]->g));
				float executiontime = (*jobs[i]->g)[*(vp.second-1)].end_time - jobs[i]->arrival_time;
				if(executiontime > jobs[i]->deadline) {
					viol_private[id] += 1.0;
				}		
				ave_time += executiontime;
			}
			cost_private[id] += moneycost;
			printf("average execution time of workflows is %f\n",ave_time/jobs.size());
			for(int i=0; i<jobs.size(); i++){
				delete jobs[i];
			}
			jobs.clear();
		}		
	}
	free(datatrace);
	for(int i=0; i<24; i++) {
		violation += viol_private[i];
		ave_cost += cost_private[i];
	}
	violation /= (float)randomsize*num_jobs;
	ave_cost /= (float)randomsize*num_jobs;
	printf("deadline meeting rate is %f, average cost is %f\n",1.0-violation,ave_cost);
	vp = vertices(*dag->g);
	for(; vp.first!=vp.second; vp.first++){
		free((*dag->g)[*vp.first].netDown);
		free((*dag->g)[*vp.first].netUp);
		free((*dag->g)[*vp.first].probestTime);
		free((*dag->g)[*vp.first].randomIO);
		free((*dag->g)[*vp.first].seqIO);
		free((*dag->g)[*vp.first].cumulativeTime);
		free((*dag->g)[*vp.first].randspot);
	}
	for(int i=0; i<24; i++) {
		violation += fail_private[i];
	}
	printf("number of fails is %f\n",violation);
}