#include "stdafx.h"
#include "Algorithms.h"
#include "PricingModel.h"
#include "ReadTrace.h"
#include <ctime>
#include <fstream>
#include <sstream>
#include <string>
#include <omp.h>

using namespace boost;
extern float lambda;
bool debug = false;
void Dyna::Simulate(){	
	std::vector<DAG*> workflows; //continuous workflow
	float arrival_time = 0;
	dag->arrival_time = 0;
	workflows.push_back(dag);
	std::pair<vertex_iter, vertex_iter> vvp;
	vvp = vertices(*dag->g);
	for(; vvp.first!=vvp.second; vvp.first++){
		printf("type: %d, %d; price: %f, %f\n",(*dag->g)[*vvp.first].configList[0],(*dag->g)[*vvp.first].configList[1],(*dag->g)[*vvp.first].prices[0],(*dag->g)[*vvp.first].prices[1]);
	}
	std::pair<vertex_iter, vertex_iter> vp; 
	std::ifstream infile;
	std::string a = "arrivaltime_integer_";
	std::string b;
	std::ostringstream strlamda;
	strlamda << lambda;
	b = strlamda.str();
	std::string c = ".txt";
	std::string fname = a + b + c;
	char time[256];
	infile.open(fname.c_str());
	if(infile==NULL){
		printf("cannot find input file!\n");
		return;
	}
	infile.getline(time,256); //jump the lamda line
	infile.getline(time,256); //jump the 0 line
	//incomming jobs
	//while(arrival_time < max_t){
	while(workflows.size()<(int)num_jobs){
		infile.getline(time,256);
		arrival_time = atof(time);

		DAG* job = new DAG(dag->deadline+arrival_time,dag->meet_dl);		
		job->g = dag->g; job->type = dag->type;
		job->arrival_time = arrival_time;
		vp = vertices(*job->g);
		workflows.push_back(job);
	}
	
	float violation = 0;
	float ave_cost = 0;
	float ioseq[types],iorand[types],net_up[types],net_down[types];
	//start simulation
	std::clock_t starttime = std::clock();
	float* datatrace = new float[35971*4];
	ReadTrace* tr = new ReadTrace();
	tr->readtomem("data2.csv",datatrace);
	tr->closefile();
	delete tr;
	float *tp, *firstfail;
	omp_set_num_threads(24);
	float viol_private[24];
	float cost_private[24];
	for(int i=0; i<24; i++) viol_private[i]=cost_private[i]=0;
	#pragma omp parallel
	{
		#pragma omp for
		for(int globaliter=0; globaliter<randomsize; globaliter++)
		{
			//each task has a distribution of execution time
			//dag.g[0].probestTime
			//assign resources in what order?
			std::vector<DAG*> jobs;
			for(int i=0; i<workflows.size(); i++){
				DAG* newdag = new DAG(*workflows[i]);
				vp = vertices(*newdag->g);
				for(int j=0; j<(*vp.second - *vp.first); j++)
					(*newdag->g)[j].sub_deadline = (*workflows[i]->g)[j].sub_deadline;
				jobs.push_back(newdag);
			}
			std::vector<VM*> VMTP[types];
			std::vector<SpotVM*> sVMTP[types];
			std::pair<vertex_iter, vertex_iter> vp;
				
			float t = 0; 
			bool condition = false; //there are at least one job not finished
			float moneycost = 0;
			////trace data
			//ReadTrace* tr = new ReadTrace();
			//tr->readfile("data.csv", tracelag);//

			float r[types];
			int totalspotfail =0;
			int totalspot = 0;
			int totalondemand = 0;
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
						//jobs.push_back(newdag);
						//printf("add new dag\n");
					}
				}

				int spotfail = 0;
				//step 0, check spotVM
				if(fmod(t, 60) == 0){ //check every minute
					int tracelag = rand()%35000;
					for(int i=0; i<types; i++)
						r[i] = datatrace[tracelag*types+i];
					//int mark = tr->readprice(r);
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
							if(!debug)
								check = function(sVMTP[i][j]->price, r[i]);//if not fail		
							else {
								//debug
								//use probability way to decide if the spot instance fails
								float lifetime = sVMTP[i][j]->life_time;								
								int rowindex = (int) (sVMTP[i][j]->price * 1000);
								int colindex = ceil(lifetime / 180.0);
								/*float succp, p1=0.0, pt=0.0;
								for(int tindex = 1; tindex <= colindex; tindex ++){
									p1 += firstfail[rowindex*440+i*110+tindex-1];
									pt += firstfail[rowindex*440+i*110+tindex-1]*tindex*180.0;
								}
								succp = 1 - p1;*/
								float ff = firstfail[rowindex*440+i*110+colindex-1];
								float rnd = (float) rand()/(RAND_MAX+1);
								if(rnd < ff) check = false;
								else check = true;
							}
							//if(check) sVMTP[i][j]->life_time = 10;
							if(!check){		//failed	
								spotfail += 1;
								totalspotfail += 1;
								if(sVMTP[i][j]->tk != NULL)	{
									sVMTP[i][j]->canAlloc = false;
									sVMTP[i][j]->tk->status = ready;
									sVMTP[i][j]->tk->readyCountdown = -1;
									//sVMTP[i][j]->tk->tasktime += t - sVMTP[i][j]->turn_on;//??

									bool isFound = false;
									while(!isFound) {
									//	sVMTP[i][j]->tk->vmID += 1;																	
										
										float pc = sVMTP[i][j]->tk->prices[sVMTP[i][j]->tk->vmID];										
										
										if(pc > 1e-12){
											isFound = true;
											int index = sVMTP[i][j]->tk->configList[sVMTP[i][j]->tk->vmID];
											sVMTP[i][j]->tk->restTime = sVMTP[i][j]->tk->probestTime[index*randomsize+globaliter];
										}
										else sVMTP[i][j]->tk->vmID += 1;																	
										
									}
									sVMTP[i][j]->tk = NULL;									
								}
							}
						}
					}
				}
				//step 1
				std::vector<taskVertex*> ready_task;
				for(int ji=0; ji<jobs.size(); ji++){
					vp = vertices((*jobs[ji]->g));
					for(int i=0; i < (*vp.second - *vp.first ); i++)
					{
						bool tag = true;
						//get parent vertices
						in_edge_iterator in_i, in_end;
						edge_descriptor e;
						boost::tie(in_i, in_end) = in_edges(i, (*jobs[ji]->g));
						if(in_i == in_end) tag = false;
						else {
							for (; in_i != in_end; ++in_i) 
							{
								e = *in_i;
								Vertex src = source(e, (*jobs[ji]->g));					
								if((*jobs[ji]->g)[src].status != finished)
								{
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

				int acqondemand = 0;
				int acqspot = 0;

				//sort according to the subdeadline, earliest deadline first
				//std::sort(ready_task.begin(),ready_task.end(), myfunction);

				for(int i=0; i<ready_task.size(); i++)//no preference on who gets resources first
				{
					taskVertex* curr_task=ready_task[i];
					if(curr_task->readyCountdown == -1)//
					{
						//find spot first
						bool suitfind = false;
						int suittp = 0;
						int suittpindex = 0;
						bool suitisSpot = true;
						
						int size = sVMTP[curr_task->configList[0]].size();
						for(int j=0; j<size && !suitfind; j++)
						{
							if(sVMTP[curr_task->configList[0]][j]->tk == NULL && sVMTP[curr_task->configList[0]][j]->life_time > curr_task->probestTime[curr_task->configList[0]*randomsize+globaliter] && sVMTP[curr_task->configList[0]][j]->canAlloc)
							{
								suittp = curr_task->configList[0];
								suittpindex = j;
								suitfind = true;									
							}
						}

						//find demand then
						if(!suitfind){
							suitisSpot = false;

							int size1 = VMTP[curr_task->configList[1]].size();
							for(int j=0; j<size1&& !suitfind; j++)
							{
								float runtime = VMTP[curr_task->configList[1]][j]->life_time;
								float availT = ceil(runtime/3600.0)*3600.0 - runtime;

								if(VMTP[curr_task->configList[1]][j]->tk == NULL)// && availT > curr_task->estTime[curr_task->configList[1]])
								{
									suittp = curr_task->configList[1];
									suittpindex = j;
									suitfind = true;										
								}
							}							
						}

						if(suitfind){
							curr_task->status = scheduled;
							curr_task->restTime = curr_task->probestTime[suittp*randomsize+globaliter] ;
							curr_task->readyCountdown = 0;

							if(suitisSpot){
								if(sVMTP[suittp][suittpindex]->has_data == curr_task->name)
									curr_task->restTime = curr_task->cpuTime[suittp] +  curr_task->netUp[suittp*randomsize+globaliter];
								else sVMTP[suittp][suittpindex]->has_data = curr_task->name;
								sVMTP[suittp][suittpindex]->tk = curr_task;
								curr_task->configList[curr_task->vmID] = suittp;
								if(curr_task->prices[curr_task->vmID]!=sVMTP[suittp][suittpindex]->price)
									printf("here");
								curr_task->prices[curr_task->vmID] = sVMTP[suittp][suittpindex]->price;
							}else{
								if(VMTP[suittp][suittpindex]->has_data == curr_task->name)
									curr_task->restTime = curr_task->cpuTime[suittp] +  curr_task->netUp[suittp*randomsize+globaliter];
								else VMTP[suittp][suittpindex]->has_data = curr_task->name;
								VMTP[suittp][suittpindex]->tk = curr_task;
								curr_task->vmID = 1; //spot+ondemand
								curr_task->configList[1] = suittp;
							}
						}else{
							bool isFound = false;
							while(!isFound){
								if(curr_task->prices[curr_task->vmID] > 1e-12){
									isFound = true;
								}else{
									curr_task->vmID += 1;
								}
							}

							int _config = curr_task->configList[curr_task->vmID];
							bool find = false;
							int vmindex = 0;
							//check VM/SpotVM list for available machine
							if(curr_task->vmID == 1)//the first is spot vm
							{
								int size = VMTP[_config].size();
								for(int j=0; j<size; j++)
								{
									if(VMTP[_config][j]->tk == NULL)
									{
										find = true;
										vmindex = j;
										VMTP[_config][j]->tk = curr_task;
										//if the needed input data is already on the instance
										curr_task->restTime = curr_task->probestTime[_config*randomsize+globaliter];
										if(VMTP[_config][j]->has_data == curr_task->name)
											curr_task->restTime = curr_task->cpuTime[_config] +  curr_task->netUp[_config*randomsize+globaliter];
										else VMTP[_config][j]->has_data = curr_task->name;
										break;
									}
								}
							}else{
								int size = sVMTP[_config].size();
								for(int j=0; j<size; j++) {
									if(sVMTP[_config][j]->tk == NULL && sVMTP[_config][j]->life_time > curr_task->probestTime[_config*randomsize+globaliter] && sVMTP[_config][j]->canAlloc)
									{
										find = true;
										sVMTP[_config][j]->tk = curr_task;
										//if the needed input data is already on the instance
										curr_task->restTime = curr_task->probestTime[_config*randomsize+globaliter];
										if(sVMTP[_config][j]->has_data == curr_task->name)
											curr_task->restTime = curr_task->cpuTime[_config] +  curr_task->netUp[_config*randomsize+globaliter];
										else sVMTP[_config][j]->has_data = curr_task->name;
										break;
									}
								}
							}
							if(find) {
								curr_task->status = scheduled;
								curr_task->taskstart = t;
								//int index = curr_task->configList[curr_task->vmID];
								//curr_task->restTime = curr_task->probestTime[index][globaliter];							
								curr_task->readyCountdown = 0;
							}
							else if(curr_task->vmID == 1) {curr_task->readyCountdown = OnDemandLag; curr_task->taskstart = t;}
							else {curr_task->readyCountdown = SpotLag; curr_task->taskstart = t;}
						}
					}
					else if(curr_task->readyCountdown == 0)
					{
						curr_task->status = scheduled;
						int index = curr_task->configList[curr_task->vmID];
						curr_task->restTime = curr_task->probestTime[index*randomsize+globaliter];

						if(curr_task->vmID == 1)//ondemand VM
						{
							VM* vm = new VM;
							vm->life_time = OnDemandLag;
							vm->tk = curr_task;
							vm->type = index;
							vm->turn_on = t;
							vm->has_data = curr_task->name;
							VMTP[index].push_back(vm);
							acqondemand += 1;
						}
						else if(curr_task->vmID == 0)
						{							
							if(r[index] > curr_task->prices[curr_task->vmID]){ //apply spotVM failed
								curr_task->status = ready;
								curr_task->vmID = curr_task->vmID + 1;

								bool isFound = false;
								while(!isFound){
									if(curr_task->prices[curr_task->vmID] > 1e-12){
										isFound = true;
									}else{
										curr_task->vmID += 1;
									}
								}

								if(curr_task->vmID == 1) curr_task->readyCountdown = OnDemandLag;
								else curr_task->readyCountdown = SpotLag;
							}else{
								SpotVM* svm = new SpotVM(curr_task->prices[curr_task->vmID]);								

								svm->tk = curr_task; //it's spot vm
								svm->type = index;
								svm->life_time = 3600; // - SpotLag
								svm->turn_on = t;
								svm->has_data = curr_task->name;
								sVMTP[index].push_back(svm);
								acqspot += 1;
							}
						}						
					}			
				}
				//delete VMs without task
				//int totalondemand = 0;
				//int totalspot = 0;
				int delondemand = 0;
				int delspot = 0;

				for(int i=0; i<types; i++)
				{
					int size1 = VMTP[i].size();
					int size2 = sVMTP[i].size();
					totalondemand += size1;
					totalspot += size2;
					
					for(int j=0; j<size1; j++)
					{
						float runtime = VMTP[i][j]->life_time;	
						if(VMTP[i][j]->tk == NULL&&(int)runtime%3600==0)
						{
							//moneycost += priceOnDemand[i]*runtime/60.0;
							moneycost += priceOnDemand[i]*ceil(runtime/3600.0);
							VM* vm = VMTP[i][j];
							delete vm;
							VMTP[i].erase(VMTP[i].begin()+j);
							j--;
							size1 --;
							delondemand += 1;
						}
					}
					
					for(int j=0; j<size2; j++) {
						float runtime = sVMTP[i][j]->life_time;
						if((sVMTP[i][j]->tk == NULL&&(int)runtime%3600==0) || (!sVMTP[i][j]->canAlloc)) {
							if(sVMTP[i][j]->tk == NULL && sVMTP[i][j]->canAlloc){
								moneycost += r[i];
							}

							SpotVM* svm = sVMTP[i][j];
							delete svm;
							sVMTP[i].erase(sVMTP[i].begin()+j);
						
							size2 --;
							j--;
							delspot += 1;
						}
					}
				}
				//step 2
				std::vector<taskVertex*> scheduled_task;
				for(int ji=0; ji<jobs.size(); ji++)
					for(int i=0; i<(*vp.second - *vp.first ); i++)
						if((*jobs[ji]->g)[i].status == scheduled)
							scheduled_task.push_back(&(*jobs[ji]->g)[i]);
				for(int i=0; i<scheduled_task.size(); i++) {
					scheduled_task[i]->restTime -= 1;////////////////////////////
					if(scheduled_task[i]->restTime <= 0) {
						scheduled_task[i]->status = finished;
						scheduled_task[i]->tasktime += t - scheduled_task[i]->taskstart;
						scheduled_task[i]->end_time = t;
						//deadline refinement and re-configuration
						//first check if the child task has been changed by other tasks
						/*out_edge_iterator out_i, out_end;
						Vertex v = scheduled_task[i]->name;
						for (boost::tie(out_i, out_end) = out_edges(v, dag.g); out_i != out_end; ++out_i) {
							//if(dag.g[target(*out_i,dag.g)].start_time < t || (scheduled_task[i]->dl - t)/scheduled_task[i]->dl > 0.2)
							//online tune
							if(dag.g[target(*out_i,dag.g)].start_time < t || (dag.g[target(*out_i,dag.g)].start_time - t)/t > 0.2)
							{
								dag.g[target(*out_i, dag.g)].start_time = t;
								dag.g[target(*out_i, dag.g)].instance_config();
							}
						}*/
						//make the vm.task = NULL
						int index = scheduled_task[i]->configList[scheduled_task[i]->vmID];
						if(scheduled_task[i]->vmID == 1) //VM type
						{
							for(int j=0; j<VMTP[index].size(); j++)
								if(VMTP[index][j]->tk == scheduled_task[i])
								{
									VMTP[index][j]->tk = NULL;
									break;
								}
						}else{
							for(int j=0; j<sVMTP[index].size(); j++)
								if(sVMTP[index][j]->tk == scheduled_task[i])
								{
									sVMTP[index][j]->tk = NULL;
									break;
								}
						}
					}	
				}
				//step 3
				for(int i=0; i<types; i++)
				{
					int size1 = VMTP[i].size();			
					
					for(int j=0; j<size1; j++)
					{
						VMTP[i][j]->life_time += 1;				
					}

					int size = sVMTP[i].size();
					for(int j=0; j<size; j++)
					{
						sVMTP[i][j]->life_time -= 1;//
						if(sVMTP[i][j]->life_time == 0){
							sVMTP[i][j]->life_time = 3600;
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
				for(int ji=0; ji<jobs.size(); ji++)
					for(int i=0; i < (*vp.second - *vp.first ); i++){
						if((*jobs[ji]->g)[i].status!= finished){
							condition = true;
							unfinishednum += 1;
					}					
				}
			//	if(!condition)
			//		printf("debug");
			}while(condition);
			//step 4 finalize
			for(int i=0; i<types; i++)
			{
				int size1 = VMTP[i].size();						
				for(int j=0; j<size1; j++)
				{
					float runtime = VMTP[i][j]->life_time;
					moneycost += (priceOnDemand[i] * ceil(runtime/3600.0));
				}

				int size = sVMTP[i].size();
				for(int j=0; j<size; j++)
				{				
					moneycost += r[i];				
				}
			}

			//tr->closefile();
			//delete tr;

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
	std::clock_t endtime = std::clock();
	float timeelapsed = (float)(endtime - starttime) / (float)CLOCKS_PER_SEC;
	printf("time elapsed for Dyna Algorithm is: %.4f\n", timeelapsed);		
}
