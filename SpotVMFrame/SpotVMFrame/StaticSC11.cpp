#include "stdafx.h"
#include "Algorithms.h"
#include "PricingModel.h"
#include <boost/random/variate_generator.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <ctime>
#include <time.h>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <string>
#include <omp.h>

using namespace boost;

extern float lambda;

void Autoscaling::Initialize(){
	
	std::pair<vertex_iter, vertex_iter> vp;
	vp = vertices(*dag->g);
	for(; vp.first!=vp.second; ++vp.first) {//edge weight for communication cost		
		Vertex v1 = *vp.first;
		(*dag->g)[v1].prefer_type = types-1;
	}
			
	dag->deadline_assign();
	//task configuration, find the prefered VM type for tasks
	vp = vertices(*dag->g);
	for(; vp.first != vp.second; ++vp.first)
		(*dag->g)[*vp.first].instance_config();
}
void Autoscaling::Simulate(){
	float arrival_time = 0;
	dag->arrival_time = 0;	
	float violation = 0;
	float ave_cost = 0;
	float ioseq[types],iorand[types],net_up[types],net_down[types];
	omp_set_num_threads(24);
	float viol_private[24];
	float cost_private[24];
	for(int i=0; i<24; i++) viol_private[i]=cost_private[i]=0;

	//prepare for runtime
	float* random_sequential_io = (float*)malloc(types*randomsize*sizeof(float));
	float* random_random_io = (float*)malloc(types*randomsize*sizeof(float));
	float* random_network_up = (float*)malloc(types*randomsize*sizeof(float));
	float* random_network_down = (float*)malloc(types*randomsize*sizeof(float));
		float* random_tmp = (float*)malloc(types*10000*sizeof(float));
	//read from file	
	FILE* rFile;
	char str[1024];
	char buf[256];
	char *ptr, *ptr2;
	rFile = fopen("randio.csv","r");
	if(rFile == NULL){
		printf("cannot open randio.csv\n");
		exit(1);
	}
	for(int i=0; i<types*10000; i++){
		if(fgets(str,1024,rFile)!=NULL)
			random_tmp[i] = atof(str);
	}	
	for(int i=0; i<types; i++){
		for(int j=0; j<randomsize; j++){
			random_random_io[i*randomsize+j] = random_tmp[i*10000+j];//10000 is fixed
		}
	}

	rFile = fopen("seqio.csv","r");
	if(rFile == NULL){
		printf("cannot open seqio.csv\n");
		exit(1);
	}
	for(int i=0; i<types*10000; i++){
		if(fgets(str,1024,rFile)!=NULL)
			random_tmp[i] = atof(str);
	}	
	for(int i=0; i<types; i++){
		for(int j=0; j<randomsize; j++){
			random_sequential_io[i*randomsize+j] = random_tmp[i*10000+j];//10000 is fixed
		}
	}
	rFile = fopen("netup.csv","r");
	if(rFile == NULL){
		printf("cannot open netup.csv\n");
		exit(1);
	}
	for(int i=0; i<types*10000; i++){
		if(fgets(str,1024,rFile)!=NULL)
			random_tmp[i] = atof(str);
	}	
	for(int i=0; i<types; i++){
		for(int j=0; j<randomsize; j++){
			random_network_up[i*randomsize+j] = random_tmp[i*10000+j];//10000 is fixed
		}
	}
	rFile = fopen("netdown.csv","r");
	if(rFile == NULL){
		printf("cannot open netdown.csv\n");
		exit(1);
	}
	for(int i=0; i<types*10000; i++){
		if(fgets(str,1024,rFile)!=NULL)
			random_tmp[i] = atof(str);
	}	
	for(int i=0; i<types; i++){
		for(int j=0; j<randomsize; j++){
			random_network_down[i*randomsize+j] = random_tmp[i*10000+j];//10000 is fixed
		}
	}
	free(random_tmp);

	std::pair<vertex_iter, vertex_iter> vp;
	property_map<Graph, edge_weight_t>::type weightmap = get(edge_weight, *dag->g);
	vp = vertices(*dag->g);

	int quantile = dag->meet_dl * randomsize;;//0.5;
	for(; vp.first != vp.second; vp.first++){
		//dag.g[*vp.first].probestTime = new float[types][randomsize];
		for(int t=0; t<types; t++){
			for(int j=0; j<randomsize; j++){
				(*dag->g)[*vp.first].netUp[t*randomsize+j] = (*dag->g)[*vp.first].trans_data * random_network_up[t*randomsize+j] / 8000;
				(*dag->g)[*vp.first].netDown[t*randomsize+j] = (*dag->g)[*vp.first].rec_data * random_network_down[t*randomsize+j] / 8000;
				(*dag->g)[*vp.first].randomIO[t*randomsize+j] = (*dag->g)[*vp.first].read_data / random_random_io[t*randomsize+j];
				(*dag->g)[*vp.first].seqIO[t*randomsize+j] = (*dag->g)[*vp.first].seq_data / random_sequential_io[t*randomsize+j];
				(*dag->g)[*vp.first].probestTime[t*randomsize+j] = (*dag->g)[*vp.first].cpuTime[t] + (*dag->g)[*vp.first].netUp[t*randomsize+j]
					+ (*dag->g)[*vp.first].netDown[t*randomsize+j] + (*dag->g)[*vp.first].randomIO[t*randomsize+j] + (*dag->g)[*vp.first].seqIO[t*randomsize+j];
			}
			//calculate the estimate time as the expected value of the proestTime
			std::sort((*dag->g)[*vp.first].probestTime+t*randomsize,(*dag->g)[*vp.first].probestTime+(t+1)*randomsize-1);
			(*dag->g)[*vp.first].estTime[t] = (*dag->g)[*vp.first].probestTime[t*randomsize+quantile];
			printf("task: %d, type: %d, time: %f\n",*vp.first,t,(*dag->g)[*vp.first].estTime[t]);
		}
		//in_edge_iterator in_i, in_end;
		//for (boost::tie(in_i, in_end) = in_edges(*vp.first, dag->g); in_i != in_end; ++in_i){
		//	edge_descriptor e = *in_i;
		//	weightmap[e] = -1* (dag->g[*vp.first].estTime[dag->g[*vp.first].config]);//deadline assign using cheapest machine
		//}
	}	
	free(random_sequential_io);
	free(random_random_io);
	free(random_network_up);
	free(random_network_down);

	//initialization, do deadline assignment and instance config
	time_t start,end;
	time(&start);
	Initialize();
	time(&end);
	printf("optimization overhead of static is %.4f\n",difftime(end,start));
	vp = vertices(*dag->g);
	for(;vp.first != vp.second; vp.first++)
		printf("task %d: %d\n",*vp.first,(*dag->g)[*vp.first ].prefer_type);
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
		for(int i=0; i<(*vp.second - *vp.first); i++)
			(*job->g)[i].sub_deadline += arrival_time;
		workflows.push_back(job);
	}
	infile.close();
	//start simulation
	std::clock_t starttime = std::clock();
	#pragma omp parallel
	{
		#pragma omp for
		for(int monte=0; monte < randomsize; monte++)
		{
			/*for(int i=0; i<types; i++){
				ioseq[i] = random_sequential_io[randomsize*i+monte];
				iorand[i] = random_random_io[randomsize*i+monte];
				net_up[i] = random_network_up[randomsize*i+monte];
				net_down[i] = random_network_down[randomsize*i+monte];
			}*/		
			std::vector<DAG*> jobs;
			for(int i=0; i<workflows.size(); i++){
				DAG* newdag = new DAG(*workflows[i]);
				vp = vertices(*newdag->g);
				for(int j=0; j<(*vp.second - *vp.first); j++)
					(*newdag->g)[j].sub_deadline = (*workflows[i]->g)[j].sub_deadline;
				jobs.push_back(newdag);
			}
			std::vector<VM*> VMTP[types];
			int need_VM[types]={0,0,0,0};

			//EDF scheduling
			double t = 0;
			bool condition = false;
			double moneycost = 0.0;
		
			do{	
				//accept workflows
				for(int i=0; i<jobs.size(); i++){
					if((int)t == (int)jobs[i]->arrival_time){
						if(dag->type == montage){
							for(int j=1; j<5; j++) {
								(*jobs[i]->g)[j].status = ready;
								(*jobs[i]->g)[j].readyCountdown = -1;
								(*jobs[i]->g)[j].restTime = 0;				
							}
							for(int j=5; j<21; j++){
								(*jobs[i]->g)[j].status = not_ready;
								(*jobs[i]->g)[j].readyCountdown = -1;
								(*jobs[i]->g)[j].restTime = 0;	
							}
							(*jobs[i]->g)[0].status = (*jobs[i]->g)[21].status = finished;
						}else if(dag->type == ligo){
							for(int j=1; j<10; j++) {
								(*jobs[i]->g)[j].status = ready;
								(*jobs[i]->g)[j].readyCountdown = -1;
								(*jobs[i]->g)[j].restTime = 0;				
							}
							for(int j=10; j<41; j++){
								(*jobs[i]->g)[j].status = not_ready;
								(*jobs[i]->g)[j].readyCountdown = -1;
								(*jobs[i]->g)[j].restTime = 0;	
							}
							(*jobs[i]->g)[0].status = (*jobs[i]->g)[41].status = finished;
						}else if(dag->type == epigenome){				
							(*jobs[i]->g)[1].status = ready;
							(*jobs[i]->g)[1].readyCountdown = -1;
							(*jobs[i]->g)[1].restTime = 0;	
							for(int j=2; j<21; j++){
								(*jobs[i]->g)[j].status = not_ready;
								(*jobs[i]->g)[j].readyCountdown = -1;
								(*jobs[i]->g)[j].restTime = 0;	
							}
							(*jobs[i]->g)[0].status = (*jobs[i]->g)[21].status = finished;
						}else{
							printf("what is the dag type?");
							exit(1);
						}
					//	jobs.push_back(newdag);
					//	printf("add new dag\n");
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
							for (; in_i != in_end; ++in_i) 
							{
								e = *in_i;
								Vertex src = source(e, *jobs[ji]->g);					
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


				std::sort(ready_task.begin(),ready_task.end(), myfunction);
				for(int i=0; i<ready_task.size(); i++)//earliest deadline first
				{
					taskVertex* curr_task=ready_task[i];
					if(curr_task->readyCountdown == -1)//
					{
						
						int _config = curr_task->prefer_type;
						bool find = false;
						//check VM/SpotVM list for available machine
						int size = VMTP[_config].size();
						for(int j=0; j<size; j++)
						{
							if(VMTP[_config][j]->tk == NULL)
							{
								find = true;
								VMTP[_config][j]->tk = curr_task;
								break;
							}
						}
						if(find) {
							curr_task->status = scheduled;
							curr_task->tasktime = t;
							curr_task->restTime =  curr_task->probestTime[curr_task->prefer_type*randomsize+monte] ;
						}
						else 			
						{
							curr_task->readyCountdown = OnDemandLag;
							curr_task->tasktime = t;
						}
					}
					else if(curr_task->readyCountdown == 0)
					{
						curr_task->status = scheduled;
						curr_task->restTime = curr_task->probestTime[curr_task->prefer_type*randomsize+monte] ;

						VM* vm = new VM; 
						vm->life_time = OnDemandLag;
						vm->tk = curr_task;
						vm->type = curr_task->prefer_type;
						VMTP[curr_task->prefer_type].push_back(vm);
						
					}			
				}
				//delete VMs without task
				for(int i=0; i<types; i++)//////////
				{
					int size1 = VMTP[i].size();
					
					for(int j=0; j<size1; j++)
					{
						if(VMTP[i][j]->tk == NULL)
						{
							double runtime = VMTP[i][j]->life_time;
							moneycost += (priceOnDemand[i] * ceil(runtime/3600.0));

							VM* vm = VMTP[i][j];
							delete vm;
							VMTP[i].erase(VMTP[i].begin()+j);
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
						scheduled_task[i]->cost = scheduled_task[i]->tasktime * priceOnDemand[scheduled_task[i]->prefer_type] /3600.0;
						//make the vm.task = NULL
						for(int j=0; j<VMTP[scheduled_task[i]->prefer_type].size(); j++)
							if(VMTP[scheduled_task[i]->prefer_type][j]->tk == scheduled_task[i])
							{
								VMTP[scheduled_task[i]->prefer_type][j]->tk = NULL;
								break;
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

			for(int i=0; i<types; i++)
			{
				int size1 = VMTP[i].size();						
				for(int j=0; j<size1; j++)
				{
					double runtime = VMTP[i][j]->life_time;
					moneycost += (priceOnDemand[i] * ceil(runtime/3600.0));
				}
			}
			printf("Money Cost: %.4f, Time: %.2f\n", moneycost, t);
			int id = omp_get_thread_num();
			printf("thread id is %d\n",id);
			float ave_time = 0.0;
			for(int i=0; i<jobs.size(); i++){
				vp = vertices((*jobs[i]->g));
				float executiontime = (*jobs[i]->g)[*(vp.second-2)].end_time - jobs[i]->arrival_time;
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
		violation += viol_private[i];
		ave_cost += cost_private[i];
	}
	violation /= (float)randomsize*num_jobs;
	ave_cost /= (float)randomsize*num_jobs;
	printf("deadline meeting rate is %f, average cost is %f\n",1.0-violation,ave_cost);
	std::clock_t endtime = std::clock();
	std::clock_t timeelapsed = (float)(endtime - starttime) / (float)CLOCKS_PER_SEC;
	printf("time elapsed for SC11 algorithm is: %.4f\n", timeelapsed);
}
