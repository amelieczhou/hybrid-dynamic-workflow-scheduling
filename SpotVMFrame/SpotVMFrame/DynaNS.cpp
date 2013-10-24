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

void DynaNS::Simulate(){	
	float arrival_time = 0;
	dag->arrival_time = 0;	
	float violation = 0;
	float ave_cost = 0;
	float ioseq[types],iorand[types],net_up[types],net_down[types];
	//start simulation
	std::clock_t starttime = std::clock();
	omp_set_num_threads(24);
	float viol_private[24];
	float cost_private[24];
	for(int i=0; i<24; i++) viol_private[i]=cost_private[i]=0;
	#pragma omp parallel
	{
		#pragma omp for
		for(int globaliter=0; globaliter<randomsize; globaliter++)
		{
			//copy dag for each thread to simulate separately
			DAG* newdag = new DAG(*dag);
			//each task has a distribution of execution time			
			//assign resources in what order?
			int tracelag = rand()%28800;
			std::vector<VM*> VMTP[types];
			std::vector<SpotVM*> sVMTP[types];
			std::pair<vertex_iter, vertex_iter> vp;

				
			if(newdag->type == montage){
				for(int i=0; i<4; i++) {
					newdag->g[i].status = ready;
					newdag->g[i].readyCountdown = -1;
					newdag->g[i].restTime = 0;				
				}
				for(int i=4; i<20; i++){
					newdag->g[i].status = not_ready;
					newdag->g[i].readyCountdown = -1;
					newdag->g[i].restTime = 0;	
				}
			}else if(newdag->type == ligo){
				for(int i=0; i<9; i++) {
					newdag->g[i].status = ready;
					newdag->g[i].readyCountdown = -1;
					newdag->g[i].restTime = 0;				
				}
				for(int i=9; i<40; i++){
					newdag->g[i].status = not_ready;
					newdag->g[i].readyCountdown = -1;
					newdag->g[i].restTime = 0;	
				}
			}else if(newdag->type == epigenome){				
				newdag->g[0].status = ready;
				newdag->g[0].readyCountdown = -1;
				newdag->g[0].restTime = 0;	
				for(int i=1; i<20; i++){
					newdag->g[i].status = not_ready;
					newdag->g[i].readyCountdown = -1;
					newdag->g[i].restTime = 0;	
				}
			}else{
				printf("what is the dag type?");
				exit(1);
			}
			float t = 0; 
			bool condition = false; //there are at least one job not finished
			float moneycost = 0;
			////trace data
			ReadTrace* tr = new ReadTrace();
			tr->readfile("data.csv", tracelag);//tracelag any use?

			float r[4];
			int totalspotfail =0;
			int totalspot = 0;
			int totalondemand = 0;
			do{
				//accept workflows
				std::vector<DAG*> jobs;
				jobs.push_back(newdag);

				//step 1
				std::vector<taskVertex*> ready_task;
				for(int ji=0; ji<jobs.size(); ji++){
					vp = vertices(jobs[ji]->g);
					for(int i=0; i < (*vp.second - *vp.first ); i++)
					{
						bool tag = true;
						//get parent vertices
						in_edge_iterator in_i, in_end;
						edge_descriptor e;
						for (boost::tie(in_i, in_end) = in_edges(i, jobs[ji]->g); in_i != in_end; ++in_i) 
						{
							e = *in_i;
							Vertex src = source(e, jobs[ji]->g);					
							if(jobs[ji]->g[src].status != finished)
							{
								tag = false;
								//break;
							}
						}
						if(jobs[ji]->g[i].status == ready || tag && jobs[ji]->g[i].status != scheduled && jobs[ji]->g[i].status != finished){
							ready_task.push_back(&jobs[ji]->g[i]);							
						}
					}
				}			
				//sort according to the subdeadline, earliest deadline first
				//std::sort(ready_task.begin(),ready_task.end(), myfunction);

				for(int i=0; i<ready_task.size(); i++)//no preference on who gets resources first
				{
					taskVertex* curr_task=ready_task[i];
					if(curr_task->readyCountdown == -1)//
					{
						int _config = curr_task->assigned_type;
						bool find = false;
						//check VM/SpotVM list for available machine
						int size = VMTP[_config].size();
						int vmindex = 0;
						for(int j=0; j<size; j++)
						{
							if(VMTP[_config][j]->tk == NULL)
							{
								find = true;
								vmindex = j;
								VMTP[_config][j]->tk = curr_task;
								break;
							}
						}
						if(find) {
							curr_task->status = scheduled;
							curr_task->taskstart = t;
							curr_task->restTime = curr_task->probestTime[_config*randomsize+globaliter];

							if(VMTP[_config][vmindex]->has_data == curr_task->name)
								curr_task->restTime = curr_task->cpuTime[_config] +  curr_task->netUp[_config*randomsize+globaliter];
							else VMTP[_config][vmindex]->has_data = curr_task->name;
						}
						else 			
						{
							curr_task->readyCountdown = OnDemandLag;
							curr_task->taskstart = t;
						}
					}
					else if(curr_task->readyCountdown == 0)
					{
						curr_task->status = scheduled;
						curr_task->restTime = curr_task->probestTime[curr_task->assigned_type*randomsize+globaliter];

						VM* vm = new VM; 
						vm->life_time = OnDemandLag;
						vm->tk = curr_task;
						vm->type = curr_task->assigned_type;
						VMTP[curr_task->config].push_back(vm);						
					}
				}			
							
				//delete VMs without task
				int delondemand = 0;
				for(int i=0; i<types; i++)
				{
					int size1 = VMTP[i].size();
					totalondemand += size1;
					
					for(int j=0; j<size1; j++)
					{
						if(VMTP[i][j]->tk == NULL)
						{
							float runtime = VMTP[i][j]->life_time;							
							moneycost += priceOnDemand[i]*ceil(runtime/3600.0);
							VM* vm = VMTP[i][j];
							delete vm;
							VMTP[i].erase(VMTP[i].begin()+j);
							j--;
							size1 --;
							delondemand += 1;
						}
					}
				}
				//step 2
				std::vector<taskVertex*> scheduled_task;
				for(int ji=0; ji<jobs.size(); ji++)
					for(int i=0; i<(*vp.second - *vp.first ); i++)
						if(jobs[ji]->g[i].status == scheduled)
							scheduled_task.push_back(&jobs[ji]->g[i]);
				for(int i=0; i<scheduled_task.size(); i++) {
					scheduled_task[i]->restTime -= 1;////////////////////////////
					if(scheduled_task[i]->restTime <= 0) {
						scheduled_task[i]->status = finished;
						scheduled_task[i]->tasktime += t - scheduled_task[i]->taskstart;
						scheduled_task[i]->end_time = t;
						//make the vm.task = NULL
						int index = scheduled_task[i]->assigned_type;
						for(int j=0; j<VMTP[index].size(); j++)
							if(VMTP[index][j]->tk == scheduled_task[i])
							{
								VMTP[index][j]->tk = NULL;
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
				for(int ji=0; ji<jobs.size(); ji++)
					for(int i=0; i < (*vp.second - *vp.first ); i++){
						if(jobs[ji]->g[i].status!= finished){
							condition = true;
							unfinishednum += 1;
					}					
				}
				//if(!condition)
				//	printf("debug");
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
			}

			tr->closefile();
			delete tr;

			printf("Money Cost: %.4f, Time: %.2f\n", moneycost, t);
			vp = vertices(newdag->g);
			float executiontime = newdag->g[*(vp.second-1)].end_time - newdag->arrival_time;
			int id = omp_get_thread_num();
			printf("thread id is %d\n",id);
			if(executiontime > newdag->deadline) {
				viol_private[id] += 1.0;
			}		
			cost_private[id] += moneycost;
			printf("average execution time of workflows is %f\n",executiontime);
		}
	}
	for(int i=0; i<24; i++) {
		violation += viol_private[i];
		ave_cost += cost_private[i];
	}

	violation /= (float)randomsize;
	ave_cost /= (float)randomsize;
	printf("deadline meeting rate is %.2f, average cost is %.2f\n",1.0-violation,ave_cost);
	std::clock_t endtime = std::clock();
	float timeelapsed = (float)(endtime - starttime) / (float)CLOCKS_PER_SEC;
	printf("time elapsed for Dyna Algorithm is: %.4f\n", timeelapsed);		
}