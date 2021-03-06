#include "stdafx.h"
#include "Algorithms.h"
#include "PricingModel.h"
#include <ctime>
#include <fstream>
#include <sstream>
#include <string>
#include <omp.h>

using namespace boost;
extern float lambda;
extern int num_jobs;
void DynaNS::Simulate(){	

	std::vector<DAG*> workflows; //continuous workflow
	float arrival_time = 0;
	dag->arrival_time = 0;
	workflows.push_back(dag);

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
		vp = vertices((*job->g));
		workflows.push_back(job);
	}

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
			//DAG* newdag = new DAG(*dag);
			//each task has a distribution of execution time			
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

			float r[4];
			int totalspotfail =0;
			int totalspot = 0;
			int totalondemand = 0;
			do{
				//accept workflows
				for(int i=0; i<jobs.size(); i++){
					if((int)t == (int)jobs[i]->arrival_time){
						if(jobs[i]->type == montage){
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
						}else if(jobs[i]->type == ligo){
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
						}else if(jobs[i]->type == epigenome){				
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
						if(in_i == in_end)
							tag = false;
						else {
							for (; in_i != in_end; ++in_i) 
							{
								e = *in_i;
								Vertex src = source(e, (*jobs[ji]->g));					
								if((*jobs[ji]->g)[src].status != finished)
								{
									tag = false;
									break;
								}
							}
						}
						if((*jobs[ji]->g)[i].status == ready || (tag && (*jobs[ji]->g)[i].status != scheduled && (*jobs[ji]->g)[i].status != finished)){
							(*jobs[ji]->g)[i].status = ready;
							ready_task.push_back(&(*jobs[ji]->g)[i]);							
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

							/*if(VMTP[_config][vmindex]->has_data == curr_task->name)
								curr_task->restTime = curr_task->cpuTime[_config] +  curr_task->netUp[_config*randomsize+globaliter];
							else VMTP[_config][vmindex]->has_data = curr_task->name;*/
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
						VMTP[curr_task->assigned_type].push_back(vm);						
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
						float runtime = VMTP[i][j]->life_time;
						if(VMTP[i][j]->tk == NULL&& ((int)ceil(runtime)%60==0))
						{
														
							moneycost += priceOnDemand[i]*ceil(runtime/60.0);
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
				for(int ji=0; ji<jobs.size(); ji++){
					vp = vertices((*jobs[ji]->g));
					for(int i=0; i<(*vp.second - *vp.first ); i++)
						if((*jobs[ji]->g)[i].status == scheduled)
							scheduled_task.push_back(&(*jobs[ji]->g)[i]);
				}
				
				for(int i=0; i<scheduled_task.size(); i++)
				{
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
				for(int ji=0; ji<jobs.size(); ji++){
					vp = vertices((*jobs[ji]->g));
					for(int i=0; i < (*vp.second - *vp.first ); i++){
						if((*jobs[ji]->g)[i].status!= finished){
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
					float runtime = VMTP[i][j]->life_time;
					moneycost += (priceOnDemand[i] * ceil(runtime/60.0));
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
	printf("deadline meeting rate is %.2f, average cost is %.2f\n",1.0-violation,ave_cost);
	std::clock_t endtime = std::clock();
	float timeelapsed = (float)(endtime - starttime) / (float)CLOCKS_PER_SEC;
	printf("time elapsed for Dyna Algorithm is: %.4f\n", timeelapsed);		
}
