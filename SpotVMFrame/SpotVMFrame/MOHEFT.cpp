#include "stdafx.h"
#include "Algorithms.h"
#include "PricingModel.h"
#include <ctime>
#include <fstream>
#include <sstream>
#include <string>
#include <omp.h>
#include <assert.h>

const int solutionsize = 10;

float cal_dist(DAG* dag, int taskid){
	float dist = 0.0;
	out_edge_iterator out_i, out_end;
	float maxdist = 0.0;
	int maxchild = -1;
	boost::tie(out_i, out_end) = out_edges(taskid, *dag->g);
	if(out_i != out_end){
		for (; out_i != out_end; ++out_i) {
			Vertex v = target(*out_i,*dag->g);
			float vdist = cal_dist(dag,v);
			if(vdist > maxdist){
				maxdist = vdist;
				maxchild = v;
			}
		}
		dist = maxdist + (*dag->g)[taskid].estTime[0];
	}else{
		dist = (*dag->g)[taskid].estTime[0];
	}
	return dist;
}

bool costfunction(taskVertex* a, taskVertex* b){
	return (a->cost > b->cost);
}

bool distfunction(individual* a, individual* b){//cannot be <=, if not strict weak, assertion failure!
	return (a->objectives.first < b->objectives.first);
}
bool crowdfunction(individual* a, individual* b){
	return (a->fitness > b->fitness);
}

void MOHEFT::evalObjectives(individual* individ, int task){
	
	//the first is cost fitness and the second is time fitness	
	for(int i=0; i<=task; i++){
		(*dag->g)[i].assigned_type = individ->solution[i];
	}
	float costobjective = estimateCost(*dag,0,task,false);
	float* exeTime = (float*)malloc(randomsize*sizeof(float));
	float timeobjective = 0.0; 
	estimateTime(*dag,0,task,exeTime);
	for(int i=0; i<randomsize; i++)
		timeobjective += exeTime[i];
	timeobjective /= randomsize;

	individ->objectives.first = costobjective;
	individ->objectives.second = timeobjective;

	return;
}

//sort the tasks according to the distance of the task to the end of the workflow
std::vector<individual* > MOHEFT::BRank(){
	//calculate the distance according to the cheapest instance type
	std::pair<vertex_iter,vertex_iter> vp = vertices(*dag->g);
	int numtasks = (*vp.second - *vp.first);
	std::vector<taskVertex*> tasks;
	for(; vp.first != vp.second; vp.first++){
		(*dag->g)[*vp.first].cost = cal_dist(dag,*vp.first);
		tasks.push_back(&(*dag->g)[*vp.first]);
	}
	std::sort(tasks.begin(),tasks.end(),costfunction);
	
	//create solutionsize empty workflow schedules
	std::vector<individual* > solutions;
	//iterate over the ranked tasks
	for(int i=0; i<tasks.size(); i++){
		std::vector<individual* > tempsolutions;
		//iterate over all resources
		for(int j=0; j<types; j++){
			for(int k=0; k<solutionsize; k++){
				//0 to i-1 instance from solution k, add type j for i-th task
				individual* tempsolution = new individual();
				tempsolution->solution = (int*)malloc(numtasks*sizeof(int));
				for(int iter=0; iter<i; iter++){
					tempsolution->solution[iter] =solutions[k]->solution[iter];
				}
				tempsolution->solution[i] = j;
				tempsolutions.push_back(tempsolution);
			}				
		}
		//sort according to crowding distance
		for(int j=0; j<tempsolutions.size(); j++){
			//calculate fitness
			evalObjectives(tempsolutions[j],i);
		}
		//select top k for task i
		std::vector<individual*> nondominatedsolutions;
		for(int j=0; j<tempsolutions.size(); j++){
			for(int k=0; k<tempsolutions.size(); k++){
				if(dominates(tempsolutions[k],tempsolutions[j]))
					tempsolutions[j]->fitness += 1.0; //dominated
			}
			if(tempsolutions[j]->fitness < 1e-12){
				nondominatedsolutions.push_back(tempsolutions[j]);
			}
		}
		if(nondominatedsolutions.size() <= solutionsize)
			solutions = nondominatedsolutions;
		else{
			//sort according to one objective
			std::sort(nondominatedsolutions.begin(),nondominatedsolutions.end(),distfunction);
			int nondominatedsize = nondominatedsolutions.size();
			for(int j=0; j<nondominatedsize; j++){
				//normalize the objectives
				float normfirst = nondominatedsolutions[0]->objectives.first;
				float normsecond = nondominatedsolutions[nondominatedsize-1]->objectives.second;
				nondominatedsolutions[j]->objectives.first /= normfirst;
				nondominatedsolutions[j]->objectives.second /= normsecond;
			}
			nondominatedsolutions[0]->fitness = nondominatedsolutions[nondominatedsize-1]->fitness = 100.0;
			for(int j=1; j<nondominatedsize-1; j++){ //the first and the last are preserved
				nondominatedsolutions[j]->fitness = calcDistance(nondominatedsolutions[j-1],nondominatedsolutions[j+1]); //the distance
			}
			//select k out of the nondominated solutions
			std::sort(nondominatedsolutions.begin(),nondominatedsolutions.end(),crowdfunction);		
			//clear
			solutions.erase(solutions.begin(),solutions.end());
			for(int k=0; k<solutionsize; k++)
				solutions.push_back(nondominatedsolutions[k]);
		}
	}

	return solutions;
}

void MOHEFT::Simulate(){

	float arrival_time = 0;
	dag->arrival_time = 0;	
	float violation = 0;
	float ave_cost = 0;
	float ioseq[types],iorand[types],net_up[types],net_down[types];
	omp_set_num_threads(24);
	float viol_private[24];
	float cost_private[24];
	for(int i=0; i<24; i++) viol_private[i]=cost_private[i]=0;

	//for the performance of each instance type
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
	vp = vertices((*dag->g));

	int quantile = dag->meet_dl * randomsize;
	for(; vp.first != vp.second; vp.first++){
		//dag->g[*vp.first].probestTime = new float[types][randomsize];
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
	}

	free(random_sequential_io);
	free(random_random_io);
	free(random_network_up);
	free(random_network_down);

	//start the evolution process to obtain solutions
	time_t start,end;
	time(&start);		
	std::vector<individual* > nondominatedsolutions = BRank();
	time(&end);
	printf("optimization overhead of moheft is %.4f\n",difftime(end,start));

	//use the  solution to simulate
	for(int i=0; i<nondominatedsolutions.size(); i++){
            printf("the %d-th solution found:\t",i);
            for(int j=0; j<(*dag->g).m_vertices.size(); j++){
				printf("%d,",nondominatedsolutions[i]->solution[j]);
            }
    }
	
	for(int nondominateindex=0; nondominateindex<solutionsize; nondominateindex++){
		vp = vertices(*dag->g);
		for(;vp.first != vp.second; vp.first++){
			(*dag->g)[*vp.first].prefer_type = nondominatedsolutions[nondominateindex]->solution[*vp.first]; //try different solutions in the nondominated set?
			printf("task %d: %d\n",*vp.first,(*dag->g)[*vp.first ].prefer_type);
		}
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
		printf("time elapsed for moheft algorithm is: %.4f\n", timeelapsed);
	}	

	return;
}
