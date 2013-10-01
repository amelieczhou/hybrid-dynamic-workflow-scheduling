#include "stdafx.h"
#include "InstanceConfig.h"
#include "ReadTrace.h"
#include "PricingModel.h"
#include <boost/random/variate_generator.hpp>
#include <boost/random/mersenne_twister.hpp>
//#include <boost/math/distributions/gamma.hpp>
//#include <boost/math/distributions/normal.hpp>
#include <ctime>
#include <cmath>
#include <utility>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <string>
#include <stack>

using namespace boost;

extern boost::normal_distribution<> r_norm_m(150.28, 49.98);//switched small and medium, to be reasonable
extern boost::normal_distribution<> r_norm_s(120.31, 22.45);
extern boost::normal_distribution<> r_norm_l(172.85, 34.77);
extern boost::normal_distribution<> r_norm_x(1034.04, 146.41);
//sequential I/O, MBytes/sec
extern boost::gamma_distribution<> seq_io_m(129.28,0.792);//switched small and medium, to be reasonable
extern boost::gamma_distribution<> seq_io_s(127.14,0.802);
extern boost::gamma_distribution<> seq_io_l(376.57,0.281);
extern boost::gamma_distribution<> seq_io_x(408.11,0.264);
//network upload and download performance distribution
//ms per 8MB data transfer
extern boost::gamma_distribution<> gamma_s_up(2.077,604.38);
extern boost::gamma_distribution<> gamma_m_up(0.812, 895.37);
extern boost::gamma_distribution<> gamma_x_up(1.12, 528.68);//switched large and xlarge, to be reasonable
extern boost::gamma_distribution<> gamma_l_up(1.509, 435.16);
extern boost::gamma_distribution<> gamma_s_down(2.361, 399.9);
extern boost::gamma_distribution<> gamma_m_down(0.727, 654.3);
extern boost::gamma_distribution<> gamma_l_down(0.316, 126.66);
extern boost::gamma_distribution<> gamma_x_down(0.179, 127.87);



//after this function, each task in the dag is assigned with an instance type
//the dag is associated with a distribution of total execution time
//dfs search all the tasks first, then increase the instance type
void SearchPrune::OfflineSP_DFS(){ 
	//DAG dag(this->dag);
	
	//for the performance of each instance type
	double random_sequential_io[types][randomsize];
	double random_random_io[types][randomsize];
	double random_network_up[types][randomsize];
	double random_network_down[types][randomsize];
	
	//generate random number for sequential io
	variate_generator<mt19937,gamma_distribution<> > generator(mt19937(time(0)), seq_io_s);
	for(int j=0; j<randomsize; j++)
		random_sequential_io[0][j] = generator();
	variate_generator<mt19937,gamma_distribution<> > generator1(mt19937(time(0)), seq_io_m);
	for(int j=0; j<randomsize; j++)
		random_sequential_io[1][j] = generator1();
	variate_generator<mt19937,gamma_distribution<> > generator2(mt19937(time(0)), seq_io_l);
	for(int j=0; j<randomsize; j++)
		random_sequential_io[2][j] = generator2();
	variate_generator<mt19937,gamma_distribution<> > generator3(mt19937(time(0)), seq_io_x);
	for(int j=0; j<randomsize; j++)
		random_sequential_io[3][j] = generator3();
	
	//generate random number for random io
	variate_generator<mt19937,normal_distribution<> > generator4(mt19937(time(0)), r_norm_s);
	for(int j=0; j<randomsize; j++)
		random_random_io[0][j] = generator4();
	variate_generator<mt19937,normal_distribution<> > generator5(mt19937(time(0)), r_norm_m);	
	for(int j=0; j<randomsize; j++)
		random_random_io[1][j] = generator5();
	variate_generator<mt19937,normal_distribution<> > generator6(mt19937(time(0)), r_norm_l);
	for(int j=0; j<randomsize; j++)
		random_random_io[2][j] = generator6();
	variate_generator<mt19937,normal_distribution<> > generator7(mt19937(time(0)), r_norm_x);
	for(int j=0; j<randomsize; j++)
		random_random_io[3][j] = generator7();

	//generate random number for upload network
	variate_generator<mt19937,gamma_distribution<> > generator8(mt19937(time(0)), gamma_s_up);
	for(int j=0; j<randomsize; j++)
		random_network_up[0][j] = generator8();
	variate_generator<mt19937,gamma_distribution<> > generator9(mt19937(time(0)), gamma_m_up);
	for(int j=0; j<randomsize; j++)
		random_network_up[1][j] = generator9();
	variate_generator<mt19937,gamma_distribution<> > generator10(mt19937(time(0)), gamma_l_up);
	for(int j=0; j<randomsize; j++)
		random_network_up[2][j] = generator10();
	variate_generator<mt19937,gamma_distribution<> > generator11(mt19937(time(0)), gamma_x_up);
	for(int j=0; j<randomsize; j++)
		random_network_up[3][j] = generator11();

	//generate random number for download network
	variate_generator<mt19937,gamma_distribution<> > generator12(mt19937(time(0)), gamma_s_down);
	for(int j=0; j<randomsize; j++)
		random_network_down[0][j] = generator12();
	variate_generator<mt19937,gamma_distribution<> > generator13(mt19937(time(0)), gamma_m_down);
	for(int j=0; j<randomsize; j++)
		random_network_down[1][j] = generator13();
	variate_generator<mt19937,gamma_distribution<> > generator14(mt19937(time(0)), gamma_l_down);
	for(int j=0; j<randomsize; j++)
		random_network_down[2][j] = generator14();
	variate_generator<mt19937,gamma_distribution<> > generator15(mt19937(time(0)), gamma_x_down);
	for(int j=0; j<randomsize; j++)
		random_network_down[3][j] = generator15();

	std::pair<vertex_iter, vertex_iter> vp;
	vp = vertices(dag.g);
	int quantile = dag.meet_dl * randomsize;
	for(; vp.first != vp.second; vp.first++){
		//dag.g[*vp.first].probestTime = new double[types][randomsize];
		for(int t=0; t<types; t++){
			for(int j=0; j<randomsize; j++){
				dag.g[*vp.first].netUp[t][j] = dag.g[*vp.first].trans_data * random_network_up[t][j] / 8000;
				dag.g[*vp.first].netDown[t][j] = dag.g[*vp.first].rec_data * random_network_down[t][j] / 8000;
				dag.g[*vp.first].randomIO[t][j] = dag.g[*vp.first].read_data / random_random_io[t][j];
				dag.g[*vp.first].seqIO[t][j] = dag.g[*vp.first].seq_data / random_sequential_io[t][j];
				dag.g[*vp.first].probestTime[t][j] = dag.g[*vp.first].cpuTime[t] + dag.g[*vp.first].netUp[t][j]
					+ dag.g[*vp.first].netDown[t][j] + dag.g[*vp.first].randomIO[t][j] + dag.g[*vp.first].seqIO[t][j];
			}
			//calculate the estimate time as the expected value of the proestTime
			std::sort(std::begin(dag.g[*vp.first].probestTime[t]),std::end(dag.g[*vp.first].probestTime[t]));
			dag.g[*vp.first].estTime[t] = dag.g[*vp.first].probestTime[t][quantile];
			printf("task: %d, type: %d, time: %f\n",*vp.first,t,dag.g[*vp.first].estTime[t]);
		}
	}	
	
	//A* Search
	//each node in the search tree is in fact a DAG, their differences are the instance type assigned to each taask
	//DAGGraph optGraph;
	
	vp = vertices(dag.g);
	for(; vp.first != vp.second; vp.first++){
		dag.g[*vp.first].assigned_type = 0;//initially all assigned to small
	}
	//firstly, find a feasible solution and use it as lower bound
	std::stack<configstack> DAGstack;
	configstack initialstate;
	for(vp=vertices(dag.g); vp.first!=vp.second; vp.first++)
		initialstate.configurations.push_back(dag.g[*vp.first].assigned_type);
	initialstate.taskno = -1;
	bool colors[types];
	colors[0] = true;
	for(int i=1; i<types; i++) colors[i] = false;
	initialstate.childcolor = colors;
	vp=vertices(dag.g);
	int numoftasks = (*vp.second - *vp.first);
	double globalBestCost = 0;
	double exeTime[randomsize];
	//start the search along the optGraph
	DAGstack.push(initialstate);
	bool continuesearch = true;
	do{
		estimateTime(dag,exeTime); ///////////////////////////////////////////start from certain task, look up for the known part
		int count = 0;
		for(int i=0; i<randomsize; i++){
			if(exeTime[i]<=dag.deadline)
				count ++;
		}
		double ratio = (double)count / (double)randomsize;
		//if satisfy users' requirement, select as the lower bound
		if(ratio >= dag.meet_dl){
			globalBestCost = estimateCost(dag,0,false);
			continuesearch = false; //out of the while loop
			printf("initial ratio in search prune: %f\n",ratio);
		}else{
			int nexttask = DAGstack.top().taskno + 1;
			if(nexttask < numoftasks){
				configstack state;
				state.configurations = DAGstack.top().configurations;
				state.taskno = nexttask;
				for(int t=0; t<types; t++){
					if(!DAGstack.top().childcolor[t]){ //type t has not been visited
						state.configurations[nexttask] = t;
						DAGstack.top().childcolor[t] = true;
						state.childcolor = new bool[types];
						state.childcolor[0] = true;
						for(int tt=1; tt<types; tt++)
							state.childcolor[tt] = false;
						DAGstack.push(state);
						break;						
					}
					if(t == types-1){ //all types have been visited
						DAGstack.pop();
						//????
					}
				}
				for(int i=0; i<numoftasks; i++)
					dag.g[i].assigned_type = DAGstack.top().configurations[i];	
			}else{ //nexttask >= numoftasks				
				DAGstack.pop();
				//?????
			}
		}
		if(DAGstack.empty()) {
			printf("cannot find one valid solution\n");
			exit(1);
		}
	}while(continuesearch);

	//lower bound is found in DAGstack.top()
	//continue search with the lower bound
	std::vector<configstack> configsolutions;
	configsolutions.push_back(DAGstack.top());
	//configstack configsolution;
	int typecount = 0;
	int itercount = 0;
	do{
		int nexttask = DAGstack.top().taskno + 1;
		if(nexttask < numoftasks){
			configstack state;
			state.configurations = DAGstack.top().configurations;
			state.taskno = nexttask;
			for(int t=0; t<types; t++){
				if(!DAGstack.top().childcolor[t]){ //type t has not been visited
					state.configurations[nexttask] = t;
					DAGstack.top().childcolor[t] = true;
					state.childcolor = new bool[types];
					state.childcolor[0] = true;
					for(int tt=1; tt<types; tt++)
						state.childcolor[tt] = false;
					DAGstack.push(state);
					//update dag configuration
					for(int i=0; i<numoftasks; i++)
						dag.g[i].assigned_type = DAGstack.top().configurations[i];
					//if select as solution
					double currentcost = estimateCost(dag,0,false);
					if(currentcost < globalBestCost){
						estimateTime(dag,exeTime);
						int count = 0;
						for(int i=0; i<randomsize; i++){
							if(exeTime[i]<=dag.deadline)
								count ++;
						}
						double ratio = (double)count / (double)randomsize;
						if(ratio >= dag.meet_dl){
							globalBestCost = currentcost;
							//preserve the current configuration
							configsolutions.push_back(DAGstack.top());
							printf("found ratio in search prune: %f\n",ratio);
							//configsolution = DAGstack.top();
						}
					}else{ //no need to continue search since the rest must be more expensive
						//based on the assumption that adopting better instances must be more expensive
						DAGstack.pop();
						//jump to the next next task?????????????????????????????????????
						/*int nexttask = DAGstack.top().taskno+2;
						if(nexttask < numoftasks){
							configstack state;
							state.configurations = DAGstack.top().configurations;
							state.taskno = nexttask;
							for(int t=0; t<types; t++){

							}
						}*/
						for(int i=t; i<types; i++)
							DAGstack.top().childcolor[i] = true;
						break;
					}						
				}
				if(t == types-1){ //all types have been visited
					DAGstack.pop();
					for(int i=0; i<numoftasks; i++)
						dag.g[i].assigned_type = DAGstack.top().configurations[i];
				}					
			}
		}else{ //nexttask >= numoftasks
			DAGstack.pop();
			for(int i=0; i<numoftasks; i++)
				dag.g[i].assigned_type = DAGstack.top().configurations[i];
		}
		typecount = 0;
		for(int i=0; i<numoftasks; i++)
			if(dag.g[i].assigned_type == types-1)
				typecount += 1;
		itercount ++;
	}while(typecount != numoftasks);//(itercount < 10000000);

	//offline solution is found in configsolution
	for(int i=0; i<numoftasks; i++)
		dag.g[i].assigned_type = configsolutions.back().configurations[i];
		//dag.g[i].assigned_type = configsolution.configurations[i];
}
void SearchPrune::SpotTune(){
	//tune the offline search result with the spot instances
	//change to spot where can be changed
	//initially no spot instances
	std::pair<vertex_iter,vertex_iter> vp = vertices(dag.g);
	for(;vp.first != vp.second; vp.first ++){
		dag.g[*vp.first].configList = new int[2];
		dag.g[*vp.first].configList[0] = -2;
		dag.g[*vp.first].configList[1] = dag.g[*vp.first].assigned_type;
		//debug
		printf("task %d, config %d\n",*vp.first,dag.g[*vp.first].assigned_type);
		dag.g[*vp.first].prices = new double[2];
		dag.g[*vp.first].vmID = 0;
	}
	vp = vertices(dag.g);
	for(; vp.first != vp.second; vp.first++){
		std::vector<int*> cslot;
		int iter = 0;
		//if can reduce cost and still satisfy the deadline meet rate
		PricingModel* pm = new PricingModel();
		pm->init();
		std::vector<double*> bidps(types);
		for(int l=0; l<types; l++) //the spot instance type
		{			
			int config[3] = {l,dag.g[*vp.first].assigned_type,-1};
			double estT[3] = {dag.g[*vp.first].estTime[l],dag.g[*vp.first].estTime[dag.g[*vp.first].assigned_type],-1};
				
			bidps[l] = new double[3];//current cost, bidding price, ondemand price
			pm->getPricing(config,estT, bidps[l]);
			//bidps.push_back(bidp);
		}
		pm->finalize();
		delete pm;

		double minp=bidps[0][0], min_index=0;
		for(int j=1; j<types; j++)
			if(bidps[j][0]<minp) min_index = j;
		dag.g[*vp.first].configList[0] = min_index; 
		dag.g[*vp.first].prices[0] = bidps[min_index][1]; //may be 0
		dag.g[*vp.first].prices[1] = priceOnDemand[dag.g[*vp.first].assigned_type];
		double tmpmincost = bidps[min_index][0];

		//compare with demandonly
		int* tmpconfig = new int[3];
		tmpconfig[0] = -2;
		tmpconfig[1] = dag.g[*vp.first].assigned_type;
		tmpconfig[2] = -1;
			
		//TaskList[i]->restTime = TaskList[i]->estimateTime[TaskList[i]->configList[2]];		
		pm = new PricingModel();
		pm->init();						
		int demandconf = tmpconfig[1];

		double estT[3] = {dag.g[*vp.first].estTime[0],dag.g[*vp.first].estTime[demandconf],-1};
		double* bidp = new double[3];
		pm->getPricing(tmpconfig, estT, bidp);

		delete pm;

		if(bidp[0] <= tmpmincost){ //cost of ondemand only is less than with spot
			dag.g[*vp.first].configList[0] = -2;
			dag.g[*vp.first].prices[0] = 0;
		}else{
			double exeTime[randomsize];
			//check whether satisfy deadline///////////////////////
			estimateTimeSpot(dag,exeTime); 
			int count = 0;
			for(int i=0; i<randomsize; i++){
				if(exeTime[i]<dag.deadline)
					count ++;
			}
			double ratio = (double)count / (double)randomsize;
			if(ratio < dag.meet_dl){//does not satisfy deadline
				dag.g[*vp.first].configList[0] = -2;
				dag.g[*vp.first].prices[0] = 0;
			}
		}

		delete[] bidp; delete[] tmpconfig;
		int count = types;
		for(int i=0; i<count; i++)
		{			 
			delete[] bidps[i];
			bidps.erase(bidps.begin() +i);
			count--;
			i--;
		}		
	}
}
void SearchPrune::OnlineSimulate(){
	
	std::vector<DAG*> workflows; //continuous workflow
	double arrival_time = 0;
	DAG* job = new DAG(this->dag.deadline,this->dag.meet_dl);
	job->g =  this->dag.g; job->type = this->dag.type;
	job->arrival_time = 0;
	workflows.push_back(job);
	std::ifstream infile;
	//read in the arrival time of workflows from files
	std::string a = "arrivaltime_integer_", b, c = ".txt";
	std::ostringstream strlamda;
	strlamda << lambda;
	b = strlamda.str();
	std::string fname = a + b + c;
	infile.open(fname.c_str());
	char timing[256];
	if(infile==NULL){
		printf("cannot find input file!\n");
		return;
	}
	infile.getline(timing,256); //jump the lamda line
	infile.getline(timing,256); //jump the 0 line
	//incomming jobs, read the first num_jobs workflows
	while(workflows.size()<num_jobs){
		infile.getline(timing,256);
		arrival_time = atof(timing);
		arrival_time *= 60; //minutes to seconds

		DAG* job = new DAG(this->dag.deadline+arrival_time,this->dag.meet_dl);
		job->g = this->dag.g; job->type = this->dag.type;
		job->arrival_time = arrival_time;
		//std::pair<vertex_iter, vertex_iter> vp = vertices(job->g);
		workflows.push_back(job);
	}

	double ioseq[types],iorand[types],net_up[types],net_down[types];
	//start simulation
	std::clock_t starttime = std::clock();
	//for(int monte=0; monte < (num_monte+1); monte++)
	
	for(int globaliter=0; globaliter<randomsize; globaliter++)
	{
		//each task has a distribution of execution time
		//dag.g[0].probestTime
		//assign resources in what order?
		int tracelag = rand()%28800;
		std::vector<VM*> VMTP[types];
		std::vector<SpotVM*> sVMTP[types];
		std::pair<vertex_iter, vertex_iter> vp;

				
		if(dag.type == montage){
			for(int i=0; i<4; i++) {
				for(int ij=0; ij<workflows.size(); ij++){
					workflows[ij]->g[i].status = ready;
					workflows[ij]->g[i].readyCountdown = -1;
					workflows[ij]->g[i].restTime = 0;
				}
			}
			for(int i=4; i<20; i++){
				for(int ij=0; ij<workflows.size(); ij++){
					workflows[ij]->g[i].status = not_ready;
					workflows[ij]->g[i].readyCountdown = -1;
					workflows[ij]->g[i].restTime = 0;
				}
			}
		}else{
			printf("what is the dag type?");
			exit(1);
		}
		double t = 0; 
		bool condition = false; //there are at least one job not finished
		double moneycost = 0;
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
			for(int i=0; i<workflows.size(); i++) {
				//printf("i is:%d ",i);
				if(workflows[i]->arrival_time <= t){
					vp = vertices(workflows[i]->g);
					if(workflows[i]->g[*(vp.second-1)].status != finished)
						jobs.push_back(workflows[i]);
				}
				else { break;}
			}

			int spotfail = 0;
			//step 0, check spotVM
			if(fmod(t, 60) == 0){ //check every minute
				int mark = tr->readprice(r);
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
										
									double pc = sVMTP[i][j]->tk->prices[sVMTP[i][j]->tk->vmID];										
										
									if(pc > 1e-12){
										isFound = true;
										int index = sVMTP[i][j]->tk->configList[sVMTP[i][j]->tk->vmID];
										sVMTP[i][j]->tk->restTime = sVMTP[i][j]->tk->probestTime[index][globaliter];
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
					for(int ti = 0; ti<types && !suitfind; ti++)
					{
						int size = sVMTP[ti].size();
						for(int j=0; j<size && !suitfind; j++)
						{
							if(sVMTP[ti][j]->tk == NULL && sVMTP[ti][j]->life_time > curr_task->estTime[ti] && sVMTP[ti][j]->canAlloc)
							{
								suittp = ti;
								suittpindex = j;
								suitfind = true;									
							}
						}
					}

					//find demand then
					if(!suitfind){
						suitisSpot = false;

						for(int ti = 0; ti<types&& !suitfind; ti++){
							int size = VMTP[ti].size();
							for(int j=0; j<size&& !suitfind; j++)
							{
								double runtime = VMTP[ti][j]->life_time;
								double availT = ceil(runtime/3600.0)*3600.0 - runtime;

								if(VMTP[ti][j]->tk == NULL && availT > curr_task->estTime[ti])
								{
									suittp = ti;
									suittpindex = j;
									suitfind = true;										
								}
							}
						}
					}

					if(suitfind){
						curr_task->status = scheduled;
						curr_task->restTime = curr_task->probestTime[suittp][globaliter] ;
						curr_task->readyCountdown = 0;

						if(suitisSpot){
							if(sVMTP[suittp][suittpindex]->has_data == curr_task->name)
								curr_task->restTime = curr_task->cpuTime[suittp] +  curr_task->netUp[suittp][globaliter];
							else sVMTP[suittp][suittpindex]->has_data = curr_task->name;
							sVMTP[suittp][suittpindex]->tk = curr_task;
							curr_task->configList[curr_task->vmID] = suittp;
							curr_task->prices[curr_task->vmID] = sVMTP[suittp][suittpindex]->price;
						}else{
							if(VMTP[suittp][suittpindex]->has_data == curr_task->name)
								curr_task->restTime = curr_task->cpuTime[suittp] +  curr_task->netUp[suittp][globaliter];
							else VMTP[suittp][suittpindex]->has_data = curr_task->name;
							VMTP[suittp][suittpindex]->tk = curr_task;
							curr_task->vmID = 1; //spot+ondemand
							curr_task->configList[1] = suittp;
						}
					}else{
						//if(curr_task->vmID == 0){ //recalculate deadline
						//	curr_task->start_time = t;
						//	curr_task->instance_config();
						//}

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
						//check VM/SpotVM list for available machine
						if(curr_task->vmID == 1)//the first is spot vm
						{
							int size = VMTP[_config].size();
							for(int j=0; j<size; j++)
							{
								if(VMTP[_config][j]->tk == NULL)
								{
									find = true;
									VMTP[_config][j]->tk = curr_task;
									//if the needed input data is already on the instance
									curr_task->restTime = curr_task->probestTime[_config][globaliter];
									if(VMTP[_config][j]->has_data == curr_task->name)
										curr_task->restTime = curr_task->cpuTime[_config] +  curr_task->netUp[_config][globaliter];
									else VMTP[_config][j]->has_data = curr_task->name;
									break;
								}
							}
						}else{
							int size = sVMTP[_config].size();
							for(int j=0; j<size; j++) {
								if(sVMTP[_config][j]->tk == NULL && sVMTP[_config][j]->life_time > curr_task->estTime[_config] && sVMTP[_config][j]->canAlloc)
								{
									find = true;
									sVMTP[_config][j]->tk = curr_task;
									//if the needed input data is already on the instance
									curr_task->restTime = curr_task->probestTime[_config][globaliter];
									if(sVMTP[_config][j]->has_data == curr_task->name)
										curr_task->restTime = curr_task->cpuTime[_config] +  curr_task->netUp[_config][globaliter];
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
					curr_task->restTime = curr_task->probestTime[index][globaliter];

					if(curr_task->vmID == 1)//ondemand VM
					{
						VM* vm = new VM;
						vm->life_time = 0; //OnDemandLag
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
					if(VMTP[i][j]->tk == NULL)
					{
						double runtime = VMTP[i][j]->life_time;							
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
					if(sVMTP[i][j]->tk == NULL || (!sVMTP[i][j]->canAlloc)) {
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
					if(jobs[ji]->g[i].status == scheduled)
						scheduled_task.push_back(&jobs[ji]->g[i]);
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
					if(jobs[ji]->g[i].status!= finished){
						condition = true;
						unfinishednum += 1;
				}					
			}
			if(!condition)
				printf("debug");
		}while(condition);
		//step 4 finalize
		for(int i=0; i<types; i++)
		{
			int size1 = VMTP[i].size();						
			for(int j=0; j<size1; j++)
			{
				double runtime = VMTP[i][j]->life_time;
				moneycost += (priceOnDemand[i] * ceil(runtime/3600.0));
			}

			int size = sVMTP[i].size();
			for(int j=0; j<size; j++)
			{				
				moneycost += r[i];				
			}
		}

		tr->closefile();
		delete tr;

		printf("Money Cost: %.4f, Time: %.2f\n", moneycost, t);
		double averagetime = 0;
		double violation = 0;
		for(int i=0; i<workflows.size(); i++){
			vp = vertices(workflows[i]->g);
			double executiontime = workflows[i]->g[*(vp.second-1)].end_time - workflows[i]->arrival_time;
			averagetime += executiontime;
			if(executiontime > dag.deadline) violation += 1;
		}
		averagetime /= workflows.size(); violation /= workflows.size();
		printf("average execution time of workflows is %f, deadline violation is %f\n",averagetime,violation);
//		for(int i=0; i<(*vp.second - *vp.first ); i++)
//			printf("%.4f\t",dag.g[i].tasktime);
		//printf("\t task cost: ");
		//for(int i=0; i<(*vp.second - *vp.first ); i++)
		//	printf("%.4f\t",dag.g[i].cost);
		//printf("spot: %d\n", totalspot);
	}
	std::clock_t endtime = std::clock();
	double timeelapsed = (double)(endtime - starttime) / (double)CLOCKS_PER_SEC;
	printf("time elapsed for Dyna Algorithm is: %.4f\n", timeelapsed);		
}
void SearchPrune::OfflineSP_BFS(){
	//DAG* dag(this->dag);
	//DAG* dag = new DAG(this->dag);
	//for the performance of each instance type
	double random_sequential_io[types][randomsize];
	double random_random_io[types][randomsize];
	double random_network_up[types][randomsize];
	double random_network_down[types][randomsize];
	
	//generate random number for sequential io
	variate_generator<mt19937,gamma_distribution<> > generator(mt19937(time(0)), seq_io_s);
	for(int j=0; j<randomsize; j++)
		random_sequential_io[0][j] = generator();
	variate_generator<mt19937,gamma_distribution<> > generator1(mt19937(time(0)), seq_io_m);
	for(int j=0; j<randomsize; j++)
		random_sequential_io[1][j] = generator1();
	variate_generator<mt19937,gamma_distribution<> > generator2(mt19937(time(0)), seq_io_l);
	for(int j=0; j<randomsize; j++)
		random_sequential_io[2][j] = generator2();
	variate_generator<mt19937,gamma_distribution<> > generator3(mt19937(time(0)), seq_io_x);
	for(int j=0; j<randomsize; j++)
		random_sequential_io[3][j] = generator3();
	
	//generate random number for random io
	variate_generator<mt19937,normal_distribution<> > generator4(mt19937(time(0)), r_norm_s);
	for(int j=0; j<randomsize; j++)
		random_random_io[0][j] = generator4();
	variate_generator<mt19937,normal_distribution<> > generator5(mt19937(time(0)), r_norm_m);	
	for(int j=0; j<randomsize; j++)
		random_random_io[1][j] = generator5();
	variate_generator<mt19937,normal_distribution<> > generator6(mt19937(time(0)), r_norm_l);
	for(int j=0; j<randomsize; j++)
		random_random_io[2][j] = generator6();
	variate_generator<mt19937,normal_distribution<> > generator7(mt19937(time(0)), r_norm_x);
	for(int j=0; j<randomsize; j++)
		random_random_io[3][j] = generator7();

	//generate random number for upload network
	variate_generator<mt19937,gamma_distribution<> > generator8(mt19937(time(0)), gamma_s_up);
	for(int j=0; j<randomsize; j++)
		random_network_up[0][j] = generator8();
	variate_generator<mt19937,gamma_distribution<> > generator9(mt19937(time(0)), gamma_m_up);
	for(int j=0; j<randomsize; j++)
		random_network_up[1][j] = generator9();
	variate_generator<mt19937,gamma_distribution<> > generator10(mt19937(time(0)), gamma_l_up);
	for(int j=0; j<randomsize; j++)
		random_network_up[2][j] = generator10();
	variate_generator<mt19937,gamma_distribution<> > generator11(mt19937(time(0)), gamma_x_up);
	for(int j=0; j<randomsize; j++)
		random_network_up[3][j] = generator11();

	//generate random number for download network
	variate_generator<mt19937,gamma_distribution<> > generator12(mt19937(time(0)), gamma_s_down);
	for(int j=0; j<randomsize; j++)
		random_network_down[0][j] = generator12();
	variate_generator<mt19937,gamma_distribution<> > generator13(mt19937(time(0)), gamma_m_down);
	for(int j=0; j<randomsize; j++)
		random_network_down[1][j] = generator13();
	variate_generator<mt19937,gamma_distribution<> > generator14(mt19937(time(0)), gamma_l_down);
	for(int j=0; j<randomsize; j++)
		random_network_down[2][j] = generator14();
	variate_generator<mt19937,gamma_distribution<> > generator15(mt19937(time(0)), gamma_x_down);
	for(int j=0; j<randomsize; j++)
		random_network_down[3][j] = generator15();

	std::pair<vertex_iter, vertex_iter> vp;
	vp = vertices(dag.g);
	int quantile = dag.meet_dl * randomsize;
	for(; vp.first != vp.second; vp.first++){
		//dag.g[*vp.first].probestTime = new double[types][randomsize];
		for(int t=0; t<types; t++){
			for(int j=0; j<randomsize; j++){
				dag.g[*vp.first].netUp[t][j] = dag.g[*vp.first].trans_data * random_network_up[t][j] / 8000;
				dag.g[*vp.first].netDown[t][j] = dag.g[*vp.first].rec_data * random_network_down[t][j] / 8000;
				dag.g[*vp.first].randomIO[t][j] = dag.g[*vp.first].read_data / random_random_io[t][j];
				dag.g[*vp.first].seqIO[t][j] = dag.g[*vp.first].seq_data / random_sequential_io[t][j];
				dag.g[*vp.first].probestTime[t][j] = dag.g[*vp.first].cpuTime[t] + dag.g[*vp.first].netUp[t][j]
					+ dag.g[*vp.first].netDown[t][j] + dag.g[*vp.first].randomIO[t][j] + dag.g[*vp.first].seqIO[t][j];
			}
			//calculate the estimate time as the expected value of the proestTime
			std::sort(std::begin(dag.g[*vp.first].probestTime[t]),std::end(dag.g[*vp.first].probestTime[t]));
			dag.g[*vp.first].estTime[t] = dag.g[*vp.first].probestTime[t][quantile];
			printf("task: %d, type: %d, time: %f\n",*vp.first,t,dag.g[*vp.first].estTime[t]);
		}
	}	
	
	//A* Search
	//each node in the search tree is in fact a DAG, their differences are the instance type assigned to each taask
	//DAGGraph optGraph;
	
	vp = vertices(dag.g);
	for(; vp.first != vp.second; vp.first++){
		dag.g[*vp.first].assigned_type = 0;//initially all assigned to small
	}
	//firstly, find a feasible solution and use it as lower bound
	std::stack<configstack> DAGstack;
	configstack initialstate;
	for(vp=vertices(dag.g); vp.first!=vp.second; vp.first++)
		initialstate.configurations.push_back(dag.g[*vp.first].assigned_type);
	initialstate.taskno = 0;
	bool colors[types];
	colors[0] = true;
	for(int i=1; i<types; i++) colors[i] = false;
	initialstate.childcolor = colors;
	vp=vertices(dag.g);
	int numoftasks = (*vp.second - *vp.first);
	double globalBestCost = 0;
	double exeTime[randomsize];
	double totalTime = 0;
	//start the search along the optGraph
	DAGstack.push(initialstate);
	bool continuesearch = true;
	do{
		estimateTime(dag,exeTime); ///////////////////////////////////////////start from certain task, look up for the known part
		int count = 0;
		for(int i=0; i<randomsize; i++){
			if(exeTime[i]<=dag.deadline)
				count ++;
		}
		double ratio = (double)count / (double)randomsize;
		//totalTime = estimateTime2(dag);
		//if(totalTime <= dag.deadline){	
		//if satisfy users' requirement, select as the lower bound
		if(ratio >= dag.meet_dl){
			globalBestCost = estimateCost(dag,0,false);
			continuesearch = false; //out of the while loop
			printf("initial ratio in search prune: %f\n",ratio);
		}else{
			int nexttype = types;
			for(int i=0; i<types; i++){
				if(!DAGstack.top().childcolor[i]){
					nexttype = i;
					break;
				}
			}			
			if(nexttype < types){
				//go to the next type
				configstack state;
				state.configurations = DAGstack.top().configurations;					
				state.taskno = DAGstack.top().taskno;
				state.configurations[state.taskno] = nexttype;
				state.childcolor = new bool[types];									
				state.childcolor[0] = true;
				for(int tt=1; tt<types; tt++)
					state.childcolor[tt] = DAGstack.top().childcolor[tt];
				state.childcolor[nexttype] = true;	
				DAGstack.push(state);	
			}else{//go for the next task
				configstack state;
				state.configurations = DAGstack.top().configurations;					
				state.taskno = DAGstack.top().taskno + 1;
				if(state.taskno == numoftasks){
					printf("cannot find a solution that satisfy deadline!\n");
					exit(1);
				}
				state.configurations[state.taskno] = 1; //search from the smallest
				state.childcolor = new bool[types];									
				state.childcolor[0] = state.childcolor[1] = true;
				for(int tt=2; tt<types; tt++)
					state.childcolor[tt] = false;					
				DAGstack.push(state);	
			}
			for(int i=0; i<numoftasks; i++)
				dag.g[i].assigned_type = DAGstack.top().configurations[i];
		}
	}while(continuesearch);
	//lower bound is found in DAGstack.top()
	//continue search with the lower bound
	std::vector<configstack> configsolutions;
	configsolutions.push_back(DAGstack.top());
	int typecount = 0;
	do{
		//continue search for better solutions
		int nexttype = types;
		for(int i=0; i<types; i++){
			if(!DAGstack.top().childcolor[i]){
				nexttype = i;
				break;
			}
		}			
		if(nexttype < types){
			//go to the next type
			configstack state;
			state.configurations = DAGstack.top().configurations;					
			state.taskno = DAGstack.top().taskno;
			state.configurations[state.taskno] = nexttype;
			state.childcolor = new bool[types];									
			state.childcolor[0] = true;
			for(int tt=1; tt<types; tt++)
				state.childcolor[tt] = DAGstack.top().childcolor[tt];
			state.childcolor[nexttype] = true;	
			DAGstack.push(state);	
		}else{//go for the next task
			configstack state;
			state.configurations = DAGstack.top().configurations;					
			state.taskno = DAGstack.top().taskno + 1;
			if(state.taskno == numoftasks){
				//printf("cannot find a solution that satisfy deadline!\n");
				//exit(1);
				break;
			}
			state.configurations[state.taskno] = 1; //search from the smallest
			state.childcolor = new bool[types];									
			state.childcolor[0] = state.childcolor[1] = true;
			for(int tt=2; tt<types; tt++)
				state.childcolor[tt] = false;					
			DAGstack.push(state);	
		}
		for(int i=0; i<numoftasks; i++)
			dag.g[i].assigned_type = DAGstack.top().configurations[i];
		double currentcost = estimateCost(dag,0,false);
		if(currentcost < globalBestCost){
			estimateTime(dag,exeTime);
			//totalTime = estimateTime2(dag);
			int count = 0;
			for(int i=0; i<randomsize; i++){
				if(exeTime[i]<=dag.deadline)
					count ++;
			}
			double ratio = (double)count / (double)randomsize;
			//if(totalTime <= dag.deadline){
			if(ratio >= dag.meet_dl){
				globalBestCost = currentcost;
				//preserve the current configuration
				configsolutions.push_back(DAGstack.top());
				printf("found ratio in search prune: %f\n",ratio);
				//configsolution = DAGstack.top();
			}
		}else{ //no need to continue search since the rest must be more expensive
			//based on the assumption that adopting better instances must be more expensive!!!!!!!
			//DAGstack.pop();
			//go to the next task
			for(int i=0; i<types; i++)
				DAGstack.top().childcolor[i] = true;
		}
		/*for(int i=0; i<numoftasks; i++)
			if(dag.g[i].assigned_type == types-1)
				typecount += 1;*/
	}while(typecount!=numoftasks);
}