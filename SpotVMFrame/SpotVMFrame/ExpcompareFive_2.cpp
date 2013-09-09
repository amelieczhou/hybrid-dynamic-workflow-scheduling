// Deadline_assign.cpp : Defines the entry point for the console application.
// Input: binary tree type DAG
//compare to previous version: add network and I/O operation time to all algorithms
//change estTime in order to make the instance configuration adapt to performance variance
//after discussion on 6th, Oct add runtime deadline refinement
#include "stdafx.h"
#include <stdlib.h>
#include <vector>
#include <utility>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <time.h>
#include "ReadTrace.h"
#include "PricingModel.h"
#include "InstanceConfig.h"
#include <boost/graph/topological_sort.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/gamma.hpp>

using namespace boost;

bool NoSpotVM = true;
double priceOnDemand[] = {0.095, 0.19, 0.38, 0.76};
double Times[4][types] = {{120,65,38,24},{90,50,30,20},{60,35,23,17},{30,20,15,13}};

//I/O performance distribution, seek/sec
math::normal r_norm_s(150.28, 49.98);
math::normal r_norm_m(120.31, 22.45);
math::normal r_norm_l(172.85, 34.77);
math::normal r_norm_x(1034.04, 146.41);
//sequential I/O, MBytes/sec
math::gamma_distribution<> seq_io_s(129.28,0.792);
math::gamma_distribution<> seq_io_m(127.14,0.802);
math::gamma_distribution<> seq_io_l(376.57,0.281);
math::gamma_distribution<> seq_io_x(408.11,0.264);
//network upload and download performance distribution
math::gamma_distribution<> gamma_s_up(2.077,604.38);
math::gamma_distribution<> gamma_m_up(0.812, 895.37);
math::gamma_distribution<> gamma_l_up(1.12, 528.68);
math::gamma_distribution<> gamma_x_up(1.509, 435.16);
math::gamma_distribution<> gamma_s_down(2.361, 399.9);
math::gamma_distribution<> gamma_m_down(0.727, 654.3);
math::gamma_distribution<> gamma_l_down(0.316, 126.66);
math::gamma_distribution<> gamma_x_down(0.179, 127.87);


int main(int argc, char** argv)
{	
	//construct the DAG
	srand( (unsigned)time( NULL ) );
	int depth = atoi(argv[2]);
	int width = atoi(argv[3]);
	double deadline = atof(argv[4]);
	double budget = atof(argv[5]);
	int num_monte = atoi(argv[6]);
	double cpu_perc = atof(argv[7]);
	double io_perc = atof(argv[8]);
	double network_perc = atof(argv[9]);
	double rand_io_ratio = atof(argv[10]);
	double exe_time = atof(argv[11]); //The parallel parameter used to modify Times[][]
	double meet_dl = atof(argv[13]); //the deadline meet rate submitted by the user

	for(int i=0; i<4; i++)
		for(int j=1; j<types; j++)
			Times[i][j]=Times[i][0]*exe_time/pow(2.0,j)+Times[i][0]*(1-exe_time);

	//I/O speed, seeks/sec
	double ioseq[4], iorand[4];
	//network speed, sec/8M
	double net_up[4], net_down[4];
	double r = (double) rand()/RAND_MAX;
	if(strcmp (argv[12], "sens") == 0) ////for senstivity test, fix the DAG type and task types
		r = 0.5;
	ioseq[0] = math::quantile(seq_io_s, r);
	iorand[0] = math::quantile(r_norm_s, r);
	net_up[0] = math::quantile(gamma_s_up, r);
	//Graph dag;
	DAG dag(deadline);
	if(strcmp(argv[1], "pipeline") == 0)
	{
		int rnds[] = {0,3,3,1,2,1,0};
		//generate DAG
		for(int i=0; i<depth; i++)
		{
			taskVertex tk;
			tk.name = i;
			tk.mark = tk.config = 0;
			tk.cost = tk.tasktime = tk.taskstart = 0;
			tk.start_time = tk.end_time = tk.dl = tk.restTime = 0;
			tk.readyCountdown = -1;
			tk.status = not_ready;
			//int lv[types]={0,0,0,0};
			//tk.LV = lv;
			tk.configList = new int[3];
			tk.prices = new double[3];
			int rnd = (double)rand() / RAND_MAX * types; //[range_min, range_max)	
			//int rnd = rn_integers(0, types-1);
			
			if(strcmp (argv[12], "sens") == 0)
				rnd = rnds[i];
			printf("rntypes:%d\n",rnd);
			tk.type = (int)rnd;
			tk.estTime = new double[types];
			tk.actTime = new double [types];
			double averT = 0;
			tk.read_data = Times[tk.type][0] /cpu_perc *io_perc *rand_io_ratio *iorand[0];
			tk.seq_data = Times[tk.type][0] /cpu_perc *io_perc *(1-rand_io_ratio) *ioseq[0];
			tk.trans_data = 8 * Times[tk.type][0] /cpu_perc *network_perc /net_up[0];
			tk.rec_data = 0;
			add_vertex(tk, dag.g);
		}
		
		//add edges to the pipeline graph
		std::pair<vertex_iter, vertex_iter> vp; 
		vp = vertices(dag.g);
		for(; vp.first != (vp.second-1); )
		{
			Vertex v1 = *vp.first;
			vp.first++;
			Vertex v2 = *vp.first;

			edge_descriptor e; 
			bool inserted;
			tie(e, inserted) = add_edge(v1, v2, dag.g);
		}				
	}
	//for parallel application
	if(strcmp(argv[1], "parallel") == 0) //(depth-3)/2 is the depth of each parallel part
	{
		//generate DAG
		for(int i=0; i<(3+(depth-3)*width); i++)
		{
			taskVertex tk;
			tk.name = i;
			tk.mark = tk.config = 0;
			tk.cost = tk.tasktime = tk.taskstart = 0;
			tk.start_time = tk.end_time = tk.dl = tk.restTime = 0;
			tk.readyCountdown = -1;
			tk.status = not_ready;
			tk.configList = new int[3];
			tk.prices = new double[3];
			int rnd = (double)rand() / RAND_MAX * types; //[range_min, range_max)	
			//int rnd = rn_integers(0, types-1);
			printf("rntypes:%d\n",rnd);
			tk.type = (int)rnd;
			tk.estTime = new double[types];
			tk.actTime = new double[types];
			double averT = 0;
			tk.read_data = Times[tk.type][types-1] /cpu_perc *io_perc *rand_io_ratio *iorand[0];
			tk.seq_data = Times[tk.type][types-1] /cpu_perc *io_perc *(1-rand_io_ratio) *ioseq[0];
			tk.trans_data = 8 * Times[tk.type][types-1] /cpu_perc *network_perc /net_up[0];
			tk.rec_data = 0;
			add_vertex(tk, dag.g);
		}

		//add edges to the parallel graph
		std::pair<vertex_iter, vertex_iter> vp; 
		vp = vertices(dag.g);

		int i=0 , j=(1+(depth-3)/2*width);
		for(int k=1; k<=width; k++)
		{
			add_edge(i, i+k, dag.g );
			add_edge(j, j+k, dag.g );
		}
		int firstdep = (depth -3)/2;
		int secdep = depth-3-firstdep;
		for(int i=1; i<=std::ceil(firstdep/2.0)*width; i++)
			add_edge(i, i+width, dag.g );
		for(int i=1+std::ceil(firstdep/2.0)*width; i<=firstdep*width; i++)
			add_edge(i, firstdep*width+1, dag.g );
		for(int i=firstdep*width+2; i<=(firstdep*width+1+std::ceil(secdep/2.0)*width); i++)
			add_edge(i, i+width, dag.g );
		for(int i=firstdep*width+1+std::ceil(secdep/2.0)*width+1; i<(vp.second -vp.first-1); i++)
			add_edge(i, vp.second-vp.first-1, dag.g );
		
	}
	
	//for hybrid application
	if(strcmp(argv[1], "hybrid") == 0) //fixed shape as indicated in the paper
	{
		int rnds[] = {0,3,3,1,2,1,0,1,2,2,3,3,0,1,2};
		//generate DAG
		for(int i=0; i<15; i++)
		{
			taskVertex tk;
			tk.name = i;
			tk.mark = tk.config = 0;
			tk.cost = tk.tasktime = tk.taskstart = 0;
			tk.start_time = tk.end_time = tk.dl = tk.restTime = 0;
			tk.readyCountdown = -1;
			tk.status = not_ready;
			tk.configList = new int[3];
			tk.prices = new double[3];
			int rnd = (double)rand() / RAND_MAX * types; //[range_min, range_max)	
			//int rnd = rn_integers(0, types-1);

			if(strcmp (argv[12], "sens") == 0)
				rnd = rnds[i];
			printf("rntypes:%d\n",rnd);
			tk.type = (int)rnd;
			tk.estTime = new double[types];
			tk.actTime = new double[types];
			double averT = 0;
			tk.read_data = Times[tk.type][0] /cpu_perc *io_perc *rand_io_ratio *iorand[0];
			tk.seq_data = Times[tk.type][0] /cpu_perc *io_perc *(1-rand_io_ratio) *ioseq[0];
			tk.trans_data = 8 * Times[tk.type][0] /cpu_perc *network_perc /net_up[0];
			tk.rec_data = 0;
			add_vertex(tk, dag.g);
		}
		
		//add edges to the hybrid graph
		std::pair<vertex_iter, vertex_iter> vp; 
		vp = vertices(dag.g);

		int p[4] = {5, 7, 10, 12};
		for(int i=0; i<4; i++)
			add_edge(p[i], p[i]+1, dag.g);
		int start[5] = {1, 3, 4, 5, 7};
		for(int i=0; i<5; i++)
			add_edge(0, start[i], dag.g);
		add_edge(6, 9, dag.g);
		add_edge(8, 9, dag.g);
		add_edge(9, 10, dag.g);
		add_edge(9, 12, dag.g);
		add_edge(1, 2, dag.g);
		add_edge(2, 14, dag.g);
		add_edge(3, 11, dag.g);
		add_edge(4, 10, dag.g);
		add_edge(11, 14, dag.g);
		add_edge(13, 14, dag.g);
	}

	//suitable for all kinds of DAG
	std::pair<vertex_iter, vertex_iter> vp; 
	vp = vertices(dag.g);
	for(; vp.first != vp.second ; ++vp.first )
	{
		in_edge_iterator in_i, in_end;
		edge_descriptor e;
		for (boost::tie(in_i, in_end) = in_edges(*vp.first , dag.g); in_i != in_end; ++in_i) 
		{
			e = *in_i;
			Vertex src = source(e, dag.g);
			dag.g[*vp.first ].rec_data += dag.g[src].trans_data;			
		}
		printf("receive data for node %d: %f\n", *vp.first, dag.g[*vp.first].rec_data);
	}

	printf("----------------------------------start Dyna algorithm------------------------------------------\n");

	double ave_ioseq[4];
	double ave_iorand[4];
	double ave_netup[4];
	double ave_netdown[4];
	double static_rnd = 0.5;
	ave_ioseq[0] = math::quantile(seq_io_s, static_rnd);
	ave_ioseq[1] = math::quantile(seq_io_m, static_rnd);
	ave_ioseq[2] = math::quantile(seq_io_l, static_rnd);
	ave_ioseq[3] = math::quantile(seq_io_x, static_rnd);
	ave_iorand[0] = math::quantile(r_norm_s, static_rnd);
	ave_iorand[1] = math::quantile(r_norm_m, static_rnd);
	ave_iorand[2] = math::quantile(r_norm_l, static_rnd);
	ave_iorand[3] = math::quantile(r_norm_x, static_rnd);
	ave_netup[0] = math::quantile(gamma_s_up, static_rnd);
	ave_netup[1] = math::quantile(gamma_m_up, static_rnd);
	ave_netup[2] = math::quantile(gamma_l_up, static_rnd);
	ave_netup[3] = math::quantile(gamma_x_up, static_rnd);
	ave_netdown[0] = math::quantile(gamma_s_down, static_rnd);
	ave_netdown[1] = math::quantile(gamma_m_down, static_rnd);
	ave_netdown[2] = math::quantile(gamma_l_down, static_rnd);
	ave_netdown[3] = math::quantile(gamma_x_down, static_rnd);

	double ioseq_req[4];
	double iorand_req[4];
	double net_up_req[4];
	double net_down_req[4];
	ioseq_req[0] = math::quantile(seq_io_s, meet_dl);
	ioseq_req[1] = math::quantile(seq_io_m, meet_dl);
	ioseq_req[2] = math::quantile(seq_io_l, meet_dl);
	ioseq_req[3] = math::quantile(seq_io_x, meet_dl);
	iorand_req[0] = math::quantile(r_norm_s, meet_dl);
	iorand_req[1] = math::quantile(r_norm_m, meet_dl);
	iorand_req[2] = math::quantile(r_norm_l, meet_dl);
	iorand_req[3] = math::quantile(r_norm_x, meet_dl);
	net_up_req[0] = math::quantile(gamma_s_up, meet_dl);
	net_up_req[1] = math::quantile(gamma_m_up, meet_dl);
	net_up_req[2] = math::quantile(gamma_l_up, meet_dl);
	net_up_req[3] = math::quantile(gamma_x_up, meet_dl);
	net_down_req[0] = math::quantile(gamma_s_down, meet_dl);
	net_down_req[1] = math::quantile(gamma_m_down, meet_dl);
	net_down_req[2] = math::quantile(gamma_l_down, meet_dl);
	net_down_req[3] = math::quantile(gamma_x_down, meet_dl);

	std::clock_t starttime = std::clock();
	// monte carlo execution
	for(int monte=0; monte < (num_monte+1); monte++)
	{		
		double rnd = (double)rand() / RAND_MAX;
		while(rnd<0.0001)
			rnd = (double)rand() / RAND_MAX;
		//double rnd = rn_01();
		printf("rnd:%f\t",rnd);

		ioseq[0] = math::quantile(seq_io_s, rnd);
		ioseq[1] = math::quantile(seq_io_m, rnd);
		ioseq[2] = math::quantile(seq_io_l, rnd);
		ioseq[3] = math::quantile(seq_io_x, rnd);
		iorand[0] = math::quantile(r_norm_s, rnd);
		iorand[1] = math::quantile(r_norm_m, rnd);
		iorand[2] = math::quantile(r_norm_l, rnd);
		iorand[3] = math::quantile(r_norm_x, rnd);
		net_up[0] = math::quantile(gamma_s_up, rnd);
		net_up[1] = math::quantile(gamma_m_up, rnd);
		net_up[2] = math::quantile(gamma_l_up, rnd);
		net_up[3] = math::quantile(gamma_x_up, rnd);
		net_down[0] = math::quantile(gamma_s_down, rnd);
		net_down[1] = math::quantile(gamma_m_down, rnd);
		net_down[2] = math::quantile(gamma_l_down, rnd);
		net_down[3] = math::quantile(gamma_x_down, rnd);

		//before deadline assign
		dag.reset();

		vp = vertices(dag.g);
		property_map<Graph, edge_weight_t>::type weightmap = get(edge_weight, dag.g);
		for(; vp.first != vp.second; ++vp.first)
		{
			Vertex v1 = *vp.first;
			for(int t=0; t<types; t++)
			{
				dag.g[v1].estTime[t] = Times[dag.g[v1].type][t] + dag.g[v1].read_data / iorand_req[t] + dag.g[v1].seq_data / ioseq_req[t] + dag.g[v1].trans_data * net_up_req[t]/8 + dag.g[v1].rec_data *net_down_req[t] /8;
				printf("t=%d : %f\t",t,dag.g[v1].estTime[t]);
				dag.g[v1].actTime[t] = Times[dag.g[v1].type][t] + dag.g[v1].read_data / iorand[t] + dag.g[v1].seq_data / ioseq[t] + dag.g[v1].trans_data * net_up[t]/8 + dag.g[v1].rec_data *net_down[t] /8;
			}
			in_edge_iterator in_i, in_end;
			for (boost::tie(in_i, in_end) = in_edges(v1, dag.g); in_i != in_end; ++in_i)
			{
				edge_descriptor e = *in_i;
				weightmap[e] = -1* (dag.g[v1].estTime[0]);//deadline assign using cheapest machine
			}
		}

		dag.deadline_assign();

		for(vp = vertices(dag.g); vp.first != vp.second; vp.first ++)
			dag.g[*vp.first].instance_config();

		int tracelag = rand()%28800;
		std::vector<VM*> VMTP[types];
		std::vector<SpotVM*> sVMTP[types];
		/*int need_VM[types]={0,0,0,0};*/

		//EDF scheduling
		vp = vertices(dag.g);
		dag.g[*vp.first ].status = ready;

		double t = 0;
		bool condition = false;
		double moneycost = 0.0;

		////trace data
		ReadTrace* tr = new ReadTrace();
		tr->readfile("data.csv", tracelag);

		float r[4];
		int totalspotfail =0;
		int totalspot = 0;
		do{		
			int spotfail = 0;
			//step 0, periodically check spotVM --> change to progress in minutes
			if(true)//(fmod(t, 10) == 0) 
			{
				int mark = tr->readprice(r);
				if(r != NULL) {			
					//cout<<fixed <<r[0]<<"\t"<<r[1]<<"\t"<<r[2]<<"\t"<<r[3]<<"\n";
				}
				else{
					std::cout<< "cannot read the price\n";	
				}
				//////////////
				for(int i=0; i<types; i++)
				{
					int size = sVMTP[i].size();
					for(int j=0; j<size; j++) {
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
								sVMTP[i][j]->tk->tasktime += t - sVMTP[i][j]->turn_on;

								bool isFound = false;
								while(!isFound) {
								//	sVMTP[i][j]->tk->vmID += 1;																	
										
									double pc = sVMTP[i][j]->tk->prices[sVMTP[i][j]->tk->vmID];										
										
									if(pc > 1e-12){
										isFound = true;
										int index = sVMTP[i][j]->tk->configList[sVMTP[i][j]->tk->vmID];
										sVMTP[i][j]->tk->restTime = sVMTP[i][j]->tk->actTime[index];
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
			vp = vertices(dag.g);
			for(int i=0; i < (*vp.second - *vp.first ); i++)
			{
				bool tag = true;
				//get parent vertices
				in_edge_iterator in_i, in_end;
				edge_descriptor e;
				for (boost::tie(in_i, in_end) = in_edges(i, dag.g); in_i != in_end; ++in_i) 
				{
					e = *in_i;
					Vertex src = source(e, dag.g);					
					if(dag.g[src].status != finished)
					{
						tag = false;
						//break;
					}
				}
				if(dag.g[i].status == ready || tag && dag.g[i].status != scheduled && dag.g[i].status != finished){
					ready_task.push_back(&dag.g[i]);							
				}
			}

			int acqondemand = 0;
			int acqspot = 0;

			std::sort(ready_task.begin(),ready_task.end(), myfunction);

			for(int i=0; i<ready_task.size(); i++)//earliest deadline first
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
								double availT = ceil(runtime/60.0)*60.0 - runtime;

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
						curr_task->restTime = curr_task->actTime[suittp] ;
						curr_task->readyCountdown = 0;

						if(suitisSpot){
							sVMTP[suittp][suittpindex]->tk = curr_task;
							curr_task->configList[curr_task->vmID] = suittp;
							curr_task->prices[curr_task->vmID] = sVMTP[suittp][suittpindex]->price;
						}else{
							VMTP[suittp][suittpindex]->tk = curr_task;
							curr_task->vmID = 2;
							curr_task->configList[2] = suittp;
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
						if(curr_task->vmID == 2)//the first 2 are spot vms
						{
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
						}else{
							int size = sVMTP[_config].size();
							for(int j=0; j<size; j++) {
								if(sVMTP[_config][j]->tk == NULL && sVMTP[_config][j]->life_time > curr_task->estTime[_config] && sVMTP[_config][j]->canAlloc)
								{
									find = true;
									sVMTP[_config][j]->tk = curr_task;
									break;
								}
							}
						}
						if(find) {
							curr_task->status = scheduled;
							curr_task->taskstart = t;
							int index = curr_task->configList[curr_task->vmID];
							curr_task->restTime = curr_task->actTime[index];
							curr_task->readyCountdown = 0;
						}
						else if(curr_task->vmID == 2) {curr_task->readyCountdown = OnDemandLag; curr_task->taskstart = t;}
						else {curr_task->readyCountdown = SpotLag; curr_task->taskstart = t;}
					}
				}
				else if(curr_task->readyCountdown == 0)
				{
					curr_task->status = scheduled;
					int index = curr_task->configList[curr_task->vmID];
					curr_task->restTime = curr_task->actTime[index];

					if(curr_task->vmID == 2)//ondemand VM
					{
						VM* vm = new VM;
						vm->life_time = 0; //OnDemandLag
						vm->tk = curr_task;
						vm->type = index;
						vm->turn_on = t;
						VMTP[index].push_back(vm);
						acqondemand += 1;
					}
					else if(curr_task->vmID == 0 || curr_task->vmID == 1)
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

							if(curr_task->vmID == 2) curr_task->readyCountdown = OnDemandLag;
							else curr_task->readyCountdown = SpotLag;
						}else{
							SpotVM* svm = new SpotVM(curr_task->prices[curr_task->vmID]);								

							svm->tk = curr_task; //it's spot vm
							svm->type = index;
							svm->life_time = 60; // - SpotLag
							svm->turn_on = t;
							sVMTP[index].push_back(svm);
							acqspot += 1;
						}
					}						
				}			
			}
			//delete VMs without task
			int totalondemand = 0;
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
						moneycost += priceOnDemand[i]*ceil(runtime/60.0);
						VM* vm = VMTP[i][j];
						delete vm;
						VMTP[i].erase(VMTP[i].begin()+j);
						j--;
						size1 --;
						delondemand += 1;
					}
				}
					
				for(int j=0; j<size2; j++)
				{
					if(sVMTP[i][j]->tk == NULL || (!sVMTP[i][j]->canAlloc))
					{
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
			for(int i=0; i<(*vp.second - *vp.first ); i++)
				if(dag.g[i].status == scheduled)
					scheduled_task.push_back(&dag.g[i]);
			for(int i=0; i<scheduled_task.size(); i++)
			{
				scheduled_task[i]->restTime -= 1;////////////////////////////
				if(scheduled_task[i]->restTime <= 0) 
				{
					scheduled_task[i]->status = finished;
					scheduled_task[i]->tasktime += t - scheduled_task[i]->taskstart;
					//deadline refinement and re-configuration
					//first check if the child task has been changed by other tasks
					out_edge_iterator out_i, out_end;
					Vertex v = scheduled_task[i]->name;
					for (boost::tie(out_i, out_end) = out_edges(v, dag.g); out_i != out_end; ++out_i) 
					{
						//if(dag.g[target(*out_i,dag.g)].start_time < t || (scheduled_task[i]->dl - t)/scheduled_task[i]->dl > 0.2)
						if(dag.g[target(*out_i,dag.g)].start_time < t || (dag.g[target(*out_i,dag.g)].start_time - t)/t > 0.2)
						{
							dag.g[target(*out_i, dag.g)].start_time = t;
							dag.g[target(*out_i, dag.g)].instance_config();
						}
					}
					//make the vm.task = NULL
					int index = scheduled_task[i]->configList[scheduled_task[i]->vmID];
					if(scheduled_task[i]->vmID == 2) //VM type
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
			for(int i=0; i < (*vp.second - *vp.first ); i++){
				if(dag.g[i].status!= finished){
					condition = true;
					unfinishednum += 1;
				}					
			}				
		}while(condition);//there is a task not finished

		//step 4 finalize
		for(int i=0; i<types; i++)
		{
			int size1 = VMTP[i].size();						
			for(int j=0; j<size1; j++)
			{
				double runtime = VMTP[i][j]->life_time;
				moneycost += (priceOnDemand[i] * ceil(runtime/60.0));
			}

			int size = sVMTP[i].size();
			for(int j=0; j<size; j++)
			{				
				moneycost += r[i];				
			}
		}

		tr->closefile();
		delete tr;

		printf("Money Cost: %.4f, Time: %.2f\t tasktime:", moneycost, t);
		for(int i=0; i<(*vp.second - *vp.first ); i++)
			printf("%.4f\t",dag.g[i].tasktime);
		printf("\t task cost: ");
		for(int i=0; i<(*vp.second - *vp.first ); i++)
			printf("%.4f\t",dag.g[i].cost);
		printf("spot: %d\n", totalspot);
	}
	std::clock_t endtime = std::clock();
	double timeelapsed = (double)(endtime - starttime) / (double)CLOCKS_PER_SEC;
	printf("time elapsed for Dyna Algorithm is: %.4f\n", timeelapsed);
	
	printf("----------------------------------start Uncertain_sc11 algorithm------------------------------------------\n");
	
	// monte carlo execution
	starttime = std::clock();
	for(int monte=0; monte < (num_monte+1); monte++)
	{
		double rnd =  (double)rand() / RAND_MAX;
		while(rnd<0.0001)
			rnd = (double)rand() / RAND_MAX;
		//double rnd = rn_01();
		printf("rnd:%f\t",rnd);

		ioseq[0] = math::quantile(seq_io_s, rnd);
		ioseq[1] = math::quantile(seq_io_m, rnd);
		ioseq[2] = math::quantile(seq_io_l, rnd);
		ioseq[3] = math::quantile(seq_io_x, rnd);
		iorand[0] = math::quantile(r_norm_s, rnd);
		iorand[1] = math::quantile(r_norm_m, rnd);
		iorand[2] = math::quantile(r_norm_l, rnd);
		iorand[3] = math::quantile(r_norm_x, rnd);
		net_up[0] = math::quantile(gamma_s_up, rnd);
		net_up[1] = math::quantile(gamma_m_up, rnd);
		net_up[2] = math::quantile(gamma_l_up, rnd);
		net_up[3] = math::quantile(gamma_x_up, rnd);
		net_down[0] = math::quantile(gamma_s_down, rnd);
		net_down[1] = math::quantile(gamma_m_down, rnd);
		net_down[2] = math::quantile(gamma_l_down, rnd);
		net_down[3] = math::quantile(gamma_x_down, rnd);

		//before deadline assign
		dag.reset();

		vp = vertices(dag.g);
		property_map<Graph, edge_weight_t>::type weightmap = get(edge_weight, dag.g);
		for(; vp.first != vp.second; ++vp.first)
		{
			Vertex v1 = *vp.first;
			for(int t=0; t<types; t++)
			{
				dag.g[v1].estTime[t] = Times[dag.g[v1].type][t] + dag.g[v1].read_data / iorand_req[t] + dag.g[v1].seq_data / ioseq_req[t] + dag.g[v1].trans_data * net_up_req[t]/8 + dag.g[v1].rec_data *net_down_req[t] /8;
				dag.g[v1].actTime[t] = Times[dag.g[v1].type][t] + dag.g[v1].read_data / iorand[t] + dag.g[v1].seq_data / ioseq[t] + dag.g[v1].trans_data * net_up[t]/8 + dag.g[v1].rec_data *net_down[t] /8;
			}
			in_edge_iterator in_i, in_end;
			for (boost::tie(in_i, in_end) = in_edges(v1, dag.g); in_i != in_end; ++in_i)
			{
				edge_descriptor e = *in_i;
				weightmap[e] = -1*(dag.g[v1].estTime[0]);//deadline assign using cheapest machine
			}
		}

		//deadline assignment
		dag.deadline_assign();

		vp = vertices(dag.g);
		for(; vp.first != vp.second; ++vp.first)
		{
			Vertex v1 = *vp.first;
			for(int t=0; t<types; t++)
			{
				if(dag.g[v1].estTime[t]+OnDemandLag < (dag.g[v1].dl - dag.g[v1].start_time))
				{
					dag.g[v1].config = t;
					break;
				}
				if(t == types-1)
					dag.g[v1].config = t;
			}
		}

		std::vector<VM*> VMTP[types];
		int need_VM[types]={0,0,0,0};

		//EDF scheduling
		vp = vertices(dag.g);
		dag.g[*vp.first ].status = ready;
		double t = 0;
		bool condition = false;
		double moneycost = 0.0;
		
		do{		
			//step 1
			std::vector<taskVertex*> ready_task;
			vp = vertices(dag.g);
			for(int i=0; i < (*vp.second - *vp.first ); i++)
			{
				bool tag = true;
				//get parent vertices
				in_edge_iterator in_i, in_end;
				edge_descriptor e;
				for (boost::tie(in_i, in_end) = in_edges(i, dag.g); in_i != in_end; ++in_i) 
				{
					e = *in_i;
					Vertex src = source(e, dag.g);					
					if(dag.g[src].status != finished)
					{
						tag = false;
						//break;
					}
				}
				if(dag.g[i].status == ready || tag && dag.g[i].status != scheduled && dag.g[i].status != finished){
					ready_task.push_back(&dag.g[i]);							
				}
			}

			std::sort(ready_task.begin(),ready_task.end(), myfunction);
			for(int i=0; i<ready_task.size(); i++)//earliest deadline first
			{
				taskVertex* curr_task=ready_task[i];
				if(curr_task->readyCountdown == -1)//
				{
						
					int _config = curr_task->config;
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
						curr_task->taskstart = t;
						curr_task->restTime = curr_task->actTime[curr_task->config];
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
					curr_task->restTime = curr_task->actTime[curr_task->config] ;

					VM* vm = new VM; 
					vm->life_time = OnDemandLag;
					vm->tk = curr_task;
					vm->type = curr_task->config;
					VMTP[curr_task->config].push_back(vm);
						
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
						moneycost += (priceOnDemand[i] * ceil(runtime/60.0));

						VMTP[i].erase(VMTP[i].begin()+j);
						size1 --;
					}
				}
			}
			//step 2
			std::vector<taskVertex*> scheduled_task;
			for(int i=0; i<(*vp.second - *vp.first ); i++)
				if(dag.g[i].status == scheduled)
					scheduled_task.push_back(&dag.g[i]);
			for(int i=0; i<scheduled_task.size(); i++)
			{
				scheduled_task[i]->restTime -= 1;////////////////////////////
				if(scheduled_task[i]->restTime <= 0) 
				{
					scheduled_task[i]->status = finished;
					scheduled_task[i]->tasktime += t - scheduled_task[i]->taskstart;
					scheduled_task[i]->cost = scheduled_task[i]->tasktime * priceOnDemand[scheduled_task[i]->config] /60.0;
					//deadline refinement and re-configuration
					//first check if the child task has been changed by other tasks
					out_edge_iterator out_i, out_end;
					Vertex v = scheduled_task[i]->name;
					for (boost::tie(out_i, out_end) = out_edges(v, dag.g); out_i != out_end; ++out_i) 
					{
						//if(dag.g[target(*out_i,dag.g)].start_time < t || (scheduled_task[i]->dl - t)/scheduled_task[i]->dl > 0.2)
						if(dag.g[target(*out_i,dag.g)].start_time < t || (dag.g[target(*out_i,dag.g)].start_time - t)/t > 0.2)
						{
							dag.g[target(*out_i, dag.g)].start_time = t;
							dag.g[target(*out_i, dag.g)].instance_config();
						}
					}

					//make the vm.task = NULL
					for(int j=0; j<VMTP[scheduled_task[i]->config].size(); j++)
						if(VMTP[scheduled_task[i]->config][j]->tk == scheduled_task[i])
						{
							VMTP[scheduled_task[i]->config][j]->tk = NULL;
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
			for(int i=0; i < (*vp.second - *vp.first ); i++){
				if(dag.g[i].status!= finished)
				{
					condition = true;
					unfinishednum += 1;
				}					
			}				
		}while(condition);//there is a task not finished

		for(int i=0; i<types; i++)
		{
			int size1 = VMTP[i].size();						
			for(int j=0; j<size1; j++)
			{
				double runtime = VMTP[i][j]->life_time;
				moneycost += (priceOnDemand[i] * ceil(runtime/60.0));
			}
		}
		printf("Money Cost: %.4f, Time: %.2f\t tasktime: ", moneycost, t);
		for(int i=0; i<(*vp.second - *vp.first ); i++)
			printf("%.4f\t",dag.g[i].tasktime);
		printf("\t task cost: ");
		for(int i=0; i<(*vp.second - *vp.first ); i++)
			printf("%.4f\t",dag.g[i].cost);
		printf("\n");
	}
	endtime = std::clock();
	timeelapsed = (double)(endtime - starttime) / (double)CLOCKS_PER_SEC;
	printf("time elapsed for uncertain_sc11 algorithm is: %.4f\n", timeelapsed);

	printf("-------------------------------------start HPDC algorithm--------------------------------------------\n");
	starttime = std::clock();

	// monte carlo execution
	for(int monte=0; monte < (num_monte+1); monte++)
	{
		double rnd =  (double)rand() / RAND_MAX;
		while(rnd<0.0001)
			rnd = (double)rand() / RAND_MAX;
		//double rnd = rn_01();
		printf("rnd:%f\t",rnd);

		ioseq[0] = math::quantile(seq_io_s, rnd);
		ioseq[1] = math::quantile(seq_io_m, rnd);
		ioseq[2] = math::quantile(seq_io_l, rnd);
		ioseq[3] = math::quantile(seq_io_x, rnd);
		iorand[0] = math::quantile(r_norm_s, rnd);
		iorand[1] = math::quantile(r_norm_m, rnd);
		iorand[2] = math::quantile(r_norm_l, rnd);
		iorand[3] = math::quantile(r_norm_x, rnd);
		net_up[0] = math::quantile(gamma_s_up, rnd);
		net_up[1] = math::quantile(gamma_m_up, rnd);
		net_up[2] = math::quantile(gamma_l_up, rnd);
		net_up[3] = math::quantile(gamma_x_up, rnd);
		net_down[0] = math::quantile(gamma_s_down, rnd);
		net_down[1] = math::quantile(gamma_m_down, rnd);
		net_down[2] = math::quantile(gamma_l_down, rnd);
		net_down[3] = math::quantile(gamma_x_down, rnd);

		//before deadline assign
		dag.reset();

		vp = vertices(dag.g);
		property_map<Graph, edge_weight_t>::type weightmap = get(edge_weight, dag.g);
		for(; vp.first != vp.second; ++vp.first)
		{
			Vertex v1 = *vp.first;
			for(int t=0; t<types; t++)
			{
				dag.g[v1].estTime[t] = Times[dag.g[v1].type][t] + dag.g[v1].read_data / ave_iorand[t] + dag.g[v1].seq_data / ave_ioseq[t]
			+ dag.g[v1].trans_data * ave_netup[t]/8 + dag.g[v1].rec_data * ave_netdown[t] /8;
				dag.g[v1].actTime[t] = Times[dag.g[v1].type][t] + dag.g[v1].read_data / iorand[t] + dag.g[v1].seq_data / ioseq[t]
			+ dag.g[v1].trans_data * net_up[t]/8 + dag.g[v1].rec_data * net_down[t] /8;
			}
			in_edge_iterator in_i, in_end;
			for (boost::tie(in_i, in_end) = in_edges(v1, dag.g); in_i != in_end; ++in_i)
			{
				edge_descriptor e = *in_i;
				weightmap[e] = -1* (dag.g[v1].estTime[0]);//deadline assign using cheapest machine
			}
		}

		dag.deadline_assign();

		for(vp = vertices(dag.g); vp.first != vp.second; vp.first ++)
			dag.g[*vp.first].instance_config();


		int tracelag = rand()%28800;
		std::vector<VM*> VMTP[types];
		std::vector<SpotVM*> sVMTP[types];


		//EDF scheduling
		vp = vertices(dag.g);
		dag.g[*vp.first ].status = ready;

		double t = 0;
		bool condition = false;
		double moneycost = 0.0;

		////trace data
		ReadTrace* tr = new ReadTrace();
		tr->readfile("data.csv", tracelag);

		float r[4];
		int totalspotfail = 0;
		int totalspot = 0;
		do{		
			int spotfail = 0;
			//step 0, periodically check spotVM --> change to progress in minutes
			if(true)//(fmod(t, 10) == 0) 
			{
				int mark = tr->readprice(r);
				if(r != NULL) {			
					//cout<<fixed <<r[0]<<"\t"<<r[1]<<"\t"<<r[2]<<"\t"<<r[3]<<"\n";
				}
				else{
					std::cout<< "cannot read the price\n";	
				}
				//////////////
				for(int i=0; i<types; i++)
				{
					int size = sVMTP[i].size();
					for(int j=0; j<size; j++) {
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
								sVMTP[i][j]->tk->tasktime += t - sVMTP[i][j]->turn_on;

								bool isFound = false;
								while(!isFound) {
								//	sVMTP[i][j]->tk->vmID += 1;																	
										
									double pc = sVMTP[i][j]->tk->prices[sVMTP[i][j]->tk->vmID];										
										
									if(pc > 1e-12){
										isFound = true;
										int index = sVMTP[i][j]->tk->configList[sVMTP[i][j]->tk->vmID];
										sVMTP[i][j]->tk->restTime = sVMTP[i][j]->tk->actTime[index] ;
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
			vp = vertices(dag.g);
			for(int i=0; i < (*vp.second - *vp.first ); i++)
			{
				bool tag = true;
				//get parent vertices
				in_edge_iterator in_i, in_end;
				edge_descriptor e;
				for (boost::tie(in_i, in_end) = in_edges(i, dag.g); in_i != in_end; ++in_i) 
				{
					e = *in_i;
					Vertex src = source(e, dag.g);					
					if(dag.g[src].status != finished)
					{
						tag = false;
						//break;
					}
				}
				if(dag.g[i].status == ready || tag && dag.g[i].status != scheduled && dag.g[i].status != finished){
					ready_task.push_back(&dag.g[i]);							
				}
			}

			int acqondemand = 0;
			int acqspot = 0;

			std::sort(ready_task.begin(),ready_task.end(), myfunction);

			for(int i=0; i<ready_task.size(); i++)//earliest deadline first
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
								double availT = ceil(runtime/60.0)*60.0 - runtime;

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
						curr_task->restTime = curr_task->actTime[suittp] ;
						curr_task->readyCountdown = 0;

						if(suitisSpot){
							sVMTP[suittp][suittpindex]->tk = curr_task;
							curr_task->configList[curr_task->vmID] = suittp;
							curr_task->prices[curr_task->vmID] = sVMTP[suittp][suittpindex]->price;
						}else{
							VMTP[suittp][suittpindex]->tk = curr_task;
							curr_task->vmID = 2;
							curr_task->configList[2] = suittp;
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
						if(curr_task->vmID == 2)//the first 2 are spot vms
						{
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
						}else{
							int size = sVMTP[_config].size();
							for(int j=0; j<size; j++) {
								if(sVMTP[_config][j]->tk == NULL && sVMTP[_config][j]->life_time > curr_task->estTime[_config] && sVMTP[_config][j]->canAlloc)
								{
									find = true;
									sVMTP[_config][j]->tk = curr_task;
									break;
								}
							}
						}
						if(find) {
							curr_task->status = scheduled;
							curr_task->taskstart = t;
							int index = curr_task->configList[curr_task->vmID];
							curr_task->restTime = curr_task->actTime[index];
							curr_task->readyCountdown = 0;
						}
						else if(curr_task->vmID == 2) {curr_task->readyCountdown = OnDemandLag; curr_task->taskstart = t;}
						else {curr_task->readyCountdown = SpotLag; curr_task->taskstart = t;}
					}
				}
				else if(curr_task->readyCountdown == 0)
				{
					curr_task->status = scheduled;
					int index = curr_task->configList[curr_task->vmID];
					curr_task->restTime = curr_task->actTime[index];

					if(curr_task->vmID == 2)//ondemand VM
					{
						VM* vm = new VM;
						vm->life_time = 0; //OnDemandLag
						vm->tk = curr_task;
						vm->type = index;
						vm->turn_on = t;
						VMTP[index].push_back(vm);
						acqondemand += 1;
					}
					else if(curr_task->vmID == 0 || curr_task->vmID == 1)
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

							if(curr_task->vmID == 2) curr_task->readyCountdown = OnDemandLag;
							else curr_task->readyCountdown = SpotLag;
						}else{
							SpotVM* svm = new SpotVM(curr_task->prices[curr_task->vmID]);								

							svm->tk = curr_task; //it's spot vm
							svm->type = index;
							svm->life_time = 60; // - SpotLag
							svm->turn_on = t;
							sVMTP[index].push_back(svm);
							acqspot += 1;
						}
					}						
				}			
			}
			//delete VMs without task
			int totalondemand = 0;
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
						moneycost += priceOnDemand[i]*ceil(runtime/60.0);
						VM* vm = VMTP[i][j];
						delete vm;
						VMTP[i].erase(VMTP[i].begin()+j);
						j--;
						size1 --;
						delondemand += 1;
					}
				}
					
				for(int j=0; j<size2; j++)
				{
					if(sVMTP[i][j]->tk == NULL || (!sVMTP[i][j]->canAlloc))
					{
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
			for(int i=0; i<(*vp.second - *vp.first ); i++)
				if(dag.g[i].status == scheduled)
					scheduled_task.push_back(&dag.g[i]);
			for(int i=0; i<scheduled_task.size(); i++)
			{
				scheduled_task[i]->restTime -= 1;////////////////////////////
				if(scheduled_task[i]->restTime <= 0) 
				{
					scheduled_task[i]->status = finished;
					scheduled_task[i]->tasktime += t - scheduled_task[i]->taskstart;
					//make the vm.task = NULL
					int index = scheduled_task[i]->configList[scheduled_task[i]->vmID];
					if(scheduled_task[i]->vmID == 2) //VM type
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
			for(int i=0; i < (*vp.second - *vp.first ); i++){
				if(dag.g[i].status!= finished){
					condition = true;
					unfinishednum += 1;
				}					
			}				
		}while(condition);//there is a task not finished

		//step 4 finalize
		for(int i=0; i<types; i++)
		{
			int size1 = VMTP[i].size();						
			for(int j=0; j<size1; j++)
			{
				double runtime = VMTP[i][j]->life_time;
				moneycost += (priceOnDemand[i] * ceil(runtime/60.0));
			}

			int size = sVMTP[i].size();
			for(int j=0; j<size; j++)
			{				
				moneycost += r[i];				
			}
		}

		tr->closefile();
		delete tr;

		printf("Money Cost: %.4f, Time: %.2f\t tasktime: ", moneycost, t);
		for(int i=0; i<(*vp.second - *vp.first ); i++)
			printf("%.4f\t",dag.g[i].tasktime);
		printf("\t task cost: ");
		for(int i=0; i<(*vp.second - *vp.first ); i++)
			printf("%.4f\t",dag.g[i].cost);
		printf("spot: %d\n",totalspot);
	}
	endtime = std::clock();
	timeelapsed = (double)(endtime - starttime) / (double)CLOCKS_PER_SEC;
	printf("time elapsed for HPDC algorithm is: %.4f\n", timeelapsed);

	printf("-------------------------------------start SC11 static algorithm--------------------------------------------\n");

	starttime = std::clock();
	for(int monte=0; monte < (num_monte+1); monte++)
	{
		double rnd =  (double)rand() / RAND_MAX;
		while(rnd<0.0001)
			rnd = (double)rand() / RAND_MAX;
		//double rnd = rn_01();
		printf("rnd:%f\t",rnd);

		ioseq[0] = math::quantile(seq_io_s, rnd);
		ioseq[1] = math::quantile(seq_io_m, rnd);
		ioseq[2] = math::quantile(seq_io_l, rnd);
		ioseq[3] = math::quantile(seq_io_x, rnd);
		iorand[0] = math::quantile(r_norm_s, rnd);
		iorand[1] = math::quantile(r_norm_m, rnd);
		iorand[2] = math::quantile(r_norm_l, rnd);
		iorand[3] = math::quantile(r_norm_x, rnd);
		net_up[0] = math::quantile(gamma_s_up, rnd);
		net_up[1] = math::quantile(gamma_m_up, rnd);
		net_up[2] = math::quantile(gamma_l_up, rnd);
		net_up[3] = math::quantile(gamma_x_up, rnd);
		net_down[0] = math::quantile(gamma_s_down, rnd);
		net_down[1] = math::quantile(gamma_m_down, rnd);
		net_down[2] = math::quantile(gamma_l_down, rnd);
		net_down[3] = math::quantile(gamma_x_down, rnd);

		//before deadline assign
		dag.reset();

		vp = vertices(dag.g);
		property_map<Graph, edge_weight_t>::type weightmap = get(edge_weight, dag.g);
		for(; vp.first != vp.second; ++vp.first)
		{
			Vertex v1 = *vp.first;
			for(int t=0; t<types; t++)
			{
				dag.g[v1].estTime[t] = Times[dag.g[v1].type][t] + dag.g[v1].read_data / ave_iorand[t] + dag.g[v1].seq_data / ave_ioseq[t]
			+ dag.g[v1].trans_data * ave_netup[t]/8 + dag.g[v1].rec_data * ave_netdown[t] /8;
				dag.g[v1].actTime[t] = Times[dag.g[v1].type][t] + dag.g[v1].read_data / iorand[t] + dag.g[v1].seq_data / ioseq[t]
			+ dag.g[v1].trans_data * net_up[t]/8 + dag.g[v1].rec_data * net_down[t] /8;
			}
			in_edge_iterator in_i, in_end;
			for (boost::tie(in_i, in_end) = in_edges(v1, dag.g); in_i != in_end; ++in_i)
			{
				edge_descriptor e = *in_i;
				weightmap[e] = -1* (dag.g[v1].estTime[dag.g[v1].config]);//deadline assign using cheapest machine
			}
		}

		dag.deadline_assign();

		vp = vertices(dag.g);
		for(; vp.first != vp.second; ++vp.first)
		{
			Vertex v1 = *vp.first;
			for(int i=0; i<types; i++)
			{
				if(dag.g[v1].estTime[i]+OnDemandLag < (dag.g[v1].dl-dag.g[v1].start_time))
				{
					dag.g[v1].config = i;
					break;
				}
				if(i == types-1)
					dag.g[v1].config = i;
			}
		}


		std::vector<VM*> VMTP[types];
		int need_VM[types]={0,0,0,0};

		//EDF scheduling
		vp = vertices(dag.g);
		dag.g[*vp.first ].status = ready;
		double t = 0;
		bool condition = false;
		double moneycost = 0.0;
		
		do{		
			//step 1
			std::vector<taskVertex*> ready_task;
			vp = vertices(dag.g);
			for(int i=0; i < (*vp.second - *vp.first ); i++)
			{
				bool tag = true;
				//get parent vertices
				in_edge_iterator in_i, in_end;
				edge_descriptor e;
				for (boost::tie(in_i, in_end) = in_edges(i, dag.g); in_i != in_end; ++in_i) 
				{
					e = *in_i;
					Vertex src = source(e, dag.g);					
					if(dag.g[src].status != finished)
					{
						tag = false;
						//break;
					}
				}
				if(dag.g[i].status == ready || tag && dag.g[i].status != scheduled && dag.g[i].status != finished){
					ready_task.push_back(&dag.g[i]);							
				}
			}

			std::sort(ready_task.begin(),ready_task.end(), myfunction);
			for(int i=0; i<ready_task.size(); i++)//earliest deadline first
			{
				taskVertex* curr_task=ready_task[i];
				if(curr_task->readyCountdown == -1)//
				{
						
					int _config = curr_task->config;
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
						curr_task->restTime =  curr_task->actTime[curr_task->config] ;
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
					curr_task->restTime = curr_task->actTime[curr_task->config] ;

					VM* vm = new VM; 
					vm->life_time = OnDemandLag;
					vm->tk = curr_task;
					vm->type = curr_task->config;
					VMTP[curr_task->config].push_back(vm);
						
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
						moneycost += (priceOnDemand[i] * ceil(runtime/60.0));

						VMTP[i].erase(VMTP[i].begin()+j);
						size1 --;
					}
				}
			}
			//step 2
			std::vector<taskVertex*> scheduled_task;
			for(int i=0; i<(*vp.second - *vp.first ); i++)
				if(dag.g[i].status == scheduled)
					scheduled_task.push_back(&dag.g[i]);
			for(int i=0; i<scheduled_task.size(); i++)
			{
				scheduled_task[i]->restTime -= 1;////////////////////////////
				if(scheduled_task[i]->restTime <= 0) 
				{
					scheduled_task[i]->status = finished;
					scheduled_task[i]->tasktime = t - scheduled_task[i]->tasktime;
					scheduled_task[i]->cost = scheduled_task[i]->tasktime * priceOnDemand[scheduled_task[i]->config] /60.0;
					//make the vm.task = NULL
					for(int j=0; j<VMTP[scheduled_task[i]->config].size(); j++)
						if(VMTP[scheduled_task[i]->config][j]->tk == scheduled_task[i])
						{
							VMTP[scheduled_task[i]->config][j]->tk = NULL;
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
			for(int i=0; i < (*vp.second - *vp.first ); i++){
				if(dag.g[i].status!= finished)
				{
					condition = true;
					unfinishednum += 1;
				}					
			}				
		}while(condition);//there is a task not finished

		for(int i=0; i<types; i++)
		{
			int size1 = VMTP[i].size();						
			for(int j=0; j<size1; j++)
			{
				double runtime = VMTP[i][j]->life_time;
				moneycost += (priceOnDemand[i] * ceil(runtime/60.0));
			}
		}
		printf("Money Cost: %.4f, Time: %.2f\t tasktime: ", moneycost, t);
		for(int i=0; i<(*vp.second - *vp.first ); i++)
				printf("%.4f\t",dag.g[i].tasktime);
		printf("\t task cost: ");
		for(int i=0; i<(*vp.second - *vp.first ); i++)
			printf("%.4f\t",dag.g[i].cost);
		printf("\n");
	}
	endtime = std::clock();
	timeelapsed = (double)(endtime - starttime) / (double)CLOCKS_PER_SEC;
	printf("time elapsed for SC11 algorithm is: %.4f\n", timeelapsed);

	printf("-------------------------------------start SpotOnly algorithm--------------------------------------------\n");
	double tp[][2] = {{0.038, 0.040}, {0.064, 0.068}, {0.152, 0.155},{0.304, 0.334}};

	int successnum = 0;
	int totaltrial = 0;
	double totalcost = 0.0;
	double maxcost = 0.0;
	double mincost = 2000.0;

	int expt = 0;
	double expp = 0.0;

	starttime = std::clock();
	for(int monte=0; monte < (num_monte+1); monte++)
	{
		double rnd =  (double)rand() / RAND_MAX;
		while(rnd<0.0001)
			rnd = (double)rand() / RAND_MAX;
		//double rnd = rn_01();
		printf("rnd:%f\t",rnd);

		ioseq[0] = math::quantile(seq_io_s, rnd);
		ioseq[1] = math::quantile(seq_io_m, rnd);
		ioseq[2] = math::quantile(seq_io_l, rnd);
		ioseq[3] = math::quantile(seq_io_x, rnd);
		iorand[0] = math::quantile(r_norm_s, rnd);
		iorand[1] = math::quantile(r_norm_m, rnd);
		iorand[2] = math::quantile(r_norm_l, rnd);
		iorand[3] = math::quantile(r_norm_x, rnd);
		net_up[0] = math::quantile(gamma_s_up, rnd);
		net_up[1] = math::quantile(gamma_m_up, rnd);
		net_up[2] = math::quantile(gamma_l_up, rnd);
		net_up[3] = math::quantile(gamma_x_up, rnd);
		net_down[0] = math::quantile(gamma_s_down, rnd);
		net_down[1] = math::quantile(gamma_m_down, rnd);
		net_down[2] = math::quantile(gamma_l_down, rnd);
		net_down[3] = math::quantile(gamma_x_down, rnd);

		dag.reset();

		////randomly select bidding price between min and max
		vp = vertices(dag.g);
		for(; vp.first != vp.second ; ++vp.first )
		{
			Vertex v1 = *vp.first;
			dag.g[v1].prices[0] = (double)rand() / RAND_MAX * (tp[0][1]- tp[0][0]) + tp[0][0];
			for(int t=0; t<types; t++)
			{
				dag.g[v1].estTime[t] = Times[dag.g[v1].type][t] + dag.g[v1].read_data / ave_iorand[t] + dag.g[v1].seq_data / ave_ioseq[t]
			+ dag.g[v1].trans_data * ave_netup[t]/8 + dag.g[v1].rec_data * ave_netdown[t] /8;
				dag.g[v1].actTime[t] = Times[dag.g[v1].type][t] + dag.g[v1].read_data / iorand[t] + dag.g[v1].seq_data / ioseq[t]
			+ dag.g[v1].trans_data * net_up[t]/8 + dag.g[v1].rec_data * net_down[t] /8;
			}
		}
		//

		int tracelag = rand()%28800;
		std::vector<SpotVM*> sVMTP[types];

		vp = vertices(dag.g);
		dag.g[*vp.first ].status = ready;

		double t = 0;
		double condition = false;
		double moneycost = 0;

		////trace data
		ReadTrace* tr = new ReadTrace();
		tr->readfile("data.csv", tracelag);

		float r[4];
		do{
			int spotfail = 0;
			//step 0, periodically check spotVM
			if(true)//(fmod(t, 10) == 0) 
			{
				int mark = tr->readprice(r);
				if(r != NULL) {			
					//cout<<fixed <<r[0]<<"\t"<<r[1]<<"\t"<<r[2]<<"\t"<<r[3]<<"\n";
				}
				else{
					std::cout<< "cannot read the price\n";	
				}
				//////////////
				for(int i=0; i<types; i++)
				{
					int size = sVMTP[i].size();
					for(int j=0; j<size; j++) {
						bool check;
						check = function(sVMTP[i][j]->price, r[i]);//if not fail							
						//if(check) sVMTP[i][j]->life_time = 10;
						if(!check){		//failed	
							spotfail += 1;
							if(sVMTP[i][j]->tk != NULL)	{
								sVMTP[i][j]->canAlloc = false;
								sVMTP[i][j]->tk->status = ready;
								sVMTP[i][j]->tk->readyCountdown = -1;
								sVMTP[i][j]->tk->tasktime += t - sVMTP[i][j]->turn_on;
								sVMTP[i][j]->tk->cost += r[i] * (int)(t - sVMTP[i][j]->turn_on)/60;
									
								sVMTP[i][j]->tk->config = (double)rand() / RAND_MAX * types;
								//put the bidding price in the first slot of prices[]
								sVMTP[i][j]->tk->prices[0] = (double)rand() / RAND_MAX 
							* (tp[sVMTP[i][j]->tk->config][1]- tp[sVMTP[i][j]->tk->config][0]) + tp[sVMTP[i][j]->tk->config][0]; 

								sVMTP[i][j]->tk->restTime = sVMTP[i][j]->tk->actTime[sVMTP[i][j]->tk->config];
															
								sVMTP[i][j]->tk = NULL;	
							}
						}								
					}	
				}
			}

			//step 1
			std::vector<taskVertex*> ready_task;
			vp = vertices(dag.g);
			for(int i=0; i < (*vp.second - *vp.first ); i++)
			{
				bool tag = true;
				//get parent vertices
				in_edge_iterator in_i, in_end;
				edge_descriptor e;
				for (boost::tie(in_i, in_end) = in_edges(i, dag.g); in_i != in_end; ++in_i) 
				{
					e = *in_i;
					Vertex src = source(e, dag.g);					
					if(dag.g[src].status != finished)
					{
						tag = false;
						//break;
					}
				}
				if(dag.g[i].status == ready || tag && dag.g[i].status != scheduled && dag.g[i].status != finished){
					ready_task.push_back(&dag.g[i]);							
				}
			}
			int acqspot = 0;

			std::sort(ready_task.begin(),ready_task.end(), myfunction);
			for(int i=0; i<ready_task.size(); i++)//earliest deadline first
			{
				taskVertex* curr_task=ready_task[i];
				if(curr_task->readyCountdown == -1)//
				{
					//find spot first
					bool suitfind = false;
					int suittp = 0;
					int suittpindex = 0;
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
					if(suitfind){//found a spot instance
						curr_task->status = scheduled;
						curr_task->taskstart = t;
						curr_task->restTime = curr_task->actTime[suittp] ;
						curr_task->readyCountdown = 0;

						sVMTP[suittp][suittpindex]->tk = curr_task;
						curr_task->config = suittp;
						curr_task->prices[0] = sVMTP[suittp][suittpindex]->price;						
					}
					else {curr_task->readyCountdown = SpotLag; curr_task->taskstart = t;}
				}
				else if (curr_task->readyCountdown == 0)
				{
					curr_task->status = scheduled;
					int index = curr_task->config;
					curr_task->restTime = curr_task->actTime[index] ;

					if(r[curr_task->config] > curr_task->prices[0]){ //apply spotVM failed
						curr_task->status = ready;						
						curr_task->prices[0] = (double)rand() / RAND_MAX * (tp[curr_task->config][1]- tp[curr_task->config][0]) + tp[curr_task->config][0];
						curr_task->readyCountdown = SpotLag;
						curr_task->tasktime += t - curr_task->taskstart;
						curr_task->taskstart = t;
					}else{
						SpotVM* svm = new SpotVM(curr_task->prices[0]);								

						svm->tk = curr_task; //it's spot vm
						svm->type = curr_task->config;
						svm->life_time = 60; // - SpotLag
						svm->turn_on = t;
						sVMTP[curr_task->config].push_back(svm);
						acqspot += 1;
					}
				}
			}

			//delete VMs without task
			int totalondemand = 0;
			int totalspot = 0;
			int delondemand = 0;
			int delspot = 0;

			for(int i=0; i<types; i++)
			{
				int size2 = sVMTP[i].size();
				totalspot += size2;
					
				for(int j=0; j<size2; j++)
				{
					if(sVMTP[i][j]->tk == NULL || (!sVMTP[i][j]->canAlloc))
					{
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
			for(int i=0; i<(*vp.second - *vp.first ); i++)
				if(dag.g[i].status == scheduled)
					scheduled_task.push_back(&dag.g[i]);
			for(int i=0; i<scheduled_task.size(); i++)
			{
				scheduled_task[i]->restTime -= 1;////////////////////////////
				if(scheduled_task[i]->restTime <= 0) 
				{
					scheduled_task[i]->status = finished;
					//make the vm.task = NULL
					//scheduled_task[i]->tasktime += t - scheduled_task[i]->taskstart;
					for(int j=0; j<sVMTP[scheduled_task[i]->config].size(); j++)
						if(sVMTP[scheduled_task[i]->config][j]->tk == scheduled_task[i])
						{
							scheduled_task[i]->tasktime += t - sVMTP[scheduled_task[i]->config][j]->turn_on;
							scheduled_task[i]->cost += r[i] *(int)(t - sVMTP[scheduled_task[i]->config][j]->turn_on)/60;
							sVMTP[scheduled_task[i]->config][j]->tk = NULL;
							break;
						}					
				}
			}

			//step 3
			for(int i=0; i<types; i++)
			{
				int size = sVMTP[i].size();
				for(int j=0; j<size; j++)
				{
					sVMTP[i][j]->life_time -= 1;//
					if(sVMTP[i][j]->life_time == 0){
						sVMTP[i][j]->life_time = 60;
						moneycost += r[i];
						//sVMTP[i][j]->tk->cost += r[i];
					}
				}
			}
			for(int i=0; i<ready_task.size(); i++)//////////////////////////////////if >0
				if(ready_task[i]->readyCountdown > 0)
					ready_task[i]->readyCountdown -= 1;
			t += 1;

			condition = false;
			int unfinishednum = 0;
			for(int i=0; i < (*vp.second - *vp.first ); i++){
				if(dag.g[i].status!= finished){
					condition = true;
					unfinishednum += 1;
				}					
			}
		}while(condition);

		//step 4 finalize
		for(int i=0; i<types; i++)
		{
			int size = sVMTP[i].size();
			for(int j=0; j<size; j++)
				moneycost += r[i];		

		}
		tr->closefile();
		delete tr;

		printf("Money Cost: %.4f, Time: %.2f\t tasktime: ", moneycost, t);
		for(int i=0; i<(*vp.second - *vp.first ); i++)
			printf("%.4f\t",dag.g[i].tasktime);
		printf("\t task cost: ");
		for(int i=0; i<(*vp.second - *vp.first ); i++)
			printf("%.4f\t",dag.g[i].cost);
		printf("\n");
	}//monte carlo ends
	
	endtime = std::clock();
	timeelapsed = (double)(endtime - starttime) / (double)CLOCKS_PER_SEC;
	printf("time elapsed for SpotOnly algorithm is: %.4f\n", timeelapsed);

	return 0;
}
