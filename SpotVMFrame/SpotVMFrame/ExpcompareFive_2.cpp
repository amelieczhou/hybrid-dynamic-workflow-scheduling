// Deadline_assign.cpp : Defines the entry point for the console application.
// Input: binary tree type DAG
//compare to previous version: add network and I/O operation time to all algorithms
//change estTime in order to make the instance configuration adapt to performance variance
//after discussion on 6th, Oct add runtime deadline refinement
#include "stdafx.h"
#include <stdlib.h>
#include <vector>
#include <iostream>
#include <time.h>
#include "ReadTrace.h"
#include "PricingModel.h"
#include "InstanceConfig.h"
#include <boost/graph/topological_sort.hpp>


using namespace boost;
//I/O performance distribution, seek/sec
double OnDemandLag = 90;//seconds//10;//0.5;
double SpotLag = 60;//seconds//1;
bool NoSpotVM = true;
double Times[4][types] = {{120,65,38,24},{90,50,30,20},{60,35,23,17},{30,20,15,13}};
double lambda;
int num_jobs;

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
	lambda = atof(argv[14]);
	num_jobs = atoi(argv[15]);

	for(int i=0; i<4; i++)
		for(int j=1; j<types; j++)
			Times[i][j]=Times[i][0]*exe_time/pow(2.0,j)+Times[i][0]*(1-exe_time);

	//I/O speed, seeks/sec
	double ioseq[4], iorand[4];
	//network speed, sec/8M
	double net_up[4], net_down[4];
	double r = (double) rand()/(RAND_MAX+1);
	if(strcmp (argv[12], "sens") == 0) ////for senstivity test, fix the DAG type and task types
		r = 0.5;
	//ioseq[0] = math::quantile(seq_io_s, r);
	//iorand[0] = math::quantile(r_norm_s, r);
	//net_up[0] = math::quantile(gamma_s_up, r);

	//input workflows
	//Graph dag;
	DAG dag(deadline,meet_dl);
	
	if(strcmp(argv[1],"montage") == 0)
	{
		dag.type = montage;
		//generate a montage DAG
		//unit time is second
		double tProjectPP[types] = {540,264,150,90};//{100,100,100,1};//{8,4,2,1};
		double tDiffFit[types] = {1200,600,339,200};//{100,100,100,1};//{8,4,2,1};
		double tConcatFit[types] = {2268,1122,636,372};//{100,100,100,1};//{8,4,2,1};
		double tBgModel[types] = {17244,8200,4633,2725};//{100,100,100,1};//{8,4,2,1};
		double tBackground[types] = {40,20,11,7};//{100,100,100,1};//{8,4,2,1};
		double tImgTbl[types] = {56,30,17,10};//{100,100,100,1};//{8,4,2,1};
		double tAdd[types] = {1080,179,101,59};//{100,100,100,1};//{8,4,2,1};
		double tShrink[types] = {324,150,85,50};//{100,100,100,1};//{8,4,2,1};
		double tJPEG[types] = {186,112,63,37};//{100,100,100,1};//{8,4,2,1};

		for(int i=0; i<20; i++){
			taskVertex tk;
			tk.name = i;
			tk.estTime = new double [types];
			tk.cpuTime = new double [types];
			tk.config = tk.mark = tk.assigned_type = tk.restTime = 0; tk.rec_data = 0; 
			tk.readyCountdown = -1; tk.status = not_ready;
            tk.actTime = new double [types];
			//if(i == 0 || i == 21) //two dummy tasks
			//{	
			//	for(int j=0; j<types; j++)
			//		tk.estTime[j] = 0;
			//}else 
			for(int i=0; i<types; i++)
				tk.estTime[i] = tk.actTime[i] = 0;
			if(i>=0 && i<4)	{
				for(int j=0; j<types; j++)
					tk.cpuTime[j] = tProjectPP[j];
				tk.read_data = 0; //random io in seeks
				tk.trans_data = 2070; //network transfer in MB
				tk.seq_data = 1050; //sequential io in MB
				tk.rec_data = 0; //network download data in MB
			}else if(i>=4 && i<10){
				for(int j=0; j<types; j++)
					tk.cpuTime[j] = tDiffFit[j];
				tk.read_data = 0; //random io in seeks
				tk.trans_data = 86.6; //network transfer in MB
				tk.seq_data = 2100; //sequential io in MB
				tk.rec_data = 2100; //network download data in MB
			}else if(i == 10) {
				for(int j=0; j<types; j++)
					tk.cpuTime[j] = tConcatFit[j];
				tk.read_data = 0; //random io in seeks
				tk.trans_data = 1; //network transfer in MB
				tk.seq_data = 519.6; //sequential io in MB
				tk.rec_data = 519.6; //network download data in MB
			}else if(i == 11) {
				for(int j=0; j<types; j++)
					tk.cpuTime[j] = tBgModel[j];
				tk.read_data = 0; //random io in seeks
				tk.trans_data = 0.1; //network transfer in MB
				tk.seq_data = 1; //sequential io in MB
				tk.rec_data = 1; //network download data in MB
			}else if(i>=12 && i<16){
				for(int j=0; j<types; j++)
					tk.cpuTime[j] = tBackground[j];
				tk.read_data = 0; //random io in seeks
				tk.trans_data = 2070; //network transfer in MB
				tk.seq_data = 2070; //sequential io in MB
				tk.rec_data = 2070; //network download data in MB
			}else if(i == 16){
				for(int j=0; j<types; j++)
					tk.cpuTime[j] = tImgTbl[j];
				tk.read_data = 0; //random io in seeks
				tk.trans_data = 0; //network transfer in MB
				tk.seq_data = 8280; //sequential io in MB
				tk.rec_data = 8280; //network download data in MB
			}else if(i == 17){
				for(int j=0; j<types; j++)
					tk.cpuTime[j] = tAdd[j];
				tk.read_data = 0; //random io in seeks
				tk.trans_data = 4950; //network transfer in MB
				tk.seq_data = 8280; //sequential io in MB
				tk.rec_data = 8280; //network download data in MB
			}else if(i == 18){
				for(int j=0; j<types; j++)
					tk.cpuTime[j] = tShrink[j];
				tk.read_data = 0; //random io in seeks
				tk.trans_data = 154.4; //network transfer in MB
				tk.seq_data = 2470; //sequential io in MB
				tk.rec_data = 2470; //network download data in MB
			}else if(i == 19) {
				for(int j=0; j<types; j++)
					tk.cpuTime[j] = tJPEG[j];
				tk.read_data = 0; //random io in seeks
				tk.trans_data = 6.8; //network transfer in MB
				tk.seq_data = 154.4; //sequential io in MB
				tk.rec_data = 154.4; //network download data in MB
			}
			//for debuging, dont forget to comment it!!!!!!!!!!
			//tk.read_data = tk.trans_data = tk.seq_data = 0;
			add_vertex(tk, dag.g);
		}
		//add edges to the graph
		//add edges to the graph
		int l1[4] = {0,1,2,3};
		int l2[6] = {4,5,6,7,8,9};
		int l3[4] = {12,13,14,15};
		add_edge(0,4,dag.g);
		add_edge(1,4,dag.g);
		add_edge(0,5,dag.g);
		add_edge(1,5,dag.g);
		add_edge(1,6,dag.g);
		add_edge(2,6,dag.g);
		add_edge(2,7,dag.g);
		add_edge(3,7,dag.g);
		add_edge(1,8,dag.g);
		add_edge(3,8,dag.g);
		add_edge(3,9,dag.g);

		for(int i=0; i<6; i++){
			add_edge(l2[i],10,dag.g);
		}
		add_edge(10,11,dag.g);
		for(int i=0; i<4;i++){
			add_edge(i,l3[i],dag.g);
			add_edge(11,l3[i],dag.g);
			add_edge(l3[i],16,dag.g);
		}
		add_edge(16,17,dag.g);
		add_edge(17,18,dag.g);
		add_edge(18,19,dag.g);
	}

	//suitable for all kinds of DAG
	if(dag.type != montage) {
		std::pair<vertex_iter, vertex_iter> vp; 
		vp = vertices(dag.g);
		for(; vp.first != vp.second ; ++vp.first ) {
			in_edge_iterator in_i, in_end;
			edge_descriptor e;
			for (boost::tie(in_i, in_end) = in_edges(*vp.first , dag.g); in_i != in_end; ++in_i) {
				e = *in_i;
				Vertex src = source(e, dag.g);
				dag.g[*vp.first ].rec_data += dag.g[src].trans_data;			
			}
			printf("receive data for node %d: %f\n", *vp.first, dag.g[*vp.first].rec_data);
		}
	}

	//offline optimization
	SearchPrune* optimizer = new SearchPrune();
	optimizer->dag = dag;
	optimizer->OfflineSP_DFS();
	optimizer->SpotTune();
	optimizer->OnlineSimulate();
	return 0;
}
