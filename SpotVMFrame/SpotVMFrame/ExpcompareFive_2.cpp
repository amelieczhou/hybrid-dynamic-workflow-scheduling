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
#include "Algorithms.h"
#include <boost/graph/topological_sort.hpp>
#include <boost/math/distributions/gamma.hpp>
#include <boost/math/distributions/normal.hpp>

using namespace boost;

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
	double r = (double) rand()/(RAND_MAX+1);
	if(strcmp (argv[12], "sens") == 0) ////for senstivity test, fix the DAG type and task types
		r = 0.5;
	ioseq[0] = math::quantile(seq_io_s, r);
	iorand[0] = math::quantile(r_norm_s, r);
	net_up[0] = math::quantile(gamma_s_up, r);

	//input workflows
	//Graph dag;
	DAG dag(deadline);
	
	if(strcmp(argv[1],"montage") == 0)
	{
		dag.type = montage;
		//generate a montage DAG
		double tProjectPP[types] = {};
		double tDiffFit[types] = {};
		double tConcatFit[types] = {};
		double tBgModel[types] = {};
		double tBackground[types] = {};
		double tImgTbl[types] = {};
		double tAdd[types] = {};
		double tShrink[types] = {};
		double tJPEG[types] = {};

		for(int i=0; i<20; i++){
			taskVertex tk;
			tk.name = i;
			tk.estTime = new double [types];
			tk.config = tk.mark = tk.assigned_type = 0; tk.rec_data = 0; 
			tk.readyCountdown = -1; tk.status = not_ready;
            tk.actTime = new double [types];
			//if(i == 0 || i == 21) //two dummy tasks
			//{	
			//	for(int j=0; j<types; j++)
			//		tk.estTime[j] = 0;
			//}else 
			if(i>=0 && i<=4)	{
				for(int j=0; j<types; j++)
					tk.estTime[j] = tProjectPP[j];
				tk.read_data = ; //random io in seeks
				tk.trans_data = ; //network transfer in MB
				tk.seq_data = ; //sequential io in seeks
			}else if(i>=5 && i<=10){
				for(int j=0; j<types; j++)
					tk.estTime[j] = tDiffFit[j];
				tk.read_data = ; //random io in seeks
				tk.trans_data = ; //network transfer in MB
				tk.seq_data = ; //sequential io in seeks
			}else if(i == 11) {
				for(int j=0; j<types; j++)
					tk.estTime[j] = tConcatFit[j];
				tk.read_data = ; //random io in seeks
				tk.trans_data = ; //network transfer in MB
				tk.seq_data = ; //sequential io in seeks
			}else if(i == 12) {
				for(int j=0; j<types; j++)
					tk.estTime[j] = tBgModel[j];
				tk.read_data = ; //random io in seeks
				tk.trans_data = ; //network transfer in MB
				tk.seq_data = ; //sequential io in seeks
			}else if(i>=13 && i<=16){
				for(int j=0; j<types; j++)
					tk.estTime[j] = tBackground[j];
				tk.read_data = ; //random io in seeks
				tk.trans_data = ; //network transfer in MB
				tk.seq_data = ; //sequential io in seeks
			}else if(i == 17){
				for(int j=0; j<types; j++)
					tk.estTime[j] = tImgTbl[j];
				tk.read_data = ; //random io in seeks
				tk.trans_data = ; //network transfer in MB
				tk.seq_data = ; //sequential io in seeks
			}else if(i == 18){
				for(int j=0; j<types; j++)
					tk.estTime[j] = tAdd[j];
				tk.read_data = ; //random io in seeks
				tk.trans_data = ; //network transfer in MB
				tk.seq_data = ; //sequential io in seeks
			}else if(i == 19){
				for(int j=0; j<types; j++)
					tk.estTime[j] = tShrink[j];
				tk.read_data = ; //random io in seeks
				tk.trans_data = ; //network transfer in MB
				tk.seq_data = ; //sequential io in seeks
			}else if(i == 20) {
				for(int j=0; j<types; j++)
					tk.estTime[j] = tJPEG[j];
				tk.read_data = ; //random io in seeks
				tk.trans_data = ; //network transfer in MB
				tk.seq_data = ; //sequential io in seeks
			}
			add_vertex(tk, dag.g);
		}
		//add edges to the graph
		//add edges to the graph
		int l1[4] = {1,2,3,4};
		int l2[6] = {5,6,7,8,9,10};
		int l3[4] = {13,14,15,16};
		for(int i=0; i<4; i++) //for the dummy entry node
			add_edge(0,l1[i],dag.g);
		add_edge(20,21,dag.g); //for the exit dummy node
		for(int i=0; i<4; i++) {
			add_edge(l1[i],l3[i],dag.g);
			add_edge(l1[i],l2[i],dag.g);
			add_edge(l1[i],l2[i]+1,dag.g);
		}
		add_edge(4,10,dag.g);
		add_edge(2,5,dag.g);
		add_edge(2,9,dag.g);
		for(int i=0; i<6; i++)
			add_edge(l2[i],11,dag.g);
		add_edge(11,12,dag.g);
		for(int i=0; i<4; i++) {
			add_edge(12,l3[i],dag.g);
			add_edge(l3[i],17,dag.g);
		}
		add_edge(17,18,dag.g);
		add_edge(18,19,dag.g);
		add_edge(19,20,dag.g);		
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

	//offline optimization
	SearchPrune* optimizer = new SearchPrune(dag);
	optimizer->OfflineSP();

	return 0;
}
