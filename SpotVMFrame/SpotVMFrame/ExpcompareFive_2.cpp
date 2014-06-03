// Deadline_assign.cpp : Defines the entry point for the console application.
// Input: binary tree type DAG
//compare to previous version: add network and I/O operation time to all algorithms
//change estTime in order to make the instance configuration adapt to performance variance
//after discussion on 6th, Oct add runtime deadline refinement
#include "stdafx.h"
#include <windows.h>
#include <stdlib.h>
#include <vector>
#include <iostream>
#include <ctime>
#include "ReadTrace.h"
#include "PricingModel.h"
#include "Algorithms.h"
#include <boost/graph/topological_sort.hpp>


using namespace boost;
//I/O performance distribution, seek/sec
float OnDemandLag = 0;//seconds//10;//0.5;
float SpotLag = 0;//seconds//1;
bool NoSpotVM = true;
float Times[4][types] = {{120,65,38,24},{90,50,30,20},{60,35,23,17},{30,20,15,13}};
float lambda;
int num_jobs;
int randomsize;

int main(int argc, char* argv[])
{	
	//float* haha = (float*)malloc(100*sizeof(float));
	//int size = sizeof(haha);
	//int sizes = sizeof(*haha);
	//PricingModel* pm = new PricingModel();
	//pm->generateFile(0);
	//construct the DAG
	srand( (unsigned)time( NULL ) );
//	int depth = atoi(argv[2]);
//	int width = atoi(argv[3]);
	float deadline = atof(argv[2]);
//	float budget = atof(argv[5]);
//	int num_monte = atoi(argv[6]);
//	float cpu_perc = atof(argv[7]);
//	float io_perc = atof(argv[8]);
	//float network_perc = atof(argv[9]);
	//float rand_io_ratio = atof(argv[10]);
	//float exe_time = atof(argv[11]); //The parallel parameter used to modify Times[][]
	float meet_dl = atof(argv[3]); //the deadline meet rate submitted by the user
	lambda = atof(argv[4]);
	num_jobs = atoi(argv[5]);
	randomsize = atoi(argv[6]);

	/*for(int i=0; i<4; i++)
		for(int j=1; j<types; j++)
			Times[i][j]=Times[i][0]*exe_time/pow(2.0,j)+Times[i][0]*(1-exe_time);*/

	//I/O speed, seeks/sec
	float ioseq[4], iorand[4];
	//network speed, sec/8M
	float net_up[4], net_down[4];
	//float r = (float) rand()/(RAND_MAX+1);
	//if(strcmp (argv[12], "sens") == 0) ////for senstivity test, fix the DAG type and task types
	//	r = 0.5;
	//ioseq[0] = math::quantile(seq_io_s, r);
	//iorand[0] = math::quantile(r_norm_s, r);
	//net_up[0] = math::quantile(gamma_s_up, r);

	//input workflows
	//Graph dag;
	DAG* dag = new DAG(deadline,meet_dl);
	dag->g = new Graph();
	if(strcmp(argv[1],"montage") == 0)
	{
		dag->type = montage;
		//generate a montage DAG
		//unit time is second
		float tProjectPP[types] = {540,264,150,90};//{100,100,100,1};//{8,4,2,1};
		float tDiffFit[types] = {1200,600,339,200};//{100,100,100,1};//{8,4,2,1};
		float tConcatFit[types] = {2268,1122,636,372};//{100,100,100,1};//{8,4,2,1};
		float tBgModel[types] = {17244,8200,4633,2725};//{100,100,100,1};//{8,4,2,1};
		float tBackground[types] = {40,20,11,7};//{100,100,100,1};//{8,4,2,1};
		float tImgTbl[types] = {56,30,17,10};//{100,100,100,1};//{8,4,2,1};
		float tAdd[types] = {1080,179,101,59};//{100,100,100,1};//{8,4,2,1};
		float tShrink[types] = {324,150,85,50};//{100,100,100,1};//{8,4,2,1};
		float tJPEG[types] = {186,112,63,37};//{100,100,100,1};//{8,4,2,1};

		int indexs[9];
		indexs[0]=0; indexs[1]=4; indexs[2]=10; 
		indexs[3]=11; indexs[4]=12; indexs[5]=16; 
		indexs[6]=17; indexs[7]=18; indexs[8]=19;
		if(strcmp(argv[7],"Autoscaling") == 0 ){
			//for(int i=0; i<9; i++)
			//	indexs[i] += 1;
			//add dummy entry and exit
			taskVertex* tk = new taskVertex();
			tk->name = 0; 
			tk->estTime = new float[types];
			for(int i=0; i<types; i++){
				tk->cpuTime[i] = 0.0;
			}
			tk->rec_data = tk->read_data= tk->seq_data  = tk->trans_data = 0;
			add_vertex(*tk,*dag->g);
		}

		for(int i=0; i<20; i++){
			taskVertex* tk = new taskVertex();
			tk->name = i;
			if(strcmp(argv[7],"Autoscaling") == 0 )
				tk->name += 1;
			//if(i == 0 || i == 21) //two dummy tasks
			//{	
			//	for(int j=0; j<types; j++)
			//		tk.estTime[j] = 0;
			//}else 
			for(int k=0; k<types; k++)
				tk->estTime[k] = tk->actTime[k] = 0;
			if(i>=indexs[0] && i<indexs[1])	{
				for(int j=0; j<types; j++)
					tk->cpuTime[j] = tProjectPP[j];
				tk->read_data = 0; //random io in seeks
				tk->trans_data = 2070; //network transfer in MB
				tk->seq_data = 1050; //sequential io in MB
				tk->rec_data = 0; //network download data in MB
			}else if(i>=indexs[1] && i<indexs[2]){
				for(int j=0; j<types; j++)
					tk->cpuTime[j] = tDiffFit[j];
				tk->read_data = 0; //random io in seeks
				tk->trans_data = 86.6; //network transfer in MB
				tk->seq_data = 2100; //sequential io in MB
				tk->rec_data = 2100; //network download data in MB
			}else if(i == indexs[2]) {
				for(int j=0; j<types; j++)
					tk->cpuTime[j] = tConcatFit[j];
				tk->read_data = 0; //random io in seeks
				tk->trans_data = 1; //network transfer in MB
				tk->seq_data = 519.6; //sequential io in MB
				tk->rec_data = 519.6; //network download data in MB
			}else if(i == indexs[3]) {
				for(int j=0; j<types; j++)
					tk->cpuTime[j] = tBgModel[j];
				tk->read_data = 0; //random io in seeks
				tk->trans_data = 0.1; //network transfer in MB
				tk->seq_data = 1; //sequential io in MB
				tk->rec_data = 1; //network download data in MB
			}else if(i>=indexs[4] && i<indexs[5]){
				for(int j=0; j<types; j++)
					tk->cpuTime[j] = tBackground[j];
				tk->read_data = 0; //random io in seeks
				tk->trans_data = 2070; //network transfer in MB
				tk->seq_data = 2070; //sequential io in MB
				tk->rec_data = 2070; //network download data in MB
			}else if(i == indexs[5]){
				for(int j=0; j<types; j++)
					tk->cpuTime[j] = tImgTbl[j];
				tk->read_data = 0; //random io in seeks
				tk->trans_data = 0; //network transfer in MB
				tk->seq_data = 8280; //sequential io in MB
				tk->rec_data = 8280; //network download data in MB
			}else if(i == indexs[6]){
				for(int j=0; j<types; j++)
					tk->cpuTime[j] = tAdd[j];
				tk->read_data = 0; //random io in seeks
				tk->trans_data = 4950; //network transfer in MB
				tk->seq_data = 8280; //sequential io in MB
				tk->rec_data = 8280; //network download data in MB
			}else if(i == indexs[7]){
				for(int j=0; j<types; j++)
					tk->cpuTime[j] = tShrink[j];
				tk->read_data = 0; //random io in seeks
				tk->trans_data = 154.4; //network transfer in MB
				tk->seq_data = 2470; //sequential io in MB
				tk->rec_data = 2470; //network download data in MB
			}else if(i == indexs[8]) {
				for(int j=0; j<types; j++)
					tk->cpuTime[j] = tJPEG[j];
				tk->read_data = 0; //random io in seeks
				tk->trans_data = 6.8; //network transfer in MB
				tk->seq_data = 154.4; //sequential io in MB
				tk->rec_data = 154.4; //network download data in MB
			}else {
				printf("there is a task not considered in montage.\n");
				exit(1);
			}
			//for debuging, dont forget to comment it!!!!!!!!!!
			//tk.read_data = tk.trans_data = tk.seq_data = 0;
			add_vertex(*tk, *dag->g);
		}
		if(strcmp(argv[7],"Autoscaling") == 0 ){
			taskVertex* tk = new taskVertex();
			tk->name = 21; 
			tk->estTime = new float[types];
			for(int i=0; i<types; i++){
				tk->cpuTime[i] = 0.0;
			}
			tk->rec_data = tk->read_data= tk->seq_data  = tk->trans_data = 0;
			add_vertex(*tk,*dag->g);
		}

		//add edges to the graph
		//add edges to the graph
		int l1[4] = {0,1,2,3};
		int l2[6] = {4,5,6,7,8,9};
		int l3[4] = {12,13,14,15};
		if(strcmp(argv[7],"Autoscaling") == 0 ){
			add_edge(1,5,*dag->g);
			add_edge(2,5,*dag->g);
			add_edge(1,6,*dag->g);
			add_edge(2,6,*dag->g);
			add_edge(2,7,*dag->g);
			add_edge(3,7,*dag->g);
			add_edge(3,8,*dag->g);
			add_edge(4,8,*dag->g);
			add_edge(2,9,*dag->g);
			add_edge(4,9,*dag->g);
			add_edge(4,10,*dag->g);
			for(int i=0; i<6; i++){
				add_edge(l2[i]+1,11,*dag->g);
			}
			add_edge(11,12,*dag->g);
			for(int i=0; i<4;i++){
				add_edge(0,l1[i]+1,*dag->g);
				add_edge(i+1,l3[i]+1,*dag->g);
				add_edge(11+1,l3[i]+1,*dag->g);
				add_edge(l3[i]+1,17,*dag->g);
			}
			add_edge(19,20,*dag->g);
			add_edge(17,18,*dag->g);
			add_edge(18,19,*dag->g);
			add_edge(20,21,*dag->g);
		}else{
			add_edge(0,4,*dag->g);
			add_edge(1,4,*dag->g);
			add_edge(0,5,*dag->g);
			add_edge(1,5,*dag->g);
			add_edge(1,6,*dag->g);
			add_edge(2,6,*dag->g);
			add_edge(2,7,*dag->g);
			add_edge(3,7,*dag->g);
			add_edge(1,8,*dag->g);
			add_edge(3,8,*dag->g);
			add_edge(3,9,*dag->g);
			for(int i=0; i<6; i++){
				add_edge(l2[i],10,*dag->g);
			}
			add_edge(10,11,*dag->g);
			for(int i=0; i<4;i++){
				add_edge(i,l3[i],*dag->g);
				add_edge(11,l3[i],*dag->g);
				add_edge(l3[i],16,*dag->g);
			}
			add_edge(16,17,*dag->g);
			add_edge(17,18,*dag->g);
			add_edge(18,19,*dag->g);
		}
		
	}else if(strcmp(argv[1],"ligo") == 0){
		dag->type = ligo;
		
		//unit time in seconds
		float TmpltBank[types] = {5296,2913,1602,881};
		float Inspiral[types] = {31092,17100,9405,5173};
		float Thinca[types] = {153,84,47,26};
		float TrigBank[types] = {13,7,4,2};

		int indexs[8];
		indexs[0]=0; indexs[1]=9; indexs[2]=18; indexs[3]=29;
		indexs[4]=38; indexs[5]=19; indexs[6]=39; indexs[7]=20;
		if(strcmp(argv[7],"Autoscaling") == 0 ){
			taskVertex* tk = new taskVertex();
			tk->name = 0; 
			tk->estTime = new float[types];
			for(int i=0; i<types; i++){
				tk->cpuTime[i] = 0.0;
			}
			tk->rec_data = tk->read_data= tk->seq_data  = tk->trans_data = 0;
			add_vertex(*tk,*dag->g);
		}
		for(int i=0; i<40; i++) //40tasks in total
		{
			taskVertex tk;
			tk.name = i;
			if(strcmp(argv[7],"Autoscaling") == 0 )
				tk.name += 1;
			for(int k=0; k<types; k++)
				tk.estTime[k] = tk.actTime[k] = 0;
			if(i>=indexs[0] && i<indexs[1]){
				//TmpltBank
				for(int j=0; j<types; j++)
					tk.cpuTime[j] = TmpltBank[j];
				tk.read_data = 0; //random io in seeks
				tk.trans_data = 0.13; //network transfer in MB
				tk.seq_data = 1780.8; //sequential io in MB
				tk.rec_data = 0; //network download data in MB
			}else if((i>=indexs[1] && i<indexs[2])||(i>=indexs[3] && i<indexs[4])){
				//Inspiral
				for(int j=0; j<types; j++)
					tk.cpuTime[j] = Inspiral[j];
				tk.read_data = 0; //random io in seeks
				tk.trans_data = 1; //network transfer in MB
				tk.seq_data = 11002; //sequential io in MB
				tk.rec_data = 11002; //network download data in MB
			}else if(i==indexs[2] || i==indexs[5] || i==indexs[4] || i==indexs[6]){
				//Thinca
				for(int j=0; j<types; j++)
					tk.cpuTime[j] = Thinca[j];
				tk.read_data = 0; //random io in seeks
				tk.trans_data = 0.7; //network transfer in MB
				tk.seq_data = 38.4; //sequential io in MB
				tk.rec_data = 38.4; //network download data in MB
			}else if(i>=indexs[7] && i<indexs[3]){
				//TrigBank
				for(int j=0; j<types; j++)
					tk.cpuTime[j] = TrigBank[j];
				tk.read_data = 0; //random io in seeks
				tk.trans_data = 0; //network transfer in MB
				tk.seq_data = 0.7; //sequential io in MB
				tk.rec_data = 0.7; //network download data in MB
			}else {
				printf("there's a task not considered in Ligo.\n");
				exit(1);
			}
			add_vertex(tk, *dag->g);
		}
		if(strcmp(argv[7],"Autoscaling") == 0 ){
			taskVertex* tk = new taskVertex();
			tk->name = 41; 
			tk->estTime = new float[types];
			for(int i=0; i<types; i++){
				tk->cpuTime[i] = 0.0;
			}
			tk->rec_data = tk->read_data= tk->seq_data  = tk->trans_data = 0;
			add_vertex(*tk,*dag->g);
		}
		//add edges
		int l1[9] = {0,1,2,3,4,5,6,7,8};
		int l2[5] = {9,10,11,12,13};
		int l3[4] = {14,15,16,17};
		if(strcmp(argv[7],"Autoscaling") == 0 ){
			for(int i=0; i<9; i++) {
				add_edge(0,l1[i]+1,*dag->g);
				add_edge(l1[i]+1,l1[i]+10,*dag->g);
				add_edge(l1[i]+21,l1[i]+30,*dag->g);
			}
			for(int i=0; i<5; i++) {
				add_edge(l2[i]+1,19,*dag->g);
				add_edge(19,l2[i]+12,*dag->g);
				add_edge(l2[i]+21,39,*dag->g);
			}
			for(int i=0; i<4; i++){
				add_edge(l3[i]+1,20,*dag->g);
				add_edge(20,l3[i]+12,*dag->g);
				add_edge(l3[i]+21,40,*dag->g);
			}
			add_edge(39,41,*dag->g);
			add_edge(40,41,*dag->g);
		}else{
			for(int i=0; i<9; i++) {
				add_edge(l1[i],l1[i]+9,*dag->g);
				add_edge(l1[i]+20,l1[i]+29,*dag->g);
			}
			for(int i=0; i<5; i++) {
				add_edge(l2[i],18,*dag->g);
				add_edge(18,l2[i]+11,*dag->g);
				add_edge(l2[i]+20,38,*dag->g);
			}
			for(int i=0; i<4; i++){
				add_edge(l3[i],19,*dag->g);
				add_edge(19,l3[i]+11,*dag->g);
				add_edge(l3[i]+20,39,*dag->g);
			}		
		}

	}else if(strcmp(argv[1],"epigenome") == 0){
		dag->type = epigenome;

		float fastQSplit[types] = {791,434,238,133};
		float filterConstams[types] = {261,144,79,43};
		float sol2sanger[types] = {53,29,16,9};
		float fastq2bfq[types] = {149,82,45,25};
		float map[types] = {21355,11745,6460,3553};
		float mapMerge[types] = {291,160,88,48};
		float maqIndex[types] = {145,80,44,24};
		float pileup[types] = {185,102,56,31};

		int indexs[8];
		indexs[0]=0; indexs[1]=1; indexs[2]=5; 
		indexs[3]=9; indexs[4]=13; indexs[5]=17; 
		indexs[6]=18; indexs[7]=19;
		if(strcmp(argv[7],"Autoscaling") == 0 ){
			taskVertex* tk = new taskVertex();
			tk->name = 0; 
			tk->estTime = new float[types];
			for(int i=0; i<types; i++){
				tk->cpuTime[i] = 0.0;
			}
			tk->rec_data = tk->read_data= tk->seq_data  = tk->trans_data = 0;
			add_vertex(*tk,*dag->g);
		}
		for(int i=0; i<20; i++){ //20 tasks in total
			taskVertex tk;
			tk.name = i;
			if(strcmp(argv[7],"Autoscaling") == 0 )
				tk.name += 1;
			for(int k=0; k<types; k++)
				tk.estTime[k] = tk.actTime[k] = 0;

			if(i==indexs[0]){
				//fastQSplit
				for(int j=0; j<types; j++)
					tk.cpuTime[j] = fastQSplit[j];
				tk.read_data = 0; //random io in seeks
				tk.trans_data = 1696; //network transfer in MB
				tk.seq_data = 1696; //sequential io in MB
				tk.rec_data = 0; //network download data in MB
			}else if(i>=indexs[1] && i<indexs[2]){
				//filterConstams
				for(int j=0; j<types; j++)
					tk.cpuTime[j] = filterConstams[j];
				tk.read_data = 0; //random io in seeks
				tk.trans_data = 424; //network transfer in MB
				tk.seq_data = 424; //sequential io in MB
				tk.rec_data = 424; //network download data in MB
			}else if(i>=indexs[2] && i<indexs[3]){
				//sol2sanger
				for(int j=0; j<types; j++)
					tk.cpuTime[j] = sol2sanger[j];
				tk.read_data = 0; //random io in seeks
				tk.trans_data = 323; //network transfer in MB
				tk.seq_data = 420.8; //sequential io in MB
				tk.rec_data = 420.8; //network download data in MB
			}else if(i>=indexs[3] && i<indexs[4]){
				//fastq2bfq
				for(int j=0; j<types; j++)
					tk.cpuTime[j] = fastq2bfq[j];
				tk.read_data = 0; //random io in seeks
				tk.trans_data = 71; //network transfer in MB
				tk.seq_data = 323; //sequential io in MB
				tk.rec_data = 323; //network download data in MB
			}else if(i>=indexs[4] && i<indexs[5]){
				//map
				for(int j=0; j<types; j++)
					tk.cpuTime[j] = map[j];
				tk.read_data = 0; //random io in seeks
				tk.trans_data = 28.8; //network transfer in MB
				tk.seq_data = 4440; //sequential io in MB
				tk.rec_data = 4440; //network download data in MB
			}else if(i==indexs[5]){
				//mapMerge
				for(int j=0; j<types; j++)
					tk.cpuTime[j] = mapMerge[j];
				tk.read_data = 0; //random io in seeks
				tk.trans_data = 213.68; //network transfer in MB
				tk.seq_data = 221.44; //sequential io in MB
				tk.rec_data = 221.44; //network download data in MB
			}else if(i==indexs[6]){
				//maqIndex
				for(int j=0; j<types; j++)
					tk.cpuTime[j] = maqIndex[j];
				tk.read_data = 0; //random io in seeks
				tk.trans_data = 107.53; //network transfer in MB
				tk.seq_data = 214.1; //sequential io in MB
				tk.rec_data = 214.1; //network download data in MB
			}else if(i==indexs[7]){
				//pileup
				for(int j=0; j<types; j++)
					tk.cpuTime[j] = pileup[j];
				tk.read_data = 0; //random io in seeks
				tk.trans_data = 84; //network transfer in MB
				tk.seq_data = 151.82; //sequential io in MB
				tk.rec_data = 151.82; //network download data in MB
			}else{
				printf("there is a task not considered in epigenome.\n");
				exit(1);
			}
			add_vertex(tk, *dag->g);
		}
		if(strcmp(argv[7],"Autoscaling") == 0 ){
			taskVertex* tk = new taskVertex();
			tk->name = 21; 
			tk->estTime = new float[types];
			for(int i=0; i<types; i++){
				tk->cpuTime[i] = 0.0;
			}
			tk->rec_data = tk->read_data= tk->seq_data  = tk->trans_data = 0;
			add_vertex(*tk,*dag->g);
		}
		if(strcmp(argv[7],"Autoscaling") == 0 ){
			//add edges
			for(int i=0; i<4; i++){
				add_edge(1,i+2,*dag->g);
				add_edge(i+2,i+6,*dag->g);
				add_edge(i+6,i+10,*dag->g);
				add_edge(i+10,i+14,*dag->g);
				add_edge(i+14,18,*dag->g);
			}
			add_edge(19,20,*dag->g);
			add_edge(18,19,*dag->g);
			add_edge(0,1,*dag->g);
			add_edge(20,21,*dag->g);
		}else{
			//add edges
			for(int i=0; i<4; i++){
				add_edge(0,i+1,*dag->g);
				add_edge(i+1,i+5,*dag->g);
				add_edge(i+5,i+9,*dag->g);
				add_edge(i+9,i+13,*dag->g);
				add_edge(i+13,17,*dag->g);
			}
			add_edge(17,18,*dag->g);
			add_edge(18,19,*dag->g);
		}
	}else if(strcmp(argv[1],"test") == 0){
		//a pipeline with two tasks
		for(int i=0; i<2; i++){
			taskVertex tk;
			tk.name = i;
			for(int k=0; k<types; k++)
				tk.estTime[k] = tk.actTime[k] = 0;
			float pileup[types] = {185,102,56,31};
			for(int j=0; j<types; j++)
				tk.cpuTime[j] = pileup[j];
			tk.read_data = 0; //random io in seeks
			tk.trans_data = 84; //network transfer in MB
			tk.seq_data = 151.82; //sequential io in MB
			tk.rec_data = 151.82; //network download data in MB

			add_vertex(tk,*dag->g);
		}
		add_edge(0,1,*dag->g);
	}else {
		printf("Please choose a workflow type from: montage, ligo or epigenome.\n");
		exit(1);
	}
	
	if(strcmp(argv[7],"Dyna")==0){
		//offline optimization
		SearchPrune* optimizer = new SearchPrune();
		optimizer->dag = *dag;
		std::clock_t starttime = std::clock();
		optimizer->OfflineAstar();
		std::clock_t endtime = std::clock();
		float timeelapsed = (float)(endtime - starttime) / (float)CLOCKS_PER_SEC;
		printf("time for Astar search is: %.4f\n", timeelapsed);	
		int flag;
		if(strcmp(argv[8],"oracle")==0)
			flag = 0;
		else if(strcmp(argv[8],"heuristic")==0)
			flag = 1;
		else {
			printf("please select between oracle and heuristic\n");
			exit(1);
		}
		time_t start,end;
        time(&start);
		optimizer->SpotTune2(flag);
		time(&end);
		printf("time for spot tune is %.4f\n",difftime(end,start));
		Dyna* dynaOptimizer = new Dyna(&optimizer->dag);
		dynaOptimizer->Simulate();
	}else if(strcmp(argv[7],"DynaNS")==0){
		//offline optimization
		SearchPrune* optimizer = new SearchPrune();
		optimizer->dag = *dag;
		std::clock_t starttime = std::clock();
		optimizer->OfflineAstar();
		std::clock_t endtime = std::clock();
		float timeelapsed = (float)(endtime - starttime) / (float)CLOCKS_PER_SEC;
		printf("time for Astar search is: %.4f\n", timeelapsed);	

		DynaNS* dnyansOptimizer = new DynaNS(&optimizer->dag);
		dnyansOptimizer->Simulate();
	}else if(strcmp(argv[7],"SpotOnly")==0){
		SearchPrune* optimizer = new SearchPrune();
		optimizer->dag = *dag;
		optimizer->OfflineAstar();
		//every assigned type is spot
		SpotOnly* spotOptimizer = new SpotOnly(&optimizer->dag);
		time_t start,end;
		time(&start);
		spotOptimizer->bidpDeter();
		time(&end);
		printf("time for spotonly determining bidding price is %.4f\n",difftime(end,start));
		spotOptimizer->Simulate();
	}else if(strcmp(argv[7],"Autoscaling")==0){
		Autoscaling* autoOptimizer = new Autoscaling(dag);
		autoOptimizer->Simulate();
	}
	return 0;
}
