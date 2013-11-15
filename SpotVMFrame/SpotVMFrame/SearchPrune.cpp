#include "stdafx.h"
#include "ReadTrace.h"
#include "PricingModel.h"
#include "InstanceConfig.h"
#include <boost/random/variate_generator.hpp>
#include <boost/random/mersenne_twister.hpp>
//#include <boost/math/distributions/gamma.hpp>
//#include <boost/math/distributions/normal.hpp>
#include <ctime>
#include <cmath>
#include <utility>
#include <algorithm>
#include <stack>

using namespace boost;

//after this function, each task in the dag is assigned with an instance type
//the dag is associated with a distribution of total execution time
//dfs search all the tasks first, then increase the instance type
void SearchPrune::SpotTune2(int flag){
	//for each bidding price, there is a distribution f(t)=success probability 
	//read in the bidding price probability from tp.dat, 800lines, 440 columns
	//float tp[352000];//800*440
	float *firstfail = new float[760*800];
	ReadTrace* tr = new ReadTrace();
	tr->readtomem("newData/ffp.dat",firstfail);
	tr->closefile();
	delete tr;

	int spotdistrsize = 400;
	float spotdistrgap = 50;//mongtage deadline 17000//ceil(dag.deadline / spotdistrsize);
	if(dag.type == ligo) spotdistrgap = 125.0;//deadline 42000
	else if(dag.type == epigenome) spotdistrgap = 40.0;//deadline 15000

	std::pair<vertex_iter,vertex_iter> vp = vertices((*dag.g));
	for(;vp.first != vp.second; vp.first ++){
		(*dag.g)[*vp.first].configList = new int[2];
		(*dag.g)[*vp.first].configList[0] = -2;
		(*dag.g)[*vp.first].configList[1] = (*dag.g)[*vp.first].assigned_type;

		(*dag.g)[*vp.first].prices = new float[2];
		(*dag.g)[*vp.first].prices[0] = 0.0;
		(*dag.g)[*vp.first].prices[1] = priceOnDemand[(*dag.g)[*vp.first].assigned_type];
		(*dag.g)[*vp.first].vmID = 0;
		
		//initialize the probabilistic distribution of execution time for each task
		for(int i=0; i<randomsize; i++){
			for(int j=0; j<spotdistrsize; j++){
				(*dag.g)[*vp.first].randspot[i*spotdistrsize+j]=0.0;
			}
			int id = (int)std::ceil((*dag.g)[*vp.first].probestTime[i+randomsize*(*dag.g)[*vp.first].assigned_type]/spotdistrgap);
			if(id >= spotdistrsize) (*dag.g)[*vp.first].randspot[i*spotdistrsize+spotdistrsize-1] = 1.0;
			else (*dag.g)[*vp.first].randspot[i*spotdistrsize+id] = 1.0;
		}
	}
	vp = vertices((*dag.g));
	for(; vp.first != vp.second; vp.first++){
		int ondemandtype = (*dag.g)[*vp.first].assigned_type;
		//define similar tasks which can share the same bidp
		if(dag.type == montage){
			if(*vp.first == 1 || *vp.first == 2 || *vp.first == 3) {
				(*dag.g)[*vp.first].configList[0] = (*dag.g)[0].configList[0];
				(*dag.g)[*vp.first].prices[0] = (*dag.g)[0].prices[0];
				continue;
			}else if(*vp.first == 5 || *vp.first == 6 || *vp.first == 7 ||*vp.first == 8 || *vp.first == 9){
				(*dag.g)[*vp.first].configList[0] = (*dag.g)[4].configList[0];
				(*dag.g)[*vp.first].prices[0] = (*dag.g)[4].prices[0];
				continue;
			}else if(*vp.first == 13 || *vp.first == 14 || *vp.first == 15){
				(*dag.g)[*vp.first].configList[0] = (*dag.g)[12].configList[0];
				(*dag.g)[*vp.first].prices[0] = (*dag.g)[12].prices[0];
				continue;
			}
		}else if(dag.type == ligo){
			if(*vp.first>0 && *vp.first <9){
				(*dag.g)[*vp.first].configList[0] = (*dag.g)[0].configList[0];
				(*dag.g)[*vp.first].prices[0] = (*dag.g)[0].prices[0];
				continue;
			}else if(*vp.first>9 && *vp.first <18){
				(*dag.g)[*vp.first].configList[0] = (*dag.g)[9].configList[0];
				(*dag.g)[*vp.first].prices[0] = (*dag.g)[9].prices[0];
				continue;
			}else if(*vp.first >20 && *vp.first<29){
				(*dag.g)[*vp.first].configList[0] = (*dag.g)[20].configList[0];
				(*dag.g)[*vp.first].prices[0] = (*dag.g)[20].prices[0];
				continue;
			}else if(*vp.first >29 && *vp.first<38){
				(*dag.g)[*vp.first].configList[0] = (*dag.g)[29].configList[0];
				(*dag.g)[*vp.first].prices[0] = (*dag.g)[29].prices[0];
				continue;
			}
		}else if(dag.type == epigenome){
			if(*vp.first == 2 || *vp.first == 3 || *vp.first == 4){
				(*dag.g)[*vp.first].configList[0] = (*dag.g)[1].configList[0];
				(*dag.g)[*vp.first].prices[0] = (*dag.g)[1].prices[0];
				continue;
			}else if(*vp.first == 6 || *vp.first == 7 || *vp.first == 8){
				(*dag.g)[*vp.first].configList[0] = (*dag.g)[5].configList[0];
				(*dag.g)[*vp.first].prices[0] = (*dag.g)[5].prices[0];
				continue;
			}else if(*vp.first == 10 || *vp.first == 11 || *vp.first == 12){
				(*dag.g)[*vp.first].configList[0] = (*dag.g)[9].configList[0];
				(*dag.g)[*vp.first].prices[0] = (*dag.g)[9].prices[0];
				continue;
			}else if(*vp.first == 14 || *vp.first == 15 || *vp.first == 16){
				(*dag.g)[*vp.first].configList[0] = (*dag.g)[13].configList[0];
				(*dag.g)[*vp.first].prices[0] = (*dag.g)[13].prices[0];
				continue;
			}
		}
		for(int spottype = ondemandtype; spottype < types; spottype ++){
			//for each spot type, search through the bidding price
			//rowindex*440+colindex,rowindex in [0,799]
			//define the bidp for Binary Search!!!!!!!!!!!!!
			float bidp = 0.0;//initial bidding price
			float averagecost = 0;
			//if the hybrid cost is lower, calculate the distribution
			//else this bidding price does not do
			int searchbidp;
			int maxbidp = (int)((*dag.g)[*vp.first].prices[1]*0.8*1000.0);
			if(flag == 0) //oracle
				searchbidp = binary_search(firstfail,1,maxbidp,*vp.first,spottype);//this probability transition correct?
			else if(flag == 1)
				searchbidp = binary_search_heuristic(firstfail,1,maxbidp,*vp.first,spottype);//this probability transition correct?
			if(searchbidp == -1){
				continue;//go to the next spottype
			}else {
				//found a bid price
				(*dag.g)[*vp.first].configList[0] = spottype;
				(*dag.g)[*vp.first].prices[0] = (float)searchbidp*0.001;
				break;
			}
		}//for each spot type
		//if there is no spot type for the task, make sure to clean the intermediate data
	}//for each task
	free(firstfail);
}
void SearchPrune::OfflineAstar(){
	std::vector<configstack*> Openset;
	std::vector<configstack*> Closeset;
	std::vector<configstack> solutions;
	std::stack<configstack*> DAGstack;

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
	vp = vertices((*dag.g));

	int quantile = dag.meet_dl * randomsize;
	for(; vp.first != vp.second; vp.first++){
		//dag.g[*vp.first].probestTime = new float[types][randomsize];
		for(int t=0; t<types; t++){
			for(int j=0; j<randomsize; j++){
				(*dag.g)[*vp.first].netUp[t*randomsize+j] = (*dag.g)[*vp.first].trans_data * random_network_up[t*randomsize+j] / 8000;
				(*dag.g)[*vp.first].netDown[t*randomsize+j] = (*dag.g)[*vp.first].rec_data * random_network_down[t*randomsize+j] / 8000;
				(*dag.g)[*vp.first].randomIO[t*randomsize+j] = (*dag.g)[*vp.first].read_data / random_random_io[t*randomsize+j];
				(*dag.g)[*vp.first].seqIO[t*randomsize+j] = (*dag.g)[*vp.first].seq_data / random_sequential_io[t*randomsize+j];
				(*dag.g)[*vp.first].probestTime[t*randomsize+j] = (*dag.g)[*vp.first].cpuTime[t] + (*dag.g)[*vp.first].netUp[t*randomsize+j]
					+ (*dag.g)[*vp.first].netDown[t*randomsize+j] + (*dag.g)[*vp.first].randomIO[t*randomsize+j] + (*dag.g)[*vp.first].seqIO[t*randomsize+j];
				(*dag.g)[*vp.first].probestTime[t*randomsize+j] /= 60.0;
			}
			//calculate the estimate time as the expected value of the proestTime
			std::sort((*dag.g)[*vp.first].probestTime+t*randomsize,(*dag.g)[*vp.first].probestTime+(t+1)*randomsize-1);
			(*dag.g)[*vp.first].estTime[t] = (*dag.g)[*vp.first].probestTime[t*randomsize+quantile];
			printf("task: %d, type: %d, time: %f\n",*vp.first,t,(*dag.g)[*vp.first].estTime[t]);
		}
	}	
	
	//A* Search
	vp = vertices((*dag.g));
	for(; vp.first != vp.second; vp.first++){
		(*dag.g)[*vp.first].assigned_type = 0;//initially all assigned to small
	}
	if(false){
		float* exeTime1 = (float*)malloc(randomsize*sizeof(float));
		estimateTime(dag,exeTime1);
		float maxtime = 0.0, mintime = 0.0;
	
		for(int i=0; i<randomsize; i++)
			maxtime += exeTime1[i];
		maxtime /= randomsize;
		vp = vertices((*dag.g));
		for(; vp.first != vp.second; vp.first++){
			(*dag.g)[*vp.first].assigned_type = types-1;//initially all assigned to small
		}
		estimateTime(dag,exeTime1);
		for(int i=0; i<randomsize; i++)
			mintime += exeTime1[i];
		mintime /= randomsize;
	}
	//firstly, find a feasible solution and use it as lower bound
	configstack* initialstate = new configstack();
	//profile section, find the cheapest initialstate according to the execution time of tasks on different instance type
	for(vp=vertices((*dag.g)); vp.first!=vp.second; vp.first++){
		float tmp = (*dag.g)[*vp.first].estTime[0];
		int initialtype = 0;
		for(int instype=0; instype<types; instype++)
			if((*dag.g)[*vp.first].estTime[instype]*pow(2.0,instype)<tmp){
				tmp = (*dag.g)[*vp.first].estTime[instype]*pow(2.0,instype);
				initialtype = instype;
			}
		initialstate->configurations.push_back(initialtype);
	}
	initialstate->taskno = -1;
	vp=vertices((*dag.g));
	int numoftasks = (*vp.second - *vp.first);
	float globalBestCost = 1000000;
	float* exeTime = (float*)malloc(randomsize*sizeof(float));
	DAGstack.push(initialstate);
    bool continuesearch = true;
	//find a lower bound first
	do{
		estimateTime(dag,exeTime); ///////////////////////////////////////////start from certain task, look up for the known part
		int count = 0;
		for(int i=0; i<randomsize; i++){
			if(exeTime[i]<=dag.deadline)
				count ++;
		}
		float ratio = (float)count / (float)randomsize;
		//if satisfy users' requirement, select as the lower bound
		if(ratio >= dag.meet_dl){
			globalBestCost = estimateCost(dag,0,false);
			DAGstack.top()->fvalue = globalBestCost;
			continuesearch = false; //out of the while loop
			printf("initial ratio in search prune: %f\n",ratio);
		}else{
			//Closeset.push_back(DAGstack.top());
			int nexttask = DAGstack.top()->taskno + 1;
			if(nexttask < numoftasks){
				configstack* state = new configstack();
				state->configurations=DAGstack.top()->configurations;
				state->taskno = nexttask;
				for(int t=0; t<types; t++){
					if(!DAGstack.top()->childcolor[t]){ //type t has not been visited
						state->configurations[nexttask] = t;
						DAGstack.top()->childcolor[t] = true;
						DAGstack.push(state);
						break;                                          
					}
					if(t == types-1){ //all types have been visited
						DAGstack.pop();
						//????
					}
				}
				for(int i=0; i<numoftasks; i++)
					(*dag.g)[i].assigned_type = DAGstack.top()->configurations[i];      
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
	solutions.push_back(*DAGstack.top());
	Openset.push_back(initialstate);
	
	int searchcount = 0;
	while(!(Openset.empty() || searchcount > 10000)){
		//no need to sort, only find the smallest one is enough
		std::nth_element(Openset.begin(),Openset.end(),Openset.end(),configsortfunction);//sort from large to small,in order to reduce the erase complexity
		configstack* headnode = Openset.back();//back has the smallest fvalue

		//check if satisfy deadline 
		for(int i=0; i<numoftasks; i++)
			(*dag.g)[i].assigned_type = headnode->configurations[i];
		estimateTime(dag,exeTime); ///////////////////////////////////////////start from certain task, look up for the known part
        int count = 0;
        for(int i=0; i<randomsize; i++){
			if(exeTime[i]<=dag.deadline)
				count ++;
        }
        float ratio = (float)count / (float)randomsize;
        //if satisfy users' requirement, select as the lower bound
        if(ratio >= dag.meet_dl){
			if(headnode->fvalue < globalBestCost){
				globalBestCost = headnode->fvalue;
				solutions.push_back(*headnode);
				//remove headnode from openset and add it to closedset
				printf("search prune find a solution with ratio: %f\n",ratio);
			}
        }
		Openset.erase(Openset.end()-1);
		Closeset.push_back(headnode);

		//for each neighboring node
		int nexttask = headnode->taskno + 1;
		if(nexttask < numoftasks){
			for(int t=headnode->configurations[nexttask]+1; t<types; t++){
				configstack* state = new configstack();
				state->taskno = nexttask;
				state->configurations = headnode->configurations;
				state->configurations[nexttask] = t;
				//is it a feasible solution?
				//dag.g[nexttask].assigned_type = t;				
				for(int ii=0; ii<numoftasks; ii++)
					(*dag.g)[ii].assigned_type = state->configurations[ii];
				float currentcost = estimateCost(dag,0,false);
				if(currentcost >= globalBestCost || std::find(Closeset.begin(),Closeset.end(),state) != Closeset.end()){//std::binary_search(Closeset.begin(),Closeset.end(),state)){
					//just ignore this configuration
				}else{
					std::vector<configstack*>::iterator iter = std::find(Openset.begin(),Openset.end(),state);
					if(iter == Openset.end()){//state is not in Openset
					//bool found = std::binary_search(Openset.begin(),Openset.end(),state)
					//if(!found){
						state->fvalue = currentcost;
						Openset.push_back(state);
					}else{
						//state is already in Openset
						printf("when would this happen?\n");
					}
				}				
			}
		}	
		searchcount ++;
	}
	for(int i=0; i<numoftasks; i++)
		(*dag.g)[i].assigned_type = solutions.back().configurations[i];
	//calculate the cumulative time distribution of the dag
	float* cumulative=(float*)malloc(randomsize*sizeof(float));
	estimateTime(dag,cumulative);
	dag.cumulativetime = cumulative;

	free(cumulative);
	free(exeTime);
	free(random_sequential_io);
	free(random_random_io);
	free(random_network_up);
	free(random_network_down);
}

int SearchPrune::binary_search( const float firstfail[], int low, int high, int taskno, int spottype)//key1 is the onDemand cost
{
	int mid = ceil((low+high)/2.0);
	float bidp = mid * 0.001;
	int ondemandtype = (*dag.g)[taskno].assigned_type;
	int columns = 200;
	int spotdistrsize = 400;
	float spotdistrgap = 50.0;
    if(dag.type == ligo) spotdistrgap = 125.0;
    else if(dag.type == epigenome) spotdistrgap = 40.0;
	if(low>high)
		return -1; //this spot type is not suitable
	else
	{
		//first calculate cost
		float spotcost = 0.0;
		float ondemandcost = 0.0;
		int rowindex = (int)(bidp * 1000)-1;
		for(int rndindex=0; rndindex < randomsize; rndindex++){
			float spottime = (*dag.g)[taskno].probestTime[spottype*randomsize+rndindex];//exeution time on spot instance
			float ondemandtime = (*dag.g)[taskno].probestTime[ondemandtype*randomsize+rndindex];
			float p1 = 0.0, pt = 0.0;	//p1 is the probability of spot instance successfully executing this long 
			int colindex = ceil(spottime / 180.0); //one prob node every 3 minutes
			for(int tindex = 1; tindex <= colindex; tindex ++){
				p1 += firstfail[rowindex*columns*4+spottype*columns+tindex-1];
				pt += firstfail[rowindex*columns*4+spottype*columns+tindex-1]*tindex*180.0;
			}
			float expectedtime = (spottime + SpotLag)*(1-p1) + pt + (ondemandtime + OnDemandLag)*p1;
			//float succps = 1 - p1;			
			//dag.g[taskno].spotTime[rndindex] = expectedtime;
			spotcost += bidp*(spottime+SpotLag)/3600.0 + p1*priceOnDemand[ondemandtype]*(ondemandtime+OnDemandLag)/3600.0;
			ondemandcost += priceOnDemand[ondemandtype]*(ondemandtime+OnDemandLag)/3600.0;
		}
		spotcost /= randomsize;
		ondemandcost /= randomsize;
		if(spotcost > ondemandcost)
			return binary_search(firstfail,low,mid-1,taskno,spottype);
		else {
			//compare the time distribution with deadline meet rate
			float finalprob = 0.0; //prob of meeting deadline
			std::vector<std::pair<float, float> > finalspotdistr;
			for(int rndindex = 0; rndindex < randomsize; rndindex++){
				float spottime = (*dag.g)[taskno].probestTime[spottype*randomsize+rndindex];//exeution time on spot instance
				float ondemandtime = (*dag.g)[taskno].probestTime[ondemandtype*randomsize+rndindex];
				//first calculate the spot random for each rndindex
				std::vector<std::pair<float, float>> spotdistr; //first is time, second is probability
				float succp = 0.0;
				int iterindex = 0, count=0; 
				float expectedtime = 0.0;
				for(int kt=1; kt<(int)spottime; kt++){
					float timee = kt + ondemandtime;
					int index = std::floor(kt /180.0);
					if(index>iterindex){
						float prob = firstfail[rowindex*columns*4+spottype*columns+iterindex];
						succp += prob;
						expectedtime /= count;
						spotdistr.push_back(std::pair<float, float>(expectedtime,prob));
						iterindex++;
						expectedtime = 0.0;
						count = 0;
					}
					expectedtime += timee;
					count ++;				
				}
				succp = 1- succp;
				spotdistr.push_back(std::pair<float,float>(ondemandtime,succp));
				//change spotrand according to spotdistr
				printf("the spot rand of task %d before spot instances\n",taskno);
				//for(int i=0; i<spotdistrsize; i++) printf("%f\n",(*dag.g)[taskno].randspot[rndindex*spotdistrsize+i]);
				(*dag.g)[taskno].randspot[rndindex*spotdistrsize+(int)std::ceil(ondemandtime/spotdistrgap)] = 0.0;
				for(int i=0; i<spotdistr.size(); i++){
					int id = (int)std::ceil(spotdistr[i].first/spotdistrgap);
					if(id >= spotdistrsize) 
						(*dag.g)[taskno].randspot[rndindex*spotdistrsize+spotdistrsize-1] += spotdistr[i].second;
					else (*dag.g)[taskno].randspot[rndindex*spotdistrsize+id] += spotdistr[i].second;
				}
				printf("new spot rand distribution\n");
				//for(int i=0; i<spotdistrsize; i++) printf("%f\n",(*dag.g)[taskno].randspot[rndindex*spotdistrsize+i]);
				//do the convolcation
				std::pair<vertex_iter,vertex_iter> vvp = vertices((*dag.g));
				int* sizes = (int*)malloc((*vvp.second-*vvp.first)*sizeof(int));//save the distribution size for each task
				for(; vvp.first != vvp.second; vvp.first++){
					(*dag.g)[*vvp.first].tag = false;
					sizes[*vvp.first] = spotdistrsize;
					for(int i=0; i<randomsize; i++)
						(*dag.g)[*vvp.first].cumulativeTime[i]=0.0;//used as an temporary container, size is irrelevant to the content
					for(int i=0; i<spotdistrsize; i++)
						(*dag.g)[*vvp.first].cumulativeTime[i]=(*dag.g)[*vvp.first].randspot[rndindex*spotdistrsize+i];
				}

				vvp = vertices((*dag.g));
				for(; vvp.first != vvp.second; vvp.first++){
					//find all its parents
					in_edge_iterator in_i, in_end;
					edge_descriptor e;
					boost::tie(in_i, in_end) = in_edges(*vvp.first, (*dag.g));
					if(in_i == in_end) {
						(*dag.g)[*vvp.first].tag = true;
						float sumtask = 0.0;
						for(int i=0; i<spotdistrsize; i++){
							(*dag.g)[*vvp.first].cumulativeTime[i] = (*dag.g)[*vvp.first].randspot[rndindex*spotdistrsize+i];
							sumtask += (*dag.g)[*vvp.first].cumulativeTime[i];
						}
						if(sumtask>1.0)
							printf("");
					}
					else{//max operation, parents may have different cumulative sizes
						int maxparentsize = 0;
						int maxiter = 0;
						for(; in_i!=in_end; ++in_i){
							e = *in_i;
							int src = source(e,(*dag.g));
							if(maxparentsize<sizes[src]){
								maxparentsize=sizes[src];
								maxiter = src;
							}
						}
						float* maxdistr = (float*)malloc(maxparentsize*sizeof(float));
						float* tmpdistr = (float*)malloc(maxparentsize*sizeof(float));
						for(int i=0; i<maxparentsize; i++){
							maxdistr[i] = (*dag.g)[maxiter].cumulativeTime[i];
							tmpdistr[i] = 0.0;
						}
						float sumtask = 0.0;
						boost::tie(in_i, in_end) = in_edges(*vvp.first, (*dag.g));
						for (; in_i != in_end; ++in_i) {
							e = *in_i;
							Vertex src = source(e,(*dag.g));
							if(src != maxiter) {
								calmaxdistr(maxdistr,(*dag.g)[src].cumulativeTime,tmpdistr,maxparentsize,sizes[src]);
								for(int kk=0; kk<maxparentsize; kk++){
									//maxdistr[kk]=maxdistr[kk]>(*dag.g)[src].cumulativeTime[kk]?maxdistr[kk]:(*dag.g)[src].cumulativeTime[kk];
									maxdistr[kk]=tmpdistr[kk];
									//sumtask += maxdistr[kk];
								}								
							}
						}
						for(int kk=0; kk<maxparentsize; kk++)
							sumtask += maxdistr[kk];
						if(sumtask>1.0)
							printf("");
						free(tmpdistr);
						int colsize = sizes[*vvp.first]+maxparentsize-1;
						float* newdistr = (float*)malloc(colsize*sizeof(float));
						conv(maxdistr,(*dag.g)[*vvp.first].cumulativeTime,newdistr,maxparentsize,sizes[*vvp.first]);
						sumtask = 0.0;
						for(int j=0; j<colsize; j++){
							(*dag.g)[*vvp.first].cumulativeTime[j]=newdistr[j];
							sumtask += newdistr[j];
						}
						if(sumtask>1.0) 
							printf("");
						sizes[*vvp.first] = colsize;
						(*dag.g)[*vvp.first].tag = true;
						free(newdistr);
						free(maxdistr);
					}
				}//for vvp.first
				FILE* wFile;
				char line[256];
				wFile = fopen("log.txt","w");
				float sum = 0.0;
				for(int j=0; j<sizes[*(vvp.second-1)]; j++){
					sprintf(line,"cumulative time %d prob is %f\n", j, (*dag.g)[*(vvp.second-1)].cumulativeTime[j]);
					fputs(line,wFile);
					sum += (*dag.g)[*(vvp.second-1)].cumulativeTime[j];
					if(j*spotdistrgap <= dag.deadline)
						finalprob += (*dag.g)[*(vvp.second-1)].cumulativeTime[j];
				}
				fclose(wFile);
			}//rndindex
			finalprob /= (float)randomsize;
			//convolcation finished	
			if(finalprob >= dag.meet_dl)
				return mid;
			else {
				for(int i=0; i<randomsize; i++){
					for(int j=0; j<spotdistrsize; j++){
						(*dag.g)[taskno].randspot[i*spotdistrsize+j]=0.0;
					}
					int id = (int)std::ceil((*dag.g)[taskno].probestTime[randomsize*(*dag.g)[taskno].assigned_type+i]/spotdistrgap);
					if(id >= spotdistrsize) (*dag.g)[taskno].randspot[i*spotdistrsize+spotdistrsize-1] = 1.0;
					else (*dag.g)[taskno].randspot[i*spotdistrsize+id] = 1.0;
				}
				return binary_search(firstfail,mid+1,high,taskno,spottype);
			}
		}
	}
}
int SearchPrune::binary_search_heuristic(const float firstfail[], int low, int high, int taskno, int spottype){
	int mid = ceil((low+high)/2.0);
	float bidp = mid * 0.001;
	int ondemandtype = (*dag.g)[taskno].assigned_type;
	int columns = 200;
	int spotdistrsize = 400;
	float spotdistrgap = 50.0;
   	if(dag.type == ligo) spotdistrgap = 125.0;
	else if(dag.type == epigenome) spotdistrgap = 40.0;
	if(low>high)
		return -1; //this spot type is not suitable
	else
	{
		//first calculate cost
		float spotcost = 0.0;
		float ondemandcost = 0.0;
		int rowindex = (int)(bidp * 1000)-1;
		for(int rndindex=0; rndindex < randomsize; rndindex++){
			float spottime = (*dag.g)[taskno].probestTime[spottype*randomsize+rndindex];//exeution time on spot instance
			float ondemandtime = (*dag.g)[taskno].probestTime[ondemandtype*randomsize+rndindex];
			float p1 = 0.0, pt = 0.0;	//p1 is the probability of spot instance successfully executing this long 
			int colindex = ceil(spottime / 180.0); //one prob node every 3 minutes
			for(int tindex = 1; tindex <= colindex; tindex ++){
				p1 += firstfail[rowindex*columns*4+spottype*columns+tindex-1];
			}
			spotcost += bidp*(spottime+SpotLag)/3600.0 + p1*priceOnDemand[ondemandtype]*(ondemandtime+OnDemandLag)/3600.0;
			ondemandcost += priceOnDemand[ondemandtype]*(ondemandtime+OnDemandLag)/3600.0;
		}
		spotcost /= randomsize;
		ondemandcost /= randomsize;
		if(spotcost > ondemandcost)
			return binary_search_heuristic(firstfail,low,mid-1,taskno,spottype);
		else {
			//compare the time distribution with dag.cumulative
			std::vector<std::pair<float, float> > spotdistr; //first is time, second is probability
			for(int rndindex = 0; rndindex < randomsize; rndindex++){
				float spottime = (*dag.g)[taskno].probestTime[spottype*randomsize+rndindex];//exeution time on spot instance
				float ondemandtime = (*dag.g)[taskno].probestTime[ondemandtype*randomsize+rndindex];
				//first calculate the spot random for each rndindex
				float succp = 0.0;
				//	for(int kt=1; kt<(int)spottime; kt++){
			//		float time = kt + ondemandtime;
			//		int index = std::floor(kt /180.0);
			//		float prob = firstfail[rowindex*columns*4+spottype*columns+index];
			//		spotdistr.push_back(std::pair<float,float>(time,prob));
			//		succp += prob;
			//	}
				int iterindex = 0, count=0; 
				float expectedtime = 0.0;
				for(int kt=1; kt<(int)spottime; kt++){
					float timee = kt + ondemandtime;
					int index = std::floor(kt /180.0);
					if(index>iterindex){
						float prob = firstfail[rowindex*columns*4+spottype*columns+iterindex];
						succp += prob;
						expectedtime /= count;
						spotdistr.push_back(std::pair<float, float>(expectedtime,prob));
						iterindex++;
						expectedtime = 0.0;
						count = 0;
					}
					expectedtime += timee;
					count ++;				
				}
				succp = 1- succp;
				spotdistr.push_back(std::pair<float,float>(ondemandtime,succp));
			}//rndindex
			float maxondemandrun=0.0;
			float* originalcdf = (float*)malloc(spotdistrsize*sizeof(float));
			float* newcdf = (float*)malloc(spotdistrsize*sizeof(float));
			for(int i=0; i<spotdistrsize; i++)
				originalcdf[i] =  newcdf[i] = 0.0;
			for(int i=0; i<randomsize; i++){
				int index = (int)std::floor((*dag.g)[taskno].probestTime[ondemandtype*randomsize+i]/spotdistrgap );
				//if(index >= columns) originalcdf[columns-1] += 1.0;
				for(int j=index; j<spotdistrsize; j++)
					originalcdf[j] += 1.0;
			}

			for(int i=0; i<spotdistr.size(); i++){
				int index = (int)std::floor(spotdistr[i].first/spotdistrgap );
				//if(index >=columns) newcdf[columns-1] += spotdistr[i].second;
				for(int j=index; j<spotdistrsize; j++)
					newcdf[j] += spotdistr[i].second;
			}
			for(int i=0; i<spotdistrsize; i++){
				originalcdf[i] /= randomsize;
				newcdf[i] /= randomsize;
				if(originalcdf[i] > newcdf[i]){
					free(originalcdf);
					free(newcdf);
					return binary_search_heuristic(firstfail,mid+1,high,taskno,spottype);
				}
			}
			free(originalcdf);
			free(newcdf);			
			return mid;
		}
	}
}
