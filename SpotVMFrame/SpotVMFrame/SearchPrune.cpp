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
		dag.g[*vp.first].prices = new float[2];
		dag.g[*vp.first].vmID = 0;
	}
	vp = vertices(dag.g);
	float* exeTime = (float*)malloc(randomsize*sizeof(float));
	for(; vp.first != vp.second; vp.first++){
		std::vector<int*> cslot;
		int iter = 0;
		//if can reduce cost and still satisfy the deadline meet rate
		PricingModel* pm = new PricingModel();
		pm->init();
		std::vector<float*> bidps(types);
		for(int l=0; l<types; l++) //the spot instance type
		{			
			int config[3] = {l,dag.g[*vp.first].assigned_type,-1};
			float estT[3] = {dag.g[*vp.first].estTime[l],dag.g[*vp.first].estTime[dag.g[*vp.first].assigned_type],-1};
				
			bidps[l] = new float[3];//current cost, bidding price, ondemand price
			pm->getPricing(config,estT, bidps[l]);
			//bidps.push_back(bidp);
		}
		pm->finalize();
		delete pm;

		float minp=bidps[0][0], min_index=0;
		for(int j=1; j<types; j++)
			if(bidps[j][0]<minp) min_index = j;
		dag.g[*vp.first].configList[0] = min_index; 
		dag.g[*vp.first].prices[0] = bidps[min_index][1]; //may be 0
		dag.g[*vp.first].prices[1] = priceOnDemand[dag.g[*vp.first].assigned_type];
		float tmpmincost = bidps[min_index][0];

		//compare with demandonly
		int* tmpconfig = new int[3];
		tmpconfig[0] = -2;
		tmpconfig[1] = dag.g[*vp.first].assigned_type;
		tmpconfig[2] = -1;
			
		//TaskList[i]->restTime = TaskList[i]->estimateTime[TaskList[i]->configList[2]];		
		pm = new PricingModel();
		pm->init();						
		int demandconf = tmpconfig[1];

		float estT[3] = {dag.g[*vp.first].estTime[0],dag.g[*vp.first].estTime[demandconf],-1};
		float* bidp = new float[3];
		pm->getPricing(tmpconfig, estT, bidp);

		delete pm;

		if(bidp[0] <= tmpmincost){ //cost of ondemand only is less than with spot
			dag.g[*vp.first].configList[0] = -2;
			dag.g[*vp.first].prices[0] = 0;
		}else{
			//check whether satisfy deadline///////////////////////
			estimateTimeSpot(dag,exeTime); 
			int count = 0;
			for(int i=0; i<randomsize; i++){
				if(exeTime[i]<=dag.deadline)
					count ++;
			}
			float ratio = (float)count / (float)randomsize;
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
void SearchPrune::SpotTune2(){
	//for each bidding price, there is a distribution f(t)=success probability 
	//read in the bidding price probability from tp.dat, 800lines, 440 columns
	//float tp[352000];//800*440
	float *tp = new float[800*440];
	ReadTrace* tr = new ReadTrace();
	tr->readtomem("tp.dat",tp);
	tr->closefile();
	delete tr;

	//the probability that task the first time fail at time t
	//float firstfail[352000];
	float *firstfail = new float[800*440];
	for(int i=0; i<800; i++){
		for(int j=0; j<types; j++){
			//k=0
			firstfail[i*440+j*110]= 1-tp[i*440+j*110];
			//printf("%f\n",firstfail[i*440+j*110]);
			for(int k=1; k<110; k++){
				//firstfail[i*440+j*110+k] = (1-tp[i*440+j*110+k])*tp[i*440+j*110+k-1];//original version
				firstfail[i*440+j*110+k] = tp[i*440+j*110+k-1] - tp[i*440+j*110+k];
				if(firstfail[i*440+j*110+k]<-1e-12)
					firstfail[i*440+j*110+k] = 0;
				//printf("%f\n",firstfail[i*440+j*110+k]);
			}
		}
	}

	std::pair<vertex_iter,vertex_iter> vp = vertices(dag.g);
	for(;vp.first != vp.second; vp.first ++){
		dag.g[*vp.first].configList = new int[2];
		dag.g[*vp.first].configList[0] = -2;
		dag.g[*vp.first].configList[1] = dag.g[*vp.first].assigned_type;

		dag.g[*vp.first].prices = new float[2];
		dag.g[*vp.first].prices[0] = 0.0;
		dag.g[*vp.first].prices[1] = priceOnDemand[dag.g[*vp.first].assigned_type];
		dag.g[*vp.first].vmID = 0;
	}
	vp = vertices(dag.g);
	for(; vp.first != vp.second; vp.first++){
		int ondemandtype = dag.g[*vp.first].assigned_type;
		for(int spottype = ondemandtype; spottype < types; spottype ++){
			//for each spot type, search through the bidding price
			//rowindex*440+colindex,rowindex in [0,799]
			//define the bidp for Binary Search!!!!!!!!!!!!!
			float bidp = 0.0;//initial bidding price
			float averagecost = 0;
			//if the hybrid cost is lower, calculate the distribution
			//else this bidding price does not do
			int searchbidp = binary_search(firstfail,1,799,*vp.first,spottype);//this probability transition correct?
			if(searchbidp == -1)
				continue;//go to the next spottype
			else {
				//found a bid price
				dag.g[*vp.first].configList[0] = spottype;
				dag.g[*vp.first].prices[0] = (float)searchbidp*0.001;
				break;
			}
		}//for each spot type
	}//for each task
	delete tp;
	delete firstfail;
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
	
	//generate random number for sequential io
	variate_generator<mt19937,gamma_distribution<> > generator(mt19937(time(0)), seq_io_s);
	for(int j=0; j<randomsize; j++)
		random_sequential_io[j] = generator();
	variate_generator<mt19937,gamma_distribution<> > generator1(mt19937(time(0)), seq_io_m);
	for(int j=0; j<randomsize; j++)
		random_sequential_io[randomsize+j] = generator1();
	variate_generator<mt19937,gamma_distribution<> > generator2(mt19937(time(0)), seq_io_l);
	for(int j=0; j<randomsize; j++)
		random_sequential_io[2*randomsize+j] = generator2();
	variate_generator<mt19937,gamma_distribution<> > generator3(mt19937(time(0)), seq_io_x);
	for(int j=0; j<randomsize; j++)
		random_sequential_io[3*randomsize+j] = generator3();
	
	//generate random number for random io
	variate_generator<mt19937,normal_distribution<> > generator4(mt19937(time(0)), r_norm_s);
	for(int j=0; j<randomsize; j++)
		random_random_io[j] = generator4();
	variate_generator<mt19937,normal_distribution<> > generator5(mt19937(time(0)), r_norm_m);	
	for(int j=0; j<randomsize; j++)
		random_random_io[randomsize+j] = generator5();
	variate_generator<mt19937,normal_distribution<> > generator6(mt19937(time(0)), r_norm_l);
	for(int j=0; j<randomsize; j++)
		random_random_io[2*randomsize+j] = generator6();
	variate_generator<mt19937,normal_distribution<> > generator7(mt19937(time(0)), r_norm_x);
	for(int j=0; j<randomsize; j++)
		random_random_io[3*randomsize+j] = generator7();

	//generate random number for upload network
	variate_generator<mt19937,gamma_distribution<> > generator8(mt19937(time(0)), gamma_s_up);
	for(int j=0; j<randomsize; j++)
		random_network_up[j] = generator8();
	variate_generator<mt19937,gamma_distribution<> > generator9(mt19937(time(0)), gamma_m_up);
	for(int j=0; j<randomsize; j++)
		random_network_up[randomsize+j] = generator9();
	variate_generator<mt19937,gamma_distribution<> > generator10(mt19937(time(0)), gamma_l_up);
	for(int j=0; j<randomsize; j++)
		random_network_up[2*randomsize+j] = generator10();
	variate_generator<mt19937,gamma_distribution<> > generator11(mt19937(time(0)), gamma_x_up);
	for(int j=0; j<randomsize; j++)
		random_network_up[3*randomsize+j] = generator11();

	//generate random number for download network
	variate_generator<mt19937,gamma_distribution<> > generator12(mt19937(time(0)), gamma_s_down);
	for(int j=0; j<randomsize; j++)
		random_network_down[j] = generator12();
	variate_generator<mt19937,gamma_distribution<> > generator13(mt19937(time(0)), gamma_m_down);
	for(int j=0; j<randomsize; j++)
		random_network_down[randomsize+j] = generator13();
	variate_generator<mt19937,gamma_distribution<> > generator14(mt19937(time(0)), gamma_l_down);
	for(int j=0; j<randomsize; j++)
		random_network_down[2*randomsize+j] = generator14();
	variate_generator<mt19937,gamma_distribution<> > generator15(mt19937(time(0)), gamma_x_down);
	for(int j=0; j<randomsize; j++)
		random_network_down[3*randomsize+j] = generator15();

	std::pair<vertex_iter, vertex_iter> vp;
	vp = vertices(dag.g);
	int quantile = dag.meet_dl * randomsize;
	for(; vp.first != vp.second; vp.first++){
		//dag.g[*vp.first].probestTime = new float[types][randomsize];
		for(int t=0; t<types; t++){
			for(int j=0; j<randomsize; j++){
				dag.g[*vp.first].netUp[t*randomsize+j] = dag.g[*vp.first].trans_data * random_network_up[t*randomsize+j] / 8000;
				dag.g[*vp.first].netDown[t*randomsize+j] = dag.g[*vp.first].rec_data * random_network_down[t*randomsize+j] / 8000;
				dag.g[*vp.first].randomIO[t*randomsize+j] = dag.g[*vp.first].read_data / random_random_io[t*randomsize+j];
				dag.g[*vp.first].seqIO[t*randomsize+j] = dag.g[*vp.first].seq_data / random_sequential_io[t*randomsize+j];
				dag.g[*vp.first].probestTime[t*randomsize+j] = dag.g[*vp.first].cpuTime[t] + dag.g[*vp.first].netUp[t*randomsize+j]
					+ dag.g[*vp.first].netDown[t*randomsize+j] + dag.g[*vp.first].randomIO[t*randomsize+j] + dag.g[*vp.first].seqIO[t*randomsize+j];
			}
			//calculate the estimate time as the expected value of the proestTime
			std::sort(dag.g[*vp.first].probestTime+t*randomsize,dag.g[*vp.first].probestTime+(t+1)*randomsize-1);
			dag.g[*vp.first].estTime[t] = dag.g[*vp.first].probestTime[t*randomsize+quantile];
			printf("task: %d, type: %d, time: %f\n",*vp.first,t,dag.g[*vp.first].estTime[t]);
		}
	}	
	
	//A* Search
	vp = vertices(dag.g);
	for(; vp.first != vp.second; vp.first++){
		dag.g[*vp.first].assigned_type = 0;//initially all assigned to small
	}
	//firstly, find a feasible solution and use it as lower bound
	configstack* initialstate = new configstack();
	//profile section, find the cheapest initialstate according to the execution time of tasks on different instance type
	for(vp=vertices(dag.g); vp.first!=vp.second; vp.first++){
		float tmp = dag.g[*vp.first].estTime[0];
		int initialtype = 0;
		for(int instype=0; instype<types; instype++)
			if(dag.g[*vp.first].estTime[instype]*pow(2.0,instype)<tmp){
				tmp = dag.g[*vp.first].estTime[instype]*pow(2.0,instype);
				initialtype = instype;
			}
		initialstate->configurations.push_back(initialtype);
	}
	initialstate->taskno = -1;
	vp=vertices(dag.g);
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
					dag.g[i].assigned_type = DAGstack.top()->configurations[i];      
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
			dag.g[i].assigned_type = headnode->configurations[i];
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
					dag.g[ii].assigned_type = state->configurations[ii];
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
		dag.g[i].assigned_type = solutions.back().configurations[i];
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
	int ondemandtype = dag.g[taskno].assigned_type;
	if(low>high)
		return -1; //this spot type is not suitable
	else
	{
		//first calculate cost
		float spotcost = 0.0;
		float ondemandcost = 0.0;
		int rowindex = (int)(bidp * 1000);
		for(int rndindex=0; rndindex < randomsize; rndindex++){
			float spottime = dag.g[taskno].probestTime[spottype*randomsize+rndindex];//exeution time on spot instance
			float ondemandtime = dag.g[taskno].probestTime[ondemandtype*randomsize+rndindex];
			float p1 = 0.0, pt = 0.0;	//p1 is the probability of spot instance successfully executing this long 
			int colindex = ceil(spottime / 180.0); //one prob node every 3 minutes
			for(int tindex = 1; tindex <= colindex; tindex ++){
				p1 += firstfail[rowindex*440+spottype*110+tindex-1];
				pt += firstfail[rowindex*440+spottype*110+tindex-1]*tindex*180.0;
			}
			float expectedtime = (spottime + SpotLag)*(1-p1) + pt + (ondemandtime + OnDemandLag)*p1;
			float succps = 1 - p1;			
			dag.g[taskno].spotTime[rndindex] = expectedtime;
			spotcost += bidp*(spottime+SpotLag)/3600.0 + (1-succps)*priceOnDemand[ondemandtype]*(ondemandtime+OnDemandLag)/3600.0;
			ondemandcost += priceOnDemand[ondemandtype]*(ondemandtime+OnDemandLag)/3600.0;
		}
		spotcost /= randomsize;
		ondemandcost /= randomsize;
		if(spotcost > ondemandcost)
			return binary_search(firstfail,low,mid-1,taskno,spottype);
		else {
			//compare the time distribution with dag.cumulative
			float* exeTime = (float*)malloc(randomsize*sizeof(float));
			estimateTimeSpot2(dag,exeTime); 
			int count = 0;
			for(int i=0; i<randomsize; i++){
				if(exeTime[i]<=dag.deadline)
					count ++;
			}
			float ratio = (float)count / (float)randomsize;
			if(ratio < dag.meet_dl)
				return binary_search(firstfail,mid+1,high,taskno,spottype);
			else
				return mid;
		}
	}
}