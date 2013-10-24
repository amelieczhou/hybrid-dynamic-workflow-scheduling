#include "stdafx.h"
#include "InstanceConfig.h"
#include "PricingModel.h"


float estimateCost(const DAG& dag, int start, bool estimate)
{   //expected cost
	float totalcost = 0;
	std::pair<vertex_iter,vertex_iter> vp;
	vp = vertices(dag.g);
	//do not round up to 1 hour
	vp.first += start;
	for(; vp.first !=vp.second; vp.first++){
		float taskcost = 0;
		if(!estimate) { //calculate according to assigned instance type
			int config = dag.g[*vp.first].assigned_type;
			for(int i=0; i<randomsize; i++)
				taskcost += (dag.g[*vp.first].probestTime[config*randomsize+i]+OnDemandLag)/3600.0*priceOnDemand[config];//ceil?
			totalcost += taskcost / randomsize;
		}else{	//estimate according to the small instance type
			float taskcost = 0;
			for(int i=0; i<randomsize; i++)
				taskcost += (dag.g[*vp.first].probestTime[i]+OnDemandLag)/3600.0*priceOnDemand[0];//ceil?
			totalcost += taskcost / randomsize;
		}
	}
	return totalcost;
}

void estimateTime(DAG& dag, float* exeTime){
	//run this function after each task is assigned with a instance type
	//in this way, calculate as distribution
	std::pair<vertex_iter,vertex_iter> vp;
	vp = vertices(dag.g);
	for(; vp.first!=vp.second; vp.first++){
		dag.g[*vp.first].tag = false; //not yet calculate its cumulative time
		for(int i=0; i<randomsize; i++)
			dag.g[*vp.first].cumulativeTime[i] = 0;
	}
	vp = vertices(dag.g);
	float *tmpTime = (float*)malloc(randomsize*sizeof(float));
	for(; vp.first!=vp.second; vp.first++){
		for(int i=0; i<randomsize; i++)
			tmpTime[i] = 0;
		//find all its parents
		int config = dag.g[*vp.first].assigned_type;
		in_edge_iterator in_i, in_end;
		edge_descriptor e;
		boost::tie(in_i, in_end) = in_edges(*vp.first, dag.g);
		if(in_i == in_end) {
			for(int i=0; i<randomsize; i++)
				dag.g[*vp.first].cumulativeTime[i] = dag.g[*vp.first].probestTime[config*randomsize+i]+OnDemandLag;
			dag.g[*vp.first].tag = true;
		}
		else{
			for (; in_i != in_end; ++in_i) {
				e = *in_i;
				Vertex src = source(e, dag.g);					
				if(!dag.g[src].tag) { //this task hasn't been assigned cumulative time
					printf("????\n");
					break;
				}else{
					for(int i=0; i<randomsize; i++){
						tmpTime[i]= tmpTime[i]>dag.g[src].cumulativeTime[i]?tmpTime[i]:dag.g[src].cumulativeTime[i];
					}
				}
			}
		}
		for(int i=0; i<randomsize; i++){
			dag.g[*vp.first].cumulativeTime[i] = tmpTime[i] + dag.g[*vp.first].probestTime[config*randomsize+i] + OnDemandLag;
		}
		dag.g[*vp.first].tag = true;
	}
	//if(dag.type == montage)
		for(int i=0; i<randomsize; i++){
			exeTime[i] = dag.g[*(vp.second-1)].cumulativeTime[i];////////////////////////////
	}
	free(tmpTime);
}
float estimateTime2(DAG& dag){
	//run this function after each task is assigned with a instance type
	//in this way, calculate use the estTime only
	std::pair<vertex_iter,vertex_iter> vp;
	vp = vertices(dag.g);
	for(; vp.first!=vp.second; vp.first++){
		dag.g[*vp.first].tag = false; //not yet calculate its cumulative time
	}
	vp = vertices(dag.g);
	for(; vp.first!=vp.second; vp.first++){
		float tmpTime = 0;
		//find all its parents
		int config = dag.g[*vp.first].assigned_type;
		in_edge_iterator in_i, in_end;
		edge_descriptor e;
		boost::tie(in_i, in_end) = in_edges(*vp.first, dag.g);
		if(in_i == in_end) {			
			dag.g[*vp.first].cumulativeTime[0] = dag.g[*vp.first].estTime[config]+OnDemandLag;
			dag.g[*vp.first].tag = true;
		}
		else{
			for (; in_i != in_end; ++in_i) {
				e = *in_i;
				Vertex src = source(e, dag.g);					
				if(!dag.g[src].tag) { //this task hasn't been assigned cumulative time
					printf("????\n");
					break;
				}else{
					tmpTime= tmpTime>dag.g[src].cumulativeTime[0]?tmpTime:dag.g[src].cumulativeTime[0];
				}
			}
		}
		dag.g[*vp.first].cumulativeTime[0] = tmpTime + dag.g[*vp.first].estTime[config] + OnDemandLag;		
		dag.g[*vp.first].tag = true;
	}
	float result;
	//if(dag.type == montage)		
		result = dag.g[*(vp.second-1)].cumulativeTime[0];////////////////////////////
	return result;
}
void estimateTimeSpot(DAG& dag, float* exeTime){
	//add the spot time with the ondemand time, to guarantee deadline
	//run this function after spot tune, when each task has a configList
	std::pair<vertex_iter,vertex_iter> vp;
	vp = vertices(dag.g);
	for(; vp.first!=vp.second; vp.first++){
		dag.g[*vp.first].tag = false; //not yet calculate its cumulative time
		for(int i=0; i<randomsize; i++)
			dag.g[*vp.first].cumulativeTime[i] = 0;
	}
	vp = vertices(dag.g);
	for(; vp.first!=vp.second; vp.first++){
		float* tmpTime = (float*)malloc(randomsize*sizeof(float));
		for(int i=0; i<randomsize; i++)
			tmpTime[i] = 0;
		//find all its parents
		in_edge_iterator in_i, in_end;
		edge_descriptor e;
		boost::tie(in_i, in_end) = in_edges(*vp.first, dag.g);
		if(in_i == in_end) {
			//consider the worst case, so add(spot time, ondemand time)
			for(int i=0; i<randomsize; i++){
				dag.g[*vp.first].cumulativeTime[i] = dag.g[*vp.first].probestTime[dag.g[*vp.first].configList[1]*randomsize+i]+OnDemandLag;
				if(dag.g[*vp.first].configList[0] != -2)
					dag.g[*vp.first].cumulativeTime[i] += dag.g[*vp.first].probestTime[dag.g[*vp.first].configList[0]*randomsize+i]+SpotLag;
			}
			dag.g[*vp.first].tag = true;
		}
		else{
			for (; in_i != in_end; ++in_i) {
				e = *in_i;
				Vertex src = source(e, dag.g);					
				if(!dag.g[src].tag) { //this task hasn't been assigned cumulative time
					printf("????\n");
					break;
				}else{
					for(int i=0; i<randomsize; i++){
						tmpTime[i]= tmpTime[i]>dag.g[src].cumulativeTime[i]?tmpTime[i]:dag.g[src].cumulativeTime[i];
					}
				}
			}
		}
		for(int i=0; i<randomsize; i++){
			dag.g[*vp.first].cumulativeTime[i] = tmpTime[i] + dag.g[*vp.first].probestTime[dag.g[*vp.first].configList[1]*randomsize+i] + OnDemandLag;
			if(dag.g[*vp.first].configList[0] != -2)
				dag.g[*vp.first].cumulativeTime[i] += dag.g[*vp.first].probestTime[dag.g[*vp.first].configList[0]*randomsize+i] + SpotLag;
		}
		dag.g[*vp.first].tag = true;
	}
	//if(dag.type == montage)
		for(int i=0; i<randomsize; i++){
			exeTime[i] = dag.g[*(vp.second-1)].cumulativeTime[i];////////////////////////////
	}
}
void estimateTimeSpot2(DAG& dag, float* exeTime){
	//add the spot time with the ondemand time, to guarantee deadline
	//run this function after spot tune, when each task has a configList
	std::pair<vertex_iter,vertex_iter> vp;
	vp = vertices(dag.g);
	for(; vp.first!=vp.second; vp.first++){
		dag.g[*vp.first].tag = false; //not yet calculate its cumulative time
		for(int i=0; i<randomsize; i++)
			dag.g[*vp.first].cumulativeTime[i] = 0;
	}
	vp = vertices(dag.g);
	for(; vp.first!=vp.second; vp.first++){
		float* tmpTime = (float*)malloc(randomsize*sizeof(float));
		for(int i=0; i<randomsize; i++)
			tmpTime[i] = 0;
		//find all its parents
		in_edge_iterator in_i, in_end;
		edge_descriptor e;
		boost::tie(in_i, in_end) = in_edges(*vp.first, dag.g);
		if(in_i == in_end) {
			for(int i=0; i<randomsize; i++){
				dag.g[*vp.first].cumulativeTime[i] = dag.g[*vp.first].spotTime[i];
				if(dag.g[*vp.first].spotTime[i] == 0)
					break;
			}
			dag.g[*vp.first].tag = true;
		}
		else{
			for (; in_i != in_end; ++in_i) {
				e = *in_i;
				Vertex src = source(e, dag.g);					
				if(!dag.g[src].tag) { //this task hasn't been assigned cumulative time
					printf("????\n");
					break;
				}else{
					for(int i=0; i<randomsize; i++){
						tmpTime[i]= tmpTime[i]>dag.g[src].cumulativeTime[i]?tmpTime[i]:dag.g[src].cumulativeTime[i];
					}
				}
			}
		}
		for(int i=0; i<randomsize; i++){
			if(dag.g[*vp.first].configList[0] != -2)
				dag.g[*vp.first].cumulativeTime[i] = tmpTime[i] + dag.g[*vp.first].spotTime[i];
			else dag.g[*vp.first].cumulativeTime[i] = tmpTime[i] + dag.g[*vp.first].probestTime[i];
		}
		dag.g[*vp.first].tag = true;
	}
	//if(dag.type == montage)
		for(int i=0; i<randomsize; i++){
			exeTime[i] = dag.g[*(vp.second-1)].cumulativeTime[i];
	}
}