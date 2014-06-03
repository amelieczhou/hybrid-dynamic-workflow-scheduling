#include "stdafx.h"
#include "InstanceConfig.h"
#include "PricingModel.h"


float estimateCost(const DAG& dag, int start, bool estimate)
{   //expected cost
	float totalcost = 0;
	float onehour = 3600.0;
	std::pair<vertex_iter,vertex_iter> vp;
	vp = vertices((*dag.g));
	//do not round up to 1 hour
	vp.first += start;
	for(; vp.first !=vp.second; vp.first++){
		float taskcost = 0;
		if(!estimate) { //calculate according to assigned instance type
			int config = (*dag.g)[*vp.first].assigned_type;
			for(int i=0; i<randomsize; i++)
				taskcost += ((*dag.g)[*vp.first].probestTime[config*randomsize+i]+OnDemandLag)/onehour*priceOnDemand[config];//ceil?
			totalcost += taskcost / randomsize;
		}else{	//estimate according to the small instance type
			float taskcost = 0;
			for(int i=0; i<randomsize; i++)
				taskcost += ((*dag.g)[*vp.first].probestTime[i]+OnDemandLag)/onehour*priceOnDemand[0];//ceil?
			totalcost += taskcost / randomsize;
		}
	}
	return totalcost;
}

void estimateTime(DAG& dag, float* exeTime){
	//run this function after each task is assigned with a instance type
	//in this way, calculate as distribution
	int cumulativesize = 4000;
	std::pair<vertex_iter,vertex_iter> vp;
	vp = vertices((*dag.g));
	for(; vp.first!=vp.second; vp.first++){
		(*dag.g)[*vp.first].tag = false; //not yet calculate its cumulative time
		for(int i=0; i<cumulativesize; i++)
			(*dag.g)[*vp.first].cumulativeTime[i] = 0;
	}
	vp = vertices((*dag.g));
	float *tmpTime = (float*)malloc(cumulativesize*sizeof(float));
	for(; vp.first!=vp.second; vp.first++){
		for(int i=0; i<cumulativesize; i++)
			tmpTime[i] = 0;
		//find all its parents
		int config = (*dag.g)[*vp.first].assigned_type;
		in_edge_iterator in_i, in_end;
		edge_descriptor e;
		boost::tie(in_i, in_end) = in_edges(*vp.first, (*dag.g));
		if(in_i == in_end) {
			for(int i=0; i<cumulativesize; i++)
				(*dag.g)[*vp.first].cumulativeTime[i] = (*dag.g)[*vp.first].probestTime[config*randomsize+i]+OnDemandLag;
			(*dag.g)[*vp.first].tag = true;
		}
		else{
			for (; in_i != in_end; ++in_i) {
				e = *in_i;
				Vertex src = source(e, (*dag.g));					
				if(!(*dag.g)[src].tag) { //this task hasn't been assigned cumulative time
					printf("????\n");
					break;
				}else{
					for(int i=0; i<randomsize; i++){
						tmpTime[i]= tmpTime[i]>(*dag.g)[src].cumulativeTime[i]?tmpTime[i]:(*dag.g)[src].cumulativeTime[i];
					}
				}
			}
		}
		for(int i=0; i<randomsize; i++){
			(*dag.g)[*vp.first].cumulativeTime[i] = tmpTime[i] + (*dag.g)[*vp.first].probestTime[config*randomsize+i] + OnDemandLag;
		}
		(*dag.g)[*vp.first].tag = true;
	}
	//if(dag.type == montage)
		for(int i=0; i<randomsize; i++){
			exeTime[i] = (*dag.g)[*(vp.second-1)].cumulativeTime[i];////////////////////////////
	}
	free(tmpTime);
}

void conv(float* array1, float* array2, float* result, int length1, int length2){
	int resultlength = length1 + length2 - 1;
	for(int index=0; index < resultlength; index++){
		float tmp = 0.0;
		for(int k=0; k<length1; k++){
			if(index >= k && (index-k)<length2)
				tmp += array1[k]*array2[index-k];
		}
		result[index]=tmp;
	}
}
void calmaxdistr(float* array1, float* array2, float* result, int length1, int length2){
	int length = length1>length2?length1:length2;
	for(int i=0; i<length; i++)
		result[i] = 0.0;
	for(int i=0; i<length1; i++){
		for(int j=0; j<length2; j++){
			if(array1[i]!=0)
				printf("");
			if(i>=j){
				result[i] += array1[i]*array2[j]; 
			}else{//i<j
				result[j] += array1[i]*array2[j];
			}
		}
	}
}
//void estimateTimeSpot2(DAG& dag, float* exeTime, int taskno, float spottime,int randindex){
//	//add the spot time with the ondemand time, to guarantee deadline
//	//run this function after spot tune, when each task has a configList
//	std::pair<vertex_iter,vertex_iter> vp;
//	vp = vertices(dag.g);
//	for(; vp.first!=vp.second; vp.first++){
//		dag.g[*vp.first].tag = false; //not yet calculate its cumulative time
//		dag.g[*vp.first].cumulativeTime[randindex] = 0;
//	}
//	vp = vertices(dag.g);
//	if(randindex != -1){
//		//calc for a certain random sample
//		for(; vp.first != vp.second; vp.first ++){
//			float tmpTime = 0.0;
//			//find all its parents
//			in_edge_iterator in_i, in_end;
//			edge_descriptor e;
//			boost::tie(in_i, in_end) = in_edges(*vp.first, dag.g);
//			if(in_i == in_end) {			
//	//			if()
////				dag.g[*vp.first].cumulativeTime[randindex] = dag.g[*vp.first].spotTime[i];
//			
//			dag.g[*vp.first].tag = true;
//		}
//		else{
//		}
//	}
//	for(; vp.first!=vp.second; vp.first++){
//		float* tmpTime = (float*)malloc(randomsize*sizeof(float));
//		for(int i=0; i<randomsize; i++)
//			tmpTime[i] = 0;
//		//find all its parents
//		in_edge_iterator in_i, in_end;
//		edge_descriptor e;
//		boost::tie(in_i, in_end) = in_edges(*vp.first, dag.g);
//		if(in_i == in_end) {
//			for(int i=0; i<randomsize; i++){
//				dag.g[*vp.first].cumulativeTime[i] = dag.g[*vp.first].spotTime[i];
//			}
//			dag.g[*vp.first].tag = true;
//		}
//		else{
//			for (; in_i != in_end; ++in_i) {
//				e = *in_i;
//				Vertex src = source(e, dag.g);					
//				if(!dag.g[src].tag) { //this task hasn't been assigned cumulative time
//					printf("????\n");
//					break;
//				}else{
//					for(int i=0; i<randomsize; i++){
//						tmpTime[i]= tmpTime[i]>dag.g[src].cumulativeTime[i]?tmpTime[i]:dag.g[src].cumulativeTime[i];
//					}
//				}
//			}
//		}
//		for(int i=0; i<randomsize; i++){
//			if(dag.g[*vp.first].configList[0] != -2)
//				dag.g[*vp.first].cumulativeTime[i] = tmpTime[i] + dag.g[*vp.first].spotTime[i];
//			else dag.g[*vp.first].cumulativeTime[i] = tmpTime[i] + dag.g[*vp.first].probestTime[i];
//		}
//		dag.g[*vp.first].tag = true;
//	}
//	//if(dag.type == montage)
//		for(int i=0; i<randomsize; i++){
//			exeTime[i] = dag.g[*(vp.second-1)].cumulativeTime[i];
//	}
//}
//}