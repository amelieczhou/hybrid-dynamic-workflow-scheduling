#include "stdafx.h"
#include "InstanceConfig.h"
#include "PricingModel.h"


double estimateCost(DAG dag, int start, bool estimate)
{   //expected cost
	double totalcost = 0;
	std::pair<vertex_iter,vertex_iter> vp;
	vp = vertices(dag.g);
	//do not round up to 1 hour
	vp.first += start;
	for(; vp.first !=vp.second; vp.first++){
		double taskcost = 0;
		if(!estimate) { //calculate according to assigned instance type
			int config = dag.g[*vp.first].assigned_type;
			for(int i=0; i<randomsize; i++)
				taskcost += dag.g[*vp.first].probestTime[config][i]/60.0*priceOnDemand[config];
			totalcost += taskcost / randomsize;
		}else{	//estimate according to the small instance type
			double taskcost = 0;
			for(int i=0; i<randomsize; i++)
				taskcost += dag.g[*vp.first].probestTime[0][i]/60.0*priceOnDemand[0];
			totalcost += taskcost / randomsize;
		}
	}
	return totalcost;
}

void estimateTime(DAG dag, double* exeTime){
	//run this function after each task is assigned with a instance type
	std::pair<vertex_iter,vertex_iter> vp;
	vp = vertices(dag.g);
	for(; vp.first!=vp.second; vp.first++)
		dag.g[*vp.first].tag = false; //not yet calculate its cumulative time
	vp = vertices(dag.g);
	for(; vp.first!=vp.second; vp.first++){
		double tmpTime[randomsize];
		for(int i=0; i<randomsize; i++)
			tmpTime[i] = 0;
		//find all its parents
		int config = dag.g[*vp.first].assigned_type;
		in_edge_iterator in_i, in_end;
		edge_descriptor e;
		boost::tie(in_i, in_end) = in_edges(*vp.first, dag.g);
		if(in_i == in_end) {
			for(int i=0; i<randomsize; i++)
				dag.g[*vp.first].cumulativeTime[config][i] = dag.g[*vp.first].probestTime[config][i];
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
						tmpTime[i]= tmpTime[i]>dag.g[src].cumulativeTime[dag.g[src].assigned_type][i]?tmpTime[i]:dag.g[src].cumulativeTime[dag.g[src].assigned_type][i];
					}
				}
			}
		}
		for(int i=0; i<randomsize; i++){
			dag.g[*vp.first].cumulativeTime[config][i] = tmpTime[i] + dag.g[*vp.first].probestTime[config][i];
		}
		dag.g[*vp.first].tag = true;
	}
	if(dag.type == montage)
		for(int i=0; i<randomsize; i++){
			exeTime[i] = dag.g[*(vp.second-1)].cumulativeTime[dag.g[*(vp.second-1)].assigned_type][i];////////////////////////////
	}
}