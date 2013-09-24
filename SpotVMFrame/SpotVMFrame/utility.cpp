#include "InstanceConfig.h"
#include "Algorithms.h"
#include "PricingModel.h"

//distance heuristic for the 

double estimateCost(DAG dag, int start, bool estimate)
{ //expected cost
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
}

void estimateTime(DAG dag, double* exeTime){
	//run this function after each task is assigned with a instance type
	double tmpTime[randomsize];

	std::pair<vertex_iter,vertex_iter> vp;
	vp = vertices(dag.g);
	
	if(dag.type == montage){
		//the first four in a row
		for(int i=0; i<randomsize; i++)
			tmpTime[i] = 0;
	}
	for(; vp.first!=vp.second; vp.first++){
		for(int i=0; i<randomsize; i++){
			
		}	
	}

}