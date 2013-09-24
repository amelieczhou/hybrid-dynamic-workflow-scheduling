#include "InstanceConfig.h"

#include <boost/random/gamma_distribution.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/graph/astar_search.hpp>
//#include <boost/math/distributions/gamma.hpp>
//#include <boost/math/distributions/normal.hpp>

//I/O performance distribution, seek/sec

boost::normal_distribution<> r_norm_s(150.28, 49.98);
boost::normal_distribution<> r_norm_m(120.31, 22.45);
boost::normal_distribution<> r_norm_l(172.85, 34.77);
boost::normal_distribution<> r_norm_x(1034.04, 146.41);
//sequential I/O, MBytes/sec
boost::gamma_distribution<> seq_io_s(129.28,0.792);
boost::gamma_distribution<> seq_io_m(127.14,0.802);
boost::gamma_distribution<> seq_io_l(376.57,0.281);
boost::gamma_distribution<> seq_io_x(408.11,0.264);
//network upload and download performance distribution
boost::gamma_distribution<> gamma_s_up(2.077,604.38);
boost::gamma_distribution<> gamma_m_up(0.812, 895.37);
boost::gamma_distribution<> gamma_l_up(1.12, 528.68);
boost::gamma_distribution<> gamma_x_up(1.509, 435.16);
boost::gamma_distribution<> gamma_s_down(2.361, 399.9);
boost::gamma_distribution<> gamma_m_down(0.727, 654.3);
boost::gamma_distribution<> gamma_l_down(0.316, 126.66);
boost::gamma_distribution<> gamma_x_down(0.179, 127.87);

class Dyna{
	
};

typedef adjacency_list<vecS, vecS, bidirectionalS, Graph, property<edge_weight_t, double> > DAGGraph;//each task is a dag

template <class DAGGraph, class CostType> 
class distance_heuristic : public astar_search<DAGGraph, CostType>{
public:
	 typedef typename graph_traits<Graph>::vertex_descriptor Vertex;

};
class SearchPrune{
	DAG dag;
	std::vector<DAG> dags;//for online simulation

public:
	//use the aStar to search
	//use the expected cost to prune
	void OfflineSP();
	void OnlineSimulate();

	SearchPrune();
	SearchPrune(DAG input){dag.g = input.g;}
};

double estimateCost(DAG dag, int start, bool estimate); //estimate or calculate the total cost of dag, starting from the start task
void estimateTime(DAG dag, double* estTime); //estimate the total execution time of dag