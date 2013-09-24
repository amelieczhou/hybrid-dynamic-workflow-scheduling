#include "Algorithms.h"
#include <boost/random/variate_generator.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <time.h>

using namespace boost;


void SearchPrune::OfflineSP(){
	DAG dag = this->dag;
	
	//for the performance of each instance type
	double random_sequential_io[types][randomsize];
	double random_random_io[types][randomsize];
	double random_network_up[types][randomsize];
	double random_network_down[types][randomsize];
	
	//generate random number for sequential io
	variate_generator<mt19937,gamma_distribution<> > generator(mt19937(time(0)), seq_io_s);
	for(int j=0; j<randomsize; j++)
		random_sequential_io[0][j] = generator();
	variate_generator<mt19937,gamma_distribution<> > generator(mt19937(time(0)), seq_io_m);
	for(int j=0; j<randomsize; j++)
		random_sequential_io[1][j] = generator();
	variate_generator<mt19937,gamma_distribution<> > generator(mt19937(time(0)), seq_io_l);
	for(int j=0; j<randomsize; j++)
		random_sequential_io[2][j] = generator();
	variate_generator<mt19937,gamma_distribution<> > generator(mt19937(time(0)), seq_io_x);
	for(int j=0; j<randomsize; j++)
		random_sequential_io[3][j] = generator();
	
	//generate random number for random io
	variate_generator<mt19937,normal_distribution<> > generator(mt19937(time(0)), r_norm_s);
	for(int j=0; j<randomsize; j++)
		random_random_io[0][j] = generator();
	variate_generator<mt19937,normal_distribution<> > generator(mt19937(time(0)), r_norm_m);	
	for(int j=0; j<randomsize; j++)
		random_random_io[1][j] = generator();
	variate_generator<mt19937,normal_distribution<> > generator(mt19937(time(0)), r_norm_l);
	for(int j=0; j<randomsize; j++)
		random_random_io[2][j] = generator();
	variate_generator<mt19937,normal_distribution<> > generator(mt19937(time(0)), r_norm_x);
	for(int j=0; j<randomsize; j++)
		random_random_io[3][j] = generator();

	//generate random number for upload network
	variate_generator<mt19937,gamma_distribution<> > generator(mt19937(time(0)), gamma_s_up);
	for(int j=0; j<randomsize; j++)
		random_network_up[0][j] = generator();
	variate_generator<mt19937,gamma_distribution<> > generator(mt19937(time(0)), gamma_m_up);
	for(int j=0; j<randomsize; j++)
		random_network_up[1][j] = generator();
	variate_generator<mt19937,gamma_distribution<> > generator(mt19937(time(0)), gamma_l_up);
	for(int j=0; j<randomsize; j++)
		random_network_up[2][j] = generator();
	variate_generator<mt19937,gamma_distribution<> > generator(mt19937(time(0)), gamma_x_up);
	for(int j=0; j<randomsize; j++)
		random_network_up[3][j] = generator();

	//generate random number for download network
	variate_generator<mt19937,gamma_distribution<> > generator(mt19937(time(0)), gamma_s_down);
	for(int j=0; j<randomsize; j++)
		random_network_down[0][j] = generator();
	variate_generator<mt19937,gamma_distribution<> > generator(mt19937(time(0)), gamma_m_down);
	for(int j=0; j<randomsize; j++)
		random_network_down[1][j] = generator();
	variate_generator<mt19937,gamma_distribution<> > generator(mt19937(time(0)), gamma_l_down);
	for(int j=0; j<randomsize; j++)
		random_network_down[2][j] = generator();
	variate_generator<mt19937,gamma_distribution<> > generator(mt19937(time(0)), gamma_x_down);
	for(int j=0; j<randomsize; j++)
		random_network_down[3][j] = generator();

	std::pair<vertex_iter, vertex_iter> vp;
	vp = vertices(dag.g);
	for(; vp.first != vp.second; vp.first++){
		for(int t=0; t<types; t++){
			for(int j=0; j<randomsize; j++){
				dag.g[*vp.first].probestTime[t][j] = dag.g[*vp.first].estTime[t] + dag.g[*vp.first].trans_data * random_network_up[t][j] / 8
					+ dag.g[*vp.first].rec_data * random_network_down[t][j] + dag.g[*vp.first].read_data / random_random_io[t][j] 
				+ dag.g[*vp.first].seq_data / random_sequential_io[t][j];
			}
		}
	}

	//A* Search
	//each node in the search tree is in fact a DAG, their differences are the instance type assigned to each taask

		
	//the dummiest way of implementation
	double globalBestCost = 0;
	//calculate the most expensive global expected cost
	double exeTime[randomsize];
	estimateCost(dag,);
}