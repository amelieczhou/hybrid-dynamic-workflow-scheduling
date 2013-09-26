//#include <vector>
#include "stdafx.h"
#include "InstanceConfig.h"
#include "PricingModel.h"
#include <iostream>
#include <cstdlib>
#include <boost/graph/dag_shortest_paths.hpp>
//#include <boost/graph/graph_traits.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp> 
#include <boost/graph/topological_sort.hpp>
#include <boost/random.hpp>
//using namespace std;
//extern bool NoSpotVM;
//#include <vector>

bool function(double bid_price, double compare_price)
{
	bool b;
	if(bid_price < compare_price)
		b = false;
	else b = true;
	return b;
};
bool myfunction(taskVertex* a, taskVertex* b)
{
	return (a->dl < b->dl);
}


void DAG::reset(){
	/*int size = tasks.size();
	for(int i=0; i<size; i++)
	{
		tasks[i]->reset();
	}*/
	std::pair<vertex_iter, vertex_iter> vp;
	vp = vertices(g);
	for(; vp.first != vp.second; ++vp.first)
	{
		Vertex v=*vp.first;
		g[v].readyCountdown = -1;
		g[v].status = not_ready;
		g[v].restTime = 0;
		g[v].dl = g[v].end_time = g[v].start_time = g[v].mark = 0;
		g[v].config = 0;
		g[v].vmID = 0;
		g[v].configList[0] = g[v].configList[2] = 0; g[v].configList[1] = 1;
		g[v].prices[0]=g[v].prices[1]=g[v].prices[2] = 0;
		g[v].cost = 0; g[v].tasktime = g[v].taskstart = 0;
	}
}

std::vector<int> DAG::find_CP()
{
	/*int size = tasks.size();
	for(int i=0; i<size; i++)
	{
		vector<task*> parents = tasks[i]->get_parent();
		int size2 = parents.size();
		double max = 0;		
		for(int j=0; j<size2; j++)
			if(parents[j]->end_time > max) max = parents[j]->end_time;
		tasks[i]->start_time = max;
		tasks[i]->end_time = max + tasks[i]->estimateTime[0];//cheapest
	}
	task* t;	
	tasks[size - 1]->dl = deadline;
	t = tasks[size - 1];*/
	std::vector<int> cp;
	/*cp.push_back(t->show_name());
	while(t->show_name()!=0)
	{
		vector<task*> parents = t->get_parent();	
		double max=0;
		int flag = 0;
		for(int i=0; i<parents.size(); i++)
			if(parents[i]->end_time > max)
			{
				max = parents[i]->end_time;
				flag = parents[i]->show_name();
			}
		cp.push_back(flag);
		tasks[flag]->dl = t->dl - t->estimateTime[0]/tasks[size-1]->end_time * deadline;
		t = tasks[flag]; 
	}*/
	property_map<Graph, vertex_distance_t>::type d_map = get(vertex_distance, this->g);
	property_map<Graph, edge_weight_t>::type w_map = get(edge_weight, g);
	typedef graph_traits<Graph>::vertex_descriptor vertex_descriptor;
	vertex_descriptor s = *(vertices(this->g).first);

#if defined(BOOST_MSVC) && BOOST_MSVC <= 1300
	std::vector<default_color_type> color(num_vertices(g));
	std::vector<std::size_t> pred(num_vertices(g));
	default_dijkstra_visitor vis;
	std::less<int> compare;
	closed_plus<int> combine;
	
	//dag_shortest_paths(this->g, s, distance_map(d_map));	
	dag_shortest_paths(g, s, d_map, w_map, &color[0], &pred[0], 
     vis, compare, combine, (std::numeric_limits<int>::max)(), 0);
#else
	std::vector<vertex_descriptor> p(num_vertices(g));
	//all vertices start out as their own parent
	typedef graph_traits<Graph>::vertices_size_type size_type;
	for (size_type i =0; i< num_vertices(g); ++i)
		p[i] = i;
    std::vector<double> d(num_vertices(g));   
	//dag_shortest_paths(g, s, distance_map(d_map));
	dag_shortest_paths(g, s, predecessor_map(&p[0]).distance_map(&d[0]));    
#endif
	//std::cout << "distances and parents:" << std::endl;     
	//graph_traits <Graph>::vertex_iterator vi, vend;   
	//for (tie(vi, vend) = vertices(g); vi != vend; ++vi) 
	//{
	//	std::cout << "distance(" << g[*vi].name << ") = " << d[*vi] << ", ";
	//	std::cout << "parent(" << g[*vi].name << ") = " << g[p[*vi]].name << std::endl;
	//}
	int k = num_vertices(g);
	//cp.push_back(0); // add the source task as the first node in cp
	/*for (size_type i = 0; i< num_vertices(g); ++i)
		if (p[i] != i)	{
			cp.push_back(i);
		}*/
	for(size_type i=(num_vertices(g)-1);i !=0;)
	{
		cp.push_back(i); // assuming one exit and one entrance node
		i = p[i];
	}
	//cp.pop_back();
	cp.push_back(0); // add the source task as the first node in cp
	return cp;//the last node is a virtual node
}

void DAG::initDLAssign(){
	/*int size = tasks.size();
	tasks[0]->start_time = 0;
	tasks[0]->end_time = tasks[0]->estimateTime[0];

	for(int i=1; i<size; i++)
	{
		vector<task*> parents = tasks[i]->get_parent();
		int size2 = parents.size();
		double max = 0;		
		for(int j=0; j<size2; j++)
			if(parents[j]->end_time > max) max = parents[j]->end_time;
		tasks[i]->start_time = max;
		tasks[i]->end_time = max + tasks[i]->estimateTime[0];//cheapest
	}
	
	for(int i=size - 2; i>=0; i--)
	{
		vector<task*> child = tasks[i]->get_child();
		int size2 = child.size();
		double max = 0;		
		for(int j=0; j<size2; j++)
			if(child[j]->start_time > max) max = child[j]->start_time;
		tasks[i]->end_time = max;
	}

	double cplength = tasks[size-1]->end_time;
	double ratio  = deadline/cplength;

	for(int i=0; i<size; i++)
	{
		tasks[i]->dl = ratio*(tasks[i]->end_time-tasks[i]->start_time);
		tasks[i]->start_time = ratio*(tasks[i]->start_time);
		tasks[i]->end_time = ratio*(tasks[i]->end_time);
	}
	
	task*t;
	t = tasks[0];
*/
}

void DAG::deadline_assign()
{
	//deadline assignment
	std::vector<int> critical_path = this->find_CP(); 	
	double longest_time = 0, cheapest_cost =0;
	for (int i=0; i<critical_path.size(); i++) {
		longest_time += (this->g[critical_path[i]].estTime[0]);
		//cheapest_cost += (this->g[critical_path[i]].estTime[0] + this->g[critical_path[i]].read_data / iospeed[0] 
		//+ this->g[critical_path[i]].trans_data * net_up[0]/8 + this->g[critical_path[i]].rec_data *net_down[0] /8)*priceOnDemand[0] / 60.0;
	}
	//proportional assign deadline for tasks on critical path
	double left_time = deadline;
	for (int i = 0; i <critical_path.size(); i++)
	{		
		this->g[critical_path[i]].dl = left_time;
		this->g[critical_path[i]].end_time = left_time;
		left_time -= (this->g[critical_path[i]].estTime[0] )/longest_time * deadline;
		if(left_time < 0) left_time = 0;
		this->g[critical_path[i]].start_time = left_time;
		this->g[critical_path[i]].mark = 1; //already assigned deadline
		this->g[critical_path[i]].configList[0] = this->g[critical_path[i]].configList[2] = 0; //initially assigned to the cheapest machine
		this->g[critical_path[i]].configList[1] = 1;
	}
	//topological sort
	typedef std::vector<Vertex> container;
	container c; //the nodes will be reversely stored accordingto their orders
	topological_sort(this->g, std::back_inserter(c)); //depth firs search
	
	for (container::iterator ii=c.begin(); ii != c.end(); ++ii)
	{
		out_edge_iterator out_i, out_end;
		edge_descriptor e;
		Vertex v = *ii;
		for (boost::tie(out_i, out_end) = out_edges(v, this->g); out_i != out_end; ++out_i) 
		{
			adja_iterator ai, ai_end;
			typedef property_map<Graph, vertex_index_t>::type id ;
			id vertex_id = get(vertex_index, this->g);
			boost::tie(ai,ai_end) = adjacent_vertices(v, this->g);
			if(ai != ai_end && this->g[v].mark != 1) {
				double earliest_Start = this->g[get(vertex_id, *ai)].start_time;
				for (; ai != ai_end; ++ai)
				{
					if(this->g[get(vertex_id, *ai)].start_time < earliest_Start)
						earliest_Start = this->g[get(vertex_id, *ai)].start_time;
				}
				this->g[v].end_time = this->g[v].dl = earliest_Start;
				this->g[v].start_time = this->g[v].dl - (this->g[v].estTime[0])*deadline/longest_time;
				this->g[v].mark = 1;
				this->g[v].configList[0] = this->g[v].configList[2] = 0; this->g[v].configList[1] = 1;
			}
		}
	}
}

void taskVertex::instance_config()
{	
	if((this->dl - this->start_time) < this->estTime[types-1]+OnDemandLag)//
	{
		//std::cout<<"fail for task config, cannot find such config!"<<std::endl;
		//std::cout<<"fail for config"<<std::endl;
		this->configList = new int[3];
		this->configList[0] = this->configList[1]=-2;
		this->configList[2] = types-1;
		this->prices = new double[3];
		this->prices[0] = this->prices[1] = 0;
		//this->prices[2] = priceOnDemand[types-1];
		PricingModel* pm = new PricingModel();
		pm->init();
		int config[4] = {-2,-2,this->configList[2],-1};
				
//		int demandconf = config[2];

		double estT[4] = {this->estTime[0],this->estTime[0],this->estTime[types-1],-1};
		double* bidp = new double[4];
		pm->getPricing(config,estT, bidp);
		this->prices[2] = bidp[3];

		delete pm; delete[] bidp;
		return ;	
		//exit(1);
	}
	std::vector<int*> cslot;
	//int count = 0;
	int iter = 0;
	
	
	int starttype = 0;
	for(int j=0; j<types; j++)
		if(this->estTime[j] < (this->dl - this->start_time ))
		{
			starttype = j;
			break;
		}
			

	for(int l=starttype; l<types; l++) //the on demand vm
	{
		for(int k=l; k<types; k++)//the 1st spotvm
		{
			double est_T = SpotLag + OnDemandLag;
			est_T += this->estTime[k] + this->estTime[l];
			if(est_T < (this->dl - this->start_time ))// + SpotLag + OnDemandLag 
				{
					int* slot = new int[3];
					slot[0] = -2; slot[1]=k; slot[2]= l;
					cslot.push_back(slot);
					iter ++;
					if(cslot.size() == 20) break;
				}

			for(int j=k+1; j<types; j++) //the 2nd spotvm
			{
				//j += 3; k += 3; //in total vms
				double est_T = 2*SpotLag + OnDemandLag;
				est_T += this->estTime[j] + this->estTime[k] + this->estTime[l];
				if(est_T < (this->dl - this->start_time ))// + SpotLag + OnDemandLag 
				{
					int* slot = new int[3];
					//cslot.resize(count+1);
					//cslot[count] = new int[3];
					//cslot[count][0]=j; cslot[count][1] = k; cslot[count][2]= l;//
					slot[0] = k; slot[1]=j; slot[2]= l;
					cslot.push_back(slot);
					//count++;
					iter ++;
					if(cslot.size() == 20) break;
				}
			}
			if(cslot.size() == 20) break;
		}
		if(cslot.size() == 20 || iter > 1000) break;			
	}
	int count = cslot.size();
	if(count > 0 ) //&& !NoSpotVM)
	{
		//TaskList[i]->configList = cslot[0];
		//set prices for configList/////////////////
		//TaskList[i]->prices.resize(3, 1.0);
		PricingModel* pm = new PricingModel();
		pm->init();
		std::vector<double*> bidps(count);
		for(int j=0; j<count; j++)
		{
			int config[4] = {cslot[j][0],cslot[j][1],cslot[j][2],-1};
			double estT[4] = {this->estTime[cslot[j][0]],this->estTime[cslot[j][1]],this->estTime[cslot[j][2]],-1};
			if(cslot[j][0] == -2) {	
				estT[0] = this->estTime[0];
			}		

			bidps[j] = new double[4];
			pm->getPricing(config,estT, bidps[j]);
			//bidps.push_back(bidp);
		}
		pm->finalize();
		delete pm;

		double minp=bidps[0][0], min_index=0;
		for(int j=1; j<count; j++)
			if(bidps[j][0]<minp) min_index = j;
		this->configList[0] = cslot[min_index][0]; //may be -2
		this->configList[1] = cslot[min_index][1];
		this->configList[2] = cslot[min_index][2];
		//TaskList[i]->restTime = TaskList[i]->estimateTime[TaskList[i]->configList[0]];
		//this->prices = new double[3];
		this->prices[0] = bidps[min_index][1]; //may be 0
		this->prices[1] = bidps[min_index][2];
		this->prices[2] = bidps[min_index][3];
		double tmpmincost = bidps[min_index][0];

		//compare with demandonly
		int* tmpconfig = new int[4];
		tmpconfig[0] = tmpconfig[1] = -2;
		tmpconfig[3] = -1;
		tmpconfig[2] = types-1;
		for(int j=0; j<types; j++)
			if(this->estTime[j] < (this->dl - this->start_time ))
			{
				tmpconfig[2] = j;
				break;
			}
		//TaskList[i]->restTime = TaskList[i]->estimateTime[TaskList[i]->configList[2]];		
		pm = new PricingModel();
		pm->init();						
		int demandconf = tmpconfig[2];

		double estT[4] = {this->estTime[0],this->estTime[0],this->estTime[demandconf],-1};
		double* bidp = new double[4];
		pm->getPricing(tmpconfig, estT, bidp);

		delete pm;

		if(bidp[0] <= tmpmincost){
			this->configList[0] = this->configList[1] = -2;
			this->configList[2] = tmpconfig[2];
			//debug
			if(tmpconfig[2]>3 || tmpconfig[2]<0)
			{
				printf("getPricing error\n");
			}
			
			this->prices[0] = this->prices[1] = 0;
			this->prices[2] = bidp[3];		
		}

		delete[] bidp; delete[] tmpconfig;
		for(int i=0; i<count; i++)
		{			 
			delete[] bidps[i];
			bidps.erase(bidps.begin() +i);
			count--;
			i--;
		}

	}
	else {
//			cout<<"cannot find any config for task "<< i <<endl;
		//this->configList = new int[3];
		this->configList[0] = this->configList[1] = -2;
		for(int j=0; j<types; j++)
			if(this->estTime[j] < (this->dl - this->start_time))
			{
				this->configList[2] = j;
				break;
			}
		//TaskList[i]->restTime = TaskList[i]->estimateTime[TaskList[i]->configList[2]];
		//this->prices = new double[3];
		this->prices[0] = this->prices[1] = 0;
		PricingModel* pm = new PricingModel();
		pm->init();
		int config[4] = {-2,-2,this->configList[2],-1};
				
		int demandconf = config[2];

		double estT[4] = {this->estTime[0],this->estTime[0],this->estTime[demandconf],-1};
		double* bidp = new double[4];
		pm->getPricing(config,estT, bidp);
		this->prices[2] = bidp[3];

		delete pm; delete[] bidp;
	}

	count = cslot.size();
	for(int i=0; i<count; i++)
	{
		delete[] cslot[i];
		cslot.erase(cslot.begin() +i);
		count--;
		i--;
	}
}

double rn_01()
{
	boost::mt19937 rng(43);
	static boost::uniform_01<boost::mt19937> zeroone(rng);
	return zeroone();
}

int rn_integers(int a, int b)
{
	boost::mt19937 rng;
	mt19937::result_type random_seed= static_cast<mt19937::result_type>(time(0));
	rng.seed(random_seed);
	rng.seed(static_cast<mt19937::result_type>(random_seed));
	static uniform_int<> type(a, b);
	variate_generator< mt19937&, uniform_int<> > gettype(rng, type);
	return gettype();
}
