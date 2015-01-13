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

bool function(float bid_price, float compare_price)
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


bool configsortfunction(configstack* a, configstack* b)
{
	return (a->fvalue > b->fvalue); //sort from large to small
}
void DAG::reset(){
	/*int size = tasks.size();
	for(int i=0; i<size; i++)
	{
		tasks[i]->reset();
	}*/
	/*std::pair<vertex_iter, vertex_iter> vp;
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
	}*/
}

std::vector<int> DAG::find_CP()
{
	std::vector<int> cp;	
	property_map<Graph, vertex_distance_t>::type d_map = get(vertex_distance, *this->g);
	property_map<Graph, edge_weight_t>::type w_map = get(edge_weight, *g);
	typedef graph_traits<Graph>::vertex_descriptor vertex_descriptor;
	vertex_descriptor s = *(vertices(*this->g).first);

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
	std::vector<vertex_descriptor> p(num_vertices(*g));
	//all vertices start out as their own parent
	typedef graph_traits<Graph>::vertices_size_type size_type;
	for (size_type i =0; i< num_vertices(*g); ++i)
		p[i] = i;
    std::vector<float> d(num_vertices(*g));   
	//dag_shortest_paths(g, s, distance_map(d_map));
	dag_shortest_paths(*g, s, predecessor_map(&p[0]).distance_map(&d[0]));    
#endif
	
	int k = num_vertices(*g);	
	for(size_type i=(num_vertices(*g)-1);i !=0;)
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
		float max = 0;		
		for(int j=0; j<size2; j++)
			if(parents[j]->end_time > max) max = parents[j]->end_time;
		tasks[i]->start_time = max;
		tasks[i]->end_time = max + tasks[i]->estimateTime[0];//cheapest
	}
	
	for(int i=size - 2; i>=0; i--)
	{
		vector<task*> child = tasks[i]->get_child();
		int size2 = child.size();
		float max = 0;		
		for(int j=0; j<size2; j++)
			if(child[j]->start_time > max) max = child[j]->start_time;
		tasks[i]->end_time = max;
	}

	float cplength = tasks[size-1]->end_time;
	float ratio  = deadline/cplength;

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

void taskVertex::instance_config(){
	for(int i=0; i<types; i++) //VM types from cheap to expensive
		if(this->estTime[i] + OnDemandLag < (this->LFT - this->EST)) //data transfer time?? + OnDemandLag
		{
			prefer_type = i;
			break;
		}
}
void DAG::deadline_assign()//each task is associated with a sub-deadline, a start time and an end time
{
	//deadline assignment, according to minimum execution time	
	//new deadline assign method
	std::pair<vertex_iter, vertex_iter> vp; 
	vp = vertices(*this->g);
	int size = num_vertices(*this->g);
	(*this->g)[0].sub_deadline = (*this->g)[0].EST = 0;
	(*this->g)[size-1].sub_deadline = (*this->g)[size-1].LFT = this->deadline;
	(*this->g)[0].mark = (*this->g)[size-1].mark = 1;//have virtual task
	for(int i=1; i<size-1; i++)
		(*this->g)[i].mark = 0;

	for(; vp.first!=vp.second; ++vp.first)
	{
		Vertex v = *vp.first;
		(*this->g)[v].EST = EST((*this->g)[v],*this);
	}
	vp = vertices(*this->g);
	for(;vp.first!=vp.second; --vp.second)
	{
		Vertex v = *vp.second-1;
		(*this->g)[v].LFT = LFT((*this->g)[v],*this);
	}
	AssignParents(&(*this->g)[size-1],this);
}
double EST(taskVertex& tk, DAG job)
{
	double max = 0;
	double tmp;
	std::vector<Vertex> parents;

	//get parent vertices
	in_edge_iterator in_i, in_end;
	edge_descriptor e;
	for (boost::tie(in_i, in_end) = in_edges(tk.name, *job.g); in_i != in_end; ++in_i) 
	{
		e = *in_i;
		Vertex src = source(e, *job.g);		
		parents.push_back(src);
		int test;
		if((*job.g)[src].EST == 0)
			test = 1;
	}
	if(parents.size() > 0)
		max = (*job.g)[parents[0]].EST + (*job.g)[parents[0]].estTime[(*job.g)[parents[0]].prefer_type];
	else max = tk.EST;
	for(int i=1; i<parents.size(); i++)
	{
		tmp = (*job.g)[parents[i]].EST + (*job.g)[parents[i]].estTime[(*job.g)[parents[i]].prefer_type];
		if(tmp > max) max = tmp;
	}
	return max;
}
double LFT(taskVertex& tk, DAG job)
{
	double min = 0;
	double tmp;
	std::vector<Vertex> children;

	//get children vertices
	out_edge_iterator out_i, out_end;
	edge_descriptor e;
	for (boost::tie(out_i, out_end) = out_edges(tk.name, (*job.g)); out_i != out_end; ++out_i) 
	{
		e = *out_i;
		Vertex tgt = target(e,(*job.g));
		children.push_back(tgt);
	}
	if(children.size() > 0)
		min = (*job.g)[children[0]].LFT - (*job.g)[children[0]].estTime[(*job.g)[children[0]].prefer_type];
	else min = tk.LFT;
	for(int i=1; i<children.size(); i++)
	{
		tmp =  (*job.g)[children[i]].LFT - (*job.g)[children[i]].estTime[(*job.g)[children[i]].prefer_type];
		if(tmp < min) min = tmp;
	}
	return min;
}
void AssignParents(taskVertex* tk, DAG* job)
{	
	while(has_unassigned_parent(tk,job))
	{
		taskVertex* ti = tk;
		std::deque<Vertex> PCP;
		while(has_unassigned_parent(ti,job))
		{
			Vertex v = CriticalParent(ti,job);
			PCP.push_front(v);
			ti = &(*job->g)[v];
		}
		AssignPath(PCP, job);
		int size = PCP.size();
		for(int iter1=0; iter1 <size; iter1++)
		{
			taskVertex tk1 = (*job->g)[PCP[iter1]];
			//update EST for unassigned successors
			//get children vertices
			out_edge_iterator out_i, out_end;
			edge_descriptor e;
			for (boost::tie(out_i, out_end) = out_edges(tk1.name, (*job->g)); out_i != out_end; ++out_i) 
			{
				e = *out_i;
				Vertex tgt = target(e,(*job->g));
				if((*job->g)[tgt].mark == 0) {
					if(tk1.EST + tk1.restTime > (*job->g)[tgt].EST)
						(*job->g)[tgt].EST = tk1.EST + tk1.restTime;

					taskVertex* task = &(*job->g)[tgt];
					update_EST(task,job);
				}
			}
			//update LFT for unassigned predecessors
			//get parent vertices
			in_edge_iterator in_i, in_end;
			edge_descriptor e1;
			for (boost::tie(in_i, in_end) = in_edges(tk1.name, (*job->g)); in_i != in_end; ++in_i) 
			{
				e1 = *in_i;
				Vertex src = source(e1, (*job->g));	
				if((*job->g)[src].mark == 0){
					if(tk1.LFT - tk1.restTime < (*job->g)[src].LFT)
						(*job->g)[src].LFT = tk1.LFT - tk1.restTime;

				taskVertex* task = &(*job->g)[src];
				update_LFT(task,job);
				}
			}
			AssignParents(&tk1,job);
		}
	}
}
Vertex CriticalParent(taskVertex* tk, DAG* job)
{
	double max = 0;
	Vertex maxP;
	//get parent vertices
	in_edge_iterator in_i, in_end;
	edge_descriptor e;
	for (boost::tie(in_i, in_end) = in_edges(tk->name, (*job->g)); in_i != in_end; ++in_i) 
	{
		e = *in_i;
		Vertex src = source(e, (*job->g));		
		if((*job->g)[src].EST + (*job->g)[src].estTime[(*job->g)[src].prefer_type] > max && (*job->g)[src].mark == 0)
		{
			maxP = src;
			max = (*job->g)[src].EST + (*job->g)[src].estTime[(*job->g)[src].prefer_type];
		}
	}
	return maxP;
}
bool has_unassigned_parent(taskVertex* tk, DAG* job)
{
	std::vector<Vertex> unassigned_parent;
	//get parent vertices
	in_edge_iterator in_i, in_end;
	edge_descriptor e;
	for (boost::tie(in_i, in_end) = in_edges(tk->name, (*job->g)); in_i != in_end; ++in_i) 
	{
		e = *in_i;
		Vertex src = source(e, (*job->g));	
		if((*job->g)[src].mark == 0)
			return true;
	}
	return false;
}
void update_EST(taskVertex* tk, DAG* job)//update the EST of tk's unassigned children
{
	out_edge_iterator out_i, out_end;
	edge_descriptor e;
	for (boost::tie(out_i, out_end) = out_edges(tk->name, (*job->g)); out_i != out_end; ++out_i) 
	{
		e = *out_i;
		Vertex tgt = target(e,(*job->g));
		if((*job->g)[tgt].mark == 0){
			if(tk->EST+tk->estTime[tk->prefer_type]>(*job->g)[tgt].EST)
				(*job->g)[tgt].EST = tk->EST+tk->estTime[tk->prefer_type];
			update_EST(&(*job->g)[tgt],job);
		}
	}
}
void update_LFT(taskVertex* tk, DAG* job)//update the LFT of tk's unassigned parents
{
	in_edge_iterator in_i, in_end;
	edge_descriptor e;
	for (boost::tie(in_i, in_end) = in_edges(tk->name, (*job->g)); in_i != in_end; ++in_i) 
	{
		e = *in_i;
		Vertex src = source(e,(*job->g));
		if((*job->g)[src].mark == 0){
			if(tk->LFT - tk->estTime[tk->prefer_type] < (*job->g)[src].LFT)
				(*job->g)[src].LFT = tk->LFT - tk->estTime[tk->prefer_type];
			update_LFT(&(*job->g)[src],job);
		}
	}
}
void AssignPath(std::deque<Vertex> PCP, DAG* job)
{
	double longest_time = (*job->g)[PCP.back()].LFT - (*job->g)[PCP.front()].EST;
	double total_execution = 0;
	int size = PCP.size();
	for(int iter=0; iter < size; iter++)
	{
		total_execution += (*job->g)[PCP[iter]].estTime[(*job->g)[PCP[iter]].prefer_type];
	}
	for(int iter=0; iter <size; iter++)
		(*job->g)[PCP[iter]].restTime = longest_time*(*job->g)[PCP[iter]].estTime[(*job->g)[PCP[iter]].prefer_type]/total_execution; //the execution time assigned
	for(int iter=0; iter <size; iter++)
	{
		if(iter == 0)
		{
			(*job->g)[PCP[iter]].sub_deadline = (*job->g)[PCP[iter]].EST + (*job->g)[PCP[iter]].restTime;
			(*job->g)[PCP[iter]].mark = 1;
		}
		else{
			int prev_iter = iter-1;
			(*job->g)[PCP[iter]].EST = (*job->g)[PCP[prev_iter]].sub_deadline;
			(*job->g)[PCP[iter]].restTime = longest_time*(*job->g)[PCP[iter]].estTime[(*job->g)[PCP[iter]].prefer_type]/total_execution;
			(*job->g)[PCP[iter]].sub_deadline =  (*job->g)[PCP[iter]].EST + (*job->g)[PCP[iter]].restTime;
			(*job->g)[PCP[size-iter-1]].LFT = (*job->g)[PCP[size-prev_iter-1]].LFT - (*job->g)[PCP[size-prev_iter-1]].restTime;
			(*job->g)[PCP[iter]].mark = 1;
		}
	}
}


float rn_01()
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

//for MOEA
int cmp_int(const void *p_i1, const void *p_i2)
/* Compares the two integers '*p_i1' and '*p_i2'.
   Returns 0 if *p_i1 == *p_i2.
   Returns 1 if *p_i1 > *p_i2.
   Returns -1 if *p_i1 < *p_i2. */
{
     int i1 = *((int *)p_i1);
     int i2 = *((int *)p_i2);
     
     if(i1 == i2)
          return (0);

     if(i1 > i2)
          return (1);
     else
          return (-1);
}