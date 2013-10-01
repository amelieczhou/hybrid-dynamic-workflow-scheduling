#include <vector>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <vector>
#include <boost/random/gamma_distribution.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/graph/astar_search.hpp>
using namespace boost;

#ifndef INSTANCECONFIG_H_
#define INSTANCECONFIG_H_
//const int depth = 5; //only odd nums
//const int width = 10;//3;
//const double deadline = 1800;
const int types = 4;//types for spotvm and vm
const int randomsize = 10;//the size of random number generated
extern double OnDemandLag;//0.5;
extern double SpotLag;//1;
extern bool NoSpotVM;
extern double Times[4][types];
extern double lambda;
extern int num_jobs;
extern int num_monte;
//const int ConfigSlot = 3;
#endif

enum Integer_vm
{
	not_ready = 0,
	ready = 1,
	scheduled = 2,
	finished = 3
};
enum DAG_type
{
	montage = 0
};


struct taskVertex
{
    int name;
	int type;
    bool mark;
    double start_time;
    double end_time;
    double dl;
    double readyCountdown;
    double* estTime; //expected time, used for deadline assign and configuration
	double* cpuTime;
	double* actTime; //actual time, vary with cloud dynamics
    double restTime;
    Integer_vm status;
    int config;
	double cost; //the monetary cost spent on this task
	double tasktime; //the execution time spent on this task
	double taskstart;
	int assigned_type; 
		
	int* configList;
	int vmID;
	double* prices;
//	int* LV; //load vector
	double read_data; //I/O data in #seeks
	double seq_data; //sequential I/O data in MBytes
	double trans_data; //network data in MBytes
	double rec_data; //network data in MBytes

	double randomIO[types][randomsize];
	double seqIO[types][randomsize];
	double netUp[types][randomsize];
	double netDown[types][randomsize];
	double probestTime[types][randomsize];
	double cumulativeTime[types][randomsize];
	bool tag;

	void instance_config();
};

class VM{
public:
	taskVertex* tk;
	int type;
	int life_time;
	double turn_on; //turn on time

	int has_data;//has data of task has_data on local disk
};

class SpotVM{
public:
	taskVertex* tk;
	bool canAlloc;
	int type;
	double price;
	double life_time;//life_time<1.0
	double turn_on; //turn on time

	int has_data;//has data of task has_data on local disk

	SpotVM(double p) { price = p; canAlloc = true;}
};

typedef adjacency_list<boost::vecS, boost::vecS, bidirectionalS, taskVertex, property<edge_weight_t, double> > Graph;
typedef graph_traits<Graph>::vertex_descriptor Vertex;
typedef graph_traits<Graph>::vertex_iterator vertex_iter;
typedef graph_traits<Graph>::out_edge_iterator out_edge_iterator;
typedef graph_traits<Graph>::in_edge_iterator in_edge_iterator;
typedef graph_traits<Graph>::edge_descriptor edge_descriptor;
typedef graph_traits<Graph>::adjacency_iterator adja_iterator;
typedef graph_traits<Graph>::edge_descriptor edge_descriptor;

class DAG
{
public:
	Graph g;
	double deadline;
	double meet_dl;
	DAG_type type;
	double arrival_time;

	DAG() {};
	DAG(double d,double d1){deadline = d; meet_dl = d1;}
	void reset();
	DAG(Graph dag){g = dag;}
	DAG(const DAG& dag) {g=dag.g;deadline=dag.deadline;meet_dl=dag.meet_dl;type=dag.type;}
	std::vector<int> find_CP();
	void initDLAssign();
	void deadline_assign();
};
class Dyna{
	
};
struct configstack{
	//this struct is used to store the searched configuration list in the stack, and maintain their parent and child info
	std::vector<int> configurations;
	//configstack* parent; //parent on the search path
	//configstack* child; //child on the search path
	int taskno; //the task modified to form this configuration 
	//int nextsearchedtype; //the type searched for the next task
	bool* childcolor; //mark which child has been searched
};
//typedef adjacency_list<vecS, vecS, bidirectionalS, DAG, property<edge_weight_t, double> > DAGGraph;//each task is a dag

//template <class DAGGraph, class CostType> 
//class distance_heuristic : public astar_search<DAGGraph, CostType>{
//public:
//	 typedef typename graph_traits<Graph>::vertex_descriptor Vertex;
//
//};
class SearchPrune{
public:
	DAG dag;
	std::vector<DAG> dags;//for online simulation

	//use the aStar to search
	//use the expected cost to prune
	void OfflineSP_BFS();//go through all the types of one task first, go through different tasks later
	void OfflineSP_DFS();//go through each task first, go through different instance types later
	void SpotTune();
	void OnlineSimulate();

//	SearchPrune();
//	SearchPrune(DAG input){dag = input;}
};

bool function(double bid_price, double compare_price);
double rn_01();
int rn_integers(int a, int b);
bool myfunction(taskVertex* a, taskVertex* b);
double estimateCost(DAG dag, int start, bool estimate); //estimate or calculate the total cost of dag, starting from the start task
void estimateTime(DAG dag, double* estTime); //estimate the total execution time of dag, used in offline search
double estimateTime2(DAG dag); //estimate the execution time with estTime
void estimateTimeSpot(DAG dag, double* estTime);//used in spot tune