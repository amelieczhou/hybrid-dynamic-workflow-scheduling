#include <vector>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/random/gamma_distribution.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/graph/astar_search.hpp>
using namespace boost;

#ifndef INSTANCECONFIG_H_
#define INSTANCECONFIG_H_
//const int depth = 5; //only odd nums
//const int width = 10;//3;
//const float deadline = 1800;
const int types = 4;//types for spotvm and vm
extern int randomsize;//the size of random number generated
extern float OnDemandLag;//0.5;
extern float SpotLag;//1;
extern bool NoSpotVM;
extern float Times[4][types];
extern float lambda;
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
	montage = 0,
	ligo = 1,
	epigenome = 2
};


class taskVertex
{
public:
    int name;
	int type;
    bool mark;
    float start_time;
    float end_time;
    float dl;
    float readyCountdown;
    float* estTime; //expected time, used for deadline assign and configuration
	float* cpuTime;
	float* actTime; //actual time, vary with cloud dynamics
    float restTime;
    Integer_vm status;
    int config;
	float cost; //the monetary cost spent on this task
	float tasktime; //the execution time spent on this task
	float taskstart;
	int assigned_type; 
		
	int* configList;
	int vmID;
	float* prices;
//	int* LV; //load vector
	float read_data; //I/O data in #seeks
	float seq_data; //sequential I/O data in MBytes
	float trans_data; //network data in MBytes
	float rec_data; //network data in MBytes

	float* randomIO; //types*randomsize
	float* seqIO;
	float* netUp;
	float* netDown;
	float* probestTime;
	float* spotTime;
	float* cumulativeTime; //randomsize
	bool tag;

	void instance_config();

	taskVertex() {
		estTime = (float*)malloc(types*sizeof(float));
		cpuTime = (float*)malloc(types*sizeof(float));
		config = 0; mark = false; assigned_type = 0; restTime = rec_data = 0.0; 
		readyCountdown = -1; status = not_ready;
        actTime = (float*)malloc(types*sizeof(float));
		netDown = (float*)malloc(types*randomsize*sizeof(float));
		netUp = (float*)malloc(types*randomsize*sizeof(float));
		probestTime = (float*)malloc(types*randomsize*sizeof(float));
		spotTime = (float*)malloc(randomsize*sizeof(float));
		randomIO = (float*)malloc(types*randomsize*sizeof(float));
		seqIO = (float*)malloc(types*randomsize*sizeof(float));
		cumulativeTime = (float*)malloc(randomsize*sizeof(float));
	}
};

class VM{
public:
	taskVertex* tk;
	int type;
	float life_time;
	float turn_on; //turn on time

	int has_data;//has data of task has_data on local disk
};

class SpotVM{
public:
	taskVertex* tk;
	bool canAlloc;
	int type;
	float price;
	float life_time;//life_time<1.0
	float turn_on; //turn on time

	int has_data;//has data of task has_data on local disk

	SpotVM(float p) { price = p; canAlloc = true;}
};

typedef adjacency_list<boost::vecS, boost::vecS, bidirectionalS, taskVertex, property<edge_weight_t, float> > Graph;
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
	float deadline;
	float meet_dl;
	DAG_type type;
	float arrival_time;
	float* cumulativetime;

	DAG() {};
	//~DAG() {printf("destroyed DAG\n");}
	DAG(float d,float d1){deadline = d; meet_dl = d1;}
	void reset();
	DAG(Graph dag){g = dag;}
	DAG(const DAG& dag) {g=dag.g;deadline=dag.deadline;meet_dl=dag.meet_dl;type=dag.type;arrival_time=dag.arrival_time;}
	std::vector<int> find_CP();
	void initDLAssign();
	void deadline_assign();
};

class configstack{
	//this struct is used to store the searched configuration list in the stack, and maintain their parent and child info
public:
	std::vector<int> configurations;
	//configstack* parent; //parent on the search path
	//configstack* child; //child on the search path
	int taskno; //the task modified to form this configuration 
	//int nextsearchedtype; //the type searched for the next task
	bool* childcolor; //mark which child has been searched
	float fvalue; //used for astar to sort the openset

	configstack(){taskno=0; fvalue=0; childcolor=new bool[types]; childcolor[0]=true; for(int i=1; i<types; i++) childcolor[i]=false;}
	configstack(const configstack& a) {configurations=a.configurations; taskno=a.taskno; fvalue=a.fvalue; childcolor=new bool[types]; for(int i=0; i<types; i++) childcolor[i]=a.childcolor[i];}
	inline bool operator==(configstack* a) {return this->configurations == a->configurations;}
};
//typedef adjacency_list<vecS, vecS, bidirectionalS, DAG, property<edge_weight_t, float> > DAGGraph;//each task is a dag

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
	//void OfflineSP_BFS();//go through all the types of one task first, go through different tasks later
	//void OfflineSP_DFS();//go through each task first, go through different instance types later
	void OfflineAstar();
	void SpotTune();
	void SpotTune2();
	//void OnlineSimulate();
	//if return -1, no suitable bid for this spottype; otherwise returned value*0.001 is the bidding price
	int binary_search(const float arr[], int low, int high, int taskno, int spottype);

//	SearchPrune();
//	SearchPrune(DAG input){dag = input;}
};

bool function(float bid_price, float compare_price);
float rn_01();
int rn_integers(int a, int b);
bool myfunction(taskVertex* a, taskVertex* b);
bool configsortfunction(configstack* a, configstack* b);
float estimateCost(const DAG& dag, int start, bool estimate); //estimate or calculate the total cost of dag, starting from the start task
void estimateTime(DAG& dag, float* estTime); //estimate the total execution time of dag, used in offline search
float estimateTime2(DAG& dag); //estimate the execution time with estTime
void estimateTimeSpot(DAG& dag, float* estTime);//used in spot tune
void estimateTimeSpot2(DAG& dag, float* estTime);//used in spot tune
