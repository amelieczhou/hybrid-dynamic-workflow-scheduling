#ifndef INSTANCECONFIG_H
#define INSTANCECONFIG_H

#include <vector>
#include <queue>
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
//extern int num_monte;
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
	epigenome = 2,
	pipeline = 3
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

	//for autoscaling deadline assign
	float sub_deadline;
	float EST;
	float LFT;
	int prefer_type;
		
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
	float* cumulativeTime; //randomsize
	bool tag;

	float* randspot; //randomsize*spotrandsize=100, probabilistic distribution

	void instance_config();

	taskVertex() {
		printf("constructor of task\n");
		estTime = (float*)malloc(types*sizeof(float));
		cpuTime = (float*)malloc(types*sizeof(float));
		config = 0; mark = false; assigned_type = 0; restTime = rec_data = 0.0;  prefer_type = 0.0; 
		readyCountdown = -1; status = not_ready; EST = LFT=0; sub_deadline=0;
        actTime = (float*)malloc(types*sizeof(float));
		netDown = (float*)malloc(types*randomsize*sizeof(float));
		netUp = (float*)malloc(types*randomsize*sizeof(float));
		probestTime = (float*)malloc(types*randomsize*sizeof(float));
		randomIO = (float*)malloc(types*randomsize*sizeof(float));
		seqIO = (float*)malloc(types*randomsize*sizeof(float));
		if(randomsize>4000)
			cumulativeTime = (float*)malloc(randomsize*sizeof(float));
		else
		cumulativeTime = (float*)malloc(4000*sizeof(float));
		randspot = (float*)malloc(400*sizeof(float));//spotdistrsize
	}
	/*taskVertex(const taskVertex& task){
		printf("copy constructor of taskVertex\n");
		config = task.config; mark = task.mark; assigned_type = task.assigned_type; prefer_type = task.prefer_type;
		restTime = task.restTime; rec_data = task.rec_data;  name = task.name; vmID = task.vmID;
		readyCountdown = task.readyCountdown; status = task.status; EST = task.EST; LFT=task.LFT; sub_deadline=task.sub_deadline;
		estTime = task.estTime;
		cpuTime = task.cpuTime;
		actTime = task.actTime;
		netDown = task.netDown;
		netUp = task.netUp;
		probestTime = task.probestTime;
		randomIO = task.randomIO;
		seqIO = task.seqIO;
		cumulativeTime = task.cumulativeTime;
		randspot = task.randspot;
	}*/
	/*taskVertex& operator=(const taskVertex& task){
		printf("= operator of task\n");
		config = task.config; mark = task.mark; assigned_type = task.assigned_type; prefer_type = task.prefer_type;
		restTime = task.restTime; rec_data = task.rec_data; name = task.name; vmID = task.vmID;
		readyCountdown = task.readyCountdown; status = task.status; EST = task.EST; LFT=task.LFT; sub_deadline=task.sub_deadline;
		estTime = task.estTime;
		cpuTime = task.cpuTime;
		actTime = task.actTime;
		netDown = task.netDown;
		netUp = task.netUp;
		probestTime = task.probestTime;
		randomIO = task.randomIO;
		seqIO = task.seqIO;
		cumulativeTime = task.cumulativeTime;
		randspot = task.randspot;
		return *this;
	}*/
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
typedef graph_traits<Graph>::edge_iterator edge_iter;
typedef graph_traits<Graph>::out_edge_iterator out_edge_iterator;
typedef graph_traits<Graph>::in_edge_iterator in_edge_iterator;
typedef graph_traits<Graph>::edge_descriptor edge_descriptor;
typedef graph_traits<Graph>::adjacency_iterator adja_iterator;

class DAG
{
public:
	Graph* g;
	float deadline;
	float meet_dl;
	DAG_type type;
	float arrival_time;
	float* cumulativetime;

	DAG() {};
	//~DAG() {printf("destroyed DAG\n");}
	DAG(float d,float d1){deadline = d; meet_dl = d1;}
	void reset();
	DAG(Graph* dag){g = dag;}
	DAG(const DAG& dag) {
		g=new Graph();//dont know the constructor of Graph, so it may cause memory leak!!!!
		Graph* graph = dag.g;
		//g->m_edges = graph->m_edges;
		//g->m_property = graph->m_property;
		std::pair<vertex_iter, vertex_iter> vp = vertices(*graph);
		for(; vp.first != vp.second; vp.first++){			
			taskVertex v = (*graph)[*vp.first];
			add_vertex(v,*g);
			in_edge_iterator in_i, in_end;
			edge_descriptor e;
			for (boost::tie(in_i, in_end) = in_edges(*vp.first, *graph); in_i != in_end; ++in_i){
				e = *in_i;
				Vertex src = source(e, *graph);
				add_edge(src,*vp.first,*g);
			}
		}

		deadline=dag.deadline;
		meet_dl=dag.meet_dl;
		type=dag.type;
		arrival_time=dag.arrival_time;
	}
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
	void SpotTune2(int i);
	//void OnlineSimulate();
	//if return -1, no suitable bid for this spottype; otherwise returned value*0.001 is the bidding price
	int binary_search(const float* arr, int low, int high, int taskno, int spottype);
	int binary_search_heuristic(const float* arr, int low, int high, int taskno, int spottype);

//	SearchPrune();
//	SearchPrune(DAG input){dag = input;}
};

class individual{//for MOEA; solution, fitness, evaluation of cost and time objectives
public:
	int index; //ID
	int* solution;
	float fitness;
	std::pair<float,float> objectives;

	individual(){fitness = 0.0; objectives.first =0.0; objectives.second = 0.0;}
	~individual() {if(solution) delete[] solution;}	
};

bool function(float bid_price, float compare_price);
float rn_01();
int rn_integers(int a, int b);
bool myfunction(taskVertex* a, taskVertex* b);
bool configsortfunction(configstack* a, configstack* b);
float estimateCost(const DAG& dag, int start, int end, bool estimate); //estimate or calculate the total cost of dag, starting from the start task to end, inclusive
void estimateTime(DAG& dag, int start, int end, float* estTime); //estimate the total execution time of dag, used in offline search
void estimateTimeSpot2(DAG& dag, float* estTime);//used in spot tune
void conv(float* array1, float* array2, float* result, int length1, int length2);
void calmaxdistr(float* array1, float* array2, float* result, int length1, int length2);
//FOR DEADLINE ASSIGN
//double MET(taskVertex& tk);
double MTT(taskVertex* tk1, taskVertex* tk2);
double EST(taskVertex& tk, DAG job);
double LFT(taskVertex& tk, DAG job);
void AssignParents(taskVertex* tk, DAG* job);
Vertex CriticalParent(taskVertex* tk, DAG* job);
void AssignPath(std::deque<Vertex> PCP,DAG* job);
bool has_unassigned_parent(taskVertex* tk, DAG* job);
bool has_unassigned_child(taskVertex* tk, DAG* job);
void update_EST(taskVertex* tk, DAG* job);
void update_LFT(taskVertex* tk, DAG* job);

bool dominates(individual*,individual*);
float calcDistance(individual*,individual*);

#endif