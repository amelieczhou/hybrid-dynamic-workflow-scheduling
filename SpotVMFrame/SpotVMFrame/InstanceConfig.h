#include <vector>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>

using namespace boost;

//const int depth = 5; //only odd nums
//const int width = 10;//3;
//const double deadline = 1800;
const int types = 4;//types for spotvm and vm

const double OnDemandLag = 5;//0.5;
const double SpotLag = 10;//1;
//const int ConfigSlot = 3;

enum Integer_vm
{
	not_ready = 0,
	ready = 1,
	scheduled = 2,
	finished = 3
};

class task //with special structure as in SIGMOD11
{
	int name;
	//double time;
	bool mark; //assigned deadline = 1

	std::vector<task*> child_TaskList;
	std::vector<task*> parent_TaskList;
public:
	double start_time;
	double end_time;
	double dl;

	double readyCountdown;
	double estimateTime[types];//depend on Machine Type, including VMs and SpotVMs
	double restTime;//time to complete
	int config;
	int* configList; //2 SpotVMs + 1 on-demand VM
	double* prices; //3D
	Integer_vm status;

	task(double* t, int n);
	void addconfig(){this->config = this->config + 1;}
	void reset(){readyCountdown = -1; status = not_ready;	config = 0; dl = end_time-start_time;}
	//void set_time(double t) {time = t;}
	//double get_time(){return time;}
	void instance_config();
	void set_mark(bool b) {mark = b;}
	bool get_mark() {return mark;}
	void set_child(std::vector<task*> child){child_TaskList = child;}
	int show_name(){return name;}
	void set_name(int n) {name = n;}
	std::vector<task*> get_child(){return child_TaskList;}
	void set_parent(std::vector<task*> parent){parent_TaskList = parent;}
	std::vector<task*> get_parent(){return parent_TaskList;}
	~task()	{child_TaskList.clear(); parent_TaskList.clear();}
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
	double* actTime; //actual time, vary with cloud dynamics
    double restTime;
    Integer_vm status;
    int config;
	double cost; //the monetary cost spent on this task
	double tasktime; //the execution time spent on this task
	double taskstart;
		
	int* configList;
	int vmID;
	double* prices;
//	int* LV; //load vector
	double read_data; //I/O data in #seeks
	double seq_data; //sequential I/O data in MBytes
	double trans_data; //network data in MBytes
	double rec_data; //network data in MBytes

	void instance_config();
};

class VM{
public:
	taskVertex* tk;
	int type;
	int life_time;
	double turn_on; //turn on time
};

class SpotVM{
public:
	taskVertex* tk;
	bool canAlloc;
	int type;
	double price;
	double life_time;//life_time<1.0
	double turn_on; //turn on time

	SpotVM(double p) { price = p; canAlloc = true;}
};

typedef adjacency_list<vecS, vecS, bidirectionalS, taskVertex, property<edge_weight_t, double> > Graph;
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

	DAG(double d){deadline = d;}
	void reset();
	DAG(Graph dag){g = dag;}
	std::vector<int> find_CP();
	void initDLAssign();
	void deadline_assign();
};


bool function(double bid_price, double compare_price);
double rn_01();
int rn_integers(int a, int b);
bool myfunction(taskVertex* a, taskVertex* b);
