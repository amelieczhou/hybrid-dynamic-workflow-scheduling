#include "stdafx.h"
#include "Algorithms.h"
#include "PricingModel.h"
#include "mutation.h"
#include <ctime>
#include <fstream>
#include <sstream>
#include <string>
#include <omp.h>
#include <assert.h>

//selection
#define PISA_MAXFLOAT 1E99 //define internal maximum double number
extern double budget;//middle for epi
const int populationsize = 10;
const int archivesize = 10;
const int maxiteration = 100;
const int tournament = 2;//binary
const int mu = 10; //number of individuals selected as parents
int** NN;
//int* copies;
std::vector<int> copies;
int* fitness_bucket;
int* fitness_bucket_mod;
//int* old_index;
std::vector<int> old_index;
float** distances;
//variation
int *current_population; // storing the IDs of the individuals


std::vector<individual* > nondominatedsolutions;
std::vector<individual* > generation; //initial size = populationsize 
std::vector<individual* > variation; //initially empty
std::vector<individual* > matingpool;


void MOEA::evalObjectives(individual* individ){
	//the first is cost fitness and the second is time fitness
	std::pair<vertex_iter,vertex_iter> vp = vertices(*dag->g);
	for(;vp.first!=vp.second; vp.first++){
		(*dag->g)[*vp.first].assigned_type = individ->solution[*vp.first];
	}
	int numoftasks = num_vertices(*dag->g);
	float costobjective = estimateCost(*dag,0,numoftasks-1,false);
	float* exeTime = (float*)malloc(randomsize*sizeof(float));
	float timeobjective = 0.0; 
	estimateTime(*dag,0,numoftasks-1,exeTime);
	for(int i=0; i<randomsize; i++)
		timeobjective += exeTime[i];
	timeobjective /= randomsize;
	costobjective /= budget;
	timeobjective /= dag->deadline;

	float penalty = 0.0;
	if(costobjective > 1.0)
		 penalty += costobjective;
	if(timeobjective > 1.0)
		penalty += timeobjective;

	costobjective += penalty;
	timeobjective += penalty;
	individ->objectives.first = costobjective;
	individ->objectives.second = timeobjective;
	return;
}

void MOEA::calFitnesses(std::vector<individual*> populationall){
	
	int popusizeall = populationall.size();
	int* strengths = (int*)malloc(popusizeall*sizeof(int));	
	//initialize fitness and strength values
	for(int i=0; i<popusizeall; i++){
		evalObjectives(populationall[i]);
		strengths[i] = 0;		
	}
	
	for(int i=0; i<popusizeall; i++){
		for(int j=0; j<popusizeall; j++){
			if(dominates(populationall[i],populationall[j])){
				strengths[i] ++;
			}
		}
	}
	for(int i=0; i<popusizeall; i++){
		int sum = 0;
		for(int j=0; j<popusizeall; j++){
			if(dominates(populationall[j],populationall[i]))
				sum +=strengths[j];
		}
		populationall[i]->fitness = sum;
		fitness_bucket[sum] ++;
		fitness_bucket_mod[(sum / popusizeall)] ++;
	}	
	free(strengths);	
	return ;
}
void calcDistances(std::vector<individual*> populationall){
	int popusize = populationall.size();
	//initialize 
	for(int i=0; i<popusize; i++) {
		copies[i] = 1;
		for(int j=0; j<popusize; j++)
			NN[i][j] = -1;
	}
	//calculate distances
	for(int i=0; i<popusize; i++){
		NN[i][0] = i;
		for(int j=i+1; j<popusize; j++){
			distances[i][j] = calcDistance(populationall[i],populationall[j]);
			distances[j][i] = distances[i][j];
			if(distances[i][j] == 0){
				NN[i][copies[i]] = j;
				NN[j][copies[j]] = i;
				copies[i]++;
				copies[j]++;
			}
		}
		distances[i][i] = 0.0;
	}
}
int getNN(int index, int k, std::vector<individual*> population)
/* lazy evaluation of the k-th nearest neighbor
   pre-condition: (k-1)-th nearest neigbor is known already */
{
    assert(index >= 0);
    assert(k >= 0);
    assert(copies[index] > 0);
    
    if (NN[index][k] < 0)
    {
	int i;
	float min_dist = PISA_MAXFLOAT;
	int min_index = -1;
	int prev_min_index = NN[index][k-1];
	double prev_min_dist = distances[index][prev_min_index];
	assert(prev_min_dist >= 0);
	
	for (i = 0; i < population.size(); i++)
	{
	    double my_dist = distances[index][i];
	    
	    if (my_dist < min_dist && index != i)
	    {
            if (my_dist > prev_min_dist || (my_dist == prev_min_dist && i > prev_min_index))
            {
    		    min_dist = my_dist;
    		    min_index = i;
    		}
	    }
	}
	
	NN[index][k] = min_index;
    }

    return (NN[index][k]);
}
double getNNd(int index, int k,std::vector<individual*> population)
/* Returns the distance to the k-th nearest neigbor
   if this individual is still in the population.
   For for already deleted individuals, returns -1 */
{
    int neighbor_index = getNN(index, k, population);
    
    if (copies[neighbor_index] == 0)
		return (-1);
    else
	return (distances[index][neighbor_index]);
}

//truncate from nondominated individuals (if too many)
void truncate_nondominated(std::vector<individual*>& populationall){
	int size = populationall.size();
	//delete all dominated individuals
	for(int i=0; i<size; i++){
		if(populationall[i]->fitness > 0){//dominated
			//free(populationall[i]->solution);
			populationall[i] = NULL;
			//populationall.erase(populationall.begin()+i); 
			copies[i] = 0;
		}
	}

	//truncate from non-dominated individuals
	while(fitness_bucket[0] > populationsize){
		int *marked;
        int max_copies = 0;
        int count = 0;
        int delete_index;
    
    	marked = (int*) malloc(populationall.size() * sizeof(int));
    
    	/* compute inds with maximal copies */
    	for (int i = 0; i < populationall.size(); i++) {
    	    if (copies[i] > max_copies) {
        		count = 0;
            	max_copies = copies[i];
    	    }
    	    if (copies[i] == max_copies)  {
                marked[count] = i;
                count++;
    	    }
    	}
	
        assert(count >= max_copies);
    
    	if (count > max_copies)	{    
    	    int *neighbor;
    	    neighbor = (int*) malloc(count * sizeof(int));
    	    for (int i = 0; i < count; i++)
        		neighbor[i] = 1;  /* pointers to next neighbor */
	    
    	    while (count > max_copies) {
        		float min_dist = PISA_MAXFLOAT;
        		int count2 = 0;
		
        		for (int i = 0; i < count; i++){
        		    double my_dist = -1;
        		    while (my_dist == -1 && neighbor[i] < populationall.size()) {
            			my_dist = getNNd(marked[i],neighbor[i],populationall);
            			neighbor[i]++;
        		    }
        		    if (my_dist < min_dist) {
            			count2 = 0;
            			min_dist = my_dist;
        		    }
        		    if (my_dist == min_dist) {
            			marked[count2] = marked[i];
            			neighbor[count2] = neighbor[i];
            			count2++;
        		    }
        		}
                count = count2;
                if (min_dist == -1) { /* all have equal distances */                
                    break;
                }
            }   
            free(neighbor);
        }
	
        /* remove individual from population */
		int index = rn_integers(0,count-1);
		while(index>=count || index<0)
			index = rn_integers(0,count-1);
		delete_index = marked[index];

		//free(populationall[delete_index]->solution);
		populationall[delete_index] = NULL;
		//populationall.erase(populationall.begin()+delete_index);
        for (int i = 0; i < count; i++) {
            if (distances[delete_index][marked[i]] == 0){
				copies[marked[i]]--;
				if(copies[marked[i]]==0)
					printf("");
			}
		}
		copies[delete_index] = 0; /* Indicates that this index is empty */
		fitness_bucket[0]--;
		fitness_bucket_mod[0]--;
		free(marked);
	}
	return;
}
void truncate_dominated(std::vector<individual*>& populationall){
	/* truncate from dominated individuals */
	int i, j;
    int size;
    int num = 0;
    
    size = populationall.size();
    
    i = -1;
	while (num < populationsize){
		i++;
		num += fitness_bucket_mod[i];
    }
    
    j = i * size;
    num = num - fitness_bucket_mod[i] + fitness_bucket[j];
    while (num < populationsize){
		j++;
		num += fitness_bucket[j];
    }
    
    if (num == populationsize){
		for (i = 0; i < size; i++){
			if (populationall[i]->fitness > j){
				//free(populationall[i]->solution);
				populationall[i] = NULL;
				//populationall.erase(populationall.begin()+i);
			}
		}
    }
    else {/* if not all fit into the next generation */
		int k;
		int free_spaces;
		int fill_level = 0;
		int *best;

		free_spaces = populationsize - (num - fitness_bucket[j]);
		best = (int*) malloc(free_spaces * sizeof(int));
		for (i = 0; i < size; i++){
			if (populationall[i]->fitness > j){
				//free(populationall[i]->solution);
				populationall[i] = NULL;
				//populationall.erase(populationall.begin()+i);
			}
			else if (populationall[i]->fitness == j){
				if (fill_level < free_spaces){
					best[fill_level] = i;
					fill_level++;
					for (k = fill_level - 1; k > 0; k--){
					int temp;
					if (getNNd(best[k], 1,populationall) <= getNNd(best[k - 1], 1,populationall)){
						break;
					}
					temp = best[k];
					best[k] = best[k-1];
					best[k-1] = temp;
					}
				}
				else{
					if (getNNd(i, 1,populationall) <= getNNd(best[free_spaces - 1], 1,populationall)){
						//free(populationall[i]->solution);
						populationall[i] = NULL;
						//populationall.erase(populationall.begin()+i);
					}
					else{
						//free(populationall[best[free_spaces - 1]]->solution);
						populationall[best[free_spaces - 1]] = NULL;
						//populationall.erase(populationall.begin()+best[free_spaces - 1]);
						best[free_spaces - 1] = i;
						for (k = fill_level - 1; k > 0; k--){
							int temp;
							if (getNNd(best[k], 1,populationall) <= getNNd(best[k - 1], 1,populationall)){
								break;
							}
							temp = best[k];
							best[k] = best[k-1];
							best[k-1] = temp;
						}
					}
				}
			}
		}
    }

	return;
}
void MOEA::Crossover(){
	return;
}
void MOEA::Mutation(std::vector<individual* >& generationall){
	/*
	int i;
    int result; // stores return values of called functions
    int min_valid, max_valid;
     
	 current_population = (int *) malloc(archivesize * sizeof(int)); 
     if (current_population == NULL){
          exit(1);
     }
	 // Removes the individuals from generation, that are not in the archive populationall
	 result = update_from_arc(populationall); 
	 result = variate(matingpool);*/
	//update_from_arc(generationall); 
	variate(matingpool);        
    // do evaluation
    for(int i = 0; i < mu; i++)   {
		evalObjectives(matingpool[i]);
    }
	return;
}

void MOEA::Evolution(){	
	//****************************************evolution process********************************************//
	//nondominated set
	int numtasks = num_vertices(*dag->g);	
	//Generate initial population randomly
	for(int i=0; i<populationsize; i++){
		individual* ind = new individual();
		ind->index = i;
		ind->solution = (int*)malloc(numtasks*sizeof(int));
		printf("solution %d:\n",i);
		for(int j=0; j<numtasks; j++){
			//randomly generate instance configuration for each task
			int config = rn_integers(0,types);
			while(config==types) config=rn_integers(0,types);
			int c=0; do{c++;}while(c<1E8);
			ind->solution[j] = config;
			printf("%d, ",config);
		}
		printf("\n");
		float fitness = 0.0;//Fitness_evaluate(solution);
		generation.push_back(ind);
	}
	int iter_counter = 0;
	while(true){
	
		//Fitness assignment		
		std::vector<individual* > generationall;
		for(int i=0; i<generation.size(); i++)
			generationall.push_back(generation[i]);
		for(int i=0; i<variation.size(); i++)
			generationall.push_back(variation[i]);
		//first create internal structures
		fitness_bucket = (int*)malloc(generationall.size()*generationall.size()*sizeof(int));
		fitness_bucket_mod = (int*)malloc(generationall.size()*sizeof(int));
		//copies = (int*)malloc(generationall.size()*sizeof(int));
		copies = std::vector<int>(generationall.size());
		old_index = std::vector<int>(generationall.size());
		//old_index = (int*)malloc(generationall.size()*sizeof(int));
		
		distances = (float**)malloc(generationall.size()*sizeof(float*));
		NN = (int**)malloc(generationall.size()*sizeof(int*));
		for(int i=0; i<generationall.size(); i++){
			NN[i] = (int*)malloc(generationall.size()*sizeof(int));
			distances[i] = (float*)malloc(generationall.size()*sizeof(float));
		}

		for(int i=0; i<generationall.size(); i++){
			for(int j=0; j<generationall.size(); j++){
				fitness_bucket[i*generationall.size()+j] = 0;
				distances[i][j] = 0.0;
				NN[i][j] = 0;
			}
		}
		for(int i=0; i<generationall.size(); i++){
			fitness_bucket_mod[i] = 0;
			copies[i] = 0;
			old_index[i] = 0;
		}
		calFitnesses(generationall);
		calcDistances(generationall);

		//environmental selection
		if(fitness_bucket[0] > populationsize){
			truncate_nondominated(generationall);
		}else if(generationall.size() > populationsize)
			truncate_dominated(generationall);
		//move remaining individuals to top of array in populationall
		
		for (int i = 0; i < generationall.size(); i++){
			int nullcount = 0;
			for(int j=0; j<i; j++){
				individual* temp_ind = generationall[j];
				if(temp_ind == NULL){
					nullcount++;
					//copies[j]=0;
				}
			}
			old_index[i-nullcount] = i;			
		}
		for (int i = 0; i < generationall.size(); i++){
			individual* temp_ind = generationall[i];
			if (temp_ind == NULL){
				//assert(copies[old_index[i]+i] > 0);    
    			generationall.erase(generationall.begin()+i);
				//copies[old_index[i]+i] = 0;
				i--;
			}
		}
		
		assert(generationall.size() <= populationsize);///////////////////////////////
		variation.clear();
		for(int i=0; i<generationall.size(); i++)
				variation.push_back(generationall[i]);
		//********to this point, variation should have the size of populationsize*******//
		/*******Termination*********/
		if(iter_counter>= maxiteration){
			for(int i=0; i<generationall.size(); i++)
				nondominatedsolutions.push_back(generationall[i]);
			free(fitness_bucket);
			free(fitness_bucket_mod);
			//free(copies);
			//free(old_index);
			for(int i=0; i<generationall.size(); i++){
				free(distances[i]);
				free(NN[i]);
			}
			free(distances);
			free(NN);
			break;
		}
		//mating selection
		//fills mating pool / offspring population variation
		matingpool.clear();
		for (int i = 0; i < mu; i++){
			int winner = rn_integers(0,generationall.size());
			while(winner == generationall.size()) winner = rn_integers(0,generationall.size());
	
			int timecounter=0; do{timecounter++;}while(timecounter<1e8);
			for (int j = 1; j < tournament; j++){
				int opponent = rn_integers(0,generationall.size());
				while(opponent == generationall.size()) opponent = rn_integers(0,generationall.size());
				int timecounter1=0; do{timecounter1++;}while(timecounter1<1e8);
				if (generationall[opponent]->fitness < generationall[winner]->fitness || winner == opponent)
				{
					winner = opponent;
				}
				else if (generationall[opponent]->fitness == generationall[winner]->fitness){
					if (distances[old_index[opponent]][getNN(old_index[opponent], 1,generationall)] >
						distances[old_index[winner]][getNN(old_index[winner], 1,generationall)]){						
							winner = opponent;
					}
				}
			}  
			matingpool.push_back(generationall[winner]);
		}
		//Variation: recombination and mutation
		Mutation(generationall);
		generation.clear();
		for(int i=0; i<matingpool.size(); i++)
			generation.push_back(matingpool[i]);
		

		free(fitness_bucket);
		free(fitness_bucket_mod);
		//free(copies);
		//free(old_index);
		for(int i=0; i<generationall.size(); i++){
			free(distances[i]);
			free(NN[i]);
		}
		free(distances);
		free(NN);	

		/*for(int i=0; i<generationall.size(); i++){
			free(generationall[i]->solution);
			generationall[i]=NULL;
		}*/
		generationall.clear();
		//std::vector<individual*>().swap(generationall);
		iter_counter ++;
	}
	
}
void MOEA::Simulate(){
	float arrival_time = 0;
	dag->arrival_time = 0;	
	float violation = 0;
	float ave_cost = 0;
	float ioseq[types],iorand[types],net_up[types],net_down[types];
	omp_set_num_threads(24);
	float viol_private[24];
	float cost_private[24];
	for(int i=0; i<24; i++) viol_private[i]=cost_private[i]=0;

	//for the performance of each instance type
	float* random_sequential_io = (float*)malloc(types*randomsize*sizeof(float));
	float* random_random_io = (float*)malloc(types*randomsize*sizeof(float));
	float* random_network_up = (float*)malloc(types*randomsize*sizeof(float));
	float* random_network_down = (float*)malloc(types*randomsize*sizeof(float));
	float* random_tmp = (float*)malloc(types*10000*sizeof(float));
	//read from file	
	FILE* rFile;
	char str[1024];
	char buf[256];
	char *ptr, *ptr2;
	rFile = fopen("randio.csv","r");
	if(rFile == NULL){
		printf("cannot open randio.csv\n");
		exit(1);
	}
	for(int i=0; i<types*10000; i++){
		if(fgets(str,1024,rFile)!=NULL)
			random_tmp[i] = atof(str);
	}	
	for(int i=0; i<types; i++){
		for(int j=0; j<randomsize; j++){
			random_random_io[i*randomsize+j] = random_tmp[i*10000+j];//10000 is fixed
		}
	}

	rFile = fopen("seqio.csv","r");
	if(rFile == NULL){
		printf("cannot open seqio.csv\n");
		exit(1);
	}
	for(int i=0; i<types*10000; i++){
		if(fgets(str,1024,rFile)!=NULL)
			random_tmp[i] = atof(str);
	}	
	for(int i=0; i<types; i++){
		for(int j=0; j<randomsize; j++){
			random_sequential_io[i*randomsize+j] = random_tmp[i*10000+j];//10000 is fixed
		}
	}
	rFile = fopen("netup.csv","r");
	if(rFile == NULL){
		printf("cannot open netup.csv\n");
		exit(1);
	}
	for(int i=0; i<types*10000; i++){
		if(fgets(str,1024,rFile)!=NULL)
			random_tmp[i] = atof(str);
	}	
	for(int i=0; i<types; i++){
		for(int j=0; j<randomsize; j++){
			random_network_up[i*randomsize+j] = random_tmp[i*10000+j];//10000 is fixed
		}
	}
	rFile = fopen("netdown.csv","r");
	if(rFile == NULL){
		printf("cannot open netdown.csv\n");
		exit(1);
	}
	for(int i=0; i<types*10000; i++){
		if(fgets(str,1024,rFile)!=NULL)
			random_tmp[i] = atof(str);
	}	
	for(int i=0; i<types; i++){
		for(int j=0; j<randomsize; j++){
			random_network_down[i*randomsize+j] = random_tmp[i*10000+j];//10000 is fixed
		}
	}
	free(random_tmp);
	std::pair<vertex_iter, vertex_iter> vp;
	vp = vertices((*dag->g));

	int quantile = dag->meet_dl * randomsize;
	for(; vp.first != vp.second; vp.first++){
		//dag->g[*vp.first].probestTime = new float[types][randomsize];
		for(int t=0; t<types; t++){
			for(int j=0; j<randomsize; j++){
				(*dag->g)[*vp.first].netUp[t*randomsize+j] = (*dag->g)[*vp.first].trans_data * random_network_up[t*randomsize+j] / 8000;
				(*dag->g)[*vp.first].netDown[t*randomsize+j] = (*dag->g)[*vp.first].rec_data * random_network_down[t*randomsize+j] / 8000;
				(*dag->g)[*vp.first].randomIO[t*randomsize+j] = (*dag->g)[*vp.first].read_data / random_random_io[t*randomsize+j];
				(*dag->g)[*vp.first].seqIO[t*randomsize+j] = (*dag->g)[*vp.first].seq_data / random_sequential_io[t*randomsize+j];
				(*dag->g)[*vp.first].probestTime[t*randomsize+j] = (*dag->g)[*vp.first].cpuTime[t] + (*dag->g)[*vp.first].netUp[t*randomsize+j]
					+ (*dag->g)[*vp.first].netDown[t*randomsize+j] + (*dag->g)[*vp.first].randomIO[t*randomsize+j] + (*dag->g)[*vp.first].seqIO[t*randomsize+j];
			}
			//calculate the estimate time as the expected value of the proestTime
			std::sort((*dag->g)[*vp.first].probestTime+t*randomsize,(*dag->g)[*vp.first].probestTime+(t+1)*randomsize-1);
			(*dag->g)[*vp.first].estTime[t] = (*dag->g)[*vp.first].probestTime[t*randomsize+quantile];
			printf("task: %d, type: %d, time: %f\n",*vp.first,t,(*dag->g)[*vp.first].estTime[t]);
		}
	}

	free(random_sequential_io);
	free(random_random_io);
	free(random_network_up);
	free(random_network_down);

	//start the evolution process to obtain solutions
	time_t start,end;
	time(&start);	
	Evolution();	
	time(&end);
	printf("optimization overhead of spea2 is %.4f\n",difftime(end,start));

	//use the nondominated solution to simulate
	int nondominateindex = -1;
	bool debug = true;
	if(!debug){	
	for(int i=0; i<nondominatedsolutions.size(); i++)	{
		if(nondominatedsolutions[i]->fitness == 0)	{
			nondominateindex = i;
			break;
		}
	}
	if(nondominateindex == -1){
		printf("no nondominated solution found, error!");
		exit(1);
	}
	}else{
		nondominateindex = 0;
		int msize = (*dag->g).m_vertices.size();
		individual* ind = new individual();		
		ind->solution = (int*)malloc(msize*sizeof(int));
		for(int i=0; i<msize; i++)
			ind->solution[i] = 1;
		nondominatedsolutions.push_back(ind);
	}

	for(int i=0; i<nondominatedsolutions.size(); i++){
            printf("the %d-th solution found:\t",i);
            for(int j=0; j<(*dag->g).m_vertices.size(); j++){
                    printf("%d,",nondominatedsolutions[i]->solution[j]);
            }
            printf(", fitness: %4.2f, objectives: %4.4f, %4.4f\n",nondominatedsolutions[i]->fitness,nondominatedsolutions[i]->objectives.first,nondominatedsolutions[i]->objectives.second);
    }
	return;

	vp = vertices(*dag->g);
	for(;vp.first != vp.second; vp.first++){
		(*dag->g)[*vp.first].prefer_type = nondominatedsolutions[nondominateindex]->solution[*vp.first]; //try different solutions in the nondominated set?
		printf("task %d: %d\n",*vp.first,(*dag->g)[*vp.first ].prefer_type);
	}
	std::vector<DAG*> workflows; //continuous workflow
	workflows.push_back(dag);

	std::ifstream infile;
	std::string a = "arrivaltime_integer_";
	std::string b;
	std::ostringstream strlamda;
	strlamda << lambda;
	b = strlamda.str();
	std::string c = ".txt";
	std::string fname = a + b + c;
	char time[256];
	infile.open(fname.c_str());
	if(infile==NULL){
		printf("cannot find input file!\n");
		return;
	}
	infile.getline(time,256); //jump the lamda line
	infile.getline(time,256); //jump the 0 line
	//incomming jobs
	//while(arrival_time < max_t){
	while(workflows.size()<(int)num_jobs){
		infile.getline(time,256);
		arrival_time = atof(time);

		DAG* job = new DAG(dag->deadline+arrival_time,dag->meet_dl);		
		job->g = dag->g; job->type = dag->type;
		job->arrival_time = arrival_time;
		vp = vertices(*job->g);
		for(int i=0; i<(*vp.second - *vp.first); i++)
			(*job->g)[i].sub_deadline += arrival_time;
		workflows.push_back(job);
	}
	infile.close();
	//start simulation
	std::clock_t starttime = std::clock();

	#pragma omp parallel
	{
		#pragma omp for
		for(int monte=0; monte < randomsize; monte++)
		{
			/*for(int i=0; i<types; i++){
				ioseq[i] = random_sequential_io[randomsize*i+monte];
				iorand[i] = random_random_io[randomsize*i+monte];
				net_up[i] = random_network_up[randomsize*i+monte];
				net_down[i] = random_network_down[randomsize*i+monte];
			}*/		
			std::vector<DAG*> jobs;
			for(int i=0; i<workflows.size(); i++){
				DAG* newdag = new DAG(*workflows[i]);
				vp = vertices(*newdag->g);
				for(int j=0; j<(*vp.second - *vp.first); j++)
					(*newdag->g)[j].sub_deadline = (*workflows[i]->g)[j].sub_deadline;
				jobs.push_back(newdag);
			}
			std::vector<VM*> VMTP[types];
			int need_VM[types]={0,0,0,0};

			//EDF scheduling
			double t = 0;
			bool condition = false;
			double moneycost = 0.0;
		
			do{	
				//accept workflows
				for(int i=0; i<jobs.size(); i++){
					if((int)t == (int)jobs[i]->arrival_time){
						if(dag->type == montage){
							for(int j=0; j<4; j++) {
								(*jobs[i]->g)[j].status = ready;
								(*jobs[i]->g)[j].readyCountdown = -1;
								(*jobs[i]->g)[j].restTime = 0;				
							}
							for(int j=4; j<20; j++){
								(*jobs[i]->g)[j].status = not_ready;
								(*jobs[i]->g)[j].readyCountdown = -1;
								(*jobs[i]->g)[j].restTime = 0;				
							}
						}else if(dag->type == ligo){
							for(int j=0; j<9; j++) {
								(*jobs[i]->g)[j].status = ready;
								(*jobs[i]->g)[j].readyCountdown = -1;
								(*jobs[i]->g)[j].restTime = 0;				
							}
							for(int j=9; j<40; j++){
								(*jobs[i]->g)[j].status = not_ready;
								(*jobs[i]->g)[j].readyCountdown = -1;
								(*jobs[i]->g)[j].restTime = 0;	
							}
						}else if(dag->type == epigenome){				
							(*jobs[i]->g)[0].status = ready;
							(*jobs[i]->g)[0].readyCountdown = -1;
							(*jobs[i]->g)[0].restTime = 0;	
							for(int j=1; j<20; j++){
								(*jobs[i]->g)[j].status = not_ready;
								(*jobs[i]->g)[j].readyCountdown = -1;
								(*jobs[i]->g)[j].restTime = 0;	
							}
						}else{
							printf("what is the dag type?");
							exit(1);
						}
						//jobs.push_back(newdag);
						//printf("add new dag\n");
					}
				}
				//step 1
				std::vector<taskVertex*> ready_task;
				for(int ji=0; ji<jobs.size(); ji++){
					vp = vertices(*jobs[ji]->g);
					for(int i=0; i < (*vp.second - *vp.first ); i++)
					{
						bool tag = true;
						//get parent vertices
						in_edge_iterator in_i, in_end;
						edge_descriptor e;
						boost::tie(in_i, in_end) = in_edges(i, *jobs[ji]->g);
						if(in_i == in_end) tag = false;
						else{
							for (; in_i != in_end; ++in_i) 
							{
								e = *in_i;
								Vertex src = source(e, *jobs[ji]->g);					
								if((*jobs[ji]->g)[src].status != finished)
								{
									tag = false;
									//break;
								}
							}
						}
						if((*jobs[ji]->g)[i].status == ready || (tag && (*jobs[ji]->g)[i].status != scheduled && (*jobs[ji]->g)[i].status != finished)){
							ready_task.push_back(&(*jobs[ji]->g)[i]);							
						}
					}
				
				}


				std::sort(ready_task.begin(),ready_task.end(), myfunction);
				for(int i=0; i<ready_task.size(); i++)//earliest deadline first
				{
					taskVertex* curr_task=ready_task[i];
					if(curr_task->readyCountdown == -1)//
					{
						
						int _config = curr_task->prefer_type;
						bool find = false;
						//check VM/SpotVM list for available machine
						int size = VMTP[_config].size();
						for(int j=0; j<size; j++)
						{
							if(VMTP[_config][j]->tk == NULL)
							{
								find = true;
								VMTP[_config][j]->tk = curr_task;
								break;
							}
						}
						if(find) {
							curr_task->status = scheduled;
							curr_task->tasktime = t;
							curr_task->restTime =  curr_task->probestTime[curr_task->prefer_type*randomsize+monte] ;
						}
						else 			
						{
							curr_task->readyCountdown = OnDemandLag;
							curr_task->tasktime = t;
						}
					}
					else if(curr_task->readyCountdown == 0)
					{
						curr_task->status = scheduled;
						curr_task->restTime = curr_task->probestTime[curr_task->prefer_type*randomsize+monte] ;

						VM* vm = new VM; 
						vm->life_time = OnDemandLag;
						vm->tk = curr_task;
						vm->type = curr_task->prefer_type;
						VMTP[curr_task->prefer_type].push_back(vm);
						
					}			
				}
				//delete VMs without task
				for(int i=0; i<types; i++)//////////
				{
					int size1 = VMTP[i].size();
					
					for(int j=0; j<size1; j++)
					{
						if(VMTP[i][j]->tk == NULL)
						{
							double runtime = VMTP[i][j]->life_time;
							moneycost += (priceOnDemand[i] * ceil(runtime/3600.0));

							VM* vm = VMTP[i][j];
							delete vm;
							VMTP[i].erase(VMTP[i].begin()+j);
							j--;
							size1--;
						}
					}
				}
				//step 2
				std::vector<taskVertex*> scheduled_task;
				for(int ji=0; ji<jobs.size(); ji++){
					vp = vertices((*jobs[ji]->g));
					for(int i=0; i<(*vp.second - *vp.first ); i++)
						if((*jobs[ji]->g)[i].status == scheduled)
							scheduled_task.push_back(&(*jobs[ji]->g)[i]);
				}
				
				for(int i=0; i<scheduled_task.size(); i++)
				{
					scheduled_task[i]->restTime -= 1;////////////////////////////
					if(scheduled_task[i]->restTime <= 0) 
					{
						scheduled_task[i]->status = finished;
						scheduled_task[i]->end_time = t;
						scheduled_task[i]->tasktime = t - scheduled_task[i]->tasktime;
						scheduled_task[i]->cost = scheduled_task[i]->tasktime * priceOnDemand[scheduled_task[i]->prefer_type] /3600.0;
						//make the vm.task = NULL
						for(int j=0; j<VMTP[scheduled_task[i]->prefer_type].size(); j++)
							if(VMTP[scheduled_task[i]->prefer_type][j]->tk == scheduled_task[i])
							{
								VMTP[scheduled_task[i]->prefer_type][j]->tk = NULL;
								break;
							}
					}
				}				
				//step 3
				for(int i=0; i<types; i++)
				{
					int size1 = VMTP[i].size();			
					
					for(int j=0; j<size1; j++)
					{
						VMTP[i][j]->life_time += 1;				
					}
				}
				for(int i=0; i<ready_task.size(); i++)//////////////////////////////////if >0
					if(ready_task[i]->readyCountdown > 0)
						ready_task[i]->readyCountdown -= 1;
				t += 1;

				condition = false;
				int unfinishednum = 0;
				for(int ji=0; ji<jobs.size(); ji++){
					vp = vertices((*jobs[ji]->g));
					for(int i=0; i < (*vp.second - *vp.first ); i++){
						if((*jobs[ji]->g)[i].status!= finished)
						{
							condition = true;
							unfinishednum += 1;
						}					
					}
				}								
			}while(condition);//there is a task not finished

			for(int i=0; i<types; i++)
			{
				int size1 = VMTP[i].size();						
				for(int j=0; j<size1; j++)
				{
					double runtime = VMTP[i][j]->life_time;
					moneycost += (priceOnDemand[i] * ceil(runtime/3600.0));
				}
			}
			printf("Money Cost: %.4f, Time: %.2f\n", moneycost, t);
			int id = omp_get_thread_num();
			printf("thread id is %d\n",id);
			float ave_time = 0.0;
			for(int i=0; i<jobs.size(); i++){
				vp = vertices((*jobs[i]->g));
				float executiontime = (*jobs[i]->g)[*(vp.second-2)].end_time - jobs[i]->arrival_time;
				if(executiontime > jobs[i]->deadline) {
					viol_private[id] += 1.0;
				}		
				ave_time += executiontime;
			}	
			cost_private[id] += moneycost;
			printf("average execution time of workflows is %f\n",ave_time/jobs.size());
			for(int i=0; i<jobs.size(); i++){
				delete jobs[i];
			}
			jobs.clear();
		}		
	}

	vp = vertices(*dag->g);
	for(; vp.first!=vp.second; vp.first++){
		free((*dag->g)[*vp.first].netDown);
		free((*dag->g)[*vp.first].netUp);
		free((*dag->g)[*vp.first].probestTime);
		free((*dag->g)[*vp.first].randomIO);
		free((*dag->g)[*vp.first].seqIO);
		free((*dag->g)[*vp.first].cumulativeTime);
		free((*dag->g)[*vp.first].randspot);
	}
	for(int i=0; i<24; i++) {
		violation += viol_private[i];
		ave_cost += cost_private[i];
	}
	violation /= (float)randomsize*num_jobs;
	ave_cost /= (float)randomsize*num_jobs;
	printf("deadline meeting rate is %f, average cost is %f\n",1.0-violation,ave_cost);
	std::clock_t endtime = std::clock();
	std::clock_t timeelapsed = (float)(endtime - starttime) / (float)CLOCKS_PER_SEC;
	printf("time elapsed for spea2 algorithm is: %.4f\n", timeelapsed);

	return;
}
