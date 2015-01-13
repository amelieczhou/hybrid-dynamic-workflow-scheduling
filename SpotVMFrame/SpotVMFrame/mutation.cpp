#include "stdafx.h"
#include "mutation.h"

extern std::vector<individual* > generation; //initial size = populationsize 
extern std::vector<individual* > variation; //initially empty
extern std::vector<individual* > matingpool;
int mu = 10;
int solutionsize = 10;//20;//for epi
float variable_swap_probability = 0.9;//crossover
float variable_mutation_probability = 0.5; //mutation

int update_from_arc(std::vector<individual*>& population)
// Removes the individuals, that are not in the archive
{
     int size, result; 
     //int *keep;	 
     int i, current;

	 size = population.size();
	 std::vector<int > keep(size);
     result = 0;
    // keep = (int *) malloc(sizeof(int)*size);

     /*if (keep == NULL)
     {
		 std::cout<<"SPEA variator: out of memory"<<std::endl;
         exit(1);
     }*/  

     for(i=0; i < size; i++)
         keep[i] = population[i]->index;

     // sort the array of indexes to keep,
     // so we can go through the array and delete all indexes,
     // that are in between
     //qsort(keep, (size_t) size, sizeof(int), cmp_int);      
    
     // delete all indexes in global_population not found in keep array
	 for(int i=0; i<generation.size(); i++){
		 if(std::find(keep.begin(),keep.end(),generation[i]->index) == keep.end()){//no need to keep this one
			generation.erase(generation.begin()+i);
			i--;
		 }
	 }
	/* int count = 0;
	 current = generation[count]->index;
     for(i = 0; i < size; i++)
     {
          while(current < keep[i])
          {
			  generation[count] = NULL;
              generation.erase(generation.begin()+count);   
			  //count--;
              current = generation[count]->index;
          }
          if (current == keep[i])
          {
               count++;
			   if(count<generation.size())
					current = generation[count]->index;
          } // this one we keep
          else  // current must be bigger than keep[i],
                // something went wrong...
          { 
              
			  std::cout<<"identity in archive is not in the global population!"<<std::endl;
              exit(1);
          }
     }*/

     // delete the last individuals at end of list
     while(generation.size() > population.size())
		 generation.erase(generation.end()-1);
  
     //free(keep);
     return (0);
}

int variate(std::vector<individual*>& selected)
// Performs the real variation of individuals
// *selected points to the selected individuals
// *result_ids will contain the offspring individuals
{
     int result, i, k;
     result = 1;

     // copying all individuals from selected to global population
	 //for(int j=0; j<selected.size(); j++)
	 //	generation.push_back(selected[j]);
     // if odd number of individuals, last one is left as is
     if((((double)mu/2) - (int)(mu/2)) != 0) k = mu - 1; 
     else k = mu;

     // do recombination
     for(i = 0; i < k; i+= 2)  {  
	  if (rn_01() <= variable_swap_probability)	  {
		  single_point_crossover(selected[i],selected[i+1]);	      
	  }
     }
    // do mutation
     for(i = 0; i < mu; i++)   {
        if (rn_01() <= variable_mutation_probability) { 
            mutation(selected[i]);
		}  
     }
     
     return (0);
}

void single_point_crossover(individual *ind1, individual *ind2)
{
	 //select two points in the vector and swap
	 int start, end, tmp;
	 start = rn_integers(0,solutionsize);
	 while(start == solutionsize) 
		 start = rn_integers(0,solutionsize);
	 end = rn_integers(0,solutionsize);
	 while(end == solutionsize) 
		 end = rn_integers(0,solutionsize);
	 if(start > end){
		tmp = start;
		start = end;
		end = tmp;
	 }
	 for(int i=start; i<=end; i++){
		 tmp = ind1->solution[i];
		 ind1->solution[i] = ind2->solution[i];
		 ind2->solution[i] = tmp;
	 }
	 
     return;
}
void mutation(individual* ind){


     if (ind == NULL) {
		exit(1);
     }
     
	 
     for (int i = 0; i < solutionsize; i++)  {
		 //replacing mutation
		 if(rn_01() <= variable_mutation_probability){
			 int mut = rn_integers(0,types);
			 while(mut == types) mut = rn_integers(0,types);
			 ind->solution[i] = mut;
		 }
		 //reordering mutation

	 }
     
     return ;
}