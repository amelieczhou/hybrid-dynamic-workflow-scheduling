#include "InstanceConfig.h"

#define STANDARD_SIZE 32200 


int cmp_int(const void *p_i1, const void *p_i2);

int update_from_arc(std::vector<individual*>& );

int variate(std::vector<individual*>&);

void single_point_crossover(individual*, individual*);

void mutation(individual*);