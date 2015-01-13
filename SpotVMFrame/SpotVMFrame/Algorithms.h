#include "InstanceConfig.h"
#include <boost/random/gamma_distribution.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/graph/astar_search.hpp>
//#include <boost/math/distributions/gamma.hpp>
//#include <boost/math/distributions/normal.hpp>

class Dyna{
public:
	DAG* dag;

	void Simulate();
	Dyna() {}
	Dyna(DAG* g) {dag = g;}

};

class DynaNPD{
public:
	DAG* dag;

	void Simulate();
	DynaNPD() {}
	DynaNPD(DAG* g) {dag=g;}
};

class DynaNS{
public:
	DAG* dag;

	void Simulate();
	DynaNS() {}
	DynaNS(DAG* g) {dag=g;}
};

class SpotOnly{
public:
	DAG* dag;
	
	void Simulate();
	void bidpDeter();
	SpotOnly() {}
	SpotOnly(DAG* g) {dag=g;}	
};

class Autoscaling{
public:
	DAG* dag;

	void Initialize();
	void Simulate();
	Autoscaling(){}
	Autoscaling(DAG* g){dag = g;}
};

class test{
public:
	void getCDF(int i);
};

class MOEA{
public:
	DAG* dag;

	void Crossover();
	void Mutation(std::vector<individual* >&);
	void evalObjectives(individual*);
	void calFitnesses(std::vector<individual*>); //evaluate the fitnesses
	void Evolution();
	
	void Simulate();
	MOEA(){}
	MOEA(DAG* g){dag=g;}

};

class MOHEFT{
public:
	DAG* dag;

	void evalObjectives(individual* indiv, int i);//evaluate till task i
	std::vector<individual*> BRank();
	void Simulate();
	MOHEFT(){}
	MOHEFT(DAG* g){dag=g;}
};