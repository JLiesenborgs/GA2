#include "multithreadedpopulationfitnesscalculation.h"
#include "vectorgenomefitness.h"
#include <vector>
#include <iostream>
#include <typeinfo>
#include <sstream>
#include <cstring>
#include <cassert>
#include <cstdio>

// TODO: for genome and fitness, a copyInto function? So that memory doesn't need to be
//       allocated/reallocated often

using namespace std;
using namespace errut;
using namespace mogal2;

class DummyGenomeFitnessCalculation : public GenomeFitnessCalculation
{
public:
	DummyGenomeFitnessCalculation() { }
	~DummyGenomeFitnessCalculation() { }
	bool_t pollCalculate(const Genome &genome, Fitness &fitness)
	{
		// TODO: move checking code to separate routine
		const FloatVectorGenome *pGenome = dynamic_cast<const FloatVectorGenome *>(&genome);
		if (!pGenome)
			return "Genome is not of expected type";
		FloatVectorFitness *pFitness = dynamic_cast<FloatVectorFitness *>(&fitness);
		if (!pFitness)
			return "Fitness is not of expected type";

		for (float &x : pFitness->getValues())
			x = pGenome->getValues()[0];

		fitness.setCalculated();
		return true;
	}
};


const int NUMTHREADS = 4;
const int LOOPS = 2;

int main(int argc, char *argv[])
{
	MultiThreadedPopulationFitnessCalculation calc;
	
	vector<shared_ptr<GenomeFitnessCalculation>> genomeCalculators;
	for (int i = 0 ; i < NUMTHREADS ; i++)
		genomeCalculators.push_back(make_shared<DummyGenomeFitnessCalculation>());

	int numFloatGenome = 16;
	int numFloatFitness = 2;
	FloatVectorGenome genome(numFloatGenome);
	FloatVectorFitness fitness(numFloatFitness);

	auto r = calc.initThreadPool(genomeCalculators);
	if (!r)
	{
	 	cerr << "Can't init MultiThreadedPopulationFitnessCalculation: " << r.getErrorString() << endl;
	 	return -1;
	}

	shared_ptr<Population> pop = make_shared<Population>();
	for (int i = 0 ; i < 32 ; i++)
		pop->append(make_shared<Individual>(
					make_shared<FloatVectorGenome>(numFloatGenome, i),
					make_shared<FloatVectorFitness>(numFloatFitness))); 

	for (int i = 0 ; i < LOOPS ; i++)
	{
		r = calc.calculatePopulationFitness({ pop });
		if (!r)
		{
			cerr << "Error in calculatePopulationFitness: " << r.getErrorString() << endl;
			MPI_Abort(MPI_COMM_WORLD, -1);
		}

		for (auto &i : pop->individuals())
			cout << i->toString() << endl;
		cout << endl;

		for (auto &i : pop->individuals())
			i->fitnessRef().setCalculated(false);
	}
	return 0;
}

