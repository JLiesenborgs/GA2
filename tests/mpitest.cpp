#include "vectorgenomefitness.h"
#include "mpipopulationfitnesscalculation.h"
#include "singlethreadedpopulationfitnesscalculation.h"
#include <vector>
#include <iostream>
#include <sstream>
#include <cstring>
#include <cassert>
#include <cstdio>

// TODO: namespace

// TODO: for genome and fitness, a copyInto function? So that memory doesn't need to be
//       allocated/reallocated often


using namespace std;
using namespace errut;

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

		int rank;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);

		for (float &x : pFitness->getValues())
			x = pGenome->getValues()[0];
		pFitness->getValues()[0] = (float)rank;

		fitness.setCalculated();
		return true;
	}
};

const int ROOT = 1;
const int LOOPS = 2;

int main_master(int argc, char *argv[])
{
	cerr << "Master" << endl;

	shared_ptr<MPIEventDistributor> mpiDist = make_shared<MPIEventDistributor>(ROOT, MPI_COMM_WORLD);
	shared_ptr<MPIPopulationFitnessCalculation> calc = make_shared<MPIPopulationFitnessCalculation>(mpiDist);

	// mpiDist->setHandler(MPIEventHandler::Calculation, calc);

	shared_ptr<DummyGenomeFitnessCalculation> genomeCalc = make_shared<DummyGenomeFitnessCalculation>();
	shared_ptr<PopulationFitnessCalculation> localCalc = make_shared<SingleThreadedPopulationFitnessCalculation>(genomeCalc);

	int numFloatGenome = 16;
	int numFloatFitness = 2;
	FloatVectorGenome genome(numFloatGenome);
	FloatVectorFitness fitness(numFloatFitness);

	auto r = calc->init(genome, fitness, localCalc, MPI_COMM_WORLD, ROOT);
	if (!r)
	{
		cerr << "Can't init MPIPopulationFitnessCalculation: " << r.getErrorString() << endl;
		MPI_Abort(MPI_COMM_WORLD, -1);
	}

	shared_ptr<Population> pop = make_shared<Population>();
	for (int i = 0 ; i < 32 ; i++)
		pop->append(make_shared<Individual>(
					make_shared<FloatVectorGenome>(numFloatGenome, i),
					make_shared<FloatVectorFitness>(numFloatFitness))); 

	for (int i = 0 ; i < LOOPS ; i++)
	{
		r = calc->calculatePopulationFitness({ pop });
		if (!r)
		{
			cerr << "Error in calculatePopulationFitness: " << r.getErrorString() << endl;
			MPI_Abort(MPI_COMM_WORLD, -1);
		}

		for (auto &i : pop->individuals())
			cout << i->toString() << endl;

		for (auto &i : pop->individuals())
			i->fitnessRef().setCalculated(false);
	}
	
	mpiDist->signal(MPIEventHandler::Done);

	cerr << "Master done\n";
	return 0;
}

int main_helper(int argc, char *argv[])
{
	cerr << "Helper" << endl;

	shared_ptr<MPIEventDistributor> mpiDist = make_shared<MPIEventDistributor>(ROOT, MPI_COMM_WORLD);
	shared_ptr<MPIPopulationFitnessCalculation> calc = make_shared<MPIPopulationFitnessCalculation>(mpiDist);

	mpiDist->setHandler(MPIEventHandler::Calculation, calc);

	shared_ptr<DummyGenomeFitnessCalculation> genomeCalc = make_shared<DummyGenomeFitnessCalculation>();
	shared_ptr<PopulationFitnessCalculation> localCalc = make_shared<SingleThreadedPopulationFitnessCalculation>(genomeCalc);
	
	FloatVectorGenome genome; // should not now exact layout here, just type
	FloatVectorFitness fitness; // same

	auto r = calc->init(genome, fitness, localCalc, MPI_COMM_WORLD, ROOT);
	if (!r)
	{
		cerr << "Can't init MPIPopulationFitnessCalculation in helper: " << r.getErrorString() << endl;
		MPI_Abort(MPI_COMM_WORLD, -1);
	}

	if (!(r = mpiDist->eventLoop()))
	{
		cerr <<	"Error in event loop: " << r.getErrorString() << endl;
		MPI_Abort(MPI_COMM_WORLD, -1);
	}

	cerr << "Helper done\n";
	return 0;
}

int main2(int argc, char *argv[])
{
	int rank, procs;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &procs);
	printf("%d/%d\n", rank, procs);

	if (rank == ROOT)
		return main_master(argc, argv);
	return main_helper(argc, argv);
}

int main(int argc, char *argv[])
{
	MPI_Init(&argc, &argv);
	int r = main2(argc, argv);
	MPI_Finalize();
	return r;
}

