#include "evolutionaryalgorithm.h"
#include "vectorgenomefitness.h"
#include "mersennerandomnumbergenerator.h"
#include "vectorgenomeuniformcrossover.h"
#include "singlethreadedpopulationfitnesscalculation.h"
#include "multithreadedpopulationfitnesscalculation.h"
#include "mpipopulationfitnesscalculation.h"
#include "simplesortedpopulation.h"
#include "rankparentselection.h"
#include "singlepopulationcrossover.h"
#include "singlebestelitism.h"
#include "valuefitness.h"
#include "vectorgenomeuniformmutation.h"
#include "vectorgenomefractionalmutation.h"
#include "remainingtargetpopulationsizeiteration.h"
#include "ndsortedpopulation.h"
#include "basicnondominatedsetcreator.h"
#include "fitnessbasedduplicateremoval.h"
#include "populationreusecreation.h"
#include "testfunctions.h"
#include <cassert>
#include <iostream>

// TODO: some kind of state that needs to be communicated when calculating fitness?
//	   perhaps something else needs to be done depending on the generation?

using namespace errut;
using namespace std;
using namespace eatk;

typedef double RealType;

class MyGA : public EvolutionaryAlgorithm
{
protected:
	bool_t onBeforeFitnessCalculation(size_t generation, const std::shared_ptr<Population> &population) override
	{
		cout << "BEFORE Fitness calculation for generation: " << generation << endl;
		population->print();
		cout << endl;
		return true;
	}

	bool_t onFitnessCalculated(size_t generation, const std::shared_ptr<Population> &population) override
	{
		cout << "Fitness calculated for generation: " << generation << endl;
		population->print();
		cout << endl;
		return true;
	}

	bool_t onAlgorithmDone(size_t generation, const std::vector<std::shared_ptr<Individual>> &bestIndividuals) override
	{
		cout << "Ending after " << generation << " generations, best are: " << endl;
		for (auto &i : bestIndividuals)
			cout << i->toString() << endl;
		return true;
	}
};

class MyCrossOver : public SinglePopulationCrossover
{
public:
	MyCrossOver(shared_ptr<RandomNumberGenerator> rng, shared_ptr<GenomeMutation> mutation)
		: SinglePopulationCrossover(0.1, false,
			// make_shared<SimpleSortedPopulation>(make_shared<VectorFitnessComparison<RealType>>()),
			//make_shared<SimpleSortedPopulation>(make_shared<ValueFitnessComparison<RealType>>()),
			make_shared<NDSortedPopulation>(
				make_shared<BasicNonDominatedSetCreator>(make_shared<ValueFitnessComparison<RealType>>(), 1),
				make_shared<FitnessBasedDuplicateRemoval>(make_shared<ValueFitnessComparison<RealType>>(), 1)
			),
			make_shared<RankParentSelection>(2.5, rng),
			make_shared<VectorGenomeUniformCrossover<RealType>>(rng, false),
			mutation,
			make_shared<SingleBestElitism>(true, mutation),
			make_shared<RemainingTargetPopulationSizeIteration>(),
			rng
		)
	{
	}

protected:
	bool_t onSelectionPopulationProcessed(size_t generation, const std::shared_ptr<SelectionPopulation> &selPop) override
	{
		const NDSortedPopulationInterface *pSortPop = dynamic_cast<const NDSortedPopulationInterface *>(selPop.get());
		if (!pSortPop)
			return "Selection population is not of expected type";
		cout << "Sorted population for generation " << generation << ":" << endl;
		for (size_t s = 0 ; s < pSortPop->getNumberOfSets() ; s++)
			for (size_t i = 0 ; i < pSortPop->getSetSize(s) ; i++)
				cout << pSortPop->getIndividual(s, i)->toString() << endl;
		return true;
	}
};

struct AdjustClass
{
	static double adjustValue(double oldVal, double hardMin, double hardMax, RandomNumberGenerator &rng)
	{
		double r = rng.getRandomDouble(0.9, 1.1);
		return oldVal * r;
	}

	static double adjustValue(double oldVal, double refScale, double hardMin, double hardMax, RandomNumberGenerator &rng)
	{
		double r = rng.getRandomDouble(-refScale/2.0, refScale/2.0);
		return oldVal + r;
	}
};

class TestFactory : public IndividualCreation
{
public:
	TestFactory(unsigned long seed)
	{
		m_rng = make_shared<MersenneRandomNumberGenerator>(seed);
		m_mutation = make_shared<VectorGenomeUniformMutation<RealType>>(0.2, 0, 100, m_rng);
		m_smallerMutation = make_shared<VectorGenomeFractionalMutation<RealType,AdjustClass,false>>(0.2, m_rng);
		m_smallerMutationRefScale = make_shared<VectorGenomeFractionalMutation<RealType,AdjustClass,true>>(0.2, 0.1, m_rng);
		m_crossover = make_shared<MyCrossOver>(m_rng, m_mutation);
		m_crossoverSmallerMut = make_shared<MyCrossOver>(m_rng, m_smallerMutation);
		m_crossoverSmallerMutRefScale = make_shared<MyCrossOver>(m_rng, m_smallerMutationRefScale);
	}

	shared_ptr<Genome> createInitializedGenome() override
	{
		auto g = make_shared<VectorGenome<RealType>>(2);
		for (auto &x : g->getValues())
			x = (RealType)(m_rng->getRandomDouble()*100.0);
		return g;
	}

	shared_ptr<Fitness> createEmptyFitness() override
	{
		return make_shared<ValueFitness<RealType>>();
		//return make_shared<VectorFitness<RealType>>(1);
	}

	shared_ptr<PopulationEvolver> getPopulationCrossover()
	{
		return m_crossover;
	}

	shared_ptr<PopulationEvolver> getPopulationCrossoverSmallerMutation()
	{
		return m_crossoverSmallerMut;
	}

	shared_ptr<PopulationEvolver> getPopulationCrossoverSmallerMutationRefScale()
	{
		return m_crossoverSmallerMutRefScale;
	}

private:
	shared_ptr<RandomNumberGenerator> m_rng;
	shared_ptr<SinglePopulationCrossover> m_crossover;
	shared_ptr<SinglePopulationCrossover> m_crossoverSmallerMut, m_crossoverSmallerMutRefScale;
	shared_ptr<GenomeMutation> m_mutation;
	shared_ptr<GenomeMutation> m_smallerMutation, m_smallerMutationRefScale;
};

class TestFitnessCalculation : public GenomeFitnessCalculation
{
public:
	TestFitnessCalculation() { }
	~TestFitnessCalculation() { }

	bool_t calculate(const Genome &g, Fitness &f)
	{
		const VectorGenome<RealType> &vg = static_cast<const VectorGenome<RealType> &>(g);
		assert(vg.getValues().size() == 2);

		//VectorFitness<RealType> &vf = static_cast<VectorFitness<RealType> &>(f);
		//assert(vf.getValues().size() == 1);

		ValueFitness<RealType> &vf = static_cast<ValueFitness<RealType> &>(f);

#if 0
		RealType x = vg.getValues()[0];
		RealType y = vg.getValues()[1];
		RealType dx = x - (RealType)30;
		RealType dy = y - (RealType)20;

		//vf.getValues()[0] = dx*dx + dy*dy;
		vf.setValue(dx*dx + dy*dy);
#else
		eatk::testfunctions::Rosenbrock tf({-10.0,10.0});
		double val = tf.calculate1(vg.getValues());
		vf.setValue(val);
#endif
		return true;
	}
};

const int CalcTypeSingle = 0;
const int CalcTypeMulti = 1;
const int CalcTypeMPI = 2;

const int calcType = CalcTypeMulti;

bool_t real_main(int argc, char *argv[], int rank)
{
	bool_t r;
	random_device rd;
	unsigned int seed = rd();
	if (argc > 1)
		seed = atoi(argv[1]);
	
	cout << "Seed: " << seed << endl;

	shared_ptr<SingleThreadedPopulationFitnessCalculation> calcSingle;
	shared_ptr<MultiThreadedPopulationFitnessCalculation> calcMulti;
#ifdef EATKCONFIG_MPISUPPORT
	shared_ptr<MPIPopulationFitnessCalculation> calcMPI;
	shared_ptr<MPIEventDistributor> mpiDist;
#endif // EATKCONFIG_MPISUPPORT

	shared_ptr<PopulationFitnessCalculation> calc;
	if (calcType == 0)
	{
		// single threaded, nothing to do
		calcSingle = make_shared<SingleThreadedPopulationFitnessCalculation>(make_shared<TestFitnessCalculation>());
		calc = calcSingle;
	}
	else if (calcType == 1)
	{		
		// Multi threaded
		calcMulti = make_shared<MultiThreadedPopulationFitnessCalculation>();
		calc = calcMulti;
	}
#ifdef EATKCONFIG_MPISUPPORT
	else if (calcType == 2)
	{
		// MPI
		mpiDist = make_shared<MPIEventDistributor>();
		calcMPI =  make_shared<MPIPopulationFitnessCalculation>(mpiDist);
		mpiDist->setHandler(MPIEventHandler::Calculation, calcMPI);
		calc = calcMPI;
	}
#endif // EATKCONFIG_MPISUPPORT
	else
		return "Unknown calculation type " + to_string(calcType);
	
	TestFactory factory(seed);

	if (calcType == 0)
	{
		// Nothing to do
	}
	else if (calcType == 1)
	{
		r = calcMulti->initThreadPool({
			make_shared<TestFitnessCalculation>(),
			// make_shared<TestFitnessCalculation>(),
			// make_shared<TestFitnessCalculation>(),
			make_shared<TestFitnessCalculation>() });
		if (!r)
			return "Couldn't init thread based fitness calculator: " + r.getErrorString();
	}
#ifdef EATKCONFIG_MPISUPPORT
	else if (calcType == 2)
	{
		auto refGenome = factory.createInitializedGenome();
		auto refFitness = factory.createEmptyFitness();
		r = calcMPI->init(*refGenome, *refFitness, 
						  make_shared<SingleThreadedPopulationFitnessCalculation>(make_shared<TestFitnessCalculation>()));
		if (!r)
			return "Couldn't init MPI fitness calculator: " + r.getErrorString();
	}
#endif // EATKCONFIG_MPISUPPORT
	else
		return "Unknown calculation type " + to_string(calcType);

	size_t numGenerations = 1000;
	size_t popSize = 128;
	if (rank == 0) // Should also work for the non MPI versions
	{
		#ifdef EATKCONFIG_MPISUPPORT

		// At this point, on the other ranks, the event loop will be waiting what
		// to do, so we should send a Done signal 
		auto cleanup = [mpiDist]()
		{
			if (mpiDist.get())
				mpiDist->signal(MPIEventHandler::Done);
		};
		#else
		auto cleanup = [](){};
		#endif // EATKCONFIG_MPISUPPORT

		shared_ptr<PopulationReuseCreation> reuseCreation;

		{
			FixedGenerationsStopCriterion stop(numGenerations);
			MyGA ga;

			r = ga.run(factory,
					   *factory.getPopulationCrossover(),
					   *calc, stop, popSize, 0, popSize*2+2);
			if (!r)
			{
				cleanup();
				return "Error running GA: " + r.getErrorString();
			}

			reuseCreation = make_shared<PopulationReuseCreation>(ga.getPopulation());
		}

		// Run a second time with smaller mutations, continuing from the previous population
		{
			FixedGenerationsStopCriterion stop(numGenerations);
			MyGA ga;

			r = ga.run(*reuseCreation,
					   *factory.getPopulationCrossoverSmallerMutationRefScale(),
					   *calc, stop, popSize, 0, popSize*2+2);
			if (!r)
			{
				cleanup();
				return "Error running second GA: " + r.getErrorString();
			}

			reuseCreation = make_shared<PopulationReuseCreation>(ga.getPopulation());
		}

		// And a third time with other smalr mutations, continuing from the previous population
		{
			FixedGenerationsStopCriterion stop(numGenerations);
			MyGA ga;

			r = ga.run(*reuseCreation,
					   *factory.getPopulationCrossoverSmallerMutation(),
					   *calc, stop, popSize, 0, popSize*2+2);
			if (!r)
			{
				cleanup();
				return "Error running second GA: " + r.getErrorString();
			}
		}

		cleanup();
	}
	else
	{
		#ifdef EATKCONFIG_MPISUPPORT
		r = mpiDist->eventLoop();
		if (!r)
			return "Error in event loop: " + r.getErrorString();
		#endif // EATKCONFIG_MPISUPPORT
	}

	return true;
}

#if 1
int main(int argc, char *argv[])
{
	bool_t r = real_main(argc, argv, 0);
	if (!r)
	{
		cerr << "Error: " << r.getErrorString() << endl;
		return -1;
	}
	return 0;
}
#else
int main(int argc, char *argv[])
{
	#ifdef EATKCONFIG_MPISUPPORT

	MPI_Init(&argc, &argv);
	int rank = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	bool_t r = real_main(argc, argv, rank);
	if (!r)
	{
		cerr << "Error: " << r.getErrorString() << endl;
		MPI_Abort(MPI_COMM_WORLD, -1);
	}
	MPI_Finalize();
	#else
	cerr << "No MPI support enabled" << endl;
	#endif // EATKCONFIG_MPISUPPORT
	return 0;
}
#endif
