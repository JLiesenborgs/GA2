#include "population.h"
#include "crossovermutation.h"
#include "vectorgenomefitness.h"
#include "mersennerandomnumbergenerator.h"
#include "uniformvectorgenomecrossover.h"
#include "singlethreadedpopulationfitnesscalculation.h"
#include "multithreadedpopulationfitnesscalculation.h"
#include "mpipopulationfitnesscalculation.h"
#include "simplesortedpopulation.h"
#include "rankparentselection.h"
#include "singlethreadedpopulationmutation.h"
#include "singlethreadedpopulationcrossover.h"
#include "gafactory.h"
#include <cassert>
#include <iostream>

using namespace errut;
using namespace std;

// TODO: elitism

template<class T>
class VectorGenomeUniformMutation : public GenomeMutation
{
public:
    VectorGenomeUniformMutation(double mutationFraction, T minValue, T maxValue, std::shared_ptr<RandomNumberGenerator> &rng)
        :  m_mutationFraction(mutationFraction), m_min(minValue), m_max(maxValue), m_rng(rng) { }
    ~VectorGenomeUniformMutation() { }

	errut::bool_t check(const Genome &genome) override
    {
        const VectorGenome<T> *pGenome = dynamic_cast<const VectorGenome<T>*>(&genome);
        if (!pGenome)
            return "Genome is not of the expected type";
        return true;
    }

	errut::bool_t mutate(Genome &genome, bool &isChanged) override
    {
        VectorGenome<T> &g = dynamic_cast<VectorGenome<T>&>(genome);
        vector<T> &v = g.getValues();
    
        for (auto &x : v)
        {
            // TODO: can we do this without as many random numbers? Floats? Generate indices first by some other means?
            if (m_rng->getRandomDouble() < m_mutationFraction)
            {
                x = m_rng->getRandomDouble(m_min, m_max);
                isChanged = true;
            }
        }
        return true;
    }
private:
    std::shared_ptr<RandomNumberGenerator> m_rng;
    double m_mutationFraction;
    double m_min, m_max;
};

class GeneticAlgorithm
{
public:
    GeneticAlgorithm();
    virtual ~GeneticAlgorithm();

    bool_t run(GAFactory &factory,
               PopulationFitnessCalculation &fitnessCalc,
               size_t popSize,
               size_t minPopulationSize = 0,
               size_t maxPopulationSize = 0);

};

GeneticAlgorithm::GeneticAlgorithm()
{
}

GeneticAlgorithm::~GeneticAlgorithm()
{
}

// Note that the population size does not need to be constant throughout the loop,
// more could arise so that their fitness is calculated. This is why the population
// size is passed on to the populationcrossover
bool_t GeneticAlgorithm::run(GAFactory &factory,
                             PopulationFitnessCalculation &fitnessCalc,
                             size_t popSize,
                             size_t minPopulationSize,
                             size_t maxPopulationSize)
{
    auto population = make_shared<Population>();
    auto newPopulation = make_shared<Population>();
    auto refFitness = factory.createEmptyFitness();
    auto popMutation = factory.getPopulationMutation();
    auto popCross = factory.getPopulationCrossover();
    bool_t r;

    if (maxPopulationSize == 0)
        maxPopulationSize = popSize;

    for (size_t i = 0 ; i < popSize ; i++)
    {
        auto g = factory.createInitializedGenome();
        auto f = refFitness->createCopy(false);
        population->m_individuals.push_back(make_shared<Individual>(g, f));
    }

    if (!(r = fitnessCalc.calculatePopulationFitness({population})))
        return "Error calculating fitness: " + r.getErrorString();    

    for (int generation = 0 ; generation < 100 ; generation++)
    {        
        cout << "Generation " << generation << ": " << endl;
        population->print();

        if (popCross.get())
        {
            if (generation == 0)
            {
                if (!(r = popCross->check(population)))
                    return "Error in population crossover check: " + r.getErrorString();
            }

            
            if (!(r = popCross->createNewPopulation(population, popSize)))
                return "Error creating new population: " + r.getErrorString();
        }

        const size_t curPopSize = population->m_individuals.size();
        if (curPopSize > maxPopulationSize)
            return "Population size (" + to_string(curPopSize) + ") exceeds maximum (" + to_string(maxPopulationSize) + ")";
        if (curPopSize < minPopulationSize)
            return "Population size (" + to_string(curPopSize) + ") is less than minimum (" + to_string(minPopulationSize) + ")";

        if (popMutation.get())
        {
            if (generation == 0)
            {
                if (!(r = popMutation->check(population)))
                    return "Error checking mutation: " + r.getErrorString();
            }

            // TODO: how best to skip mutation on introduced elitist solutions?
            if (!(r = popMutation->mutate(population)))
                return "Error in mutation: " + r.getErrorString();
        }

        if (!(r = fitnessCalc.calculatePopulationFitness({population})))
            return "Error calculating fitness: " + r.getErrorString();

        // TODO: elitism -> also needs what's best, perhaps the selectionpopulation is
        //       the way to track this?

        // TODO: somehow check that all fitness values have been calculated??
        //       some calculation postprocessor? Might also be a good way to
        //       store the scale factor from own GA into the genome
    }

    population->print();

    cout << "Best are: " << endl;
    for (auto &i : popCross->getBestIndividuals())
        cout << i->toString() << endl;

    return true;
}

class TestGAFactory : public GAFactory
{
public:
    TestGAFactory(unsigned long seed)
    {
        m_rng = make_shared<MersenneRandomNumberGenerator>(seed);
    }

    shared_ptr<Genome> createInitializedGenome() override
    {
        auto g = make_shared<FloatVectorGenome>(2);
        for (auto &x : g->getValues())
            x = m_rng->getRandomFloat();
        return g;
    }

    shared_ptr<Fitness> createEmptyFitness() override
    {
        return make_shared<FloatVectorFitness>(1);
    }

    shared_ptr<PopulationMutation> getPopulationMutation()
    {
        // TODO: this creates a new one, return a previously created one
        return make_shared<SingleThreadedPopulationMutation>(
            make_shared<VectorGenomeUniformMutation<float>>(0.2, 0, 1, m_rng));
    }

    shared_ptr<PopulationCrossover> getPopulationCrossover()
    {
        // TODO: this creates a new one, return a previously created one
        return make_shared<SingleThreadedPopulationCrossover>(
            0.1,
            make_shared<SimpleSortedPopulation>(make_shared<VectorFitnessComparison<float>>()),
            make_shared<RankParentSelection>(2.5, m_rng),
            make_shared<UniformVectorGenomeCrossover<float>>(m_rng, false),
            m_rng
        );
    }
private:
    shared_ptr<RandomNumberGenerator> m_rng;
};

class TestFitnessCalculation : public GenomeFitnessCalculation
{
public:
    TestFitnessCalculation() { }
    ~TestFitnessCalculation() { }

    bool_t calculate(const Genome &g, Fitness &f)
    {
        const FloatVectorGenome &vg = static_cast<const FloatVectorGenome &>(g);
        FloatVectorFitness &vf = static_cast<FloatVectorFitness &>(f);

        assert(vg.getValues().size() == 2);
        assert(vf.getValues().size() == 1);
        
        float x = vg.getValues()[0];
        float y = vg.getValues()[1];
        float dx = x - 0.3f;
        float dy = y - 0.2f;

        vf.getValues()[0] = dx*dx + dy*dy;
        return true;
    }
};

const int CalcTypeSingle = 0;
const int CalcTypeMulti = 1;
const int CalcTypeMPI = 2;

const int calcType = CalcTypeMPI;

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
    shared_ptr<MPIPopulationFitnessCalculation> calcMPI;
    shared_ptr<MPIEventDistributor> mpiDist;

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
    else if (calcType == 2)
    {
        // MPI
        mpiDist = make_shared<MPIEventDistributor>();
        calcMPI =  make_shared<MPIPopulationFitnessCalculation>(mpiDist);
        mpiDist->setHandler(MPIEventHandler::Calculation, calcMPI);
        calc = calcMPI;
    }
    else
        return "Unknown calculation type " + to_string(calcType);
    
    TestGAFactory factory(seed);

    if (calcType == 1)
    {
        r = calcMulti->initThreadPool({
            make_shared<TestFitnessCalculation>(),
            make_shared<TestFitnessCalculation>(),
            make_shared<TestFitnessCalculation>(),
            make_shared<TestFitnessCalculation>() });
        if (!r)
            return "Couldn't init thread based fitness calculator: " + r.getErrorString();
    }
    else if (calcType == 2)
    {
        auto refGenome = factory.createInitializedGenome();
        auto refFitness = factory.createEmptyFitness();
        r = calcMPI->init(*refGenome, *refFitness, 
                          make_shared<SingleThreadedPopulationFitnessCalculation>(make_shared<TestFitnessCalculation>()));
        if (!r)
            return "Couldn't init MPI fitness calculator: " + r.getErrorString();
    }

    if (rank == 0) // Should also work for the non MPI versions
    {
        // At this point, on the other ranks, the event loop will be waiting what
        // to do, so we should send a Done signal 
        auto cleanup = [mpiDist]()
        {
            if (mpiDist.get())
                mpiDist->signal(MPIEventHandler::Done);
        };

        GeneticAlgorithm ga;

        //r = ga.run(factory, calc, 16, 0, 32);
        r = ga.run(factory, *calc, 16);
        if (!r)
        {
            cleanup();
            return "Error running GA: " + r.getErrorString();
        }

        cleanup();
    }
    else
    {
        r = mpiDist->eventLoop();
        if (!r)
            return "Error in event loop: " + r.getErrorString();
    }

    return true;
}

#if 0
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
    return 0;
}
#endif