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
#include "singlebestelitism.h"
#include "valuefitness.h"
#include "vectorgenomeuniformmutation.h"
#include "stopcriterion.h"
#include "remainingtargetpopulationsizeiteration.h"
#include "probesystem.h"
#include <cassert>
#include <iostream>

// TODO: some kind of state that needs to be communicated when calculating fitness?
//       perhaps something else needs to be done depending on the generation?

using namespace errut;
using namespace std;

// TODO: do we need a separate class? Just a function
// TODO: feedback/probes?
class GeneticAlgorithm
{
public:
    GeneticAlgorithm(std::shared_ptr<ProbeSystem> probes = nullptr);
    virtual ~GeneticAlgorithm();

    bool_t run(GenomeFitnessCreation &gfc,
               shared_ptr<PopulationCrossover> crossover,
               shared_ptr<PopulationMutation> mutation,
               PopulationFitnessCalculation &fitnessCalc,
               StopCriterion &stopCriterion,
               size_t popSize,
               size_t minPopulationSize = 0,
               size_t maxPopulationSize = 0);
private:
    std::shared_ptr<ProbeSystem> m_probes;
};

GeneticAlgorithm::GeneticAlgorithm(std::shared_ptr<ProbeSystem> probes) : m_probes(probes)
{
}

GeneticAlgorithm::~GeneticAlgorithm()
{
}

// Note that the population size does not need to be constant throughout the loop,
// more could arise so that their fitness is calculated. This is why the population
// size is passed on to the populationcrossover
bool_t GeneticAlgorithm::run(GenomeFitnessCreation &gfc,
                             shared_ptr<PopulationCrossover> popCross,
                             shared_ptr<PopulationMutation> popMutation,
                             PopulationFitnessCalculation &fitnessCalc,
                             StopCriterion &stopCriterion,
                             size_t popSize,
                             size_t minPopulationSize,
                             size_t maxPopulationSize)
{
    bool_t r;
    auto population = make_shared<Population>();
    auto newPopulation = make_shared<Population>();
    auto refFitness = gfc.createEmptyFitness();
    if (!refFitness.get())
        return "Unable to create a reference fitness object";

    if (popSize == 0)
        return "No population size specified";
    
    if (maxPopulationSize == 0)
        maxPopulationSize = popSize;

    for (size_t i = 0 ; i < popSize ; i++)
    {
        auto g = gfc.createInitializedGenome();
        if (!g.get())
            return "Unable to create an inialized genome";

        auto f = refFitness->createCopy(false);
        population->append(make_shared<Individual>(g, f));
    }

    size_t generation = 0;

    if (!(r = fitnessCalc.calculatePopulationFitness({population})))
        return "Error calculating fitness: " + r.getErrorString();

    auto reportOnFitness = [&population, &generation, this]() -> bool_t
    {
        if (m_probes.get())
        {
            m_probes->setGeneration(generation);
            bool_t r = m_probes->inspect(ProbeSystem::FitnessCalculated, population);
            if (!r)
                return "Error inspecting population after fitness calculation: " + r;
        }
        return true;
    };

    if (!(r = reportOnFitness()))
        return r;

    while (true)
    {        
        population->setGenomesToSkipMutation(0);
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

        const size_t curPopSize = population->size();
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
            //       -> population now has member to keep track of this count
            if (!(r = popMutation->mutate(population)))
                return "Error in mutation: " + r.getErrorString();
        }

        generation++; // At this point we can call it the new generation

        if (!(r = fitnessCalc.calculatePopulationFitness({population})))
            return "Error calculating fitness: " + r.getErrorString();

        if (!(r = reportOnFitness()))
            return r;

        bool shouldStop = false;
        if (!(r = stopCriterion.analyze(popCross->getBestIndividuals(), generation, shouldStop)))
            return "Error in termination check: " + r.getErrorString();
        if (shouldStop)
            break;
    }

    if (m_probes.get())
    {
        if (!(r = m_probes->inspect(ProbeSystem::AlgorithmDone, popCross->getBestIndividuals())))
            return "Error inspecting best individuals upon algorithm end: " + r.getErrorString();
    }

    return true;
}

typedef int RealType;

class MyProbes : public ProbeSystem
{
public:
    MyProbes() { }
    ~MyProbes() { }

    bool_t inspect(EventType eventType, shared_ptr<Population> &population) override
    {
        cout << "Generation: " << getGeneration() << endl;
        cout << "EventType: " << eventType << endl;
        population->print();
        cout << endl;

        // TODO: check if all fitnesses are calculated?
        return true;
    }

    bool_t inspect(EventType eventType, shared_ptr<SelectionPopulation> &selPop) override
    {
        return true;
    }
    
    bool_t inspect(EventType eventType, const vector<shared_ptr<Individual>> &bestIndividuals) override
    {
        cout << "Best are: " << endl;
        for (auto &i : bestIndividuals)
            cout << i->toString() << endl;

        return true;
    }
};

class TestFactory : public GenomeFitnessCreation
{
public:
    TestFactory(unsigned long seed)
    {
        m_rng = make_shared<MersenneRandomNumberGenerator>(seed);
        m_mutation = make_shared<SingleThreadedPopulationMutation>(
            make_shared<VectorGenomeUniformMutation<RealType>>(0.2, 0, 100, m_rng));
        m_crossover = make_shared<SingleThreadedPopulationCrossover>(
            0.1,
            // make_shared<SimpleSortedPopulation>(make_shared<VectorFitnessComparison<RealType>>()),
            make_shared<SimpleSortedPopulation>(make_shared<ValueFitnessComparison<RealType>>()),
            make_shared<RankParentSelection>(2.5, m_rng),
            make_shared<UniformVectorGenomeCrossover<RealType>>(m_rng, false),
            make_shared<SingleBestElitism>(true, true),
            make_shared<RemainingTargetPopulationSizeIteration>(),
            m_rng
        );
        m_probes = make_shared<MyProbes>();
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

    shared_ptr<PopulationMutation> getPopulationMutation()
    {
        return m_mutation;
    }

    shared_ptr<PopulationCrossover> getPopulationCrossover()
    {
        return m_crossover;
    }

    shared_ptr<ProbeSystem> getProbeSystem() { return m_probes; }
private:
    shared_ptr<RandomNumberGenerator> m_rng;
    shared_ptr<SingleThreadedPopulationMutation> m_mutation;
    shared_ptr<SingleThreadedPopulationCrossover> m_crossover;
    shared_ptr<MyProbes> m_probes;
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

        RealType x = vg.getValues()[0];
        RealType y = vg.getValues()[1];
        RealType dx = x - (RealType)30;
        RealType dy = y - (RealType)20;

        //vf.getValues()[0] = dx*dx + dy*dy;
        vf.setValue(dx*dx + dy*dy);
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
    
    TestFactory factory(seed);

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

        FixedGenerationsStopCriterion stop(100);
        GeneticAlgorithm ga { factory.getProbeSystem() };

        //r = ga.run(factory, calc, 16, 0, 32);
        r = ga.run(factory,
                   factory.getPopulationCrossover(),
                   factory.getPopulationMutation(),
                   *calc, stop, 16);
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