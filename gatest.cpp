#include "population.h"
#include "vectorgenomefitness.h"
#include "singlethreadedpopulationfitnesscalculation.h"
#include <random>
#include <cassert>
#include <iostream>

using namespace errut;
using namespace std;

class RandomNumberGenerator
{
public:
    RandomNumberGenerator(unsigned int seed) : m_rng(seed) { }
    virtual ~RandomNumberGenerator() { }

    double getRandomDouble(double min = 0.0, double max = 1.0)
    {
        uniform_real_distribution<double> u(min, max);
        return u(m_rng);
    }
    
    float getRandomFloat(float min = 0.0, float max = 1.0)
    {
        uniform_real_distribution<float> u(min, max);
        return u(m_rng);
    }

private:
    std::mt19937 m_rng;
};

class GeneticAlgorithm
{
public:
    GeneticAlgorithm(shared_ptr<RandomNumberGenerator> rng,
        shared_ptr<PopulationFitnessCalculation> calc);
    virtual ~GeneticAlgorithm();

    bool_t run(size_t popSize);
protected:
    shared_ptr<Genome> createInitialGenome();
    shared_ptr<Fitness> createEmptyFitness();

    shared_ptr<RandomNumberGenerator> m_rng;
    shared_ptr<PopulationFitnessCalculation> m_fitnessCalc;
};

GeneticAlgorithm::GeneticAlgorithm(shared_ptr<RandomNumberGenerator> rng,
    shared_ptr<PopulationFitnessCalculation> calc)
    : m_rng(rng), m_fitnessCalc(calc)
{
}

GeneticAlgorithm::~GeneticAlgorithm()
{
}

bool_t GeneticAlgorithm::run(size_t popSize)
{
    auto population = make_shared<Population>();
    auto refFitness = createEmptyFitness();

    for (size_t i = 0 ; i < popSize ; i++)
    {
        auto g = createInitialGenome();
        auto f = refFitness->createCopy(false);
        population->m_individuals.push_back(make_shared<Individual>(g, f));
    }

    bool_t r = m_fitnessCalc->calculatePopulationFitness({population});
    if (!r)
        return "Error calculating fitness: " + r.getErrorString();
    
    for (auto &i : population->m_individuals)
        cout << i->m_genome->toString() << ": " << i->m_fitness->toString() << endl;
    return true;
}

shared_ptr<Genome> GeneticAlgorithm::createInitialGenome()
{
    auto g = make_shared<FloatVectorGenome>(2);
    for (auto &x : g->getValues())
        x = m_rng->getRandomFloat();
    return g;
}

shared_ptr<Fitness> GeneticAlgorithm::createEmptyFitness()
{
    return make_shared<FloatVectorFitness>(1);
}

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

int main(int argc, char const *argv[])
{
    random_device rd;
    shared_ptr<RandomNumberGenerator> rng = make_shared<RandomNumberGenerator>(rd());

    shared_ptr<SingleThreadedPopulationFitnessCalculation> calc = make_shared<SingleThreadedPopulationFitnessCalculation>(
        make_shared<TestFitnessCalculation>()
    );
    
    GeneticAlgorithm ga { rng, calc };
    auto r = ga.run(16);
    if (!r)
    {
        cerr << "Error running GA: " << r.getErrorString() << endl;
        return -1;
    }
    return 0;
}

