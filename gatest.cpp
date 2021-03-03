#include "population.h"
#include "vectorgenomefitness.h"
#include "singlethreadedpopulationfitnesscalculation.h"
#include <random>
#include <cassert>
#include <iostream>
#include <algorithm>

using namespace errut;
using namespace std;

// For a one
class SimpleSortedPopulation
{
public:
    SimpleSortedPopulation(std::shared_ptr<FitnessComparison> fitComp,
        int objectiveNumber = 0) : m_objectiveNumber(objectiveNumber), m_fitnessComp(fitComp) { }
    virtual ~SimpleSortedPopulation() { }

    void setObjectiveNumber(int objectiveNumber) { m_objectiveNumber = objectiveNumber; }
    
    // In place!
    bool_t sortPopulation(Population &population);
private:
    int m_objectiveNumber;
    std::shared_ptr<FitnessComparison> m_fitnessComp;
};

bool_t SimpleSortedPopulation::sortPopulation(Population &population)
{
    FitnessComparison &cmp = *m_fitnessComp;
    const int N = m_objectiveNumber;
    auto comp = [&cmp, N](auto &i1, auto &i2)
    {
        return cmp.isFitterThan(*i1->m_fitness, *i2->m_fitness, N);
    };

    sort(population.m_individuals.begin(), population.m_individuals.end(), comp);
    return true;
}

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
    shared_ptr<FitnessComparison> getFitnessComparison();

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

void printPopulation(const Population &population)
{
    cout << "Population: " << endl;
    for (auto &i : population.m_individuals)
        cout << i->m_genome->toString() << ": " << i->m_fitness->toString() << endl;
    cout << endl;
}

bool_t GeneticAlgorithm::run(size_t popSize)
{
    auto population = make_shared<Population>();
    auto refFitness = createEmptyFitness();
    auto fitnessComparison = getFitnessComparison();

    bool_t r = fitnessComparison->check(*refFitness);
    if (!r)
        return "Error checking fitness comparison: " + r.getErrorString();

    for (size_t i = 0 ; i < popSize ; i++)
    {
        auto g = createInitialGenome();
        auto f = refFitness->createCopy(false);
        population->m_individuals.push_back(make_shared<Individual>(g, f));
    }
    
    if (!(r = m_fitnessCalc->calculatePopulationFitness({population})))
        return "Error calculating fitness: " + r.getErrorString();
        
    SimpleSortedPopulation sorter(fitnessComparison);

    printPopulation(*population);

    if (!(r = sorter.sortPopulation(*population)))
        return "Error in sort: " + r.getErrorString();

    printPopulation(*population);

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

shared_ptr<FitnessComparison> GeneticAlgorithm::getFitnessComparison()
{
    return make_shared<VectorFitnessComparison<float>>();
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

