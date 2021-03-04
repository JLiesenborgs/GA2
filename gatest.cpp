#include "population.h"
#include "crossovermutation.h"
#include "vectorgenomefitness.h"
#include "mersennerandomnumbergenerator.h"
#include "uniformvectorgenomecrossover.h"
#include "singlethreadedpopulationfitnesscalculation.h"
#include <random>
#include <cmath>
#include <cassert>
#include <iostream>
#include <algorithm>

using namespace errut;
using namespace std;

// For a one-objective scenario
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

class RankParentSelection : public ParentSelection
{
public:
    RankParentSelection(double beta, shared_ptr<RandomNumberGenerator> rng) : m_beta(beta), m_rng(rng) { }
    ~RankParentSelection() { }

    // TODO: what to use as input here? Just a population? Sorted population is needed
    //       what with non-dominated sort?
    //       we don't need sorting for all kinds of selection mechanisms
    //       seems that an other kind of input is needed, and a check will be needed then
    bool_t selectParents(const Population &pop, std::vector<std::shared_ptr<Genome>> &parents) override
    {
        int popSize = (int)pop.m_individuals.size();
        if (popSize < 2)
            return "Not enough population members";
        
        auto getIndex = [this, popSize]()
        {
            double x = m_rng->getRandomDouble();
            double y = 1.0-std::pow(x, 1.0/(1.0+m_beta));
            int idx = (int)(y*(double)popSize);
            if (idx >= popSize)
                return popSize-1; // Population size is at least one
            return idx;
        };

        int idx1 = getIndex();
        int idx2 = getIndex();

        if (idx1 == idx2)
            idx2 = (idx1+1)%popSize;

        parents.resize(2);
        parents[0] = pop.m_individuals[idx1]->m_genome;
        parents[1] = pop.m_individuals[idx2]->m_genome;

        return true;
    }
private:
    const double m_beta;
    shared_ptr<RandomNumberGenerator> m_rng;
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
    auto newPopulation = make_shared<Population>();
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
    RankParentSelection selection(2.5, m_rng);
    UniformVectorGenomeCrossover<float> crossOver(m_rng, false);
    vector<shared_ptr<Genome>> parents, offspring;

    for (int generation = 0 ; generation < 100 ; generation++)
    {        
        cout << "Generation " << generation << ": " << endl;
        printPopulation(*population);

        if (!(r = sorter.sortPopulation(*population)))
            return "Error in sort: " + r.getErrorString();

        cout << "Generation " << generation << " (sorted): " << endl;
        printPopulation(*population);

        newPopulation->m_individuals.clear();
        for (size_t i = 0 ; i < popSize ; i++)
        {
            if (!(r = selection.selectParents(*population, parents)))
                return "Error in parent selection: " + r.getErrorString();

            if (generation == 0)
            {
                if (!(r = crossOver.check(parents))) // TODO: only for first generation
                    return "Error in crossover check: " + r.getErrorString();
            }

            if (!(r = crossOver.generateOffspring(parents, offspring)))
                return "Error generating offspring: " + r.getErrorString();

            for (auto &g : offspring)
            {
                auto f = refFitness->createCopy(false);
                newPopulation->m_individuals.push_back(make_shared<Individual>(g, f));
            }
        }

        std::swap(newPopulation, population);
        if (!(r = m_fitnessCalc->calculatePopulationFitness({population})))
            return "Error calculating fitness: " + r.getErrorString();

        // TODO: somehow check that all fitness values have been calculated
    }

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
    shared_ptr<RandomNumberGenerator> rng = make_shared<MersenneRandomNumberGenerator>(rd());

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

