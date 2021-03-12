#include "geneticalgorithm.h"
#include "mersennerandomnumbergenerator.h"
#include "valuefitness.h"
#include "vectorgenomefitness.h"
#include "singlethreadedpopulationcrossover.h"
#include "rankparentselection.h"
#include "simplesortedpopulation.h"
#include "singlebestelitism.h"
#include "remainingtargetpopulationsizeiteration.h"
#include "singlethreadedpopulationfitnesscalculation.h"
#include <iostream>

using namespace errut;
using namespace std;
using namespace mogal2;

class PermutationCrossover : public GenomeCrossover
{
public:
    PermutationCrossover(shared_ptr<RandomNumberGenerator> rng) : m_rng(rng) { }

	bool_t check(const std::vector<std::shared_ptr<Genome>> &parents) override
    {
        // TODO
        return true;
    }

	bool_t generateOffspring(const std::vector<std::shared_ptr<Genome>> &parents,
	                         std::vector<std::shared_ptr<Genome>> &generatedOffspring) override
    {
        // TODO
        return "TODO";
    };

private:
    shared_ptr<RandomNumberGenerator> m_rng;
};

class PermutationMutation : public GenomeMutation
{
public:
    PermutationMutation(const shared_ptr<RandomNumberGenerator> &rng, double prob) : m_rng(rng), m_prob(prob) { }

	bool_t check(const Genome &genome)
    {
        VectorGenome<int> *vg = dynamic_cast<VectorGenome<int>*>(genome);
        if (!vg)
            return "Wrong genome type";
        return true;
    }

	bool_t mutate(Genome &genome, bool &isChanged)
    {
        VectorGenome<int> &vg = static_cast<VectorGenome<int>&>(genome);
        vector<int> &values = vg.getValues();

        for (size_t i = 0 ; i < values.size() ; i++)
        {
            if (m_rng->getRandomDouble() < m_prob)
            {
                uint32_t idx1 = m_rng->getRandomUint32()%(uint32_t)values.size();
                uint32_t idx2 = m_rng->getRandomUint32()%(uint32_t)values.size();
                int x = values[idx1];
                int y = values[idx2];

                values[idx1] = y;
                values[idx2] = x;
                isChanged = true;
            }
        }
        return true;
    }
private:
    shared_ptr<RandomNumberGenerator> m_rng;
    double m_prob;
};

class Creation : public GenomeFitnessCreation
{
public:
    Creation(const shared_ptr<RandomNumberGenerator> rng, size_t numCities)
     : m_rng(rng), m_numCities(numCities) { }

    shared_ptr<Genome> createInitializedGenome() override
    {
        auto g = make_shared<VectorGenome<int>>(m_numCities);
        vector<int> indices(m_numCities);
        for (size_t i = 0 ; i < m_numCities ; i++)
            indices[i] = i;
        
        for (size_t i = 0 ; i < m_numCities ; i++)
        {
            uint32_t idx = m_rng->getRandomUint32() % (uint32_t)indices.size();
            g->getValues()[i] = indices[idx];

            indices[idx] = indices[indices.size()-1];
            indices.resize(indices.size()-1);
        }
        return g;
    }
    
    shared_ptr<Fitness> createEmptyFitness() override
    {
        return make_shared<ValueFitness<double>>();
    }
private:
    shared_ptr<RandomNumberGenerator> m_rng;
    size_t m_numCities;
};

class City
{
public:
    City(double X, double Y) : x(X), y(Y) { }
    double x, y;
};

class MyStop : public FixedGenerationsStopCriterion
{
public:
    MyStop(size_t n) : FixedGenerationsStopCriterion(n) { }
    bool_t analyze(const vector<shared_ptr<Individual>> &currentBest, size_t generationNumber, bool &shouldStop) override
    {
        if (currentBest.size() > 0)
        {
            assert(currentBest.size() == 1);
            cout << generationNumber << "| " << currentBest[0]->toString() << endl;
        }
        return FixedGenerationsStopCriterion::analyze(currentBest, generationNumber, shouldStop);
    }
};

class MyGA : public GeneticAlgorithm
{
protected:
    bool_t onAlgorithmDone(size_t generation, const std::vector<std::shared_ptr<Individual>> &bestIndividuals) override
    {
        cout << "Ending after " << generation << " generations, best are: " << endl;
        for (auto &i : bestIndividuals)
            cout << i->toString() << endl;
        return true;
    }
};

class TSPFitnessCalculation : public GenomeFitnessCalculation
{
public:
    TSPFitnessCalculation(const vector<City> &cities) : m_cities(cities) { }
    bool_t calculate(const Genome &genome, Fitness &fitness) override
    {
        const VectorGenome<int> &vg = static_cast<const VectorGenome<int> &>(genome);
        const vector<int> &indices = vg.getValues();
        ValueFitness<double> &vf = static_cast<ValueFitness<double> &>(fitness);

        double totalDist = 0;
        assert(indices.size() == m_cities.size());
        for (size_t i = 0 ; i < m_cities.size() ; i++)
        {
            size_t j = (i+1)%m_cities.size();
            City c1 = m_cities[indices[i]];
            City c2 = m_cities[indices[j]];
            double dx = c1.x - c2.x;
            double dy = c1.y - c2.y;
            totalDist += sqrt(dx*dx + dy*dy);
        }
        vf.setValue(totalDist);
        return true;
    }
private:
    vector<City> m_cities;
};

int main(int argc, char const *argv[])
{
    size_t numCities = 28;
    vector<City> cities;
    for (size_t i = 0 ; i < numCities ; i++)
    {
        double theta = (double)i/(double)numCities * 2.0*M_PI;
        cities.push_back(City { cos(theta), sin(theta) } );
    }

    random_device rd;
    unsigned int seed = rd();
    auto rng = make_shared<MersenneRandomNumberGenerator>(seed);
    auto creation = make_shared<Creation>(rng, cities.size());

    auto calcSingle = make_shared<SingleThreadedPopulationFitnessCalculation>(make_shared<TSPFitnessCalculation>(cities));
    auto mutation = make_shared<PermutationMutation>(rng, 0.5/cities.size());
    auto cross = make_shared<SingleThreadedPopulationCrossover>(1.0, false,
            make_shared<SimpleSortedPopulation>(make_shared<ValueFitnessComparison<double>>()),
            make_shared<RankParentSelection>(2.5, rng),
            make_shared<PermutationCrossover>(rng),
            mutation,
            make_shared<SingleBestElitism>(true, mutation),
            make_shared<RemainingTargetPopulationSizeIteration>(),
            rng
        );
    
    MyStop stop(1000);
    MyGA ga;

    auto r = ga.run(*creation,
                    *cross,
                    *calcSingle, stop, 256); //, 0, 34);
    if (!r)
    {
        cerr << "Error: " << r.getErrorString() << endl;
        return -1;
    }
    cout << "Done" << endl;
    return 0;
}
