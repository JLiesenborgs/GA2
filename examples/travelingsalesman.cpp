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

class PermutationOrderCrossover : public GenomeCrossover
{
public:
    PermutationOrderCrossover(shared_ptr<RandomNumberGenerator> rng, bool twoOffspring) : m_rng(rng), m_twoOffspring(twoOffspring) { }

	bool_t check(const std::vector<std::shared_ptr<Genome>> &parents) override
    {
        if (parents.size() != 2)
            return "Need two parents";
        for (auto g : parents)
        {
            if (!dynamic_cast<VectorGenome<int>*>(g.get()))
                return "Genome is of wrong type";
        }
        return true;
    }

	bool_t generateOffspring(const std::vector<std::shared_ptr<Genome>> &parents,
	                         std::vector<std::shared_ptr<Genome>> &offspring) override
    {
        const VectorGenome<int> *pParents[2];
        for (int i = 0 ; i < 2 ; i++)
        {
            pParents[i] = static_cast<const VectorGenome<int>*>(parents[i].get());
            assert(pParents[i]);
        }

        uint32_t num = pParents[0]->getValues().size();
        uint32_t idx1 = m_rng->getRandomUint32()%num;
        uint32_t idx2 = m_rng->getRandomUint32()%num;
        if (idx1 > idx2)
            swap(idx1, idx2);

        m_hasIndex.resize(num);
        offspring.resize((m_twoOffspring)?2:1);
        for (auto &o0 : offspring)
        {
            auto &v1 = pParents[0]->getValues();
            auto &v2 = pParents[1]->getValues();

            o0 = pParents[0]->createCopy(false);
            auto o = static_cast<VectorGenome<int> *>(o0.get());
            auto &ov = o->getValues();

            for (size_t i = 0 ; i < m_hasIndex.size() ; i++)
                m_hasIndex[i] = false;

            for (uint32_t i = idx1 ; i <= idx2 ; i++)
            {
                ov[i] = v1[i];
                m_hasIndex[ov[i]] = true;
            }

            uint32_t next = idx2+1;
            for (uint32_t j = 0 ; j < v1.size() ; j++)
            {
                uint32_t i = (j+idx2+1)%v1.size(); // start from next index
                int val = v2[i];
                if (!m_hasIndex[val])
                {
                    ov[next%v1.size()] = val;
                    next++;
                }
            }

            // swap the role of the parents for the second child
            swap(pParents[0], pParents[1]);
        }

        // cout << "parent 1 " << pParents[0]->toString() << endl;
        // cout << "parent 2 " << pParents[1]->toString() << endl;
        // cout << "idx = " << idx1 << " " << idx2 << endl;
        // cout << "offspring 1" << offspring[0]->toString() << endl;
        // if (m_twoOffspring)
        //     cout << "offspring 2" << offspring[1]->toString() << endl;
        // cout << endl;
        // return "TODO";

        return true;
    };
private:
    shared_ptr<RandomNumberGenerator> m_rng;
    vector<bool> m_hasIndex;
    bool m_twoOffspring;
};

class PermutationSwapMutation : public GenomeMutation
{
public:
    PermutationSwapMutation(const shared_ptr<RandomNumberGenerator> &rng, double prob) : m_rng(rng), m_prob(prob) { }

	bool_t check(const Genome &genome)
    {
        auto vg = dynamic_cast<const VectorGenome<int>*>(&genome);
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
    auto mutation = make_shared<PermutationSwapMutation>(rng, 0.5/cities.size());
    auto cross = make_shared<SingleThreadedPopulationCrossover>(0.1, false,
            make_shared<SimpleSortedPopulation>(make_shared<ValueFitnessComparison<double>>()),
            make_shared<RankParentSelection>(2.5, rng),
            make_shared<PermutationOrderCrossover>(rng, true),
            mutation,
            make_shared<SingleBestElitism>(true, mutation),
            make_shared<RemainingTargetPopulationSizeIteration>(),
            rng
        );
    
    MyStop stop(1000);
    MyGA ga;

    auto r = ga.run(*creation,
                    *cross,
                    *calcSingle, stop, 128, 0, 258);
    if (!r)
    {
        cerr << "Error: " << r.getErrorString() << endl;
        return -1;
    }
    cout << "Done" << endl;
    return 0;
}
