#include "population.h"
#include "mersennerandomnumbergenerator.h"
#include "vectorgenomefitness.h"
#include <iostream>
#include <algorithm>

using namespace std;
using namespace errut;
using namespace mogal2;

class BasicNDSetCreator
{
public:
    BasicNDSetCreator(const shared_ptr<FitnessComparison> &fitCmp, size_t numObjectives)
     : m_cmp(fitCmp), m_numObjectives(numObjectives) { }
    ~BasicNDSetCreator() { }

    bool_t calculateAllNDSets(const vector<shared_ptr<Individual>> &individuals);
    bool_t calculateNonDomitatedSet(const vector<shared_ptr<Individual>> &individuals,
        vector<shared_ptr<Individual>> &ndSet,
        vector<shared_ptr<Individual>> &remaining);

    size_t getNumberOfSets() const { return m_sets.size(); }
    const vector<shared_ptr<Individual>> &getSet(size_t idx) const { return m_sets[idx]; }
private:
    size_t m_numObjectives;
    shared_ptr<FitnessComparison> m_cmp;
    
    vector<vector<shared_ptr<Individual>>> m_sets;
    vector<shared_ptr<Individual>> m_tmpND, m_tmpRem[2];
};

bool_t BasicNDSetCreator::calculateAllNDSets(const vector<shared_ptr<Individual>> &individuals)
{
    m_sets.clear();

    const vector<shared_ptr<Individual>> *pIn = &individuals;
    vector<shared_ptr<Individual>> *pND = &m_tmpND;
    vector<shared_ptr<Individual>> *pRem = &m_tmpRem[0];
    size_t rIdx = 0;
    bool_t r;

    while (true)
    {
        if (pIn->size() == 0)
            break;

        pND->clear();
        pRem->clear();

        if (!(r = calculateNonDomitatedSet(*pIn, *pND, *pRem)))
            return "Error calulating one of the sets: " + r.getErrorString();

        size_t newIdx = m_sets.size();
        m_sets.resize(newIdx+1);
        swap(m_sets[newIdx], *pND);

        pIn = pRem;
        rIdx++;
        pRem = &m_tmpRem[rIdx%2];
    }
    return true;
}

bool_t BasicNDSetCreator::calculateNonDomitatedSet(const vector<shared_ptr<Individual>> &individuals,
                                                   vector<shared_ptr<Individual>> &ndSet,
                                                   vector<shared_ptr<Individual>> &remaining)
{
    ndSet.clear();
    remaining.clear();

    const size_t N = m_numObjectives;
    FitnessComparison &cmp = *m_cmp;

    auto isDominated = [&cmp, N](auto &i, auto &j)
    {
        size_t betterCount = 0;
        size_t betterEqualCount = 0;

        for (size_t k = 0 ; k < N ; k++)
        {
            const Fitness &fi = i->fitnessRef();
            const Fitness &fj = j->fitnessRef();
            if (cmp.isFitterThan(fj, fi, k))
            {
                betterCount++;
                betterEqualCount++;
            }
            else
            {
                if (!cmp.isFitterThan(fi, fj, k))
                    betterEqualCount++;
            }
        }

        if (betterEqualCount == N && betterCount > 0)
            return true;
        return false;
    };

    for (auto &i : individuals)
    {
        bool found = false;

        for (auto &j : individuals)
        {
            if (i.get() == j.get())
                continue;

            if (isDominated(i, j))
            {
                found = true;
                break;
            }
        }

        if (!found)
            ndSet.push_back(i);
        else
            remaining.push_back(i);
    }
    return true;
}

int main(int argc, char const *argv[])
{
    random_device rd;
    MersenneRandomNumberGenerator rng(rd());

    Population pop;
    size_t N = 50;
    size_t dim = 2;
    
    for (size_t i = 0 ; i < N ; i++)
    {
        VectorGenome<float> dummy;
        VectorFitness<float> f { dim };
        f.setCalculated();

        for (size_t j = 0 ; j < dim ; j++)
            f.getValues()[j] = (int)(rng.getRandomDouble()*100);

        pop.append(make_shared<Individual>(dummy.createCopy(), f.createCopy()));
    }

    vector<shared_ptr<Individual>> ndSet;
    vector<shared_ptr<Individual>> remainder;

    BasicNDSetCreator bndsc(make_shared<VectorFitnessComparison<float>>(), dim);
    bool_t r;
    
    auto printSet = [dim](auto &v)
    {
        for (auto &i : v)
        {
            const VectorFitness<float> &f = dynamic_cast<const VectorFitness<float> &>(i->fitnessRef());
            for (size_t k = 0 ; k < dim ; k++)
                cout << f.getValues()[k] << " ";
            cout << endl;
        };
    };

    // if (!(r = bndsc.calculateNonDomitatedSet(pop.individuals(), ndSet, remainder)))
    //     cerr << "Error: " << r.getErrorString() << endl;

    printSet(pop.individuals());
    cout << endl << endl;
    // printSet(ndSet);
    // cout << endl << endl;
    // printSet(remainder);
    // cout << endl << endl;

    if (!(r = bndsc.calculateAllNDSets(pop.individuals())))
    {
        cerr << "Error: " << r.getErrorString() << endl;
        return -1;
    }

    cout << "# Found " << bndsc.getNumberOfSets() <<  " sets" << endl;
    for (size_t i = 0 ; i < bndsc.getNumberOfSets() ; i++)
    {
        auto s = bndsc.getSet(i);
        sort(s.begin(), s.end(), [](auto &e1, auto &e2) {
            const VectorFitness<float> &f1 = dynamic_cast<const VectorFitness<float> &>(e1->fitnessRef());
            const VectorFitness<float> &f2 = dynamic_cast<const VectorFitness<float> &>(e2->fitnessRef());
            return f1.getValues()[0] < f2.getValues()[0];
        });
        printSet(s);
        cout << endl << endl;
    }
    return 0;
}
