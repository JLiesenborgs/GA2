#include "permutationordercrossover.h"
#include "vectorgenomefitness.h"

using namespace std;
using namespace errut;

namespace mogal2
{

PermutationOrderCrossover::PermutationOrderCrossover(const shared_ptr<RandomNumberGenerator> &rng, bool twoOffspring)
 : m_rng(rng), m_twoOffspring(twoOffspring)
{
    
}

bool_t PermutationOrderCrossover::check(const vector<shared_ptr<Genome>> &parents)
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

bool_t PermutationOrderCrossover::generateOffspring(const vector<shared_ptr<Genome>> &parents,
                            vector<shared_ptr<Genome>> &offspring)
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
}

}