#include "permutationswapmutation.h"
#include "vectorgenomefitness.h"

using namespace std;
using namespace errut;

namespace mogal2
{

PermutationSwapMutation::PermutationSwapMutation(const shared_ptr<RandomNumberGenerator> &rng, double prob)
    : m_rng(rng), m_prob(prob)
{
}

bool_t PermutationSwapMutation::check(const Genome &genome)
{
    auto vg = dynamic_cast<const VectorGenome<int>*>(&genome);
    if (!vg)
        return "Wrong genome type";
    return true;
}

bool_t PermutationSwapMutation::mutate(Genome &genome, bool &isChanged)
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

}