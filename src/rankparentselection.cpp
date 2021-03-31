#include "rankparentselection.h"
#include "ndsortedpopulationinterface.h"
#include <cassert>
#include <cmath>

using namespace std;
using namespace errut;

namespace eatk
{

RankParentSelection::RankParentSelection(double beta, shared_ptr<RandomNumberGenerator> rng)
    : m_beta(beta), m_rng(rng)
{
}

RankParentSelection::~RankParentSelection()
{
}

bool_t RankParentSelection::check(const SelectionPopulation &pop)
{
    const NDSortedPopulationInterface *pPop = dynamic_cast<const NDSortedPopulationInterface *>(&pop);
    if (!pPop)
        return "Expecting an ND sorted population";

    return true;
}

bool_t RankParentSelection::selectParents(const Population &population,
                                          const SelectionPopulation &selPop,
                                          std::vector<std::shared_ptr<Individual>> &parents)
{
    const NDSortedPopulationInterface &ndPop = static_cast<const NDSortedPopulationInterface &>(selPop);
    int numSets = ndPop.getNumberOfSets();
    assert(numSets > 0);
    
    auto getSetIndex = [this, numSets]()
    {
        double x = m_rng->getRandomDouble();
        double y = 1.0-std::pow(x, 1.0/(1.0+m_beta));
        int idx = (int)(y*(double)numSets);
        if (idx >= numSets)
            idx = numSets-1; // Population size is at least one
        return idx;
    };

    auto getGenomeIndex = [this, &ndPop](int setIdx)
    {
        int setSize = ndPop.getSetSize(setIdx);
        assert(setSize > 0);
        if (setSize <= 1)
            return 0;
        
        int idx = (int)(m_rng->getRandomDouble()*(double)setSize);
        if (idx >= setSize)
            idx = setSize-1; // Set size should be at least one
        return idx;
    };

    for (auto &p : parents)
    {
        int setIdx = getSetIndex();
        int genomeIdx = getGenomeIndex(setIdx);

        p = ndPop.getIndividual(setIdx, genomeIdx);
    }

    // TODO: something to prevent inbreeding?

    return true;
}

}
