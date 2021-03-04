#include "rankparentselection.h"
#include "ndsortedpopulation.h"
#include <cassert>
#include <cmath>

using namespace std;
using namespace errut;

RankParentSelection::RankParentSelection(double beta, shared_ptr<RandomNumberGenerator> rng)
    : m_beta(beta), m_rng(rng)
{
}

RankParentSelection::~RankParentSelection()
{
}

bool_t RankParentSelection::check(const SelectionPopulation &pop)
{
    const NDSortedPopulation *pPop = dynamic_cast<const NDSortedPopulation *>(&pop);
    if (!pPop)
        return "Expecting an ND sorted population";

    return true;
}

bool_t RankParentSelection::selectParents(const SelectionPopulation &pop, std::vector<std::shared_ptr<Genome>> &parents)
{
    const NDSortedPopulation &ndPop = dynamic_cast<const NDSortedPopulation &>(pop);
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

    int idx1 = getSetIndex();
    int idx2 = getSetIndex();
    int genomeIdx1 = getGenomeIndex(idx1);
    int genomeIdx2 = getGenomeIndex(idx2);

    parents.resize(2);
    parents[0] = ndPop.getGenome(idx1, genomeIdx1);
    parents[1] = ndPop.getGenome(idx1, genomeIdx2);

    return true;
}
