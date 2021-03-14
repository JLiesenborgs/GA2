#include "tournamentparentselection.h"
#include "simplesortedpopulation.h"
#include <algorithm>

using namespace std;
using namespace errut;

namespace mogal2
{

TournamentParentSelection::TournamentParentSelection(const shared_ptr<RandomNumberGenerator> &rng,
                                                     size_t tournamentSize,
                                                     const std::shared_ptr<FitnessComparison> &fitCmp,
                                                     size_t objectiveNumber)
    : m_rng(rng), m_cmp(fitCmp), m_objectiveNumber(objectiveNumber)
{
    m_tournament.resize(tournamentSize);
}

TournamentParentSelection::~TournamentParentSelection()
{
}

bool_t TournamentParentSelection::check(const SelectionPopulation &pop)
{
    const SimpleSortedPopulation *pPop = dynamic_cast<const SimpleSortedPopulation *>(&pop);
    if (!pPop)
        return "Selection preprocessor must be a simple sorted population";

    if (m_tournament.size() < 2)
        return "Tournament size must be at least two";

    return true;
}

bool_t TournamentParentSelection::selectParents(const SelectionPopulation &selPop, vector<shared_ptr<Individual>> &parents)
{
    const SimpleSortedPopulation &p = static_cast<const SimpleSortedPopulation&>(selPop);
    auto &pop = *p.getSortedPopulation();

    FitnessComparison &cmp = *m_cmp;
    const size_t N = m_objectiveNumber;

    auto comp = [&cmp, N](auto &i1, auto &i2)
    {
        return cmp.isFitterThan(i1->fitnessRef(), i2->fitnessRef(), N);
    };

    auto pickParent = [this, &pop, &comp]()
    {
        // Fill tournament
        // TODO: avoid doubles?
        for (auto &i : m_tournament)
        {
            uint32_t idx = m_rng->getRandomUint32() % (uint32_t)pop.size();
            i = pop.individual(idx);
        }

        sort(m_tournament.begin(), m_tournament.end(), comp);
        return m_tournament[0];
    };

    for (auto &p : parents)
        p = pickParent();
    return true;
}

}