#include "ndsortedpopulation.h"

using namespace errut;
using namespace std;

namespace mogal2
{

NDSortedPopulation::NDSortedPopulation(const std::shared_ptr<NonDominatedSetCreator> &setCreator,
                       const std::shared_ptr<DuplicateIndividualRemoval> &dupRemoval,
                       bool bestOnlyDuplicateRemoval)
    : m_ndSetCreator(setCreator), m_dupRemoval(dupRemoval), 
      m_removeFromBestOnly(bestOnlyDuplicateRemoval)
{
}

NDSortedPopulation::~NDSortedPopulation()
{
}

bool_t NDSortedPopulation::check(const Population &population)
{
    if (!m_ndSetCreator.get())
        return "No non-dominated set creator was set";

    if (m_dupRemoval.get())
    {
        bool_t r;
        if (!(r = m_dupRemoval->check(population.individuals())))
            return "Error checking duplicate removal: " + r.getErrorString();
    }
    return true;
}

bool_t NDSortedPopulation::processPopulation(const std::shared_ptr<Population> &population, size_t targetPopulationSize)
{
    bool_t r;

    assert(population->size() > 0);

    // TODO: can't do anything right now with target population size

    // Create the sets
    if (!(r = m_ndSetCreator->calculateAllNDSets(population->individuals())))
        return "Error calculating non-dominated sets: " + r.getErrorString();
    
    assert(m_ndSetCreator->getNumberOfSets() > 0);

    if (m_dupRemoval.get() && !m_removeFromBestOnly) // Remove doubles from all sets
    {
        for (size_t s = 0 ; s < m_ndSetCreator->getNumberOfSets() ; s++)
        {
            if (!(r = m_dupRemoval->removeDuplicates(m_ndSetCreator->getSet(s))))
                return "Error removing duplicats from set " + to_string(s) + ": " + r.getErrorString();
        }
    }

    // Update the best ones
    if (!(r = m_ndSetCreator->mergeNDSets(m_best, m_ndSetCreator->getSet(0))))
        return "Error merging new non-dominated front into best genome set: " + r.getErrorString();

    if (m_dupRemoval.get())
    {
        if (!(r = m_dupRemoval->removeDuplicates(m_best)))
            return "Error removing duplicates from best set: " + r.getErrorString();
    }

    return true;
}

}