#include "trackbestonlyselectionpopulation.h"
#include <cassert>

using namespace std;
using namespace errut;

namespace mogal2
{

TrackBestOnlySelectionPopulation::TrackBestOnlySelectionPopulation(std::shared_ptr<FitnessComparison> fitComp, size_t objectiveNumber)
    : m_fitnessComp(fitComp), m_objectiveNumber(objectiveNumber)
{
}

TrackBestOnlySelectionPopulation::~TrackBestOnlySelectionPopulation()
{
}

bool_t TrackBestOnlySelectionPopulation::check(const Population &population)
{
    return true;
}

bool_t TrackBestOnlySelectionPopulation::processPopulation(std::shared_ptr<Population> &population, size_t targetPopulationSize)
{
    assert(population->size() > 0);

    if (m_best.size() == 0)
        m_best.push_back(population->individual(0)->createCopy());
    
    for (auto &i : population->individuals())
    {
        if (m_fitnessComp->isFitterThan(i->fitnessRef(), m_best[0]->fitnessRef(), m_objectiveNumber))
            m_best[0] = i->createCopy();
    }
    return true;
}

}