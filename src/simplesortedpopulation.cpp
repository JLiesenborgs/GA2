#include "simplesortedpopulation.h"
#include <algorithm>

using namespace std;
using namespace errut;

namespace mogal2
{

SimpleSortedPopulation::SimpleSortedPopulation(shared_ptr<FitnessComparison> fitComp, size_t objectiveNumber)
    : m_objectiveNumber(objectiveNumber), m_fitnessComp(fitComp) 
{
}

SimpleSortedPopulation::~SimpleSortedPopulation()
{
}

bool_t SimpleSortedPopulation::check(const Population &population)
{
    FitnessComparison &cmp = *m_fitnessComp;
    for (auto &i : population.individuals())
    {
        bool_t r = cmp.check(i->fitnessRef());
        if (!r)
            return "Error in fitness comparison check: " + r.getErrorString();
    }
    return true; 
}

bool_t SimpleSortedPopulation::processPopulation(shared_ptr<Population> &population, size_t targetPopulationSize)
{
    m_lastPopulation = population;

    FitnessComparison &cmp = *m_fitnessComp;
    const size_t N = m_objectiveNumber;

    auto comp = [&cmp, N](auto &i1, auto &i2)
    {
        return cmp.isFitterThan(i1->fitnessRef(), i2->fitnessRef(), N);
    };

    sort(population->individuals().begin(), population->individuals().end(), comp);

    // TODO: allow more complex pruning?
    population->resize(targetPopulationSize);

    assert(population->size() > 0);
    auto &i = population->individuals().front();
    if (m_best.size() > 0)
    {
        auto &best = m_best[0];
        if (cmp.isFitterThan(i->fitnessRef(), best->fitnessRef(), N))
            m_best[0] = i->createCopy(); // TODO: is it safe to not make a copy?
    }
    else
        m_best.push_back(i->createCopy()); // TODO: is it safe to not make a copy?

    return true;
}

}
