#include "simplesortedpopulation.h"
#include <algorithm>

using namespace std;
using namespace errut;

SimpleSortedPopulation::SimpleSortedPopulation(shared_ptr<FitnessComparison> fitComp, int objectiveNumber)
    : m_objectiveNumber(objectiveNumber), m_fitnessComp(fitComp) 
{
}

SimpleSortedPopulation::~SimpleSortedPopulation()
{
}

bool_t SimpleSortedPopulation::check(const Population &population)
{
    FitnessComparison &cmp = *m_fitnessComp;
    for (auto &i : population.m_individuals)
    {
        bool_t r = cmp.check(*i->m_fitness);
        if (!r)
            return "Error in fitness comparison check: " + r.getErrorString();
    }
    return true; 
}

bool_t SimpleSortedPopulation::processPopulation(shared_ptr<Population> &population)
{
    m_lastPopulation = population;

    FitnessComparison &cmp = *m_fitnessComp;
    const int N = m_objectiveNumber;
    auto comp = [&cmp, N](auto &i1, auto &i2)
    {
        return cmp.isFitterThan(*i1->m_fitness, *i2->m_fitness, N);
    };

    sort(population->m_individuals.begin(), population->m_individuals.end(), comp);
    return true;
}
