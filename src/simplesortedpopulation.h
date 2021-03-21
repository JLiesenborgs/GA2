#pragma once

#include "mogal2config.h"
#include "ndsortedpopulation.h"
#include <cassert>

namespace mogal2
{

class SimpleSortedPopulation : public NDSortedPopulation
{
public:
    // TODO: add some kind of pruning operator?
    SimpleSortedPopulation(std::shared_ptr<FitnessComparison> fitComp, size_t objectiveNumber = 0);
    ~SimpleSortedPopulation();

    void setObjectiveNumber(size_t objectiveNumber) { m_objectiveNumber = objectiveNumber; }

    errut::bool_t check(const Population &population) override;
    errut::bool_t processPopulation(const std::shared_ptr<Population> &population, size_t targetPopulationSize) override;
    std::shared_ptr<Population> getSortedPopulation() const { return m_lastPopulation; }

    int getNumberOfSets() const override { return (int)m_lastPopulation->size(); }
    int getSetSize(int s) const override { return 1; }
    std::shared_ptr<Individual> getIndividual(int s, int i) const override
    {
        assert(m_lastPopulation.get());
        assert(s >= 0 && s < (int)m_lastPopulation->size());
        return m_lastPopulation->individual(s);
    };

    const std::vector<std::shared_ptr<Individual>> &getBestIndividuals() const { return m_best; }
private:
    size_t m_objectiveNumber;
    std::shared_ptr<FitnessComparison> m_fitnessComp;
    std::shared_ptr<Population> m_lastPopulation;
    std::vector<std::shared_ptr<Individual>> m_best;
};

}
