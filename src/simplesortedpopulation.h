#pragma once

#include "eatkconfig.h"
#include "ndsortedpopulationinterface.h"
#include <cassert>

namespace eatk
{

class SimpleSortedPopulation : public NDSortedPopulationInterface
{
public:
    // TODO: add some kind of pruning operator?
    SimpleSortedPopulation(std::shared_ptr<FitnessComparison> fitComp, size_t objectiveNumber = 0);
    ~SimpleSortedPopulation();

    void setObjectiveNumber(size_t objectiveNumber) { m_objectiveNumber = objectiveNumber; }

    errut::bool_t check(const Population &population) override;
    errut::bool_t processPopulation(const std::shared_ptr<Population> &population, size_t targetPopulationSize) override;
    std::shared_ptr<Population> getSortedPopulation() const { return m_lastPopulation; }

    size_t getNumberOfSets() const override { return m_lastPopulation->size(); }
    size_t getSetSize(size_t s) const override { return 1; }
    const std::shared_ptr<Individual> &getIndividual(size_t s, size_t i) const override
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
