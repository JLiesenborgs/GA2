#pragma once

#include "ndsortedpopulation.h"
#include <cassert>

class SimpleSortedPopulation : public NDSortedPopulation
{
public:
    // TODO: add some kind of pruning operator?
    SimpleSortedPopulation(std::shared_ptr<FitnessComparison> fitComp, int objectiveNumber = 0);
    ~SimpleSortedPopulation();

    void setObjectiveNumber(int objectiveNumber) { m_objectiveNumber = objectiveNumber; }

    errut::bool_t check(const Population &population) override;
    errut::bool_t processPopulation(std::shared_ptr<Population> &population, size_t targetPopulationSize) override;
    std::shared_ptr<Population> getSortedPopulation() const { return m_lastPopulation; }

    int getNumberOfSets() const override { return (int)m_lastPopulation->m_individuals.size(); }
    int getSetSize(int s) const override { return 1; }
    std::shared_ptr<Genome> getGenome(int s, int i) const override
    {
        assert(m_lastPopulation.get());
        assert(s >= 0 && s < (int)m_lastPopulation->m_individuals.size());
        return m_lastPopulation->m_individuals[s]->m_genome;
    };

    const std::vector<std::shared_ptr<Individual>> &getBestIndividuals() const { return m_best; }
private:
    int m_objectiveNumber;
    std::shared_ptr<FitnessComparison> m_fitnessComp;
    std::shared_ptr<Population> m_lastPopulation;
    std::vector<std::shared_ptr<Individual>> m_best;
};
