#pragma once

#include "eatkconfig.h"
#include "crossovermutation.h"
#include "randomnumbergenerator.h"

namespace eatk
{

class SingleThreadedPopulationCrossover : public PopulationCrossover
{
public:
    SingleThreadedPopulationCrossover(double cloneFraction, bool keepExistingPopulation,
                                      const std::shared_ptr<SelectionPopulation> &selectionPop,
                                      const std::shared_ptr<ParentSelection> &parentSelection,
                                      const std::shared_ptr<GenomeCrossover> &genomeCrossover,
                                      const std::shared_ptr<GenomeMutation> &genomeMutation,
                                      const std::shared_ptr<Elitism> &elitism,
                                      const std::shared_ptr<PopulationCrossoverIteration> &popIteration,
                                      const std::shared_ptr<RandomNumberGenerator> &rng);
    ~SingleThreadedPopulationCrossover();
    
    errut::bool_t check(const std::vector<std::shared_ptr<Population>> &populations) override;
	errut::bool_t createNewPopulation(size_t generation, std::vector<std::shared_ptr<Population>> &populations, size_t targetPopulationSize) override;

    const std::vector<std::shared_ptr<Individual>> &getBestIndividuals() const override
    {
        // Here, the SelectionPopulation is the one that inspects the population after
        // The fitnesses have all been calculated, so it should be the one to track
        // the best
        return m_selectionPop->getBestIndividuals();
    }
protected:
    virtual errut::bool_t onSelectionPopulationProcessed(size_t generation, const std::shared_ptr<SelectionPopulation> &selPop) { return true; }
private:
    double m_cloneFraction;
    bool m_keepExistingPopulation;
    std::shared_ptr<SelectionPopulation> m_selectionPop;
    std::shared_ptr<ParentSelection> m_parentSelection;
    std::shared_ptr<GenomeCrossover> m_genomeCrossover;
    std::shared_ptr<GenomeMutation> m_genomeMutation;
    std::shared_ptr<Elitism> m_elitism;
    std::shared_ptr<PopulationCrossoverIteration> m_popIteration;
    std::shared_ptr<RandomNumberGenerator> m_rng;
};

}