#pragma once

#include "crossovermutation.h"
#include "randomnumbergenerator.h"

class SingleThreadedPopulationCrossover : public PopulationCrossover
{
public:
    SingleThreadedPopulationCrossover(double cloneFraction, bool keepExistingPopulation,
                                      std::shared_ptr<SelectionPopulation> selectionPop,
                                      std::shared_ptr<ParentSelection> parentSelection,
                                      std::shared_ptr<GenomeCrossover> genomeCrossover,
                                      std::shared_ptr<GenomeMutation> genomeMutation,
                                      std::shared_ptr<Elitism> elitism,
                                      std::shared_ptr<PopulationCrossoverIteration> popIteration,
                                      std::shared_ptr<RandomNumberGenerator> rng);
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
    virtual errut::bool_t onSelectionPopulationProcessed(size_t generation, std::shared_ptr<SelectionPopulation> &selPop) { return true; }
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