#pragma once

#include "crossovermutation.h"
#include "randomnumbergenerator.h"

class SingleThreadedPopulationCrossover : public PopulationCrossover
{
public:
    SingleThreadedPopulationCrossover(double cloneFraction,
                                      std::shared_ptr<SelectionPopulation> selectionPop,
                                      std::shared_ptr<ParentSelection> parentSelection,
                                      std::shared_ptr<GenomeCrossover> genomeCrossover,
                                      std::shared_ptr<RandomNumberGenerator> rng);
    ~SingleThreadedPopulationCrossover();
    
    errut::bool_t check(const std::vector<std::shared_ptr<Population>> &populations) override;
	errut::bool_t createNewPopulation(std::vector<std::shared_ptr<Population>> &populations, int targetPopulationSize) override;
private:
    double m_cloneFraction;
    std::shared_ptr<SelectionPopulation> m_selectionPop;
    std::shared_ptr<ParentSelection> m_parentSelection;
    std::shared_ptr<GenomeCrossover> m_genomeCrossover;
    std::shared_ptr<RandomNumberGenerator> m_rng;
};