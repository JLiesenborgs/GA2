#pragma once

#include "crossovermutation.h"

namespace eatk
{

class TrackBestOnlySelectionPopulation : public SelectionPopulation
{
public:
	TrackBestOnlySelectionPopulation(std::shared_ptr<FitnessComparison> fitComp, size_t objectiveNumber = 0);
	~TrackBestOnlySelectionPopulation();

	errut::bool_t check(const Population &population) override;
	errut::bool_t processPopulation(const std::shared_ptr<Population> &population, size_t targetPopulationSize) override;
	const std::vector<std::shared_ptr<Individual>> &getBestIndividuals() const override { return m_best; }
private:
	size_t m_objectiveNumber;
	std::shared_ptr<FitnessComparison> m_fitnessComp;
	std::vector<std::shared_ptr<Individual>> m_best;
};

}
