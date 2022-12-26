#pragma once

#include "eatkconfig.h"
#include "populationevolver.h"
#include "fasternondominatedsetcreator.h"

namespace eatk
{

class MultiPopulationEvolver : public PopulationEvolver
{
public:
	MultiPopulationEvolver(const std::vector<std::shared_ptr<PopulationEvolver>> &singlePopulationEvolvers,
	                       const std::shared_ptr<FitnessComparison> &fitCmp, size_t numObjectives);
	~MultiPopulationEvolver();

	errut::bool_t check(const std::vector<std::shared_ptr<Population>> &populations) override;
	errut::bool_t createNewPopulations(size_t generation, std::vector<std::shared_ptr<Population>> &populations, const std::vector<size_t> &targetPopulationSize) override;

	const std::vector<std::shared_ptr<Individual>> &getBestIndividuals() const override;
private:
	std::vector<std::shared_ptr<PopulationEvolver>> m_singlePopEvolvers;
	std::unique_ptr<FasterNonDominatedSetCreator> m_ndSetCreator;
	std::shared_ptr<FitnessComparison> m_fitCmp;
	size_t m_numObjectives;
	mutable std::vector<std::shared_ptr<Individual>> m_best;
};

}
