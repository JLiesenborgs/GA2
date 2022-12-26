#pragma once

#include "eatkconfig.h"
#include "population.h"

namespace eatk
{

// allow in-place? some children that overwrite older parents?
class PopulationEvolver
{
public:
	PopulationEvolver() { }
	virtual ~PopulationEvolver() { }

	virtual errut::bool_t check(const std::shared_ptr<Population> &population) { return "Not implemented in base class"; }
	virtual errut::bool_t check(const std::vector<std::shared_ptr<Population>> &populations) { return "Not implemented in base class"; }
	// The populations are overwritten, if the old one is still needed it should
	// be stored externally
	virtual errut::bool_t createNewPopulations(size_t generation, std::vector<std::shared_ptr<Population>> &populations, const std::vector<size_t> &targetPopulationSize) { return "Not implemented in base class"; }
	virtual errut::bool_t createNewPopulation(size_t generation, std::shared_ptr<Population> &population, size_t targetPopulationSize) { return "Not implemented in base class"; }

	// This is the main class that creates new populations/individuals, this
	// should also be the access point to keep track of the best ones
	virtual const std::vector<std::shared_ptr<Individual>> &getBestIndividuals() const { return m_emptyBest; }
private:
	const std::vector<std::shared_ptr<Individual>> m_emptyBest;
};

} 
