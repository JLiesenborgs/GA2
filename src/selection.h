#pragma once

#include "eatkconfig.h"
#include "population.h"

namespace eatk
{

// Something to give ParentSelection as input, e.g. a simple sorted population,
// a non-dominated sorted population
class SelectionPopulation
{
public:
	SelectionPopulation() { }
	virtual ~SelectionPopulation() { }

	virtual errut::bool_t check(const Population &population) { return "Not implemented in base class"; }
	
	// May change population! (e.g sort it immediately)
	// target population size is to allow pruning
	virtual errut::bool_t processPopulation(const std::shared_ptr<Population> &population, size_t targetPopulationSize) { return "Not implemented in base class"; }

	virtual const std::vector<std::shared_ptr<Individual>> &getBestIndividuals() const { return m_emptyBest; }
private:
	const std::vector<std::shared_ptr<Individual>> m_emptyBest;
};

class ParentSelection // Don't think we need some kind of check, shouldn't depend on type of genome
{
public:
	ParentSelection() { }
	virtual ~ParentSelection() { }
	
	virtual errut::bool_t check(const SelectionPopulation &pop) { return "Not implemented in base class"; }
	// The length of 'parents' describes the number of parents that should be
	// TODO: just a single parent? Letting the amount and possibly inbreeding be
	//	   controlled outside?
	// Note that population might have been changed (e.g. the order) by the SelectionPopulation
	// step that precedes it
	virtual errut::bool_t selectParents(const Population &population, const SelectionPopulation &selPop, std::vector<std::shared_ptr<Individual>> &parents) { return "Not implemented in base class"; }
};

}
