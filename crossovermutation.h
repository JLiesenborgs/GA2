#pragma once

#include "genomefitness.h"
#include "population.h"
#include <vector>

class GenomeCrossover
{
public:
	GenomeCrossover(size_t numParents = 2) : m_numParents(numParents) { }
	virtual ~GenomeCrossover() { }

	void setNumberOfParents(size_t n) { m_numParents = n; }
	size_t getNumberOfParents() const { return m_numParents; }
	// This function is intended to check number and type of parents, so that
	// e.g. no dynamic cast is needed in generateOffspring itself, the number
	// of parents and their type is assumed to be checked
	virtual errut::bool_t check(const std::vector<std::shared_ptr<Genome>> &parents) { return "Not implemented in base class"; }

	virtual errut::bool_t generateOffspring(const std::vector<std::shared_ptr<Genome>> &parents,
	                                        std::vector<std::shared_ptr<Genome>> &generatedOffspring) { return "Not implemented in base class"; }
private:
	int m_numParents;
};

class GenomeMutation
{
public:
	GenomeMutation() { }
	virtual ~GenomeMutation() { }
	// This function is intended to check type of genome, so that
	// e.g. no dynamic cast is needed in next function call
	virtual errut::bool_t check(const Genome &genome) { return "Not implemented in base class"; }
	virtual errut::bool_t mutate(Genome &genome, bool &isChanged) { return "Not implemented in base class"; }
};

// Idea is to have e.g a single threaded one and a multithreaded one
// Should we allow multithreaded? What to do with random number generation?
//  -> Mersenne twister has 'discard' in c++11
class PopulationMutation
{
public:
	PopulationMutation() { m_tmp.resize(1); }
	virtual ~PopulationMutation() { }

	// Convenience functions
	virtual errut::bool_t check(const std::shared_ptr<Population> &population)
	{
		m_tmp[0] = population;
		return check(m_tmp);
	}
	virtual errut::bool_t mutate(std::shared_ptr<Population> &population)
	{
		m_tmp[0] = population;
		return mutate(m_tmp);
	}

	virtual errut::bool_t check(const std::vector<std::shared_ptr<Population>> &populations) { return "Not implemented in base class"; }
	// This is in-place, must reset isCalculated flag if changed
	// Idea is to apply a GenomeMutation operator to every individual
	virtual errut::bool_t mutate(const std::vector<std::shared_ptr<Population>> &populations) { return "Not implemented in base class"; }
private:
	std::vector<std::shared_ptr<Population>> m_tmp;
};

// TODO: PopulationCrossover ?
// allow in-place? some children that overwrite older parents?
class PopulationCrossover
{
public:
	PopulationCrossover() { m_tmp.resize(1); }
	virtual ~PopulationCrossover() { }

	virtual errut::bool_t check(const std::shared_ptr<Population> &population)
	{
		m_tmp[0] = population;
	 	return check(m_tmp);
	}
	virtual errut::bool_t check(const std::vector<std::shared_ptr<Population>> &populations){ return "Not implemented in base class"; }
	// The populations are overwritten, if the old one is still needed it should
	// be stored externally
	// This is in-place, must reset isCalculated flag if changed
	// Idea is to apply a GenomeMutation operator to every individual
	virtual errut::bool_t createNewPopulation(std::vector<std::shared_ptr<Population>> &populations, int targetPopulationSize) { return "Not implemented in base class"; }
	virtual errut::bool_t createNewPopulation(std::shared_ptr<Population> &population, int targetPopulationSize)
	{
		m_tmp[0] = population;
		auto r = createNewPopulation(m_tmp, targetPopulationSize);
		std::swap(m_tmp[0], population);
		return r;
	}
private:
	std::vector<std::shared_ptr<Population>> m_tmp;
};

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
	virtual errut::bool_t processPopulation(std::shared_ptr<Population> &population, int targetPopulationSize) { return "Not implemented in base class"; }
};

class ParentSelection // Don't think we need some kind of check, shouldn't depend on type of genome
{
public:
	ParentSelection() { }
	virtual ~ParentSelection() { }
	
	virtual errut::bool_t check(const SelectionPopulation &pop) { return "Not implemented in base class"; }
	// The length of 'parents' describes the number of parents that should be
	virtual errut::bool_t selectParents(const SelectionPopulation &pop, std::vector<std::shared_ptr<Genome>> &parents) { return "Not implemented in base class"; }
};

