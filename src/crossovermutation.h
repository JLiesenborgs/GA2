#pragma once

#include "mogal2config.h"
#include "genomefitness.h"
#include "population.h"
#include <vector>

namespace mogal2
{

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
	size_t m_numParents;
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

// allow in-place? some children that overwrite older parents?
class PopulationCrossover
{
public:
	PopulationCrossover() { m_tmp.resize(1); }
	virtual ~PopulationCrossover() { }

	virtual errut::bool_t check(const std::shared_ptr<Population> &population)
	{
		m_tmp[0] = population;
		auto r = check(m_tmp);
		m_tmp[0] = nullptr; // Don't keep a reference to this object
		return r;
	}
	virtual errut::bool_t check(const std::vector<std::shared_ptr<Population>> &populations){ return "Not implemented in base class"; }
	// The populations are overwritten, if the old one is still needed it should
	// be stored externally
	virtual errut::bool_t createNewPopulation(size_t generation, std::vector<std::shared_ptr<Population>> &populations, size_t targetPopulationSize) { return "Not implemented in base class"; }
	virtual errut::bool_t createNewPopulation(size_t generation, std::shared_ptr<Population> &population, size_t targetPopulationSize)
	{
		m_tmp[0] = population;
		auto r = createNewPopulation(generation, m_tmp, targetPopulationSize);
		std::swap(m_tmp[0], population);
		m_tmp[0] = nullptr; // Don't keep a reference to this object 
		return r;
	}

	// This is the main class that creates new populations/individuals, this
	// should also be the access point to keep track of the best ones
	virtual const std::vector<std::shared_ptr<Individual>> &getBestIndividuals() const { return m_emptyBest; }
private:
	std::vector<std::shared_ptr<Population>> m_tmp;
	const std::vector<std::shared_ptr<Individual>> m_emptyBest;
};

class PopulationCrossoverIteration
{
public:
	PopulationCrossoverIteration() { }
	virtual ~PopulationCrossoverIteration() { }

	// new population can already have some individuals from elitims for example
	virtual void startNewIteration(std::shared_ptr<Population> &newPopulation, size_t targetPopulationSize) { }
	// TODO: do we need other arguments here?
	virtual bool iterate(std::shared_ptr<Population> &newPopulation) { return false; }
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
	virtual errut::bool_t processPopulation(std::shared_ptr<Population> &population, size_t targetPopulationSize) { return "Not implemented in base class"; }

	virtual const std::vector<std::shared_ptr<Individual>> &getBestIndividuals() const { return m_emptyBest; }
private:
	const std::vector<std::shared_ptr<Individual>> m_emptyBest;
};

class Elitism
{
public:
	Elitism() { }
	virtual ~Elitism() { }

	virtual errut::bool_t check(const std::shared_ptr<SelectionPopulation> &selPop) { return "Not implemented in base class"; }

	// Population should be empty, number of genomes to skip should be set in this function
	virtual errut::bool_t introduceElites(size_t generation, const std::shared_ptr<SelectionPopulation> &selPop,
										  std::shared_ptr<Population> &population,
										  size_t targetPopulationSize) { return "Not implemented in base class"; }
};

class ParentSelection // Don't think we need some kind of check, shouldn't depend on type of genome
{
public:
	ParentSelection() { }
	virtual ~ParentSelection() { }
	
	virtual errut::bool_t check(const SelectionPopulation &pop) { return "Not implemented in base class"; }
	// The length of 'parents' describes the number of parents that should be
	virtual errut::bool_t selectParents(const SelectionPopulation &pop, std::vector<std::shared_ptr<Individual>> &parents) { return "Not implemented in base class"; }
};

}
