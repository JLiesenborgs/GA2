#pragma once

#include "genomefitness.h"
#include <vector>

class GenomeCrossover
{
public:
	GenomeCrossover() { }
	virtual ~GenomeCrossover() { }

	// This function is intended to check number and type of parents, so that
	// e.g. no dynamic cast is needed in generateOffspring itself, the number
	// of parents and their type is assumed to be checked
	virtual errut::bool_t check(const std::vector<const Genome*> &parents) { return "Not implemented in base class"; }
	// This needs to be implemented, results stored in m_offspringBuffer
	virtual errut::bool_t generateOffspring(const std::vector<const Genome*> &parents) { return "Not implemented in base class"; }
	const std::vector<std::shared_ptr<Genome>> &getGeneratedOffpring() const { return m_offspringBuffer; }
protected:
	std::vector<std::shared_ptr<Genome>> m_offspringBuffer;
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
	PopulationMutation() { }
	virtual ~PopulationMutation() { }
	virtual errut::bool_t check(const std::vector<std::shared_ptr<Population>> &populations) { return "Not implemented in base class"; }
	// This is in-place, must reset isCalculated flag if changed
	// Idea is to apply a GenomeMutation operator to every individual
	virtual errut::bool_t mutate(const std::vector<std::shared_ptr<Population>> &populations) { return "Not implemented in base class"; }
};

// TODO: PopulationCrossover ?
// allow in-place? some children that overwrite older parents?
class PopulationCrossover
{
public:
	PopulationCrossover() { }
	virtual ~PopulationCrossover() { }

	virtual errut::bool_t check(const std::vector<std::shared_ptr<Population>> &populations) { return "Not implemented in base class"; }
	// The populations are overwritten, if the old one is still needed it should
	// be stored externally
	// This is in-place, must reset isCalculated flag if changed
	// Idea is to apply a GenomeMutation operator to every individual
	virtual errut::bool_t createNewPopulations(std::vector<std::shared_ptr<Population>> &populations) { return "Not implemented in base class"; }
};

class ParentSelection // Don't think we need some kind of check, shouldn't depend on type of genome
{
public:
	ParentSelection() { }
	virtual ~ParentSelection() { }
	// TODO: should we consider multiple populations instead?
	virtual errut::bool_t selectParents(const Population &pop, std::vector<const Genome *> &parents) { return "Not implemented in base class"; }
};

