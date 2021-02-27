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
	virtual errut::bool_t mutate(Genome &genome) { return "Not implemented in base class"; }
};
