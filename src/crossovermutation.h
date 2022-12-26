#pragma once

#include "eatkconfig.h"
#include "genomefitness.h"
#include "population.h"
#include <vector>

namespace eatk
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

}
