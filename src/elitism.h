#pragma once

#include "eatkconfig.h"
#include "population.h"
#include "selection.h"

namespace eatk
{

class Elitism
{
public:
	Elitism() { }
	virtual ~Elitism() { }

	virtual errut::bool_t check(const std::shared_ptr<SelectionPopulation> &selPop) { return "Not implemented in base class"; }

	// Population should be empty, number of genomes to skip should be set in this function
	virtual errut::bool_t introduceElites(size_t generation, const std::shared_ptr<SelectionPopulation> &selPop,
										  const std::shared_ptr<Population> &population,
										  size_t targetPopulationSize) { return "Not implemented in base class"; }
};

}
