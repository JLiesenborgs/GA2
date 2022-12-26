#pragma once

#include "eatkconfig.h"
#include "selection.h"

namespace eatk
{

class NDSortedPopulationInterface : public SelectionPopulation
{
public:
	NDSortedPopulationInterface() { }
	~NDSortedPopulationInterface() { }

	virtual size_t getNumberOfSets() const = 0;
	virtual size_t getSetSize(size_t s) const = 0;
	virtual const std::shared_ptr<Individual> &getIndividual(size_t s, size_t i) const = 0;
};

}
