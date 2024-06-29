#pragma once

#include "eatkconfig.h"
#include "population.h"
#include <limits>

namespace eatk
{

class NonDominatedSetCreator
{
public:
	NonDominatedSetCreator() { }
	virtual ~NonDominatedSetCreator() { }

	// Idea behind requestStopSize is to stop looking for other
	// non-dominated sets once this threshold has been reached or
	// exceeded. The non-dominated set itself should be complete,
	// and it is only a request, no obligation to take this into
	// account - all nd sets may still be calculated
	virtual errut::bool_t calculateAllNDSets(const std::vector<std::shared_ptr<Individual>> &individuals,
	                                         size_t requestStopSize = std::numeric_limits<size_t>::max()) { return "Not implemented"; }
	virtual errut::bool_t calculateNonDomitatedSet(const std::vector<std::shared_ptr<Individual>> &individuals,
		std::vector<std::shared_ptr<Individual>> &ndSet,
		std::vector<std::shared_ptr<Individual>> &remaining) { return "Not implemented"; }

	// Assumes that both are in fact ND sets
	virtual errut::bool_t mergeNDSets(std::vector<std::shared_ptr<Individual>> &inOut, const std::vector<std::shared_ptr<Individual>> &added) { return "Not implemented"; }

	virtual size_t getNumberOfSets() const { return 0; }
	virtual const std::vector<std::shared_ptr<Individual>> &getSet(size_t idx) const { return m_empty; }
	virtual std::vector<std::shared_ptr<Individual>> &getSet(size_t idx) { return m_empty; }
private:
	std::vector<std::shared_ptr<Individual>> m_empty;
};

};
