#pragma once

#include "eatkconfig.h"
#include "ndsortedpopulationinterface.h"
#include "nondominatedsetcreator.h"
#include "duplicateindividualremoval.h"

namespace eatk
{

class NDSortedPopulation : public NDSortedPopulationInterface
{
public:
	NDSortedPopulation(const std::shared_ptr<NonDominatedSetCreator> &setCreator,
					   const std::shared_ptr<DuplicateIndividualRemoval> &dupRemoval = nullptr,
					   bool bestOnlyDuplicateRemoval = true);
	~NDSortedPopulation();

	void clearBest() { m_best.resize(0); }

	errut::bool_t check(const Population &population) override;
	errut::bool_t processPopulation(const std::shared_ptr<Population> &population, size_t targetPopulationSize) override;
	const std::vector<std::shared_ptr<Individual>> &getBestIndividuals() const override { return m_best; }

	size_t getNumberOfSets() const override { return m_ndSetCreator->getNumberOfSets(); }
	size_t getSetSize(size_t s) const override { return m_ndSetCreator->getSet(s).size(); }
	const std::shared_ptr<Individual> &getIndividual(size_t s, size_t i) const override { return m_ndSetCreator->getSet(s)[i]; }
private:
	std::vector<std::shared_ptr<Individual>> m_best;
	std::shared_ptr<NonDominatedSetCreator> m_ndSetCreator;
	std::shared_ptr<DuplicateIndividualRemoval> m_dupRemoval;
	bool m_removeFromBestOnly;
};

}
