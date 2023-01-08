#pragma once

#include "eatkconfig.h"
#include "populationevolver.h"
#include "nondominatedsetcreator.h"
#include "duplicateindividualremoval.h"

namespace eatk
{

class BestIndividualMerger
{
public:
	BestIndividualMerger() { }
	virtual ~BestIndividualMerger() { }

	virtual errut::bool_t check(const std::vector<std::shared_ptr<Population>> &populations,
	                            const std::vector<std::shared_ptr<PopulationEvolver>> &evolvers) { return "Not implemented in base class"; }
	virtual const std::vector<std::shared_ptr<Individual>> &mergeBestIndividuals(const std::vector<std::shared_ptr<PopulationEvolver>> &evolvers) { return m_dummy; }
public:
	std::vector<std::shared_ptr<Individual>> m_dummy;
};

class SingleObjectiveBestIndividualMerger : public BestIndividualMerger
{
public:
	SingleObjectiveBestIndividualMerger(const std::shared_ptr<FitnessComparison> &fitCmp);
	~SingleObjectiveBestIndividualMerger();

	errut::bool_t check(const std::vector<std::shared_ptr<Population>> &populations,
	                    const std::vector<std::shared_ptr<PopulationEvolver>> &evolvers) override;
	const std::vector<std::shared_ptr<Individual>> &mergeBestIndividuals(const std::vector<std::shared_ptr<PopulationEvolver>> &evolvers) override;
private:
	std::shared_ptr<FitnessComparison> m_fitCmp;
	std::vector<std::shared_ptr<Individual>> m_best;
};

class MultiObjectiveBestIndividualMerger : public BestIndividualMerger
{
public:
	MultiObjectiveBestIndividualMerger(const std::shared_ptr<NonDominatedSetCreator> &ndCreator,
	                                   const std::shared_ptr<DuplicateIndividualRemoval> &dupRemoval = nullptr);
	~MultiObjectiveBestIndividualMerger();

	errut::bool_t check(const std::vector<std::shared_ptr<Population>> &populations,
	                    const std::vector<std::shared_ptr<PopulationEvolver>> &evolvers) override;
	const std::vector<std::shared_ptr<Individual>> &mergeBestIndividuals(const std::vector<std::shared_ptr<PopulationEvolver>> &evolvers) override;
public:
	std::shared_ptr<NonDominatedSetCreator> m_ndSetCreator;
	std::shared_ptr<DuplicateIndividualRemoval> m_dupRemoval;
	std::vector<std::shared_ptr<Individual>> m_best;
};

class MultiPopulationEvolver : public PopulationEvolver
{
public:
	MultiPopulationEvolver(const std::vector<std::shared_ptr<PopulationEvolver>> &singlePopulationEvolvers,
	                       const std::shared_ptr<BestIndividualMerger> &merger);
	~MultiPopulationEvolver();

	errut::bool_t check(const std::vector<std::shared_ptr<Population>> &populations) override;
	errut::bool_t createNewPopulations(size_t generation, std::vector<std::shared_ptr<Population>> &populations, const std::vector<size_t> &targetPopulationSize) override;

	const std::vector<std::shared_ptr<PopulationEvolver>> &getSinglePopulationEvolvers() const { return m_singlePopEvolvers; }
	const std::vector<std::shared_ptr<Individual>> &getBestIndividuals() const override;
private:
	std::vector<std::shared_ptr<PopulationEvolver>> m_singlePopEvolvers;
	size_t m_lastCreatedGeneration;
	mutable std::shared_ptr<BestIndividualMerger> m_merger;
	mutable size_t m_lastMergeGeneration;
	mutable std::vector<std::shared_ptr<Individual>> m_best;
};

}
