#pragma once

#include "eatkconfig.h"
#include "populationevolver.h"
#include "crossovermutation.h"
#include "randomnumbergenerator.h"
#include "selection.h"
#include "singlepopulationcrossover.h"
#include "nondominatedsetcreator.h"
#include <memory>
#include <array>
#include <cmath>

namespace eatk
{

class NSGA2IndividualWrapper : public Individual
{
public:
	NSGA2IndividualWrapper(const std::shared_ptr<Genome> &genome, const std::shared_ptr<Fitness> &fitness,
			size_t introducedInGeneration, size_t originalPosition);

	std::shared_ptr<Individual> createNew(const std::shared_ptr<Genome> &genome, const std::shared_ptr<Fitness> &fitness,
			size_t introducedInGeneration = std::numeric_limits<size_t>::max()) const override;

	std::string toString() const override;

	static double getCrowdingValue(const std::shared_ptr<Individual> &ind);

	size_t m_originalPosition;
	std::vector<double> m_fitnessDistances;
};

class NSGA2FitnessWrapper : public Fitness
{
public:
	NSGA2FitnessWrapper(const std::shared_ptr<Fitness> &origFitness);
	std::shared_ptr<Fitness> createCopy(bool copyContents = true) const override;
	std::string toString() const override;

	bool hasRealValues() const override { return m_origFitness->hasRealValues(); }
	double getRealValue(size_t objectiveNumber) const override { return m_origFitness->getRealValue(objectiveNumber); }

	std::shared_ptr<Fitness> m_origFitness;
	double m_totalFitnessDistance;
	size_t m_ndSetIndex;
};

inline double NSGA2IndividualWrapper::getCrowdingValue(const std::shared_ptr<Individual> &ind)
{
	assert(dynamic_cast<const NSGA2IndividualWrapper*>(ind.get()));
	const Fitness &f = static_cast<const NSGA2IndividualWrapper*>(ind.get())->fitnessRef();

	assert(dynamic_cast<const NSGA2FitnessWrapper*>(&f));
	const NSGA2FitnessWrapper &fw = static_cast<const NSGA2FitnessWrapper&>(f);

	return fw.m_totalFitnessDistance;
}

class NSGA2FitnessWrapperOriginalComparison : public FitnessComparison
{
public:
	NSGA2FitnessWrapperOriginalComparison(const std::shared_ptr<FitnessComparison> &origCmp);
	errut::bool_t check(const Fitness &f) const;
	bool isFitterThan(const Fitness &first0, const Fitness &second0, size_t objectiveNumber) const;
private:
	std::shared_ptr<FitnessComparison> m_origCmp;
};

class NSGA2FitWrapperNDSetCrowdingComparison : public FitnessComparison
{
public:
	errut::bool_t check(const Fitness &f) const;
	bool isFitterThan(const Fitness &first0, const Fitness &second0, size_t objectiveNumber) const;
};

class NSGA2Evolver : public PopulationEvolver
{
public:
	NSGA2Evolver(
		const std::shared_ptr<RandomNumberGenerator> &rng,
		const std::shared_ptr<GenomeCrossover> &genomeCrossover,
		const std::shared_ptr<GenomeMutation> &genomeMutation,
		const std::shared_ptr<FitnessComparison> &fitComp, size_t numObjectives,
		// In the default usage of the algorithm, we'll still have the wrapper
		// individuals from the previous generation and can use these. If the
		// population is being changed in between, the wrappers will need to
		// be rebuilt
		bool alwaysRebuildWrapperPopulation = false
		);
	~NSGA2Evolver();

	errut::bool_t check(const std::shared_ptr<Population> &population) override;
	errut::bool_t createNewPopulation(size_t generation, std::shared_ptr<Population> &population, size_t targetPopulationSize) override;
	
	const std::vector<std::shared_ptr<Individual>> &getBestIndividuals() const override { return m_best; }

	static void buildWrapperPopulation(const Population &pop, Population &wrapperPop);
	static void unwrapPopulation(const Individual &refInd, const Population &wrapperPop, Population &pop);
	// Should be a vector wrapper individuals
	static void calculateCrowdingDistances(const std::vector<std::shared_ptr<Individual>> &ndset, size_t numObjectives);
protected:
	virtual std::shared_ptr<NonDominatedSetCreator> allocateNDSetCreator(const std::shared_ptr<FitnessComparison> &fitCmp, size_t numObjectives);
private:
	errut::bool_t createNewPopulation_Single(size_t generation, std::shared_ptr<Population> &population, size_t targetPopulationSize);
	errut::bool_t createNewPopulation_Multi(size_t generation, std::shared_ptr<Population> &population, size_t targetPopulationSize);

	std::shared_ptr<FitnessComparison> m_fitComp;
	std::unique_ptr<SinglePopulationCrossover> m_crossover;
	std::shared_ptr<Population> m_wrapperPop;
	const size_t m_numObjectives;

	std::vector<std::shared_ptr<Individual>> m_best;
	std::shared_ptr<NSGA2FitnessWrapperOriginalComparison> m_fitOrigComp;
	std::shared_ptr<NonDominatedSetCreator> m_ndSetCreator;

	bool m_alwaysRebuildWrapper;
};

inline void NSGA2Evolver::buildWrapperPopulation(const Population &pop, Population &wrapperPop)
{
	wrapperPop.clear();
	for (size_t i = 0 ; i < pop.size() ; i++)
	{
		auto &ind = pop.individual(i);
		wrapperPop.append(std::make_shared<NSGA2IndividualWrapper>(ind->genome(),
		                  std::make_shared<NSGA2FitnessWrapper>(ind->fitness()),
		                  ind->getIntroducedInGeneration(), i));
	}
}

inline void NSGA2Evolver::unwrapPopulation(const Individual &refInd, const Population &wrapperPop, Population &population)
{
	population.clear();

	for (size_t i = 0 ; i < wrapperPop.size() ; i++)
	{
		const std::shared_ptr<Individual> &wrapperInd = wrapperPop.individual(i);
		const std::shared_ptr<Fitness> &wrapperFit = wrapperInd->fitness();
		assert(dynamic_cast<NSGA2FitnessWrapper*>(wrapperFit.get()));

		const NSGA2FitnessWrapper *pWrapperFit = static_cast<const NSGA2FitnessWrapper*>(wrapperFit.get());

		auto newInd = refInd.createNew(wrapperInd->genome(), pWrapperFit->m_origFitness, wrapperInd->getIntroducedInGeneration());
		newInd->setLastMutationGeneration(wrapperInd->getLastMutationGeneration());
		population.append(newInd);
	}
}

}
