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
	NSGA2IndividualWrapper(size_t numObjectives,
			const std::shared_ptr<Genome> &genome, const std::shared_ptr<Fitness> &fitness,
			size_t introducedInGeneration, size_t originalPosition);

	std::shared_ptr<Individual> createNew(const std::shared_ptr<Genome> &genome, const std::shared_ptr<Fitness> &fitness,
			size_t introducedInGeneration = std::numeric_limits<size_t>::max()) const override;

	std::string toString() const override;

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
protected:
	virtual std::shared_ptr<NonDominatedSetCreator> allocatedNDSetCreator(const std::shared_ptr<FitnessComparison> &fitCmp, size_t numObjectives);
private:
	void buildWrapperPopulation(const Population &population);
	void calculateCrowdingDistances(const std::vector<std::shared_ptr<Individual>> &ndset) const;

	std::shared_ptr<FitnessComparison> m_fitComp;
	std::unique_ptr<SinglePopulationCrossover> m_crossover;
	std::shared_ptr<Population> m_tmpPop;
	const size_t m_numObjectives;

	std::vector<std::shared_ptr<Individual>> m_popWrapper;
	std::vector<std::shared_ptr<Individual>> m_best;
	std::shared_ptr<NSGA2FitnessWrapperOriginalComparison> m_fitOrigComp;
	std::shared_ptr<NonDominatedSetCreator> m_ndSetCreator;

	bool m_alwaysRebuildWrapper;
};

}
