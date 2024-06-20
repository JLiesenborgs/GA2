#pragma once

#include "eatkconfig.h"
#include "crossovermutation.h"
#include "randomnumbergenerator.h"
#include "populationevolver.h"
#include "nondominatedsetcreator.h"

namespace eatk
{

class DifferentialEvolutionMutation
{
public:
	DifferentialEvolutionMutation() { }
	virtual ~DifferentialEvolutionMutation() { }

	virtual errut::bool_t check(const Genome &g) { return "Not implemented in base class"; }
	virtual std::shared_ptr<Genome> mutate(const std::vector<const Genome*> &genomes, const std::vector<double> &weights) { return nullptr; }
};

class DifferentialEvolutionCrossover
{
public:
	DifferentialEvolutionCrossover() { }
	virtual ~DifferentialEvolutionCrossover() { }

	virtual errut::bool_t check(const Genome &g) { return "Not implemented in base class"; }
	virtual errut::bool_t crossover(double CR, Genome &mutantDest, const Genome &origVector) { return "Not implemented in base class"; }
};

class DifferentialEvolutionEvolver : public PopulationEvolver
{
public:
	DifferentialEvolutionEvolver(
		const std::shared_ptr<RandomNumberGenerator> &rng,
		const std::shared_ptr<DifferentialEvolutionMutation> &mut,
		double F,
		const std::shared_ptr<DifferentialEvolutionCrossover> &cross,
		double CR,
		const std::shared_ptr<FitnessComparison> &fitComp, int objectiveNumber = 0,	// negative means multi-objective
		size_t numObjectives = 1,
		const std::shared_ptr<NonDominatedSetCreator> &ndCreator = nullptr,
		bool needStrictlyBetter = true);
	~DifferentialEvolutionEvolver();

	errut::bool_t check(const std::shared_ptr<Population> &population) override;
	errut::bool_t createNewPopulation(size_t generation, std::shared_ptr<Population> &population, size_t targetPopulationSize) override;
	const std::vector<std::shared_ptr<Individual>> &getBestIndividuals() const override { return m_bestIndividuals; }
private:
	std::shared_ptr<RandomNumberGenerator> m_rng;
	std::shared_ptr<DifferentialEvolutionMutation> m_mut;
	std::shared_ptr<DifferentialEvolutionCrossover> m_cross;
	std::shared_ptr<FitnessComparison> m_fitComp;
	const int m_objectiveNumber;
	const size_t m_numObjectives;

	std::vector<const Genome*> m_mutationGenomes;
	std::vector<double> m_mutationFactors;
	double m_CR;
	bool m_needStrictlyBetter;

	std::vector<std::shared_ptr<Individual>> m_bestIndividuals;
	std::shared_ptr<NonDominatedSetCreator> m_ndCreator;
};

inline std::function<bool(const Fitness &f1, const Fitness &f2)> getShouldAcceptFunction(const std::shared_ptr<FitnessComparison> &fitComp,
		                                                                                 int objectiveNumber, size_t numObjectives,
																						 bool needStrictlyBetter)
{
	auto isFitterThan = FitnessComparison::getDominatesFunction(fitComp, objectiveNumber, numObjectives);
	std::function<bool(const Fitness &f1, const Fitness &f2)> shouldAcceptFunction;

	if (needStrictlyBetter)
		shouldAcceptFunction = isFitterThan; // [isFitterThan](const Fitness &f1, const Fitness &f2) { return isFitterThan(f1, f2); };
	else
		shouldAcceptFunction = [isFitterThan](const Fitness &f1, const Fitness &f2) { return !isFitterThan(f2, f1); };
	return shouldAcceptFunction;
}

}
