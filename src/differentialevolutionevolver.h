#pragma once

#include "eatkconfig.h"
#include "crossovermutation.h"
#include "randomnumbergenerator.h"
#include "populationevolver.h"
#include "nondominatedsetcreator.h"
#include <functional>

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
		const std::shared_ptr<FitnessComparison> &fitComp, int objectiveNumber = 0, size_t numObjectives = 1,
		const std::shared_ptr<NonDominatedSetCreator> &ndCreator = nullptr); // negative means multi-objective
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

	std::vector<std::shared_ptr<Individual>> m_bestIndividuals;
	std::shared_ptr<NonDominatedSetCreator> m_ndCreator;
};

inline std::function<bool(const Fitness &f1, const Fitness &f2)> getDominatesFunction(const std::shared_ptr<FitnessComparison> &fitComp, int objectiveNumber, size_t numObjectives)
{
	std::function<bool(const Fitness &f1, const Fitness &f2)> isFitterThan;

	if (numObjectives <= 1) // single objective
	{
		isFitterThan = [fitComp, objectiveNumber](const Fitness &f1, const Fitness &f2)
		{
			return fitComp->isFitterThan(f1, f2, objectiveNumber);
		};
	}
	else
	{
		// multi-objective, use dominance
		isFitterThan = [fitComp, numObjectives](const Fitness &f1, const Fitness &f2)
		{
			size_t betterOrEqualCount = 0;
			size_t betterCount = 0;
			for (size_t i = 0 ; i < numObjectives ; i++)
			{
				if (fitComp->isFitterThan(f1, f2, i))
				{
					betterCount++;
					betterOrEqualCount++;
				}
				else // f1 not strictly better than f2 for i
				{
					if (!fitComp->isFitterThan(f2, f1, i)) // then they must have equal fitness
					{
						betterOrEqualCount++;
					}
					else
					{
						// We can never get betterOrEqualCount == m_numObjectives
						return false;
					}
				}
			}
			// if we got here, then betterOrEqualCount == m_numObjectives
			if (betterCount > 0)
				return true;
			return false;
		};
	}
	return isFitterThan;
}

}
