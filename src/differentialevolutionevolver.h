#pragma once

#include "eatkconfig.h"
#include "crossovermutation.h"
#include "randomnumbergenerator.h"

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
		const std::shared_ptr<FitnessComparison> &fitComp, size_t objectiveNumber = 0);
	~DifferentialEvolutionEvolver();

	errut::bool_t check(const std::shared_ptr<Population> &population) override;
	errut::bool_t createNewPopulation(size_t generation, std::vector<std::shared_ptr<Population>> &populations, const std::vector<size_t> &targetPopulationSize) override;
	const std::vector<std::shared_ptr<Individual>> &getBestIndividuals() const override { return m_bestIndividual; }
private:
	std::shared_ptr<RandomNumberGenerator> m_rng;
	std::shared_ptr<DifferentialEvolutionMutation> m_mut;
	std::shared_ptr<DifferentialEvolutionCrossover> m_cross;
	std::shared_ptr<FitnessComparison> m_fitComp;
	const size_t m_objectiveNumber;

	std::vector<const Genome*> m_mutationGenomes;
	std::vector<double> m_mutationFactors;
	double m_CR;

	std::vector<std::shared_ptr<Individual>> m_bestIndividual;
};

}
