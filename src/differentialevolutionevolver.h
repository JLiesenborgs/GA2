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
	virtual std::shared_ptr<Genome> mutate(const Genome &r1, const Genome &r2, const Genome &r3) { return nullptr; }
};

class DifferentialEvolutionCrossover
{
public:
	DifferentialEvolutionCrossover() { }
	virtual ~DifferentialEvolutionCrossover() { }

	virtual errut::bool_t check(const Genome &g) { return "Not implemented in base class"; }
	virtual errut::bool_t crossover(Genome &mutantDest, const Genome &origVector) { return "Not implemented in base class"; }
};

class DifferentialEvolutionEvolver : public PopulationEvolver
{
public:
	DifferentialEvolutionEvolver(
		const std::shared_ptr<RandomNumberGenerator> &rng,
		const std::shared_ptr<DifferentialEvolutionMutation> &mut,
		const std::shared_ptr<DifferentialEvolutionCrossover> &cross,
		const std::shared_ptr<FitnessComparison> &fitComp, size_t objectiveNumber = 0);
	~DifferentialEvolutionEvolver();

	errut::bool_t check(const std::shared_ptr<Population> &population) override;
	errut::bool_t createNewPopulation(size_t generation, std::vector<std::shared_ptr<Population>> &populations, size_t targetPopulationSize) override;
	const std::vector<std::shared_ptr<Individual>> &getBestIndividuals() const override { return m_bestIndividual; }
private:
	std::shared_ptr<RandomNumberGenerator> m_rng;
	std::shared_ptr<DifferentialEvolutionMutation> m_mut;
	std::shared_ptr<DifferentialEvolutionCrossover> m_cross;
	std::shared_ptr<FitnessComparison> m_fitComp;
	size_t m_objectiveNumber;

	std::vector<std::shared_ptr<Individual>> m_bestIndividual;
};

}