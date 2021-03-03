#pragma once

#include "genomefitness.h"
#include <vector>

// TODO: record parents?
// TODO: generation of creation? To allow age?
class Individual
{
public:
	Individual(std::shared_ptr<Genome> genome, std::shared_ptr<Fitness> fitness)
		: m_genome(genome), m_fitness(fitness) { }
//private:
	std::shared_ptr<Genome> m_genome;
	std::shared_ptr<Fitness> m_fitness;
};

// Do we need this? Just a typedef perhaps?
class Population
{
public:
	Population() { }
	~Population() { }

	std::vector<std::shared_ptr<Individual>> m_individuals;
};

class PopulationFitnessCalculation
{
public:
	PopulationFitnessCalculation() { }
	virtual ~PopulationFitnessCalculation() { }

	// This function can be called eg on first generation, to check that the types etc
	// are in order. This can then skip some checks in calculatePopulationFitness
	virtual errut::bool_t check(const std::vector<std::shared_ptr<Population>> &populations) { return "Not implemented in base class"; }

	// TODO: all populations should have exactly the same genomes! (ie same number of floats)
	virtual errut::bool_t calculatePopulationFitness(const std::vector<std::shared_ptr<Population>> &populations) { return "Not implemented in base class"; }
};

