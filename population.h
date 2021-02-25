#pragma once

#include "genomefitness.h"
#include <vector>

// TODO: record parents?
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

	virtual errut::bool_t calculatePopulationFitness(const std::vector<std::shared_ptr<Population>> &populations) { return "Not implemented in base class"; }
};

