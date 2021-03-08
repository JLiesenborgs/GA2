#pragma once

#include "genomefitness.h"
#include <vector>
#include <iostream>

// TODO: record parents?
// TODO: generation of creation? To allow age?
class Individual
{
public:
	Individual(std::shared_ptr<Genome> genome, std::shared_ptr<Fitness> fitness)
		: m_genome(genome), m_fitness(fitness) { }

	std::string toString() const
	{
		return m_genome->toString() + ": " + m_fitness->toString();
	}
	
//private:
	std::shared_ptr<Genome> m_genome;
	std::shared_ptr<Fitness> m_fitness;
};

// Do we need this? Just a typedef perhaps?
class Population
{
public:
	Population() : m_skipMutation(0) { }
	~Population() { }

	size_t size() const { return m_individuals.size(); }

	void setGenomesToSkipMutation(size_t n) { m_skipMutation = n; }
	size_t getGenomesToSkipMutation() const { return m_skipMutation; }

	void print() const
	{
		using namespace std;

		cout << "Population: " << endl;
	    for (auto &i : m_individuals)
    	    cout << i->toString() << endl;
    	cout << endl;
	}

	std::vector<std::shared_ptr<Individual>> m_individuals;
	size_t m_skipMutation;
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

