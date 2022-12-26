#pragma once

#include "eatkconfig.h"
#include "individual.h"
#include <vector>
#include <iostream>
#include <cassert>
#include <limits>

namespace eatk
{

// Do we need this? Just a typedef perhaps?
class Population
{
public:
	Population() { }
	~Population() { }

	void clear() { m_individuals.clear(); }
	size_t size() const { return m_individuals.size(); }
	void append(const std::shared_ptr<Individual> &i)
	{
		assert(i.get());
		m_individuals.push_back(i);
	}
	void resize(size_t n) { m_individuals.resize(n); }

	const std::vector<std::shared_ptr<Individual>> &individuals() const { return m_individuals; }
	std::vector<std::shared_ptr<Individual>> &individuals() { return m_individuals; }
	
	std::shared_ptr<Individual> &individual(size_t n)
	{
		assert(n < m_individuals.size());
		return m_individuals[n];
	}

	const std::shared_ptr<Individual> &individual(size_t n) const
	{
		assert(n < m_individuals.size());
		return m_individuals[n];
	}

	void print() const
	{
		using namespace std;

		cout << "Population: " << endl;
		for (auto &i : m_individuals)
			cout << i->toString() << endl;
		cout << endl;
	}
private:
	std::vector<std::shared_ptr<Individual>> m_individuals;
};

class PopulationCrossoverIteration
{
public:
	PopulationCrossoverIteration() { }
	virtual ~PopulationCrossoverIteration() { }

	// new population can already have some individuals from elitims for example
	virtual void startNewIteration(const Population &newPopulation, size_t targetPopulationSize) { }
	// TODO: do we need other arguments here?
	virtual bool iterate(const Population &newPopulation) { return false; }
};

}
