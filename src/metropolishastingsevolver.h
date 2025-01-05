#pragma once

#include "eatkconfig.h"
#include "populationevolver.h"
#include "randomnumbergenerator.h"

namespace eatk
{

// TODO: use common interface with GoodmanWeareEvolver

class MetropolisHastingsEvolver : public PopulationEvolver
{
public:
	enum ProbType { Regular, Log, NegativeLog };

	MetropolisHastingsEvolver(const std::shared_ptr<RandomNumberGenerator> &rng,
			            const std::vector<double> &stepScales, // should be as many components as fitness
						ProbType t = Regular);
	~MetropolisHastingsEvolver();

	errut::bool_t check(const std::shared_ptr<Population> &population) override;

	// Members must be vector genomes, either double or float, fitness is a value
	// fitness, double or float
	
	// Every individual does an independent Metropolis-Hastings walk
	errut::bool_t createNewPopulation(size_t generation, std::shared_ptr<Population> &population, size_t targetPopulationSize) override;

	// The idea here is to keep track of the individual with the highest
	// logprob/prob (fitness)
	const std::vector<std::shared_ptr<Individual>> &getBestIndividuals() const override { return m_bestIndividual; }

	// TODO: use common code ?
	void setAnnealingExponent(double alpha) { m_alpha = alpha; } // to use prob^alpha
protected:
	// These are the actual individuals, not copies, for efficiency
	// TODO: change this?
	// TODO: use common base class as GoodmanWeareEvolver
	virtual void onSamples(const std::vector<std::shared_ptr<Individual>> &samples) { }
private:
	std::shared_ptr<RandomNumberGenerator> m_rng;
	std::vector<double> m_stepScales;
	ProbType m_probType = Regular;
	bool m_doubleGenomes = false;
	double m_alpha = 1.0;

	std::vector<std::shared_ptr<Individual>> m_bestIndividual;
	std::vector<std::shared_ptr<Individual>> m_samples;
};

} // end namespace

