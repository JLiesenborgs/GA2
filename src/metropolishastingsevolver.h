#pragma once

#include "eatkconfig.h"
#include "samplingevolver.h"
#include "randomnumbergenerator.h"

namespace eatk
{

class MetropolisHastingsEvolver : public SamplingEvolver
{
public:
	MetropolisHastingsEvolver(const std::shared_ptr<RandomNumberGenerator> &rng,
			            const std::vector<double> &stepScales, // should be as many components as fitness
						ProbType t = Regular);
	~MetropolisHastingsEvolver();

	errut::bool_t check(const std::shared_ptr<Population> &population) override;

	// Members must be vector genomes, either double or float, fitness is a value
	// fitness, double or float
	
	// Every individual does an independent Metropolis-Hastings walk
	errut::bool_t createNewPopulation(size_t generation, std::shared_ptr<Population> &population, size_t targetPopulationSize) override;
private:
	std::vector<double> m_stepScales;
	std::vector<std::shared_ptr<Individual>> m_samples;
};

} // end namespace

