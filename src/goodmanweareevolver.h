#pragma once

#include "eatkconfig.h"
#include "populationevolver.h"
#include "randomnumbergenerator.h"

namespace eatk
{

// based on https://research.edm.uhasselt.be/~jori/page/Misc/MeanWalker.html and
// on parallelization info in "emcee: The MCMC Hammer" (https://doi.org/10.1086/670067)

class GoodmanWeareEvolver : public PopulationEvolver
{
public:
	GoodmanWeareEvolver(const std::shared_ptr<RandomNumberGenerator> &rng, bool isLog = true, double a = 2.0);
	~GoodmanWeareEvolver();

	errut::bool_t check(const std::shared_ptr<Population> &population) override;

	// Members must be vector genomes, either double or float, fitness is a value
	// fitness, double or float
	
	// In the first step, N = targetPopulationSize is the number of individuals,
	// for which the probabilty/logprob (fitness) is calculated.
	// The newly added members are temporarily appended, so their logprobs can
	// be calculated. On the next iterations there are targetPopulationSize*2
	// individuals, and the new ones can be compared to the old ones to see
	// which need to be replaced. The first N are always the actual samples
	// (how best to report this?)
	errut::bool_t createNewPopulation(size_t generation, std::shared_ptr<Population> &population, size_t targetPopulationSize) override;

	// The idea here is to keep track of the individual with the highest
	// logprob/prob (fitness)
	const std::vector<std::shared_ptr<Individual>> &getBestIndividuals() const override { return m_bestIndividual; }
protected:
	// These are the actual individuals, not copies, for efficiency
	// TODO: change this?
	virtual void onSamples(const std::vector<std::shared_ptr<Individual>> &samples) { }
private:
	std::shared_ptr<RandomNumberGenerator> m_rng;
	double m_a = 0, m_aScale, m_aOffset;
	bool m_isLog = true;
	bool m_doubleGenomes = false;
	std::vector<double> m_zValues;

	const size_t m_objectiveNumber = 0;

	std::vector<std::shared_ptr<Individual>> m_bestIndividual;
	std::vector<std::shared_ptr<Individual>> m_samples;
};

} // end namespace

