#pragma once

#include "eatkconfig.h"
#include "samplingevolver.h"
#include "randomnumbergenerator.h"

namespace eatk
{

// based on https://research.edm.uhasselt.be/~jori/page/Misc/MeanWalker.html and
// on parallelization info in "emcee: The MCMC Hammer" (https://doi.org/10.1086/670067)

class GoodmanWeareEvolver : public SamplingEvolver
{
public:
	GoodmanWeareEvolver(const std::shared_ptr<RandomNumberGenerator> &rng, ProbType t = Regular, double a = 2.0);
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
private:
	double m_a = 0, m_aScale, m_aOffset;
	std::vector<double> m_zValues;
	std::vector<std::shared_ptr<Individual>> m_samples;
};

} // end namespace

