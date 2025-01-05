#pragma once

#include "eatkconfig.h"
#include "populationevolver.h"
#include "randomnumbergenerator.h"
#include "vectorgenomefitness.h"
#include <cassert>

namespace eatk
{

class SamplingEvolver : public PopulationEvolver
{
public:
	enum ProbType { Regular, Log, NegativeLog };

	SamplingEvolver(const std::shared_ptr<RandomNumberGenerator> &rng, ProbType t);
	~SamplingEvolver();

	ProbType getProbabilityType() const { return m_probType; }

	void setAnnealingExponent(double alpha) { m_alpha = alpha; } // to use prob^alpha
	double getAnnealingExponent() const { return m_alpha; }

	errut::bool_t check(const std::shared_ptr<Population> &population) override;

	// The idea here is to keep track of the individual with the highest
	// logprob/prob (fitness)
	const std::vector<std::shared_ptr<Individual>> &getBestIndividuals() const override { return m_bestIndividual; }
protected:
	void checkBest(Population &pop, size_t startIdx, size_t stopIdx);

	bool haveDoubleGenomes() const { return m_doubleGenomes; }
	size_t getDimension(const Genome &refGenome) const
	{
		return (m_doubleGenomes)?getNumberOfVariables<double>(refGenome):getNumberOfVariables<float>(refGenome);
	}

	// These are the actual individuals, not copies, for efficiency
	// TODO: change this?
	virtual void onSamples(const std::vector<std::shared_ptr<Individual>> &samples) { }

	template <class T>
	static inline size_t getNumberOfVariables(const Genome &g0)
	{
		assert(dynamic_cast<const VectorGenome<T> *>(&g0));
		return VectorGenome<T>::getSize(g0);
	}

	std::shared_ptr<RandomNumberGenerator> m_rng;
private:
	ProbType m_probType = Regular;
	bool m_doubleGenomes = false;
	double m_alpha = 1.0;

	std::vector<std::shared_ptr<Individual>> m_bestIndividual;
};

} // end namespace

