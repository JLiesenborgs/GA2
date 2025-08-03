#include "samplingevolver.h"
#include <cmath>

using namespace errut;
using namespace std;

namespace eatk
{

SamplingEvolver::SamplingEvolver(const std::shared_ptr<RandomNumberGenerator> &rng, ProbType t)
	: m_rng(rng), m_probType(t)
{
}

SamplingEvolver::~SamplingEvolver()
{
}

bool_t SamplingEvolver::check(const shared_ptr<Population> &population)
{
	if (population->size() == 0)
		return "Empty population";

	// check genome and fitness type, float or double
	for (auto &i : population->individuals())
	{
		if (dynamic_cast<const FloatVectorGenome *>(i->genomePtr()))
			m_doubleGenomes = false;
		else
		{
			if (dynamic_cast<const DoubleVectorGenome *>(i->genomePtr()))
				m_doubleGenomes = true;
			else
				return "Genome type should be either a FloatVectorGenome or DoubleVectorGenome, but is " + string(typeid(i->genomeRef()).name());
		}

		if (!i->fitnessRef().hasRealValues())
			return "Fitness should be interpretable as a real value (is prob/logprob)";
	}
	return true;
}

void SamplingEvolver::checkBest(Population &pop, size_t startIdx, size_t stopIdx)
{
	double bestFitness = (m_probType == NegativeLog)?numeric_limits<double>::infinity():-numeric_limits<double>::infinity();
#ifndef _WIN32
	auto comp = (m_probType == NegativeLog)?[](double a, double b) { return a < b; }:[](double a, double b) { return a > b; };
#else
	auto comp = [this](double a, double b) { if (m_probType == NegativeLog) return a < b; return a > b; };
#endif // _WIN32
	if (m_bestIndividual.size() > 0)
	{
		assert(m_bestIndividual.size() == 1);
		bestFitness = m_bestIndividual[0]->fitness()->getRealValue(0); // TODO: make objective configurable?
	}

	size_t bestIdx = numeric_limits<size_t>::max();
	for (size_t i = startIdx ; i < stopIdx ; i++)
	{
		assert(i < pop.size());
		double fit = pop.individual(i)->fitness()->getRealValue(0); // TODO: make objective configurable
		if (comp(fit, bestFitness))
		{
			bestIdx = i;
			bestFitness = fit;
		}
	}

	if (bestIdx < stopIdx) // We found a better one
	{
		m_bestIndividual.clear();
		m_bestIndividual.push_back(pop.individual(bestIdx)->createCopy());
	}
}

} // end namespace
