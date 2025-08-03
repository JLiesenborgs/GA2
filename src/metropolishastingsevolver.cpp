#define _USE_MATH_DEFINES
#include "metropolishastingsevolver.h"
#include "vectorgenomefitness.h"
#include <cmath>

using namespace errut;
using namespace std;

namespace eatk
{

MetropolisHastingsEvolver::MetropolisHastingsEvolver(const shared_ptr<RandomNumberGenerator> &rng,
										 const std::vector<double> &stepScales,
                                         ProbType t)
	: SamplingEvolver(rng, t), m_stepScales(stepScales)
{
}

MetropolisHastingsEvolver::~MetropolisHastingsEvolver()
{
}

bool_t MetropolisHastingsEvolver::check(const shared_ptr<Population> &population)
{
	bool_t r = SamplingEvolver::check(population);
	if (!r)
		return r;

	// check m_stepScales size against individuals
	for (const auto &ind : population->individuals())
	{
		size_t dim = getDimension(ind->genomeRef());
		if (dim != m_stepScales.size())
			return "Number of entries in the step scales (" + to_string(m_stepScales.size()) + ") does not equal number of parameters in individual (" + to_string(dim) + ")";
	}

	return true;
}

template<class T>
void appendNewIndividuals(RandomNumberGenerator &rng, Population &pop, const vector<double> &stepScales, size_t N, size_t dim, size_t generation)
{
	assert(pop.size() == N);
	assert(stepScales.size() == dim);

	const auto refInd = pop.individual(0); // Don't get a reference! We'll be resizing the internal array so this may become invalid
	const auto &refFit = refInd->fitnessRef();

	for (size_t i = 0 ; i < N ; i++)
	{
		const VectorGenome<T> &X = static_cast<const VectorGenome<T> &>(pop.individual(i)->genomeRef());
		const auto &Xv = X.getValues();

		// This will be the new genome
		auto Ygen = X.createCopy();
		auto Yfit = refFit.createCopy(false);
		VectorGenome<T> &Y = static_cast<VectorGenome<T> &>(*Ygen);
		auto &Yv = Y.getValues();

		assert(Xv.size() == dim && Yv.size() == dim);
		
		for (size_t j = 0 ; j < dim ; j++)
		{
			//Yv[j] = Xv[j] + (T)rng.getRandomDouble(-stepScales[j]*0.5, stepScales[j]*0.5);
			double z = rng.getRandomDouble();
			if (z != 0)
				Yv[j] = Xv[j] + (T)(-stepScales[j]/std::tan(z*M_PI)); // from meanwalker code, should be lorenzian
			else
				Yv[j] = Xv[j] + (T)((z-0.5)*stepScales[j]);
		}

		pop.append(refInd->createNew(Ygen, Yfit, generation));
	}
}

inline void checkAccept(RandomNumberGenerator &rng, Population &pop,
                        MetropolisHastingsEvolver::ProbType probType, size_t N, double alpha)
{
	assert(pop.size() == N*2);

	auto acceptTestLog = [&rng](double pNew, double pOld, double alpha)
	{
		if (pOld == -numeric_limits<double>::infinity())
			return true; // try to get away from bad region
		
		if (pNew == -numeric_limits<double>::infinity())
			return false; // try to stay away from bad region
		
		double logQ = alpha*(pNew-pOld);
		if (logQ >= 0) // means prob >= 1
			return true;

		double r = rng.getRandomDouble();
		if (r <= std::exp(logQ))
			return true;

		return false;
	};

	auto acceptTestNegLog = [&acceptTestLog](double pNew, double pOld, double alpha)
	{
		return acceptTestLog(-pNew, -pOld, alpha);
	};

	auto acceptTestNoLog = [&rng](double pNew, double pOld, double alpha)
	{
		if (pOld == 0)
			return true; // try to move away from bad region
		if (pNew == 0)
			return false; // try to stay away from bad region

		double q = std::pow(pNew/pOld, alpha);
		if (q >= 1.0)
			return true;

		double r = rng.getRandomDouble();
		if (r <= q)
			return true;
		return false;
	};

	auto checkAcceptLoop = [&pop, N, alpha](auto &shouldAcceptFunction)
	{
		for (size_t i = 0, j = N ; i < N ; i++, j++)
		{
			assert(i < pop.size() && j < pop.size());
			double pNew = pop.individual(j)->fitness()->getRealValue(0);
			double pOld = pop.individual(i)->fitness()->getRealValue(0);

			if (shouldAcceptFunction(pNew, pOld, alpha))
				pop.individual(i) = pop.individual(j); // accept, replace old with new
		}

		pop.individuals().resize(N); // Remove the ones that were added for comparison
	};

	if (probType == MetropolisHastingsEvolver::Regular)
		checkAcceptLoop(acceptTestNoLog);
	else if (probType == MetropolisHastingsEvolver::Log)
		checkAcceptLoop(acceptTestLog);
	else if (probType == MetropolisHastingsEvolver::NegativeLog)
		checkAcceptLoop(acceptTestNegLog);
	else
		throw runtime_error("Internal error: Invalid probability type");
}

bool_t MetropolisHastingsEvolver::createNewPopulation(size_t generation, shared_ptr<Population> &population, size_t targetPopulationSize)
{
	assert(population.get());
	size_t N = targetPopulationSize;

	Population &pop = *population;
	size_t curPopSize = pop.size();
	RandomNumberGenerator &rng = *m_rng;

	size_t dim = getDimension(pop.individual(0)->genomeRef());

	auto appendNew = [this, &pop, N, dim, generation]()
	{
		if (haveDoubleGenomes())
			appendNewIndividuals<double>(*m_rng, pop, m_stepScales, N, dim, generation);
		else
			appendNewIndividuals<float>(*m_rng, pop, m_stepScales, N, dim, generation);
	};

	if (generation == 0)
	{
		// On generation 0, pop size should be targetPopulationSize, new individuals should be
		// added based on step sizes, these first N can be reported as samples

		if (pop.size() != N)
			return "For the initial generation the number of indiviuals is expected to be " + to_string(N) + ", but is " + to_string(pop.size()); 

		checkBest(pop, 0, N);

		appendNew();
	}
	else
	{
		// For other generations, first individuals from first half should be compared to those of
		// second half, metropolis hastings should be used to check if the old or new is used as sample,
		// then the population should be resized to N. Again, based on step sizes new trial samples
		// can be added, the first N can still be reported as samples

		if (pop.size() != 2*N)
			return "For this generation the number of individuals is expected to be 2*" + to_string(N) + ", but is " + to_string(pop.size());

		checkBest(pop, N, 2*N);

		checkAccept(*m_rng, pop, getProbabilityType(), N, getAnnealingExponent());

		appendNew();
	}

	// Finally, report the samples

	assert(pop.size() == targetPopulationSize*2);

	m_samples.clear();
	for (size_t i = 0 ; i < N ; i++)
		m_samples.push_back(pop.individual(i));
	
	onSamples(m_samples);

	return true;
}

} // end namespace
