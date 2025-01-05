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
	: m_rng(rng), m_stepScales(stepScales), m_probType(t)
{
}

MetropolisHastingsEvolver::~MetropolisHastingsEvolver()
{
}

// TODO: move this to common base class!
bool_t MetropolisHastingsEvolver::check(const shared_ptr<Population> &population)
{
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

	// TODO: check m_stepScales size against individuals!

	return true;
}

template <class T>
inline size_t getNumVars(const Genome &g0)
{
	assert(dynamic_cast<const VectorGenome<T> *>(&g0));
	return VectorGenome<T>::getSize(g0);
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

	// TODO: use alpha!

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

	const Genome &refGenome = pop.individual(0)->genomeRef();
	size_t dim = (m_doubleGenomes)?getNumVars<double>(refGenome):getNumVars<float>(refGenome);

	auto appendNew = [this, &pop, N, dim, generation]()
	{
		if (m_doubleGenomes)
			appendNewIndividuals<double>(*m_rng, pop, m_stepScales, N, dim, generation);
		else
			appendNewIndividuals<float>(*m_rng, pop, m_stepScales, N, dim, generation);
	};

	// TODO: copied from GoodmanWeareEvolver, use common code!
	auto checkBest = [this, &pop](size_t startIdx, size_t stopIdx)
	{
		double bestFitness = (m_probType == NegativeLog)?numeric_limits<double>::infinity():-numeric_limits<double>::infinity();
		auto comp = (m_probType == NegativeLog)?[](double a, double b) { return a < b; }:[](double a, double b) { return a > b; };

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
	};

	if (generation == 0)
	{
		// On generation 0, pop size should be targetPopulationSize, new individuals should be
		// added based on step sizes, these first N can be reported as samples

		if (pop.size() != N)
			return "For the initial generation the number of indiviuals is expected to be " + to_string(N) + ", but is " + to_string(pop.size()); 

		checkBest(0, N);

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

		checkBest(N, 2*N);

		checkAccept(*m_rng, pop, m_probType, N, m_alpha);

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
