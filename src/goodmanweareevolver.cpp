#include "goodmanweareevolver.h"
#include "vectorgenomefitness.h"
#include <cmath>

using namespace errut;
using namespace std;

namespace eatk
{

GoodmanWeareEvolver::GoodmanWeareEvolver(const shared_ptr<RandomNumberGenerator> &rng,
                                         ProbType t, double a)
	: m_rng(rng), m_probType(t), m_a(a)
{
	m_aScale = std::sqrt(m_a) - 1.0/std::sqrt(m_a);
	m_aOffset = 1.0/std::sqrt(m_a);
}

GoodmanWeareEvolver::~GoodmanWeareEvolver()
{
}

bool_t GoodmanWeareEvolver::check(const shared_ptr<Population> &population)
{
	if (m_a <= 1)
		return "The 'a' parameter must be larger than one";

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

template <class T>
inline size_t getNumVars(const Genome &g0)
{
	assert(dynamic_cast<const VectorGenome<T> *>(&g0));
	return VectorGenome<T>::getSize(g0);
}

inline double pickZ(RandomNumberGenerator &rng, double aScale, double aOffset)
{
	double F = rng.getRandomDouble();
	double Z = F * aScale + aOffset;
	return Z*Z;
}

// Walks half of the walkers with the stretch move, based on the other half
template<class T>
inline void calculateNewPositions(size_t generation, Population &pop, size_t N, bool walkFirstHalf, RandomNumberGenerator &rng,
                                  double aScale, double aOffset, vector<double> &zValues)
{
	assert(pop.size() == N);
	assert(N%2 == 0);
	size_t N2 = N/2;

	size_t rndOffset = (walkFirstHalf)?N2:0;
	size_t kOffset = (walkFirstHalf)?0:N2;

	const auto refInd = pop.individual(0); // Don't get a reference! We'll be resizing the internal array so this may become invalid
	const auto &refFit = refInd->fitnessRef();

	zValues.clear();
	
	for (size_t i = 0, k = kOffset ; i < N2 ; i++, k++)
	{
		size_t j = ((size_t)rng.getRandomUint32())%N2 + rndOffset;
		assert(k < N && j < N);

		const VectorGenome<T> &Xj = static_cast<const VectorGenome<T> &>(pop.individual(j)->genomeRef());
		const VectorGenome<T> &Xk = static_cast<const VectorGenome<T> &>(pop.individual(k)->genomeRef());

		// Initialize a new genome
		auto Ygen = Xk.createCopy();
		auto Yfit = refFit.createCopy(false);
		VectorGenome<T> &Y = static_cast<VectorGenome<T> &>(*Ygen);

		const auto &Xjv = Xj.getValues();
		const auto &Xkv = Xk.getValues();
		auto &Yv = Y.getValues();

		double z = pickZ(rng, aScale, aOffset);

		for (size_t u = 0 ; u < Yv.size() ; u++)
			Yv[u] = Xjv[u] + ((T)z)*(Xkv[u] - Xjv[u]);

		pop.append(refInd->createNew(Ygen, Yfit, generation));

		// We'll need to calculate z^(dim-1) * probY / probXk to see if we can accept this
		// step (or log), so let's store this z
		zValues.push_back(z);
	}

	assert(zValues.size() == N2);
	assert(pop.size() == N2*3);
}

inline void checkAccept(RandomNumberGenerator &rng, Population &pop, bool compareFirstHalf,
                        GoodmanWeareEvolver::ProbType probType, size_t N, size_t dim, const vector<double> &zValues,
						double alpha)
{
	assert(pop.size() == (N*3)/2);
	assert(N%2 == 0);
	size_t N2 = N/2;
	size_t c = (compareFirstHalf)?0:N2;

	assert(zValues.size() == N2);

	auto acceptTestLog = [&rng, dim](double pNew, double pOld, double z, double alpha)
	{
		if (pOld == -numeric_limits<double>::infinity())
			return true; // try to get away from bad region
		
		if (pNew == -numeric_limits<double>::infinity())
			return false; // try to stay away from bad region
		
		double logQ = (double)(dim-1) * std::log(z) + alpha*(pNew - pOld);
		if (logQ >= 0) // means prob >= 1
			return true;

		double r = rng.getRandomDouble();
		if (r <= std::exp(logQ))
			return true;

		return false;
	};

	auto acceptTestNegLog = [&acceptTestLog](double pNew, double pOld, double z, double alpha)
	{
		return acceptTestLog(-pNew, -pOld, z, alpha);
	};

	auto acceptTestNoLog = [&rng, dim](double pNew, double pOld, double z, double alpha)
	{
		if (pOld == 0)
			return true; // try to move away from bad region
		if (pNew == 0)
			return false; // try to stay away from bad region

		double q = std::pow(z, dim-1) * std::pow(pNew/pOld, alpha); // accept prob
		
		if (q >= 1.0)
			return true;

		double r = rng.getRandomDouble();
		if (r <= q)
			return true;
		return false;
	};

	auto checkAcceptLoop = [&pop, N, &c, &zValues, alpha](auto &shouldAcceptFunction)
	{
		for (size_t n = N, zIdx = 0 ; n < pop.size() ; n++, c++, zIdx++)
		{
			assert(n < pop.size() && c < pop.size() && zIdx < zValues.size());
			double pNew = pop.individual(n)->fitness()->getRealValue(0);
			double pOld = pop.individual(c)->fitness()->getRealValue(0);
			double z = zValues[zIdx];

			if (shouldAcceptFunction(pNew, pOld, z, alpha))
				pop.individual(c) = pop.individual(n); // accept, replace old with new
		}

		pop.individuals().resize(N); // Remove the ones that were added for comparison
	};

	if (probType == GoodmanWeareEvolver::Regular)
		checkAcceptLoop(acceptTestNoLog);
	else if (probType == GoodmanWeareEvolver::Log)
		checkAcceptLoop(acceptTestLog);
	else if (probType == GoodmanWeareEvolver::NegativeLog)
		checkAcceptLoop(acceptTestNegLog);
	else
		throw runtime_error("Internal error: Invalid probability type");
}

bool_t GoodmanWeareEvolver::createNewPopulation(size_t generation, shared_ptr<Population> &population, size_t targetPopulationSize)
{
	assert(population.get());
	size_t N = targetPopulationSize;

	if (N%2 != 0)
		return "Target population size must be a multiple of two";

	Population &pop = *population;
	size_t curPopSize = pop.size();
	RandomNumberGenerator &rng = *m_rng;

	auto calculateNewPositions_General = [this, &pop, N, &rng, generation](bool walkFirstHalf) 
	{
		if (m_doubleGenomes)
			calculateNewPositions<double>(generation, pop, N, walkFirstHalf, rng, m_aScale, m_aOffset, m_zValues);
		else
			calculateNewPositions<float>(generation, pop, N, walkFirstHalf, rng, m_aScale, m_aOffset, m_zValues);
	};

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

	const Genome &refGenome = pop.individual(0)->genomeRef();
	size_t dim = (m_doubleGenomes)?getNumVars<double>(refGenome):getNumVars<float>(refGenome);

	if (curPopSize == targetPopulationSize)
	{
		if (generation != 0)
			return "Expecting generation to be zero, but is " + to_string(generation);

		if (N <= dim)
			return "The number of walkers must be at least one more as the dimension of the problem";

		checkBest(0, N);
		calculateNewPositions_General(true);

		// TODO: check all feasible?
	}
	else if (curPopSize == (targetPopulationSize*3)/2)
	{
		assert(generation > 0);

		// We've just calculated fitnesses in the upper part, see if we need
		// to update the best one
		checkBest(N, curPopSize);

		// On generation 1 we need to compare the new fitnesses to the first half
		bool compareFirstHalf = (generation%2 == 1);
		bool walkFirstHalf = (generation%2 == 0);

		// Check if we should accept the new genomes, then resize pop to N
		checkAccept(rng, pop, compareFirstHalf, m_probType, N, dim, m_zValues, m_alpha);

		// Calculate new walk		
		calculateNewPositions_General(walkFirstHalf);
	}
	else
		return "Unexpected population size: should be the target (" + to_string(targetPopulationSize) + ") or 3/2 that, but is " + to_string(curPopSize);

	// After every two generations we have a new set of samples, report those
	if (generation%2 == 0)
	{
		m_samples.clear();
		for (size_t i = 0 ; i < N ; i++)
			m_samples.push_back(pop.individual(i));
		
		onSamples(m_samples);
	}

	return true;
}

} // end namespace
