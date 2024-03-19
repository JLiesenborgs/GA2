#include "jadeevolver.h"
#include <limits>
#include <random>
#include <algorithm>

using namespace std;
using namespace errut;

namespace eatk
{

JADEEvolver::JADEEvolver(const shared_ptr<RandomNumberGenerator> &rng,
		const shared_ptr<DifferentialEvolutionMutation> &mut,
		const shared_ptr<DifferentialEvolutionCrossover> &cross,
		const shared_ptr<FitnessComparison> &fitComp, int objectiveNumber,
		double p, double c,
		bool useArchive,
		double initMuF,
		double initMuCR,
		size_t numObjectives,
		const shared_ptr<NonDominatedSetCreator> &ndCreator
		)
	: m_rng(rng),
	  m_mut(mut),
	  m_cross(cross),
	  m_fitComp(fitComp),
	  m_objectiveNumber(objectiveNumber),
	  m_p(p), m_c(c), m_useArchive(useArchive),
	  m_initMuF(initMuF),m_initMuCR(initMuCR),
	  m_numObjectives(numObjectives),
	  m_ndCreator(ndCreator)
{
	m_mutationFactors = { 1.0, 0, 0, 0, 0 };
	m_mutationGenomes = { nullptr, nullptr, nullptr, nullptr, nullptr };
}

JADEEvolver::~JADEEvolver()
{
}

bool_t JADEEvolver::check(const shared_ptr<Population> &population)
{
	bool_t r;
	for (auto &ind : population->individuals())
	{
		if (!(r = m_mut->check(ind->genomeRef())))
			return "Can't handle genome type in DE mutation: " + r.getErrorString();
		if (!(r = m_cross->check(ind->genomeRef())))
			return "Can't handle genome type in DE crossover: " + r.getErrorString();
		if (!(r = m_fitComp->check(ind->fitnessRef())))
			return "Can't handle fitness type: " + r.getErrorString();
	}

	if (m_c < 0 || m_c > 1)
		return "Invalid value for 'c', should lie between 0 and 1 but is " + to_string(m_c);
	if (m_p < 0 || m_p > 1)
		return "Invalid value for 'p', should lie between 0 and 1 but is " + to_string(m_p);

	if (m_numObjectives == 0)
		return "Number of objectives must be at least one";
	if (m_objectiveNumber >= 0 && (size_t)m_objectiveNumber >= m_numObjectives)
		return "Objective number is not compatible with number of objectives";
	if (m_objectiveNumber < 0)
	{
		if (!m_ndCreator.get())
			return "Multi-objective mode was specified, but no non-dominated set creator was specified";
	}

	return true;
}

inline double meanA(const std::vector<double> &values)
{
	if (values.size() == 0)
		return 0;

	double sum = 0;
	for (double x: values)
		sum += x;
	sum /= values.size();
	return sum;
}

inline double meanL(const std::vector<double> &values)
{
	if (values.size() == 0)
		return 0;

	double sum = 0, sum2 = 0;
	for (double x: values)
	{
		sum2 += x*x;
		sum += x;
	}

	if (sum == 0) // Guard against NaN
		return 0;

	sum2 /= sum;

#ifndef NDEBUG
	if (isnan(sum2))
	{
		cerr << "DEBUG: meanL [ ";
		for (auto x : values)
			cerr << " " << x;
		cerr << " ]\n";
		cerr << "  sum = " << sum << endl;
		cerr << "  sum2 = " << sum2 << endl;
		assert(!isnan(sum2));
	}
#endif // NDEBUG

	return sum2;
}

template<class Distribution>
inline double chooseTruncatedDistribution(RandomNumberGenerator &rng, double mu, double sigma, double x0, double x1)
{
	struct G
	{
		typedef uint64_t result_type;

		G(RandomNumberGenerator &_rng) : rng(_rng) { }
		uint64_t operator()() { return (uint64_t)rng.getRandomUint32(); }
		constexpr static uint64_t min() { return 0; }
		constexpr static uint64_t max() { return 0x100000000; }

		RandomNumberGenerator &rng;
	};

	G gen(rng);
	Distribution dist(mu, sigma);
	double x = dist(gen);
	if (x < x0)
		x = x0;
	if (x > x1)
		x = x1;
	return x;
}

bool_t JADEEvolver::createNewPopulation(size_t generation, shared_ptr<Population> &population,
										size_t targetPopulationSize)
{
#ifndef NDEBUG
	auto printVector = [this](const vector<double> &vals, const string &name)
	{
		cerr << "DEBUG: " << name << " =";
		for (auto x : vals)
			cerr << " " << x;
		cerr << endl;
	};

	auto printVectors = [this, &printVector]()
	{
		printVector(m_mutationFactors, "m_mutationFactors");
		printVector(m_CRi, "m_CRi");
		printVector(m_Fi, "m_Fi");
		printVector(m_SCR, "m_SCR");
		printVector(m_SF, "m_SF");
	};

	//cerr << "DEBUG: start vectors" << endl;
	//printVectors();
#endif // NDEBUG

	Population &pop = *population;

	auto isFitterThan = getDominatesFunction(m_fitComp, m_objectiveNumber, m_numObjectives);

	if (pop.size() < 4)
		return "Population size must be at least 4";

	if (pop.size() == targetPopulationSize)
	{
		// This is on the first iteration
		if (generation != 0)
			return "Expecting generation to be zero, but is " + to_string(generation);

		m_muF = m_initMuF;
		m_muCR = m_initMuCR;
		m_Fi.assign(targetPopulationSize, numeric_limits<double>::quiet_NaN());
		m_CRi.assign(targetPopulationSize, numeric_limits<double>::quiet_NaN());
		m_archive.clear();
	}
	else if (pop.size() == targetPopulationSize*2)
	{
		// The population is twice the size, which means we need to
		// compare each individual in the first half to the corresponding
		// one in the second half and keep the best one

		m_SCR.clear();
		m_SF.clear();

		assert(m_CRi.size() == targetPopulationSize);
		assert(m_Fi.size() == targetPopulationSize);

		for (size_t i = 0 ; i < targetPopulationSize ; i++)
		{
			Individual &indBase = *(pop.individual(i));
			Individual &indNew = *(pop.individual(i+targetPopulationSize));
			
			if (isFitterThan(indNew.fitnessRef(), indBase.fitnessRef()))
			{
				if (m_useArchive)
					m_archive.push_back(indBase.genomeRef().createCopy()); // add to archive before overwriting location
				pop.individual(i) = pop.individual(i+targetPopulationSize);

				m_SF.push_back(m_Fi[i]);
				m_SCR.push_back(m_CRi[i]);
			}
		}

		// Keep only relevant part of population
		pop.resize(targetPopulationSize);

		// Make sure the archive doesn't grow to big
		trimArchive(targetPopulationSize);

		// Update muCR and muF values
		m_muCR = (1.0-m_c)*m_muCR + m_c*meanA(m_SCR);
		m_muF = (1.0-m_c)*m_muF + m_c*meanL(m_SF);
	}
	else
		return "Unexpected population size: should be the target (" + to_string(targetPopulationSize) + ") or twice that, but is " + to_string(pop.size());

	onMutationCrossoverSettings(m_muF, m_muCR);

	if (m_numObjectives == 1) // single objective
	{
		// Sort the population, we need this to select the p fraction of best
		auto comp = [this](auto &i1, auto &i2)
		{
			return m_fitComp->isFitterThan(i1->fitnessRef(), i2->fitnessRef(), m_objectiveNumber);
		};
		sort(pop.individuals().begin(), pop.individuals().end(), comp);

		// Keep the best
		const auto &firstAfterSort = pop.individuals()[0];
		if (m_bestIndividual.size() == 0)
			m_bestIndividual.push_back(firstAfterSort->createCopy());
		else
		{
			if (comp(firstAfterSort, m_bestIndividual[0]))
				m_bestIndividual[0] = firstAfterSort->createCopy();
		}
	}
	else // multi-objective
	{
		// Use an ND sort
		bool_t r;
		r = m_ndCreator->calculateAllNDSets(pop.individuals());
		if (!r)
			return "Can't create ND sets: " + r.getErrorString();

		size_t oldPopsize = pop.size();
		size_t numSets = m_ndCreator->getNumberOfSets();

		// Clear the population and insert the individuals again according to the ND sets
		pop.clear();
		for (size_t i = 0 ; i < numSets ; i++)
		{
			for (auto &ind : m_ndCreator->getSet(i))
				pop.append(ind);
		}

		if (oldPopsize != pop.size())
			return "Internal error: population size differs from expected";

		// Keep the best set
		m_bestIndividual.clear();
		for (auto &ind : m_ndCreator->getSet(0))
			m_bestIndividual.push_back(ind->createCopy());
	}

	// Do mutation/crossover

	assert(pop.size() == targetPopulationSize);

	size_t archiveSize = m_archive.size();
	assert(archiveSize <= targetPopulationSize);
	
	auto pickIndex = [this](size_t num) { return ((size_t)m_rng->getRandomUint32())%(num); };
	auto adjust = [](size_t &idx, size_t refIdx) { if (idx >= refIdx) idx++; };

	auto pickRandomIndices = [targetPopulationSize,archiveSize,pickIndex,adjust](size_t i, size_t &r1, size_t &r2)
	{
		r1 = pickIndex(targetPopulationSize-1);
		r2 = pickIndex(targetPopulationSize+archiveSize-2);

		adjust(r2, r1);

		adjust(r1, i);
		adjust(r2, i);

		assert(r1 != i);
		assert(r2 != i && r2 != r1);
	};

	Individual &refIndividual = *(pop.individual(0));
	Fitness &refFitness = refIndividual.fitnessRef();

	for (size_t i = 0 ; i < targetPopulationSize ; i++) // Note: not using pop.size() since we'll be changing the size!
	{
		m_CRi[i] = chooseTruncatedDistribution<normal_distribution<>>(*m_rng, m_muCR, 0.1, 0, 1);
		m_Fi[i] = chooseTruncatedDistribution<cauchy_distribution<>>(*m_rng, m_muF, 0.1, 0, 1);

		size_t bestp = (size_t)(m_rng->getRandomDouble(0, m_p) * targetPopulationSize);
		assert(bestp < targetPopulationSize);
		if (bestp >= targetPopulationSize)
			bestp = targetPopulationSize-1;

		m_mutationFactors[0] = 1.0;
		m_mutationFactors[1] = m_Fi[i];
		m_mutationFactors[2] = -m_Fi[i];
		m_mutationFactors[3] = m_Fi[i];
		m_mutationFactors[4] = -m_Fi[i];

		size_t r1, r2;
		pickRandomIndices(i, r1, r2);

		m_mutationGenomes[0] = pop.individual(i)->genomePtr();
		m_mutationGenomes[1] = pop.individual(bestp)->genomePtr();
		m_mutationGenomes[2] = m_mutationGenomes[0];
		m_mutationGenomes[3] = pop.individual(r1)->genomePtr();
		if (r2 < targetPopulationSize)
			m_mutationGenomes[4] = pop.individual(r2)->genomePtr();
		else // from archive
		{
			size_t archIdx = r2-targetPopulationSize;
			assert(archIdx < archiveSize);
			m_mutationGenomes[4] = m_archive[archIdx].get();
			// cout << "Choosing " << archIdx << " from archive" << endl;
		}

		//printVector(m_mutationFactors, "m_mutationFactors");

		shared_ptr<Genome> newGenome = m_mut->mutate(m_mutationGenomes, m_mutationFactors);
		if (!newGenome.get())
			return "Unable to create new genome from three reference genomes";

		// Crossover part

		bool_t r = m_cross->crossover(m_CRi[i], *newGenome, pop.individual(i)->genomeRef());
		if (!r)
			return "Unable to crossover genomes: " + r.getErrorString();

		shared_ptr<Fitness> newFitness = refFitness.createCopy();
		newFitness->setCalculated(false);

		pop.append(refIndividual.createNew(newGenome, newFitness, generation)); // TODO: is generation+1 better here?
	}

	//cerr << "DEBUG: end vectors" << endl;
	//printVectors();

	return true;
}

void JADEEvolver::trimArchive(size_t targetPopulationSize)
{
	while (m_archive.size() > targetPopulationSize)
	{
		size_t idx = (size_t)m_rng->getRandomUint32() % (m_archive.size());
		
		// Remove that element
		size_t lastIdx = m_archive.size()-1;
		if (lastIdx != idx)
			swap(m_archive[idx], m_archive[lastIdx]); // move selected element to end of array
		m_archive.resize(lastIdx);
	}
}

}
