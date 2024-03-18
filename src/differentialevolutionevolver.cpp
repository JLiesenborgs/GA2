#include "differentialevolutionevolver.h"
#include "nondominatedsetcreator.h"
#include <functional>

using namespace errut;
using namespace std;

namespace eatk
{

DifferentialEvolutionEvolver::DifferentialEvolutionEvolver(
	const shared_ptr<RandomNumberGenerator> &rng,
	const shared_ptr<DifferentialEvolutionMutation> &mut,
	double F,
	const shared_ptr<DifferentialEvolutionCrossover> &cross,
	double CR,
	const shared_ptr<FitnessComparison> &fitComp,
	int objectiveNumber, size_t numObjectives,
	const shared_ptr<NonDominatedSetCreator> &ndCreator)
	: m_rng(rng),
	  m_mut(mut),
	  m_cross(cross),
	  m_CR(CR),
	  m_fitComp(fitComp),
	  m_objectiveNumber(objectiveNumber),
	  m_numObjectives(numObjectives),
	  m_ndCreator(ndCreator)
{
	m_mutationFactors = { 1.0, F, -F };
	m_mutationGenomes = { nullptr, nullptr, nullptr };
}

DifferentialEvolutionEvolver::~DifferentialEvolutionEvolver()
{
}

bool_t DifferentialEvolutionEvolver::check(const std::shared_ptr<Population> &population)
{
	if (m_numObjectives == 0)
		return "Number of objectives must be at least one";
	if (m_objectiveNumber >= 0 && (size_t)m_objectiveNumber >= m_numObjectives)
		return "Objective number is not compatible with number of objectives";
	if (m_objectiveNumber < 0)
	{
		if (!m_ndCreator.get())
			return "Multi-objective mode was specified, but no non-dominated set creator was specified";
	}

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
	return true;
}

bool_t DifferentialEvolutionEvolver::createNewPopulation(size_t generation, std::shared_ptr<Population> &population,
                                                         size_t targetPopulationSize)
{
	Population &pop = *population;

	std::function<bool(const Fitness &f1, const Fitness &f2)> isFitterThan;
	if (m_objectiveNumber >= 0) // single objective
	{
		isFitterThan = [this](const Fitness &f1, const Fitness &f2)
		{
			return m_fitComp->isFitterThan(f1, f2, m_objectiveNumber);
		};
	}
	else
	{
		// multi-objective, use dominance
		isFitterThan = [this](const Fitness &f1, const Fitness &f2)
		{
			size_t betterOrEqualCount = 0;
			size_t betterCount = 0;
			for (size_t i = 0 ; i < m_numObjectives ; i++)
			{
				if (m_fitComp->isFitterThan(f1, f2, i))
				{
					betterCount++;
					betterOrEqualCount++;
				}
				else // f1 not strictly better than f2 for i
				{
					if (!m_fitComp->isFitterThan(f2, f1, i)) // then they must have equal fitness
					{
						betterOrEqualCount++;
					}
					else
					{
						// We can never get betterOrEqualCount == m_numObjectives
						return false;
					}
				}
			}
			// if we got here, then betterOrEqualCount == m_numObjectives
			if (betterCount > 0)
				return true;
			return false;
		};
	}

	if (pop.size() < 4)
		return "Population size must be at least 4";

	if (pop.size() == targetPopulationSize)
	{
		// Nothing to do, this is on the first iteration, we've
		// just created a population of this size
		if (generation != 0)
			return "Expecting generation to be zero, but is " + to_string(generation);
	}
	else if (pop.size() == targetPopulationSize*2)
	{
		// The population is twice the size, which means we need to
		// compare each individual in the first half to the corresponding
		// one in the second half and keep the best one

		for (size_t i = 0 ; i < targetPopulationSize ; i++)
		{
			Individual &indBase = *(pop.individual(i));
			Individual &indNew = *(pop.individual(i+targetPopulationSize));
			
			if (isFitterThan(indNew.fitnessRef(), indBase.fitnessRef()))
				pop.individual(i) = pop.individual(i+targetPopulationSize);
		}

		pop.resize(targetPopulationSize);
	}
	else
		return "Unexpected population size: should be the target (" + to_string(targetPopulationSize) + ") or twice that, but is " + to_string(pop.size());

	// Check best

	if (m_objectiveNumber >= 0) // single-objective
	{
		for (auto &ind : pop.individuals())
		{
			if (m_bestIndividuals.size() == 0)
				m_bestIndividuals.push_back(ind->createCopy());
			else
			{
				if (m_fitComp->isFitterThan(ind->fitnessRef(), m_bestIndividuals[0]->fitnessRef(), m_objectiveNumber))
					m_bestIndividuals[0] = ind->createCopy();
			}	
		}
	}
	else // multi-objective
	{
		bool_t r;

		m_bestIndividuals.clear();
		vector<shared_ptr<Individual>> remaining_notused;
		
		if (!(r = m_ndCreator->calculateNonDomitatedSet(pop.individuals(), m_bestIndividuals, remaining_notused)))
			return "Can't create non-dominated set: " + r.getErrorString();

		// Actually copy the individuals
		for (auto &x : m_bestIndividuals)
			x = x->createCopy();
	}

	// Do mutation/crossover

	assert(pop.size() == targetPopulationSize);

	auto pickIndex = [this](size_t num) { return ((size_t)m_rng->getRandomUint32())%(num); };
	auto adjust = [](size_t &idx, size_t refIdx) { if (idx >= refIdx) idx++; };

	auto pickRandomIndices = [targetPopulationSize,pickIndex,adjust](size_t i, size_t &r1, size_t &r2, size_t &r3)
	{
		r1 = pickIndex(targetPopulationSize-1);
		r2 = pickIndex(targetPopulationSize-2);
		r3 = pickIndex(targetPopulationSize-3);
		adjust(r3, r2);
		
		adjust(r2, r1);
		adjust(r3, r1);

		adjust(r1, i);
		adjust(r2, i);
		adjust(r3, i);
	
		assert(r1 != i && r2 != i && r3 != i);
		assert(r2 != r1 && r3 != r1);
		assert(r3 != r2);
	};

	Individual &refIndividual = *(pop.individual(0));
	Fitness &refFitness = refIndividual.fitnessRef();

	for (size_t i = 0 ; i < targetPopulationSize ; i++) // Note: not using pop.size() since we'll be changing the size!
	{
		size_t r1, r2, r3;
		pickRandomIndices(i, r1, r2, r3);

		m_mutationGenomes[0] = pop.individual(r1)->genomePtr();
		m_mutationGenomes[1] = pop.individual(r2)->genomePtr();
		m_mutationGenomes[2] = pop.individual(r3)->genomePtr();

		shared_ptr<Genome> newGenome = m_mut->mutate(m_mutationGenomes, m_mutationFactors);
		if (!newGenome.get())
			return "Unable to create new genome from three reference genomes";

		//cerr << "Orig (" << i << "): " << pop.individual(i)->genome()->toString() << endl;
		//cerr << "  r1 (" << r1 << "): " << pop.individual(r1)->genome()->toString() << endl;
		//cerr << "  r2 (" << r2 << "): " << pop.individual(r2)->genome()->toString() << endl;
		//cerr << "  r3 (" << r3 << "): " << pop.individual(r3)->genome()->toString() << endl;
		//cerr << "  before crossover: " << newGenome->toString() << endl;

		bool_t r = m_cross->crossover(m_CR, *newGenome, pop.individual(i)->genomeRef());
		if (!r)
			return "Unable to crossover genomes: " + r.getErrorString();

		//cerr << "  after crossover: " << newGenome->toString() << endl;

		shared_ptr<Fitness> newFitness = refFitness.createCopy();
		newFitness->setCalculated(false);

		pop.append(refIndividual.createNew(newGenome, newFitness, generation)); // TODO: is generation+1 better here?
	}

	return true;
}

}
