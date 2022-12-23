#include "differentialevolutionevolver.h"

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
	const shared_ptr<FitnessComparison> &fitComp, size_t objectiveNumber)
	: m_rng(rng),
	  m_mut(mut),
	  m_cross(cross),
	  m_CR(CR),
	  m_fitComp(fitComp),
	  m_objectiveNumber(objectiveNumber)
{
	m_mutationFactors = { 1.0, F, -F };
	m_mutationGenomes = { nullptr, nullptr, nullptr };
}

DifferentialEvolutionEvolver::~DifferentialEvolutionEvolver()
{

}

bool_t DifferentialEvolutionEvolver::check(const std::shared_ptr<Population> &population)
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
	return true;
}

bool_t DifferentialEvolutionEvolver::createNewPopulation(size_t generation, std::vector<std::shared_ptr<Population>> &populations,
                                                         const vector<size_t> &targetPopulationSizes)
{
	if (populations.size() != 1)
		return "DE evolver only works with one population";
	
	Population &pop = *(populations[0]);

	if (pop.size() < 4)
		return "Population size must be at least 4";

	if (targetPopulationSizes.size() != 1)
		return "Exactly one target population size should be mentioned";
	size_t targetPopulationSize = targetPopulationSizes[0];

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
			
			if (m_fitComp->isFitterThan(indNew.fitnessRef(), indBase.fitnessRef(), m_objectiveNumber))
				pop.individual(i) = pop.individual(i+targetPopulationSize);
		}

		pop.resize(targetPopulationSize);
	}
	else
		return "Unexpected population size: should be the target (" + to_string(targetPopulationSize) + ") or twice that, but is " + to_string(pop.size());

	// Check best

	for (auto &ind : pop.individuals())
	{
		if (m_bestIndividual.size() == 0)
			m_bestIndividual.push_back(ind->createCopy());
		else
		{
			if (m_fitComp->isFitterThan(ind->fitnessRef(), m_bestIndividual[0]->fitnessRef(), m_objectiveNumber))
				m_bestIndividual[0] = ind->createCopy();
		}	
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
