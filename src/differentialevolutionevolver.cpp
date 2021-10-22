#include "differentialevolutionevolver.h"

using namespace errut;
using namespace std;

namespace eatk
{

DifferentialEvolutionEvolver::DifferentialEvolutionEvolver(
	const shared_ptr<RandomNumberGenerator> &rng,
	const shared_ptr<DifferentialEvolutionMutation> &mut,
	const shared_ptr<DifferentialEvolutionCrossover> &cross,
	const shared_ptr<FitnessComparison> &fitComp, size_t objectiveNumber)
	: m_rng(rng),
	  m_mut(mut),
	  m_cross(cross),
	  m_fitComp(fitComp),
	  m_objectiveNumber(objectiveNumber)
{

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

bool_t DifferentialEvolutionEvolver::createNewPopulation(size_t generation, std::vector<std::shared_ptr<Population>> &populations, size_t targetPopulationSize)
{
	if (populations.size() != 1)
		return "DE evolver only works with one population";
	
	Population &pop = *(populations[0]);

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

		shared_ptr<Genome> newGenome = m_mut->mutate(pop.individual(r1)->genomeRef(),
		                                 pop.individual(r2)->genomeRef(),
										 pop.individual(r3)->genomeRef());
		if (!newGenome.get())
			return "Unable to create new genome from three reference genomes";

		bool_t r = m_cross->crossover(*newGenome, pop.individual(i)->genomeRef());
		if (!r)
			return "Unable to crossover genomes: " + r.getErrorString();

		shared_ptr<Fitness> newFitness = refFitness.createCopy();
		newFitness->setCalculated(false);

		pop.append(refIndividual.createNew(newGenome, newFitness, generation)); // TODO: is generation+1 better here?
	}

	return true;
}

}