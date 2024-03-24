#include "evolutionaryalgorithm.h"

using namespace std;
using namespace errut;

namespace eatk
{

EvolutionaryAlgorithm::EvolutionaryAlgorithm()
	: m_generation(0)
{
}

EvolutionaryAlgorithm::~EvolutionaryAlgorithm()
{
}

inline bool_t initializePopulation(size_t popSize, shared_ptr<Population> &population, IndividualCreation &gfc,
                                   shared_ptr<Individual> &refIndividual, shared_ptr<Fitness> &refFitness,
                                   size_t generation)
{
	for (size_t i = 0 ; i < popSize ; i++)
	{
		auto g = gfc.createInitializedGenome();
		if (!g.get())
			return "Unable to create an inialized genome";

		auto f = refFitness->createCopy(false);
		population->append(refIndividual->createNew(g, f, generation));
	}
	return true;
}

// Note that the population size does not need to be constant throughout the loop,
// more could arise so that their fitness is calculated. This is why the population
// size is passed on to the populationcrossover
bool_t EvolutionaryAlgorithm::run(IndividualCreation &gfc,
							 PopulationEvolver &evolver, // We really do need this, it keeps track of the best
							 PopulationFitnessCalculation &fitnessCalc,
							 StopCriterion &stopCriterion,
							 size_t popSize,
							 size_t minPopulationSize,
							 size_t maxPopulationSize)
{
	bool_t r;
	m_population = make_shared<Population>();
	m_populations = { m_population };

	auto newPopulation = make_shared<Population>();
	auto refFitness = gfc.createEmptyFitness();
	auto refIndividual = gfc.createReferenceIndividual();

	if (!refFitness.get())
		return "Unable to create a reference fitness object";

	if (!refIndividual.get())
		return "Unable to create reference individual";

	if (popSize == 0)
		return "No population size specified";
	
	if (maxPopulationSize == 0)
		maxPopulationSize = popSize;

	m_generation = 0;

	if (!(r = initializePopulation(popSize, m_population, gfc, refIndividual, refFitness, m_generation)))
		return r;

	auto beforeFitnessCalculatedCallback = [this]() -> bool_t
	{
		bool_t r;
		if (!(r = onBeforeFitnessCalculation(m_generation, m_population)))
			return "Error inspecting population before fitness calculation in generation " + to_string(m_generation) + ": " + r.getErrorString();
		return true;
	};

	if (!(r = beforeFitnessCalculatedCallback()))
		return r;

	if (!(r = fitnessCalc.calculatePopulationFitness({m_population})))
		return "Error calculating fitness: " + r.getErrorString();

	auto fitnessCalculatedCallback = [this]() -> bool_t
	{
		bool_t r;
		if (!(r = onFitnessCalculated(m_generation, m_population)))
			return "Error inspecting population after fitness calculation in generation " + to_string(m_generation) + ": " + r.getErrorString();
		return true;
	};

	if (!(r = fitnessCalculatedCallback()))
		return r;

	while (true)
	{		
		if (m_generation == 0)
		{
			if (!(r = evolver.check(m_population)))
				return "Error in population evolver check: " + r.getErrorString();
		}
		if (!(r = evolver.createNewPopulation(m_generation, m_population, popSize)))
			return "Error creating new population: " + r.getErrorString();

		const size_t curPopSize = m_population->size();
		if (curPopSize > maxPopulationSize)
			return "Population size (" + to_string(curPopSize) + ") exceeds maximum (" + to_string(maxPopulationSize) + ")";
		if (curPopSize < minPopulationSize)
			return "Population size (" + to_string(curPopSize) + ") is less than minimum (" + to_string(minPopulationSize) + ")";

		m_generation++; // At this point we can call it the new generation

		if (!(r = beforeFitnessCalculatedCallback()))
			return r;

		if (!(r = fitnessCalc.calculatePopulationFitness({m_population})))
			return "Error calculating fitness: " + r.getErrorString();

		if (!(r = fitnessCalculatedCallback()))
			return r;

		bool shouldStop = false;

		if (!(r = stopCriterion.analyze(evolver, m_generation, shouldStop)))
			return "Error in termination check: " + r.getErrorString();
		if (shouldStop)
			break;
	}

	if (!(r = onAlgorithmDone(m_generation, evolver.getBestIndividuals())))
		return "Error inspecting best individuals upon algorithm end: " + r.getErrorString();

	return true;
}

bool_t EvolutionaryAlgorithm::run(IndividualCreation &gfc,
							 PopulationEvolver &evolver, // We really do need this, it keeps track of the best
							 PopulationFitnessCalculation &fitnessCalc,
							 StopCriterion &stopCriterion,
			                 MigrationStrategy &migrationStrategy,
			                 const vector<size_t> &popSizes,
			                 const vector<size_t> &minPopSizes,
			                 const vector<size_t> &maxPopSizes)
{
	bool_t r;

	auto refFitness = gfc.createEmptyFitness();
	auto refIndividual = gfc.createReferenceIndividual();

	if (!refFitness.get())
		return "Unable to create a reference fitness object";

	if (!refIndividual.get())
		return "Unable to create reference individual";

	if (popSizes.size() < 2)
		return "At least two population sizes should be used in this version, but " + to_string(popSizes.size()) + " was/were specified";
	for (auto popSize : popSizes)
		if (popSize == 0)
			return "At least one of the population sizes is zero";
	
	vector<size_t> minPopulationSizes, maxPopulationSizes;
	if (minPopSizes.empty())
		for (auto _ : popSizes)
			minPopulationSizes.push_back(0);

	if (maxPopSizes.empty())
		for (auto popSize : popSizes)
			maxPopulationSizes.push_back(popSize);
			
	if (minPopulationSizes.size() != popSizes.size())
		return "The number of minimum population sizes does not match the number of populations";
	if (maxPopulationSizes.size() != popSizes.size())
		return "The number of maximum population sizes does not match the number of populations";

	m_population = nullptr; // Is there something sensible to do? 
	m_populations.resize(popSizes.size());
	vector<shared_ptr<Population>> newPopulations(popSizes.size());

	m_generation = 0;

	for (size_t p = 0 ; p < popSizes.size() ; p++)
	{
		m_populations[p] = make_shared<Population>();
		newPopulations[p] = make_shared<Population>();
		if (!(r = initializePopulation(popSizes[p], m_populations[p], gfc, refIndividual, refFitness, m_generation)))
			return r;
	}

	auto beforeFitnessCalculatedCallback = [this]() -> bool_t
	{
		bool_t r;
		if (!(r = onBeforeFitnessCalculation(m_generation, m_populations)))
			return "Error inspecting population before fitness calculation in generation " + to_string(m_generation) + ": " + r.getErrorString();
		return true;
	};

	if (!(r = beforeFitnessCalculatedCallback()))
		return r;

	if (!(r = fitnessCalc.calculatePopulationFitness(m_populations)))
		return "Error calculating fitness: " + r.getErrorString();

	auto fitnessCalculatedCallback = [this]() -> bool_t
	{
		bool_t r;
		if (!(r = onFitnessCalculated(m_generation, m_populations)))
			return "Error inspecting population after fitness calculation in generation " + to_string(m_generation) + ": " + r.getErrorString();
		return true;
	};

	if (!(r = fitnessCalculatedCallback()))
		return r;
	
	auto popSizeCheck = [this, &popSizes, &maxPopulationSizes, &minPopulationSizes]() -> bool_t
	{
		if (m_populations.size() != maxPopulationSizes.size())
			return "Number of populations changed! Was " + to_string(maxPopulationSizes.size()) + " and is now " + to_string(m_populations.size());

		for (size_t p = 0 ; p < m_populations.size() ; p++)
		{
			auto &curPop = m_populations[p];
			size_t maxPopulationSize = maxPopulationSizes[p];
			size_t minPopulationSize = minPopulationSizes[p];

			const size_t curPopSize = curPop->size();
			if (curPopSize > maxPopulationSize)
				return "Population size (" + to_string(curPopSize) + ") for population (" + to_string(p) + ") exceeds maximum (" + to_string(maxPopulationSize) + ")";
			if (curPopSize < minPopulationSize)
				return "Population size (" + to_string(curPopSize) + ") for population (" + to_string(p) + ") is less than minimum (" + to_string(minPopulationSize) + ")";
		}
		return true;
	};

	// TODO: check that number of populations hasn't changed?

	while (true)
	{		
		if (m_generation == 0)
		{
			if (!(r = evolver.check(m_populations)))
				return "Error in population evolver check: " + r.getErrorString();
		}

		if (!(r = evolver.createNewPopulations(m_generation, m_populations, popSizes)))
			return "Error creating new population: " + r.getErrorString();

		if (!(r = popSizeCheck()))
			return r;

		m_generation++; // At this point we can call it the new generation

		if (!(r = beforeFitnessCalculatedCallback()))
			return r;

		if (!(r = fitnessCalc.calculatePopulationFitness(m_populations)))
			return "Error calculating fitness: " + r.getErrorString();

		if (!(r = fitnessCalculatedCallback()))
			return r;

		if (m_generation == 1) // Have already increased it, this is the first time
		{
			if (!(r = migrationStrategy.check(m_populations)))
				return "Error in migration strategy check: " + r.getErrorString();
		}
		
		if (!(r = migrationStrategy.migrate(m_generation, m_populations)))
			return "Error in migration step: " + r.getErrorString();

		if (!(r = popSizeCheck())) // Perhaps the amount changed, check it again!
			return r;

		bool shouldStop = false;

		if (!(r = stopCriterion.analyze(evolver, m_generation, shouldStop)))
			return "Error in termination check: " + r.getErrorString();
		if (shouldStop)
			break;
	}

	if (!(r = onAlgorithmDone(m_generation, evolver.getBestIndividuals())))
		return "Error inspecting best individuals upon algorithm end: " + r.getErrorString();

	return true;
}

}
