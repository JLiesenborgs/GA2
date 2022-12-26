#include "evolutionaryalgorithm.h"

using namespace std;
using namespace errut;

namespace eatk
{

EvolutionaryAlgorithm::EvolutionaryAlgorithm()
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
	auto population = make_shared<Population>();
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

	size_t generation = 0;

	if (!(r = initializePopulation(popSize, population, gfc, refIndividual, refFitness, generation)))
		return r;

	auto beforeFitnessCalculatedCallback = [&generation, &population, this]() -> bool_t
	{
		bool_t r;
		if (!(r = onBeforeFitnessCalculation(generation, population)))
			return "Error inspecting population before fitness calculation in generation " + to_string(generation) + ": " + r.getErrorString();
		return true;
	};

	if (!(r = beforeFitnessCalculatedCallback()))
		return r;

	if (!(r = fitnessCalc.calculatePopulationFitness({population})))
		return "Error calculating fitness: " + r.getErrorString();

	auto fitnessCalculatedCallback = [&generation, &population, this]() -> bool_t
	{
		bool_t r;
		if (!(r = onFitnessCalculated(generation, population)))
			return "Error inspecting population after fitness calculation in generation " + to_string(generation) + ": " + r.getErrorString();
		return true;
	};

	if (!(r = fitnessCalculatedCallback()))
		return r;

	while (true)
	{		
		if (generation == 0)
		{
			if (!(r = evolver.check(population)))
				return "Error in population evolver check: " + r.getErrorString();
		}
		if (!(r = evolver.createNewPopulation(generation, population, popSize)))
			return "Error creating new population: " + r.getErrorString();

		const size_t curPopSize = population->size();
		if (curPopSize > maxPopulationSize)
			return "Population size (" + to_string(curPopSize) + ") exceeds maximum (" + to_string(maxPopulationSize) + ")";
		if (curPopSize < minPopulationSize)
			return "Population size (" + to_string(curPopSize) + ") is less than minimum (" + to_string(minPopulationSize) + ")";

		generation++; // At this point we can call it the new generation

		if (!(r = beforeFitnessCalculatedCallback()))
			return r;

		if (!(r = fitnessCalc.calculatePopulationFitness({population})))
			return "Error calculating fitness: " + r.getErrorString();

		if (!(r = fitnessCalculatedCallback()))
			return r;

		bool shouldStop = false;

		if (!(r = stopCriterion.analyze(evolver.getBestIndividuals(), generation, shouldStop)))
			return "Error in termination check: " + r.getErrorString();
		if (shouldStop)
			break;
	}

	if (!(r = onAlgorithmDone(generation, evolver.getBestIndividuals())))
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

	vector<shared_ptr<Population>> populations(popSizes.size());
	vector<shared_ptr<Population>> newPopulations(popSizes.size());

	size_t generation = 0;

	for (size_t p = 0 ; p < popSizes.size() ; p++)
	{
		populations[p] = make_shared<Population>();
		newPopulations[p] = make_shared<Population>();
		if (!(r = initializePopulation(popSizes[p], populations[p], gfc, refIndividual, refFitness, generation)))
			return r;
	}

	auto beforeFitnessCalculatedCallback = [&generation, &populations, this]() -> bool_t
	{
		bool_t r;
		if (!(r = onBeforeFitnessCalculation(generation, populations)))
			return "Error inspecting population before fitness calculation in generation " + to_string(generation) + ": " + r.getErrorString();
		return true;
	};

	if (!(r = beforeFitnessCalculatedCallback()))
		return r;

	if (!(r = fitnessCalc.calculatePopulationFitness(populations)))
		return "Error calculating fitness: " + r.getErrorString();

	auto fitnessCalculatedCallback = [&generation, &populations, this]() -> bool_t
	{
		bool_t r;
		if (!(r = onFitnessCalculated(generation, populations)))
			return "Error inspecting population after fitness calculation in generation " + to_string(generation) + ": " + r.getErrorString();
		return true;
	};

	if (!(r = fitnessCalculatedCallback()))
		return r;
	
	auto popSizeCheck = [&populations, &popSizes, &maxPopulationSizes, &minPopulationSizes]() -> bool_t
	{
		for (size_t p = 0 ; p < populations.size() ; p++)
		{
			auto &population = populations[p];
			size_t maxPopulationSize = maxPopulationSizes[p];
			size_t minPopulationSize = minPopulationSizes[p];

			const size_t curPopSize = population->size();
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
		if (generation == 0)
		{
			if (!(r = evolver.check(populations)))
				return "Error in population evolver check: " + r.getErrorString();
		}

		if (!(r = evolver.createNewPopulations(generation, populations, popSizes)))
			return "Error creating new population: " + r.getErrorString();

		if (!(r = popSizeCheck()))
			return r;

		generation++; // At this point we can call it the new generation

		if (!(r = beforeFitnessCalculatedCallback()))
			return r;

		if (!(r = fitnessCalc.calculatePopulationFitness(populations)))
			return "Error calculating fitness: " + r.getErrorString();

		if (!(r = fitnessCalculatedCallback()))
			return r;

		if (generation == 1) // Have already increased it, this is the first time
		{
			if (!(r = migrationStrategy.check(populations)))
				return "Error in migration strategy check: " + r.getErrorString();
		}
		
		if (!(r = migrationStrategy.migrate(generation, populations)))
			return "Error in migration step: " + r.getErrorString();

		if (!(r = popSizeCheck())) // Perhaps the amount changed, check it again!
			return r;

		bool shouldStop = false;

		if (!(r = stopCriterion.analyze(evolver.getBestIndividuals(), generation, shouldStop)))
			return "Error in termination check: " + r.getErrorString();
		if (shouldStop)
			break;
	}

	if (!(r = onAlgorithmDone(generation, evolver.getBestIndividuals())))
		return "Error inspecting best individuals upon algorithm end: " + r.getErrorString();

	return true;
}

}
