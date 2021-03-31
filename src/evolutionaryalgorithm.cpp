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

	for (size_t i = 0 ; i < popSize ; i++)
	{
		auto g = gfc.createInitializedGenome();
		if (!g.get())
			return "Unable to create an inialized genome";

		auto f = refFitness->createCopy(false);
		population->append(refIndividual->createNew(g, f, generation));
	}

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

}
