#pragma once

#include "eatkconfig.h"
#include "genomefitness.h"
#include "population.h"

namespace eatk
{

class GenomeFitnessCalculation
{
public:
	GenomeFitnessCalculation() { }
	virtual ~GenomeFitnessCalculation() { }

	// TODO: do we need a check function here?

	// genomesForPopulationCalculator is e.g. the number that needs to be calculated across threads
	// The MPI implementation itself uses a different local calculator, so doesn't call this itself
	// This means that it's mainly a multi-thread thing
	//
	// 'iteration' is meant to be a count for each time calculatePopulationFitness is called. This will
	// mostly match the generation number, but if the calculator is used for several different EAs,
	// then this may deviate. The main idea is to have some kind of identifier to be able to check if
	// some initialization needs to be done when a new iteration starts
	virtual errut::bool_t onNewCalculationStart(size_t iteration, size_t genomesForThisCalculator, size_t genomesForPopulationCalculator)  { return true; }
	virtual errut::bool_t onCalculationStarted() { return true; }
	virtual errut::bool_t onCalculationEnded() { return true; }

	// These are to allow a more async version, but by default the sync version is called
	virtual errut::bool_t startNewCalculation(const Genome &genome) { return true; }
	
	// Return type says if something went wrong; fitness.isCalculated marks when done
	// May need multiple calls, to make async behaviour possible
	virtual errut::bool_t pollCalculate(const Genome &genome, Fitness &fitness)
	{
		auto r = calculate(genome, fitness);
		fitness.setCalculated();
		return r;
	}

	// Convenience function for sync operation
	virtual errut::bool_t calculate(const Genome &genome, Fitness &fitness) { return "Not implemented"; }
};

class PopulationFitnessCalculation
{
public:
	PopulationFitnessCalculation() { }
	virtual ~PopulationFitnessCalculation() { }

	// This function can be called eg on first generation, to check that the types etc
	// are in order. This can then skip some checks in calculatePopulationFitness
	virtual errut::bool_t check(const std::vector<std::shared_ptr<Population>> &populations) { return "Not implemented in base class"; }

	// TODO: all populations should have exactly the same genomes! (ie same number of floats)
	virtual errut::bool_t calculatePopulationFitness(const std::vector<std::shared_ptr<Population>> &populations) { return "Not implemented in base class"; }
};

}
