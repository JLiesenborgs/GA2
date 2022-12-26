#pragma once

#include "eatkconfig.h"
#include "crossovermutation.h"
#include "stopcriterion.h"
#include "populationevolver.h"
#include "migrationstrategy.h"
#include "calculation.h"

namespace eatk
{

class EvolutionaryAlgorithm
{
public:
	EvolutionaryAlgorithm();
	virtual ~EvolutionaryAlgorithm();

	errut::bool_t run(IndividualCreation &gfc,
			   PopulationEvolver &evolver,
			   PopulationFitnessCalculation &fitnessCalc,
			   StopCriterion &stopCriterion,
			   size_t popSize,
			   size_t minPopulationSize = 0,
			   size_t maxPopulationSize = 0);

	errut::bool_t run(IndividualCreation &gfc,
			   PopulationEvolver &evolver,
			   PopulationFitnessCalculation &fitnessCalc,
			   StopCriterion &stopCriterion,
			   MigrationStrategy &migrationStrategy,
			   const std::vector<size_t> &popSizes,
			   const std::vector<size_t> &minPopulationSizes = std::vector<size_t>(),
			   const std::vector<size_t> &maxPopulationSizes = std::vector<size_t>());
protected:
	virtual errut::bool_t onBeforeFitnessCalculation(size_t generation, const std::shared_ptr<Population> &population) { return true; }
	virtual errut::bool_t onFitnessCalculated(size_t generation, const std::shared_ptr<Population> &population) { return true; }
	virtual errut::bool_t onBeforeFitnessCalculation(size_t generation, const std::vector<std::shared_ptr<Population>> &populations) { return true; }
	virtual errut::bool_t onFitnessCalculated(size_t generation, const std::vector<std::shared_ptr<Population>> &populations) { return true; }
	virtual errut::bool_t onAlgorithmDone(size_t generation, const std::vector<std::shared_ptr<Individual>> &bestIndividuals) { return true; }
};

}
