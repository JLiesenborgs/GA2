#pragma once

#include "eatkconfig.h"
#include "crossovermutation.h"
#include "stopcriterion.h"

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
protected:
    virtual errut::bool_t onBeforeFitnessCalculation(size_t generation, std::shared_ptr<Population> &population) { return true; }
    virtual errut::bool_t onFitnessCalculated(size_t generation, std::shared_ptr<Population> &population) { return true; }
    virtual errut::bool_t onAlgorithmDone(size_t generation, const std::vector<std::shared_ptr<Individual>> &bestIndividuals) { return true; }
};

}