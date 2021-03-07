#pragma once

#include "genomefitness.h"
#include "crossovermutation.h"
#include <memory>

class GAFactory
{
public:
    GAFactory() { }
    virtual ~GAFactory() { }

    // TODO: how to signal error?
    virtual std::shared_ptr<Genome> createInitializedGenome() = 0;
    virtual std::shared_ptr<Fitness> createEmptyFitness() = 0;
    virtual std::shared_ptr<PopulationMutation> getPopulationMutation() = 0;
    virtual std::shared_ptr<PopulationCrossover> getPopulationCrossover() = 0;
};
