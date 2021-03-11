#pragma once

#include "mogal2config.h"
#include "population.h"

namespace mogal2
{

class MultiThreadedPopulationFitnessCalculation : public PopulationFitnessCalculation
{
public:
	MultiThreadedPopulationFitnessCalculation();
	~MultiThreadedPopulationFitnessCalculation();

	errut::bool_t initThreadPool(const std::vector<std::shared_ptr<GenomeFitnessCalculation>> &threadGenomeCalculations);

	errut::bool_t check(const std::vector<std::shared_ptr<Population>> &populations) override;
	// TODO: all populations should have exactly the same genomes! (ie same number of floats)
	errut::bool_t calculatePopulationFitness(const std::vector<std::shared_ptr<Population>> &populations);
private:
	std::vector<std::shared_ptr<GenomeFitnessCalculation>> m_threadGenomeCalculations;
};

}
