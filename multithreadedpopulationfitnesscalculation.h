#pragma once

#include "population.h"

class MultiThreadedPopulationFitnessCalculation : public PopulationFitnessCalculation
{
public:
	MultiThreadedPopulationFitnessCalculation();
	~MultiThreadedPopulationFitnessCalculation();

	errut::bool_t initThreadPool(const std::vector<std::shared_ptr<GenomeFitnessCalculation>> &threadGenomeCalculations);

	// TODO: all populations should have exactly the same genomes! (ie same number of floats)
	errut::bool_t calculatePopulationFitness(const std::vector<std::shared_ptr<Population>> &populations);
private:
	std::vector<std::shared_ptr<GenomeFitnessCalculation>> m_threadGenomeCalculations;
};

