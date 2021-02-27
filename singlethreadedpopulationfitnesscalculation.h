#pragma once

#include "population.h"

class SingleThreadedPopulationFitnessCalculation : public PopulationFitnessCalculation
{
public:
	SingleThreadedPopulationFitnessCalculation(std::shared_ptr<GenomeFitnessCalculation> genomeFitCalc);
	~SingleThreadedPopulationFitnessCalculation();

	errut::bool_t check(const std::vector<std::shared_ptr<Population>> &populations) override;
	errut::bool_t calculatePopulationFitness(const std::vector<std::shared_ptr<Population>> &populations) override;
private:
	std::shared_ptr<GenomeFitnessCalculation> m_genomeFitnessCalculation;
};
