#pragma once

#include "eatkconfig.h"
#include "population.h"
#include "calculation.h"

namespace eatk
{

class SingleThreadedPopulationFitnessCalculation : public PopulationFitnessCalculation
{
public:
	SingleThreadedPopulationFitnessCalculation(const std::shared_ptr<GenomeFitnessCalculation> &genomeFitCalc);
	~SingleThreadedPopulationFitnessCalculation();

	errut::bool_t check(const std::vector<std::shared_ptr<Population>> &populations) override;
	errut::bool_t calculatePopulationFitness(const std::vector<std::shared_ptr<Population>> &populations) override;
private:
	std::shared_ptr<GenomeFitnessCalculation> m_genomeFitnessCalculation;
	std::vector<Individual*> m_tmpIndividuals;
	size_t m_iteration;
};

}
