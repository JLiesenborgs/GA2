#include "singlethreadedpopulationfitnesscalculation.h"

using namespace std;
using namespace errut;

SingleThreadedPopulationFitnessCalculation::SingleThreadedPopulationFitnessCalculation(shared_ptr<GenomeFitnessCalculation> genomeFitCalc)
    : m_genomeFitnessCalculation(genomeFitCalc)
{ 
}

SingleThreadedPopulationFitnessCalculation::~SingleThreadedPopulationFitnessCalculation()
{
}

bool_t SingleThreadedPopulationFitnessCalculation::calculatePopulationFitness(const vector<shared_ptr<Population>> &populations)
{
	if (!m_genomeFitnessCalculation.get())
		return "No genome fitness calculation has been set";

	// First, initialize the calculations
	for (auto &pop : populations)
		for (auto &i : pop->m_individuals)
			if (!i->m_fitness->isCalculated())
				m_genomeFitnessCalculation->startNewCalculation(*i->m_genome);
	
	// Calculate until all is done
	bool allCalculated;
	do
	{
		allCalculated = true;
		for (auto &pop : populations)
		{
			for (auto &i : pop->m_individuals)
			{
				Fitness &f = *i->m_fitness;
				if (!f.isCalculated())
				{
					m_genomeFitnessCalculation->pollCalculate(*i->m_genome, f);
					if (!f.isCalculated())
						allCalculated = false;
				}
			}
		}
	} while (!allCalculated);

	return true;
}
