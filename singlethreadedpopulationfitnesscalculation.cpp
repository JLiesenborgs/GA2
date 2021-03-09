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

bool_t SingleThreadedPopulationFitnessCalculation::check(const vector<shared_ptr<Population>> &populations)
{
	if (!m_genomeFitnessCalculation.get())
		return "No genome fitness calculation has been set";
	return true;
}

bool_t SingleThreadedPopulationFitnessCalculation::calculatePopulationFitness(const vector<shared_ptr<Population>> &populations)
{
	if (!m_genomeFitnessCalculation.get())
		return "No genome fitness calculation has been set";

	// First, initialize the calculations
	for (auto &pop : populations)
	{
		for (auto &i : pop->m_individuals)
		{
			if (!i->fitnessRef().isCalculated())
			{		
				auto r = m_genomeFitnessCalculation->startNewCalculation(i->genomeRef());
				if (!r)
					return "Error starting genome calculation: " + r.getErrorString();
			}
		}
	}
	
	// Calculate until all is done
	bool allCalculated;
	do
	{
		allCalculated = true;
		for (auto &pop : populations)
		{
			for (auto &i : pop->m_individuals)
			{
				Fitness &f = i->fitnessRef();
				if (!f.isCalculated())
				{
					auto r = m_genomeFitnessCalculation->pollCalculate(i->genomeRef(), f);
					if (!r)
						return "Error calculating fitness for genome: " + r.getErrorString();

					if (!f.isCalculated())
						allCalculated = false;
				}
			}
		}
	} while (!allCalculated);

	return true;
}
