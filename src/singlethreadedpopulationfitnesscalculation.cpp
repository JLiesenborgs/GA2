#include "singlethreadedpopulationfitnesscalculation.h"

using namespace std;
using namespace errut;

namespace eatk
{

SingleThreadedPopulationFitnessCalculation::SingleThreadedPopulationFitnessCalculation(const shared_ptr<GenomeFitnessCalculation> &genomeFitCalc)
	: m_genomeFitnessCalculation(genomeFitCalc),
	  m_iteration(0)
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
	bool_t r;

	if (!m_genomeFitnessCalculation.get())
		return "No genome fitness calculation has been set";

	m_tmpIndividuals.clear();

	// First, initialize the calculations
	for (auto &pop : populations)
		for (auto &i : pop->individuals())
			if (!i->fitnessRef().isCalculated())
				m_tmpIndividuals.push_back(i.get());

	if (!(r = m_genomeFitnessCalculation->onNewCalculationStart(m_iteration, m_tmpIndividuals.size(), m_tmpIndividuals.size())))
		return "Can't signal new calculation start: " + r.getErrorString();

	m_iteration++;

	for (auto i : m_tmpIndividuals)
	{
		auto r = m_genomeFitnessCalculation->startNewCalculation(i->genomeRef());
		if (!r)
			return "Error starting genome calculation: " + r.getErrorString();
	}

	if (!(r = m_genomeFitnessCalculation->onCalculationStarted()))
		return "Error signalling calculation started: " + r.getErrorString();
	
	// Calculate until all is done
	bool allCalculated;
	do
	{
		allCalculated = true;
		for (auto i : m_tmpIndividuals)
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
	} while (!allCalculated);

	if (!(r = m_genomeFitnessCalculation->onCalculationEnded()))
		return "Error signalling calculation ended: " + r.getErrorString();

	return true;
}

}
