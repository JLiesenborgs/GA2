#include "multipopulationevolver.h"
#include <algorithm>

using namespace std;
using namespace errut;

namespace eatk
{

MultiPopulationEvolver::MultiPopulationEvolver(const vector<shared_ptr<PopulationEvolver>> &singlePopulationEvolvers,
                                               const std::shared_ptr<FitnessComparison> &fitCmp, size_t numObjectives)
	: m_singlePopEvolvers(singlePopulationEvolvers),
	  m_fitCmp(fitCmp),
	  m_numObjectives(numObjectives)
{
	if (numObjectives > 1)
		m_ndSetCreator = make_unique<FasterNonDominatedSetCreator>(fitCmp, numObjectives);
}

MultiPopulationEvolver::~MultiPopulationEvolver()
{
}

bool_t MultiPopulationEvolver::check(const vector<shared_ptr<Population>> &populations)
{
	if (m_numObjectives < 1)
		return "No objectives were specified";

	if (populations.size() < 1)
		return "No populations are present";

	if (populations.size() == 1)
		return "Only one population is present, no need to use this class";

	if (populations[0]->size() < 1)
		return "No individuals in first population";

	if (m_singlePopEvolvers.size() != populations.size())
		return "A different number of populations and population evolvers were specified";

	bool_t r;

	for (size_t i = 0 ; i < m_singlePopEvolvers.size() ; i++)
	{
		if (!(r = m_singlePopEvolvers[i]->check(populations[i])))
			return "Error checking evolver for population " + to_string(i) + ": " + r.getErrorString();
	}

	return true;
}

// TODO: multi-threaded version?

bool_t MultiPopulationEvolver::createNewPopulations(size_t generation, vector<shared_ptr<Population>> &populations, const vector<size_t> &targetPopulationSizes)
{
	if (m_singlePopEvolvers.size() != populations.size())
		return "A different number of populations and population evolvers were specified";

	if (targetPopulationSizes.size() != populations.size())
		return "A different number of target population sizes and populations were specified";

	bool_t r;

	for (size_t i = 0 ; i < m_singlePopEvolvers.size() ; i++)
	{
		if (!(r = m_singlePopEvolvers[i]->createNewPopulation(generation, populations[i], targetPopulationSizes[i])))
			return "Error creating new population " + to_string(i) + ": " + r.getErrorString();
	}
	return true;
}

const vector<shared_ptr<Individual>> &MultiPopulationEvolver::getBestIndividuals() const
{
	vector<shared_ptr<Individual>> best;

	if (m_numObjectives == 1)
	{
		size_t count = 0;
		size_t bestPop = 0;
		for (auto &ev : m_singlePopEvolvers)
		{
			auto &curBest = ev->getBestIndividuals();
			if (curBest.size() > 2)
			{
				cerr << "ERROR: got more than one best individual with only one objective!" << endl;
				m_best.clear();
				return m_best;
			}

			if (curBest.size() == 0)
				continue;

			if (best.size() == 0)
			{
				bestPop = count;
				best = curBest;
			}
			else
			{
				auto &bestInd = best[0];
				auto &curBestInd = curBest[0];

				if (m_fitCmp->isFitterThan(curBestInd->fitnessRef(), bestInd->fitnessRef(), 0))
				{
					best[0] = curBestInd;
					bestPop = count;
				}
			}

			count++;
		}

		cerr << "Using best from population " << bestPop << endl;
	}
	else // use ndset
	{
		// TODO
	}

	std::swap(best, m_best);
	return m_best;
}

}
