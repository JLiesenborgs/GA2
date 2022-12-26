#include "multipopulationevolver.h"
#include <algorithm>

using namespace std;
using namespace errut;

namespace eatk
{

MultiPopulationEvolver::MultiPopulationEvolver(const vector<shared_ptr<PopulationEvolver>> &singlePopulationEvolvers,
                                               const shared_ptr<BestIndividualMerger> &merger)
	: m_singlePopEvolvers(singlePopulationEvolvers),
	  m_merger(merger),
	  m_lastCreatedGeneration(numeric_limits<size_t>::max()),
	  m_lastMergeGeneration(numeric_limits<size_t>::max())
{
}

MultiPopulationEvolver::~MultiPopulationEvolver()
{
}

bool_t MultiPopulationEvolver::check(const vector<shared_ptr<Population>> &populations)
{
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

	if (!(r = m_merger->check(populations, m_singlePopEvolvers)))
		return "Error checking best individual merger: " + r.getErrorString();

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

	m_lastCreatedGeneration = generation;
	return true;
}

const vector<shared_ptr<Individual>> &MultiPopulationEvolver::getBestIndividuals() const
{
	if (m_lastMergeGeneration == m_lastCreatedGeneration) // already calculated it, just return it again
		return m_best;

	assert(m_merger.get());

	m_lastMergeGeneration = m_lastCreatedGeneration;
	m_best = m_merger->mergeBestIndividuals(m_singlePopEvolvers);
	return m_best;
}

SingleObjectiveBestIndividualMerger::SingleObjectiveBestIndividualMerger(const shared_ptr<FitnessComparison> &fitCmp)
	: m_fitCmp(fitCmp)
{
}

SingleObjectiveBestIndividualMerger::~SingleObjectiveBestIndividualMerger()
{
}

bool_t SingleObjectiveBestIndividualMerger::check(const vector<shared_ptr<Population>> &populations,
                                                  const vector<shared_ptr<PopulationEvolver>> &evolvers)
{
	if (populations.size() != evolvers.size())
		return "Number of populations doesn't match the number of population evolvers";
	return true;
}

const vector<shared_ptr<Individual>> &SingleObjectiveBestIndividualMerger::mergeBestIndividuals(const vector<shared_ptr<PopulationEvolver>> &evolvers)
{
	vector<shared_ptr<Individual>> best;

	size_t count = 0;
	size_t bestPop = 0;
	for (auto &ev : evolvers)
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

	//cerr << "Using best from population " << bestPop << endl;

	swap(best, m_best);
	return m_best;
}

MultiObjectiveBestIndividualMerger::MultiObjectiveBestIndividualMerger(const shared_ptr<NonDominatedSetCreator> &ndCreator,
								   const shared_ptr<DuplicateIndividualRemoval> &dupRemoval)
	: m_ndSetCreator(ndCreator),
	  m_dupRemoval(dupRemoval)
{
}

MultiObjectiveBestIndividualMerger::~MultiObjectiveBestIndividualMerger()
{
}

bool_t MultiObjectiveBestIndividualMerger::check(const vector<shared_ptr<Population>> &populations,
					const vector<shared_ptr<PopulationEvolver>> &evolvers)
{
	if (populations.size() != evolvers.size())
		return "Number of populations doesn't match the number of population evolvers";

	bool_t r;

	if (m_dupRemoval.get())
	{
		for (auto &pop : populations)
		{
			if (!(r = m_dupRemoval->check(pop->individuals())))
				return "Error checking duplication removal: " + r.getErrorString();
		}
	}

	return true;
}

const vector<shared_ptr<Individual>> &MultiObjectiveBestIndividualMerger::mergeBestIndividuals(const vector<shared_ptr<PopulationEvolver>> &evolvers)
{
	bool_t r;

	m_best.clear();
	for (auto &ev : evolvers)
	{
		if (!(r = m_ndSetCreator->mergeNDSets(m_best, ev->getBestIndividuals())))
		{
			cerr << "ERROR: can't merge ND sets: " << r.getErrorString() << endl;
			m_best.clear();
			return m_best;
		}
	}

	if (m_dupRemoval.get())
	{
		if (!(r = m_dupRemoval->removeDuplicates(m_best)))
		{
			cerr << "ERROR: can't remove duplicates: " + r.getErrorString();
			m_best.clear();
			return m_best;
		}
	}

	return m_best;
}

}
