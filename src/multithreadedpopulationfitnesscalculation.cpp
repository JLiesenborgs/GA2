#include "multithreadedpopulationfitnesscalculation.h"
#include <thread>
#include <mutex>
#include <condition_variable>
#include <sstream>
#include <set>

using namespace std;
using namespace errut;

namespace eatk
{

MultiThreadedPopulationFitnessCalculation::MultiThreadedPopulationFitnessCalculation()
 : m_totalThreads(0)
{
}

MultiThreadedPopulationFitnessCalculation::~MultiThreadedPopulationFitnessCalculation()
{
	destroyThreadPool();
}

bool_t MultiThreadedPopulationFitnessCalculation::initThreadPool(const std::vector<std::shared_ptr<GenomeFitnessCalculation>> &threadGenomeCalculations)
{
	if (m_workers.size() > 0)
		return "Thread pool already exists";

	if (threadGenomeCalculations.size() < 1)
		return "Need at least one thread";

	m_totalThreads = threadGenomeCalculations.size();
	m_threadGenomeCalculations = threadGenomeCalculations;
	m_workers.resize(m_totalThreads-1);

	m_allThreadsWaiter = make_unique<ThreadsReadyWaiter>(m_workers.size());
	m_individualWaiters = make_unique<IndividualWaiters>(m_workers.size());

	m_helperGenomes.resize(m_totalThreads);

	for (size_t i = 0 ; i < m_workers.size() ; i++)
		m_workers[i] = thread(staticWorkerThread, this, i);

	if (m_workers.size() > 0)
		m_allThreadsWaiter->wait();

	m_errors.resize(m_workers.size());
	m_errorStrings.resize(m_workers.size());

	return true;
}

void MultiThreadedPopulationFitnessCalculation::staticWorkerThread(MultiThreadedPopulationFitnessCalculation *pInstance, size_t idx)
{
	pInstance->workerThread(idx);
}

void MultiThreadedPopulationFitnessCalculation::workerThread(size_t idx)
{
	bool done = false;

	while (!done)
	{
		// cerr << "Thread " + to_string(idx) + " waiting to start\n";
		m_allThreadsWaiter->signal(idx);
		
		// Wait till we get some work
		m_individualWaiters->wait(idx);

		// cerr << "In thread " + to_string(idx) + ", m_pPopulations = " + to_string((long int)((void*)m_pPopulations)) + "\n";
		if (m_helperGenomes.empty())
			done = true;
		else
		{
			m_errors[idx] = false;
			m_errorStrings[idx].clear();
			bool_t r = workerCalculatePopulationFitness(idx);
			if (!r)
			{
				m_errors[idx] = true;
				m_errorStrings[idx] = r.getErrorString();
			}
		}
	}	

	// cerr << "Thread " + to_string(idx) +  " exiting\n";
}

void MultiThreadedPopulationFitnessCalculation::destroyThreadPool()
{
	if (m_helperGenomes.size() == 0)
		return;

	m_helperGenomes.clear();
	
	if (m_workers.size() > 0)
		m_individualWaiters->signalAll();

	for (auto &t : m_workers)
		t.join();

	m_workers.resize(0);
}

bool_t MultiThreadedPopulationFitnessCalculation::workerCalculatePopulationFitness(size_t workerIdx)
{
	auto pGenomeCalc = m_threadGenomeCalculations[workerIdx].get();
	bool_t r;

	if (!(r = pGenomeCalc->onNewCalculationStart(m_helperGenomes[workerIdx].size(), m_genomesToCalculateInThisIteration)))
		return "Error signalling new calculation iteration start: " + r.getErrorString();

	for (auto &gf : m_helperGenomes[workerIdx])
	{
		if (!(r = pGenomeCalc->startNewCalculation(*gf.first)))
			return "Error starting new calculation for a genome: " + r.getErrorString();
	}

	if (!(r = pGenomeCalc->onCalculationStarted()))
		return "Error signalling calculation started: " + r.getErrorString();

	bool allCalculated = false;
	while (!allCalculated)
	{
		allCalculated = true;

		for (auto &gf : m_helperGenomes[workerIdx])
		{
			Genome *pGenome = gf.first;
			Fitness *pFitness = gf.second;
			if (!pFitness->isCalculated())
			{
				if (!(r = pGenomeCalc->pollCalculate(*pGenome, *pFitness)))
					return "Error calculating genome fitness: " + r.getErrorString();
				if (!pFitness->isCalculated()) // Still not done
					allCalculated = false;
			}
		}
	}

	if (!(r = pGenomeCalc->onCalculationEnded()))
		return "Error signalling calculation end: " + r.getErrorString();

	return true;
}

errut::bool_t MultiThreadedPopulationFitnessCalculation::check(const std::vector<std::shared_ptr<Population>> &populations)
{
	if (m_totalThreads < 1)
		return "Need at least one thread";
	return true;
}

bool_t MultiThreadedPopulationFitnessCalculation::calculatePopulationFitness(const vector<shared_ptr<Population>> &populations)
{
	if (m_totalThreads < 1)
		return "Need at least one thread";

	for (auto &h : m_helperGenomes)
		h.clear();

	size_t nextHelper = 0;
	size_t genomesToCalculate = 0;

	// Split the work over the helpers
	for (auto pPop : populations)
	{		
		for (auto &i : pPop->individuals())
		{
			if (!i->fitnessRef().isCalculated())
			{
				m_helperGenomes[nextHelper].push_back({i->genomePtr(), i->fitnessPtr()});
				nextHelper = (nextHelper + 1)%m_totalThreads;
				genomesToCalculate++;
			}
		}
	}

	m_genomesToCalculateInThisIteration = genomesToCalculate;

	if (m_totalThreads > 1)
		m_individualWaiters->signalAll();

	// Calculate part ourselves
	bool_t r = workerCalculatePopulationFitness(m_totalThreads-1);

	if (m_totalThreads > 1)
		m_allThreadsWaiter->wait();

	// Check our own error
	set<string> errors;
	if (!r)
		errors.insert(r.getErrorString());

	// Check other threads' errors
	for (size_t i = 0 ; i < m_errors.size() ; i++)
	{
		if (m_errors[i])
			errors.insert(m_errorStrings[i]);
	}

	if (errors.size() > 0)
	{
		stringstream ss;
		if (errors.size() > 1)
			ss << "(";

		auto it = errors.begin();
		ss << *it;
		it++;

		while (it != errors.end())
		{
			ss << "," << *it;
			it++;
		}
		if (errors.size() > 1)
			ss << ")";
		return ss.str();
	}

	return true;
}

}
