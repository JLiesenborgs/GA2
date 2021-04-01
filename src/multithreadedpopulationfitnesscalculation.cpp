#include "multithreadedpopulationfitnesscalculation.h"
#include <thread>
#include <mutex>
#include <condition_variable>

using namespace std;
using namespace errut;

namespace eatk
{

MultiThreadedPopulationFitnessCalculation::MultiThreadedPopulationFitnessCalculation()
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

	m_threadGenomeCalculations = threadGenomeCalculations;
	m_workers.resize(m_threadGenomeCalculations.size());

	m_allThreadsWaiter = make_unique<ThreadsReadyWaiter>(m_workers.size());
	m_individualWaiters = make_unique<IndividualWaiters>(m_workers.size());

	for (size_t i = 0 ; i < m_workers.size() ; i++)
		m_workers[i] = thread(staticWorkerThread, this, i);

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
	auto performActions = [this, idx, &done]()
	{
		// cerr << "In thread " + to_string(idx) + ", m_pPopulations = " + to_string((long int)((void*)m_pPopulations)) + "\n";
		if (!m_pPopulations)
			done = true;
		else
		{
			m_errors[idx] = false;
			m_errorStrings[idx].clear();
			bool_t r = workerCalculatePopulationFitness(*m_pPopulations, idx);
			if (!r)
			{
				m_errors[idx] = true;
				m_errorStrings[idx] = r.getErrorString();
			}
		}
	};

	while (!done)
	{
		// cerr << "Thread " + to_string(idx) + " waiting to start\n";
		m_allThreadsWaiter->signal(idx);
		
		// Wait till we get some work
		m_individualWaiters->wait(idx);

		performActions();
	}	

	// cerr << "Thread " + to_string(idx) +  " exiting\n";
}

void MultiThreadedPopulationFitnessCalculation::destroyThreadPool()
{
	if (m_workers.size() == 0)
		return;

	m_pPopulations = nullptr; // Signals that we're done
	
	m_individualWaiters->signalAll();

	for (auto &t : m_workers)
		t.join();

	m_workers.resize(0);
}

bool_t MultiThreadedPopulationFitnessCalculation::workerCalculatePopulationFitness(const vector<shared_ptr<Population>> &populations, size_t workerIdx)
{
	auto performIndividualActionForWorker = [this, workerIdx, &populations](auto action) -> bool_t
	{
		size_t count = 0;
		bool_t r;
		for (auto &pop : populations)
		{
			for (auto &ind : pop->individuals())
			{
				if (count % m_workers.size() == workerIdx)
				{
					if (!(r = action(ind)))
						return r;
				}
				count++;
			}
		}
		return true;
	};			

	auto pGenomeCalc = m_threadGenomeCalculations[workerIdx].get();

	// Start new calculation
	bool_t r = performIndividualActionForWorker([pGenomeCalc](auto ind) -> bool_t {
		Genome *pGenome = ind->genomePtr();
		Fitness *pFitness = ind->fitnessPtr();

		if (!pFitness->isCalculated())
			return pGenomeCalc->startNewCalculation(*pGenome);
		return true;
	});

	if (!r)
		return "Error starting new calculation for a genome: " + r.getErrorString();

	bool allCalculated = false;
	while (!allCalculated)
	{
		allCalculated = true;
		r = performIndividualActionForWorker([&allCalculated, pGenomeCalc](auto ind) -> bool_t {
			Genome *pGenome = ind->genomePtr();
			Fitness *pFitness = ind->fitnessPtr();

			if (!pFitness->isCalculated())
			{
				bool_t r = pGenomeCalc->pollCalculate(*pGenome, *pFitness);
				if (!r)
					return r;
				if (!pFitness->isCalculated()) // Still not done
					allCalculated = false;
			}

			return true;
		});

		if (!r)
			return "Error in calculation for a genome: " + r.getErrorString();
	}

	return true;
}

errut::bool_t MultiThreadedPopulationFitnessCalculation::check(const std::vector<std::shared_ptr<Population>> &populations)
{
	if (m_threadGenomeCalculations.size() < 1)
		return "Need at least one thread";
	return true;
}

bool_t MultiThreadedPopulationFitnessCalculation::calculatePopulationFitness(const vector<shared_ptr<Population>> &populations)
{
	if (m_threadGenomeCalculations.size() < 1)
		return "Need at least one thread";

	m_pPopulations = &populations;

	m_individualWaiters->signalAll();

	m_allThreadsWaiter->wait();

	for (size_t i = 0 ; i < m_errors.size() ; i++)
	{
		if (m_errors[i])
			return m_errorStrings[i];
	}

	return true;
}

}
