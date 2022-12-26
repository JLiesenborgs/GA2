#pragma once

#include "eatkconfig.h"
#include "population.h"
#include "calculation.h"
#include <thread>
#include <mutex>
#include <condition_variable>

namespace eatk
{

class MultiThreadedPopulationFitnessCalculation : public PopulationFitnessCalculation
{
public:
	MultiThreadedPopulationFitnessCalculation();
	~MultiThreadedPopulationFitnessCalculation();

	errut::bool_t initThreadPool(const std::vector<std::shared_ptr<GenomeFitnessCalculation>> &threadGenomeCalculations);
	void destroyThreadPool();

	errut::bool_t check(const std::vector<std::shared_ptr<Population>> &populations) override;
	// TODO: all populations should have exactly the same genomes! (ie same number of floats)
	errut::bool_t calculatePopulationFitness(const std::vector<std::shared_ptr<Population>> &populations) override;
private:
	static void staticWorkerThread(MultiThreadedPopulationFitnessCalculation *pInstance, size_t idx);
	void workerThread(size_t idx);
	errut::bool_t workerCalculatePopulationFitness(size_t workerIdx);

	size_t m_totalThreads;
	std::vector<std::shared_ptr<GenomeFitnessCalculation>> m_threadGenomeCalculations;
	std::vector<std::thread> m_workers;
	std::vector<bool> m_errors;
	std::vector<std::string> m_errorStrings;

	std::vector<std::vector<std::pair<Genome *, Fitness *>>> m_helperGenomes;
	size_t m_genomesToCalculateInThisIteration;


	class ThreadsReadyWaiter
	{
	public:
		ThreadsReadyWaiter(size_t N) : m_N(N), m_count(0) { }
		void signal(size_t idx)
		{
			std::unique_lock<std::mutex> lk(m_mut);
			m_count++;
			lk.unlock();

			m_cv.notify_one();
		}

		void wait()
		{
			std::unique_lock<std::mutex> lk(m_mut);
			m_cv.wait(lk, [this]() { 
				return m_count == m_N;
			});

			m_count = 0; // Reset count again
		}
	private:
		const size_t m_N;
		size_t m_count;
		std::mutex m_mut;
		std::condition_variable m_cv;
	};

	class IndividualWaiters
	{
	public:
		IndividualWaiters(size_t N) : m_mut(N), m_cv(N), m_flags(N, 0) { }
		void wait(size_t idx)
		{
			std::unique_lock<std::mutex> lk(m_mut[idx]);
			m_cv[idx].wait(lk, [this, idx]() { return m_flags[idx]; });
			m_flags[idx] = 0;
		}

		void signal(size_t idx)
		{
			std::unique_lock<std::mutex> lk(m_mut[idx]);
			m_flags[idx] = 1;
			lk.unlock();

			m_cv[idx].notify_one();
		}

		void signalAll()
		{
			for (size_t idx = 0 ; idx < m_cv.size() ; idx++)
				signal(idx);
		}
	private:
		std::vector<std::mutex> m_mut;
		std::vector<std::condition_variable> m_cv;
		std::vector<int> m_flags;
	};

	std::unique_ptr<ThreadsReadyWaiter> m_allThreadsWaiter;
	std::unique_ptr<IndividualWaiters> m_individualWaiters;
};

}
