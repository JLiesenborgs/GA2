#include "multithreadedpopulationfitnesscalculation.h"
#include <thread>
#include <mutex>
#include <condition_variable>

using namespace std;
using namespace errut;

MultiThreadedPopulationFitnessCalculation::MultiThreadedPopulationFitnessCalculation()
{
}

MultiThreadedPopulationFitnessCalculation::~MultiThreadedPopulationFitnessCalculation()
{
}

bool_t MultiThreadedPopulationFitnessCalculation::initThreadPool(const std::vector<std::shared_ptr<GenomeFitnessCalculation>> &threadGenomeCalculations)
{
	m_threadGenomeCalculations = threadGenomeCalculations;
	if (m_threadGenomeCalculations.size() < 1)
		return "Need at least one thread";

	// TODO: for now we're not using a real thread pool, for simplicity we'll just
	//       spawn threads when needed
	return true;
}

class Countdown
{
public:
	Countdown(size_t n) : m_num (n) { }
	~Countdown() { }
    void wait()
	{
		unique_lock<mutex> l(m_mut);
		if (m_num == 0)
			return;
		m_cond.wait(l);
	}

	void decrement()
	{
		lock_guard<mutex> l(m_mut);
		if (m_num == 0)
			return;
		
		m_num--;
		if (m_num == 0)
			m_cond.notify_all();
	}
private:
	size_t m_num;
    condition_variable m_cond;
    mutex m_mut;
};

errut::bool_t MultiThreadedPopulationFitnessCalculation::check(const std::vector<std::shared_ptr<Population>> &populations)
{
	if (m_threadGenomeCalculations.size() < 1)
		return "Need at least one thread";
	return true;
}

bool_t MultiThreadedPopulationFitnessCalculation::calculatePopulationFitness(const vector<shared_ptr<Population>> &populations)
{
	// TODO: use actual thread pool
	if (m_threadGenomeCalculations.size() < 1)
		return "Need at least one thread";

	// A very simple approach
	auto helperThread = [populations](size_t id, size_t numThreads, GenomeFitnessCalculation *pGenomeCalc,
									  int *pError, string *pErrStr, Countdown *pBarrier)
	{
		vector<pair<Genome *, Fitness *>> work;

		int count = 0;
		for (auto &pop : populations)
		{
			for (auto &ind : pop->m_individuals)
			{
				if (!ind->fitnessRef().isCalculated())
				{
					if (count%numThreads == id)
					{
						Genome *pGenome = ind->genomePtr();
						Fitness *pFitness = ind->fitnessPtr();

						bool_t r = pGenomeCalc->startNewCalculation(*pGenome);
						if (!r)
						{
							*pError = 1;
							*pErrStr = "Error starting genome calculation: " + r.getErrorString();
							return;
						}
						// fprintf(stderr, "Picking up work for %d on thread %d\n", count, id);
						work.push_back({ pGenome, pFitness });
					}
					// else
					// 	fprintf(stderr, "Skipping work for %d on thread %d\n", count, id);

					count++;
				}
			}
			// fprintf(stderr, "End count on thread %d is %d\n", id, count);
		}

		if (count == 0) // Nothing needed to be done, will be detected by all threads
		{
			// cerr << "Nothing to do in thread, bailing early\n";
			return;
		}

		// Need sync here, as we're changing the isCalculated flags
		pBarrier->decrement();
		pBarrier->wait();

		// Calculate until all is done
		bool allCalculated;
		do
		{
			allCalculated = true;
			for (auto &gf : work)
			{
				Fitness &f = *gf.second;
				if (!f.isCalculated())
				{
					Genome &g = *gf.first;
					auto r = pGenomeCalc->pollCalculate(g, f);
					if (!r)
					{
						*pError = 1;
						*pErrStr = "Error calculating fitness for genome: " + r.getErrorString();
						return;
					}

					if (!f.isCalculated())
						allCalculated = false;
				}
			}
		} while (!allCalculated);
	};

	vector<thread> helpers(m_threadGenomeCalculations.size());
	vector<int> errors(helpers.size(), 0);
	vector<string> errorMessages(helpers.size());
	Countdown barrier(helpers.size());

	for (size_t id = 0 ; id < helpers.size() ; id++)
		helpers[id] = thread(helperThread, id, helpers.size(), m_threadGenomeCalculations[id].get(),
		                     &errors[id], &errorMessages[id], &barrier);
	
	for (auto &t : helpers)
		t.join();

	for (size_t i = 0 ; i < errors.size() ; i++)
	{
		if (errors[i])
			return errorMessages[i];
	}

	return true;
}
