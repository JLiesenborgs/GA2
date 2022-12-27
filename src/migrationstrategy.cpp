#include "migrationstrategy.h"
#include <numeric>

using namespace std;
using namespace errut;

namespace eatk
{

BasicMigrationStrategy::BasicMigrationStrategy(const shared_ptr<MigrationGenerationCheck> &check,
                                               const shared_ptr<MigrationIndividualExchange> &exchange)
	: m_check(check),
	  m_exchange(exchange)
{
}

BasicMigrationStrategy::~BasicMigrationStrategy()
{
}

bool_t BasicMigrationStrategy::check(const vector<shared_ptr<Population>> &populations)
{
	return true;
}

bool_t BasicMigrationStrategy::migrate(size_t generation, vector<shared_ptr<Population>> &populations)
{
	bool_t r;

	assert(m_check.get());
	assert(m_exchange.get());

	bool migrate;
	if (!(r = m_check->shouldMigrateThisGeneration(generation, migrate)))
		return "Error checking migration requirement: " + r.getErrorString();

	if (migrate)
	{
		if (!(r = m_exchange->exchange(generation, populations)))
			return "Error exchanging individuals: " + r.getErrorString();
	}

	return true;
}

UniformProbabilityMigrationCheck::UniformProbabilityMigrationCheck(const shared_ptr<RandomNumberGenerator> &rng, float fractionOfGenerations,
	                                 size_t generationGracePeriod)
	: m_rng(rng),
	  m_fraction(fractionOfGenerations),
	  m_skipGenerations(generationGracePeriod)
{
}

UniformProbabilityMigrationCheck::~UniformProbabilityMigrationCheck()
{
}

bool_t UniformProbabilityMigrationCheck::shouldMigrateThisGeneration(size_t generation, bool &migrate)
{
	migrate = false;

	if (generation <= m_skipGenerations)
		return true;

	float x = m_rng->getRandomFloat();
	if (x < m_fraction)
		migrate = true;

	return true;
}

SequentialRandomIndividualExchange::SequentialRandomIndividualExchange(const shared_ptr<RandomNumberGenerator> &rng, size_t iterations)
	: m_rng(rng),
	  m_iterations(iterations)
{
}

SequentialRandomIndividualExchange::~SequentialRandomIndividualExchange()
{
}

bool_t SequentialRandomIndividualExchange::exchangeIteration(size_t generation, vector<shared_ptr<Population>> &populations)
{
	if (populations.size() < 2)
		return "Not enough populations (" + to_string(populations.size()) + "), cannot exchange individuals";

	// Note: not using a member for this, because we don't want to hang
	// on to the shared_ptr instances too long
	vector<pair<shared_ptr<Individual>, uint32_t>> srcIndividuals(populations.size());

	// Get the source individuals
	for (size_t p = 0 ; p < populations.size() ; p++)
	{
		auto &population = populations[p];
		uint32_t num = (uint32_t)population->size();
		if (population->size() == 0)
			return "Population " + to_string(p) + " size is zero!";

		uint32_t idx = m_rng->getRandomUint32()%num;
		srcIndividuals[p] = { population->individual(idx), idx };
	}

	// Get the population ids to migrate to
	

	vector<uint32_t> destPop(populations.size());

	/*
	auto printDestPop = [&destPop]()
	{
		cout << "Destination populations:";
		for (auto x : destPop)
			cout << " " << x;
		cout << endl;
	};
	*/

	if (populations.size() == 2)
	{
		destPop[0] = 1;
		destPop[1] = 0;
	}
	else
	{
		std::iota(destPop.begin(), destPop.end(), 0);

		// printDestPop();

		// Fisher Yates shuffle for 'derangement'
		for (uint32_t p = 0 ; p < (uint32_t)populations.size()-1 ; p++)
		{
			uint32_t i = (uint32_t)populations.size()-1-p;
			uint32_t j = m_rng->getRandomUint32() % i;

			std::swap(destPop[i], destPop[j]);
		}
	}

#ifndef NDEBUG
	vector<int> srcPopCounts(populations.size(), 0);
	vector<int> dstPopCounts(populations.size(), 0);
#endif
	
	//printDestPop();

	// Migrate!
	for (size_t p = 0 ; p < populations.size() ; p++)
	{
		auto &srcInfo = srcIndividuals[p];
		uint32_t srcPopIdx = srcInfo.second;

		size_t srcPop = p;
		size_t dstPop = destPop[p];

		assert(srcPop != dstPop);

		// See who we're replacing, get that position in the array
		
		uint32_t dstPopIdx = srcIndividuals[dstPop].second;
		populations[dstPop]->individual(dstPopIdx) = srcInfo.first;

		onExchange(generation, srcPop, srcPopIdx, dstPop, dstPopIdx);

		srcPopCounts[srcPop]++;
		dstPopCounts[dstPop]++;
	}

#ifndef NDEBUG
	for (size_t i = 0 ; i < populations.size() ; i++)
	{
		assert(srcPopCounts[i] == 1);
		assert(dstPopCounts[i] == 1);
	}
#endif 

	return true;
}

bool_t SequentialRandomIndividualExchange::exchange(size_t generation, vector<shared_ptr<Population>> &populations)
{
	bool_t r;
	for (size_t it = 0 ; it < m_iterations ; it++)
		if (!(r = exchangeIteration(generation, populations)))
			return r;

	return true;
}

}
