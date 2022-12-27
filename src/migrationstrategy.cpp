#include "migrationstrategy.h"

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
		if (!(r = m_exchange->exchange(populations)))
			return "Error exchanging individuals: " + r.getErrorString();
	}

	return true;
}

UniformProbabilityMigrationCheck::UniformProbabilityMigrationCheck(const std::shared_ptr<RandomNumberGenerator> &rng, float fractionOfGenerations,
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

}
