#pragma once

#include "eatkconfig.h"
#include "population.h"
#include "randomnumbergenerator.h"

namespace eatk
{

class MigrationStrategy
{
public:
	MigrationStrategy() { }
	virtual ~MigrationStrategy() { }

	virtual errut::bool_t check(const std::vector<std::shared_ptr<Population>> &populations) { return "Not implemented in base class"; } 
	// Note that this is run after the generation has been updated (so generation will
	// be at least 1), and when
	// all fitness values have been calculated, in case the migration strategy
	// wants to take this into account. The individuals are not sorted in any
	// way at this point
	virtual errut::bool_t migrate(size_t generation, std::vector<std::shared_ptr<Population>> &populations) { return "Not implemented in base class"; }
};

class MigrationGenerationCheck
{
public:
	MigrationGenerationCheck() { }
	virtual ~MigrationGenerationCheck() { }

	virtual errut::bool_t shouldMigrateThisGeneration(size_t generation, bool &migrate) { return "Not implemented in base class"; }
};

class UniformProbabilityMigrationCheck : public MigrationGenerationCheck
{
public:
	UniformProbabilityMigrationCheck(const std::shared_ptr<RandomNumberGenerator> &rng, float fractionOfGenerations,
	                                 size_t generationGracePeriod);
	~UniformProbabilityMigrationCheck();

	errut::bool_t shouldMigrateThisGeneration(size_t generation, bool &migrate);
private:
	std::shared_ptr<RandomNumberGenerator> m_rng;
	size_t m_skipGenerations;
	float m_fraction;
};

class MigrationIndividualExchange
{
public:
	MigrationIndividualExchange() { }
	virtual ~MigrationIndividualExchange() { }

	virtual errut::bool_t exchange(std::vector<std::shared_ptr<Population>> &populations) { return "Not implemented in base class"; }
};

class BasicMigrationStrategy : public MigrationStrategy
{
public:
	BasicMigrationStrategy(const std::shared_ptr<MigrationGenerationCheck> &check,
	                       const std::shared_ptr<MigrationIndividualExchange> &exchange);
	~BasicMigrationStrategy();

	errut::bool_t check(const std::vector<std::shared_ptr<Population>> &populations) override;
	errut::bool_t migrate(size_t generation, std::vector<std::shared_ptr<Population>> &populations) override;
private:
	std::shared_ptr<MigrationGenerationCheck> m_check;
	std::shared_ptr<MigrationIndividualExchange> m_exchange;
};

}
