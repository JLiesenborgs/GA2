#pragma once

#include "population.h"

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

}
