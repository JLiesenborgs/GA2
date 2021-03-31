#pragma once

#include "eatkconfig.h"
#include "population.h"

namespace eatk
{

class DuplicateIndividualRemoval
{
public:
    DuplicateIndividualRemoval() { }
    virtual ~DuplicateIndividualRemoval() { }

    virtual errut::bool_t check(const std::vector<std::shared_ptr<Individual>> &individuals) { return "Not implemented in base class"; }
    virtual errut::bool_t removeDuplicates(std::vector<std::shared_ptr<Individual>> &individuals) { return "Not implemented in base class"; }
};

}
