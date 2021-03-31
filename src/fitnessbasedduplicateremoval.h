#pragma once

#include "eatkconfig.h"
#include "duplicateindividualremoval.h"
#include "genomefitness.h"

namespace eatk
{

class FitnessBasedDuplicateRemoval : public DuplicateIndividualRemoval
{
public:
    FitnessBasedDuplicateRemoval(const std::shared_ptr<FitnessComparison> &fitComp, size_t numObjectives);
    ~FitnessBasedDuplicateRemoval();

    errut::bool_t check(const std::vector<std::shared_ptr<Individual>> &individuals) override;
    errut::bool_t removeDuplicates(std::vector<std::shared_ptr<Individual>> &individuals) override;
private:
    std::shared_ptr<FitnessComparison> m_cmp;
    size_t m_numObjectives;
};

}