#pragma once

#include "eatkconfig.h"
#include "crossovermutation.h"
#include "randomnumbergenerator.h"

namespace eatk
{

class PermutationSwapMutation : public GenomeMutation
{
public:
    PermutationSwapMutation(const std::shared_ptr<RandomNumberGenerator> &rng, double prob);

	errut::bool_t check(const Genome &genome) override;
	errut::bool_t mutate(Genome &genome, bool &isChanged) override;
private:
    std::shared_ptr<RandomNumberGenerator> m_rng;
    double m_prob;
};

}
