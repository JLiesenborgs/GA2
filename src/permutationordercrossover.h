#pragma once

#include "mogal2config.h"
#include "crossovermutation.h"
#include "randomnumbergenerator.h"

namespace mogal2
{

class PermutationOrderCrossover : public GenomeCrossover
{
public:
    PermutationOrderCrossover(const std::shared_ptr<RandomNumberGenerator> &rng, bool twoOffspring);

	errut::bool_t check(const std::vector<std::shared_ptr<Genome>> &parents) override;
	errut::bool_t generateOffspring(const std::vector<std::shared_ptr<Genome>> &parents,
	                         std::vector<std::shared_ptr<Genome>> &offspring) override;
private:
    std::shared_ptr<RandomNumberGenerator> m_rng;
    std::vector<bool> m_hasIndex;
    bool m_twoOffspring;
};

}
