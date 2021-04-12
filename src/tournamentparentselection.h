#pragma once

#include "eatkconfig.h"
#include "crossovermutation.h"
#include "randomnumbergenerator.h"

namespace eatk
{

class TournamentParentSelection : public ParentSelection
{
public:
	TournamentParentSelection(const std::shared_ptr<RandomNumberGenerator> &rng,
							  size_t tournamentSize,
							  const std::shared_ptr<FitnessComparison> &fitCmp,
							  size_t objectiveNumber = 0);
	~TournamentParentSelection();

	errut::bool_t check(const SelectionPopulation &pop) override;
	errut::bool_t selectParents(const Population &population, const SelectionPopulation &pop, std::vector<std::shared_ptr<Individual>> &parents) override;
private:
	size_t m_tournamentSize;
	std::shared_ptr<RandomNumberGenerator> m_rng;
	std::shared_ptr<FitnessComparison> m_cmp;
	size_t m_objectiveNumber;
};

}
