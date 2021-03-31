#include "tournamentparentselection.h"
#include "simplesortedpopulation.h"
#include <algorithm>

using namespace std;
using namespace errut;

namespace eatk
{

TournamentParentSelection::TournamentParentSelection(const shared_ptr<RandomNumberGenerator> &rng,
													 size_t tournamentSize,
													 const std::shared_ptr<FitnessComparison> &fitCmp,
													 size_t objectiveNumber)
	: m_rng(rng), m_cmp(fitCmp), m_objectiveNumber(objectiveNumber), m_tournamentSize(tournamentSize)
{
}

TournamentParentSelection::~TournamentParentSelection()
{
}

bool_t TournamentParentSelection::check(const SelectionPopulation &pop)
{
	if (m_tournamentSize < 1)
		return "Tournament size must be at least one";

	return true;
}

bool_t TournamentParentSelection::selectParents(const Population &pop,
												const SelectionPopulation &selPop,
												vector<shared_ptr<Individual>> &parents)
{
	FitnessComparison &cmp = *m_cmp;
	const size_t N = m_objectiveNumber;

	auto pickParent = [this, &pop, &cmp, N]()
	{
		uint32_t idx = m_rng->getRandomUint32() % (uint32_t)pop.size();
		auto winner = pop.individual(idx);

		for (size_t i = 1 ; i < m_tournamentSize ; i++)
		{
			idx = m_rng->getRandomUint32() % (uint32_t)pop.size();
			auto &ind = pop.individual(idx);

			if (cmp.isFitterThan(ind->fitnessRef(), winner->fitnessRef(), N))
				winner = ind;
		}

		return winner;
	};

	for (auto &p : parents)
		p = pickParent();
	return true;
}

}