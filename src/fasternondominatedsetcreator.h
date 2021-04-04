#pragma once

#include "eatkconfig.h"
#include "nondominatedsetcreator.h"
#include "genomefitness.h"

namespace eatk
{

class FasterNonDominatedSetCreator : public NonDominatedSetCreator
{
public:
	FasterNonDominatedSetCreator(const std::shared_ptr<FitnessComparison> &fitCmp, size_t numObjectives);
	~FasterNonDominatedSetCreator();

	errut::bool_t calculateAllNDSets(const std::vector<std::shared_ptr<Individual>> &individuals) override;
	errut::bool_t calculateNonDomitatedSet(const std::vector<std::shared_ptr<Individual>> &individuals,
		std::vector<std::shared_ptr<Individual>> &ndSet,
		std::vector<std::shared_ptr<Individual>> &remaining) override;

	// Assumes that both are in fact ND sets
	errut::bool_t mergeNDSets(std::vector<std::shared_ptr<Individual>> &inOut, const std::vector<std::shared_ptr<Individual>> &added) override;

	size_t getNumberOfSets() const override { return m_sets.size(); }
	const std::vector<std::shared_ptr<Individual>> &getSet(size_t idx) const override { return m_sets[idx]; }
	std::vector<std::shared_ptr<Individual>> &getSet(size_t idx) override { return m_sets[idx]; }
private:
    template<typename Update> void buildCounts(const std::vector<std::shared_ptr<Individual>> &individuals, Update update);
    void checkDominated(const Fitness &fi, const Fitness &fj, bool &jDominatedByI, bool &iDominatedByJ);

	size_t m_numObjectives;
	std::shared_ptr<FitnessComparison> m_cmp;
	
	std::vector<std::vector<std::shared_ptr<Individual>>> m_sets;

    std::vector<int> m_dominatedCount;
    std::vector<int> m_dominatedCount2;
    std::vector<std::vector<size_t>> m_dominatesList;

    std::vector<std::shared_ptr<Individual>> m_tmpInd;
    std::vector<size_t> m_tmpSelected;
};

inline void FasterNonDominatedSetCreator::checkDominated(const Fitness &fi, const Fitness &fj,
                                    bool &jDominatedByI, bool &iDominatedByJ)
{
    FitnessComparison &cmp = *m_cmp;
    size_t N = m_numObjectives;

    size_t betterCount = 0;
    size_t betterCount2 = 0;
    size_t betterEqualCount = 0;
    size_t betterEqualCount2 = 0;

    for (size_t k = 0 ; k < N ; k++)
    {
        if (cmp.isFitterThan(fj, fi, k)) // j is fitter than i
        {
            betterCount++;
            betterEqualCount++;
        }
        else // j is not fitter than i
        {
            if (!cmp.isFitterThan(fi, fj, k)) // i is not fitter than j
            {
                betterEqualCount++;
                betterEqualCount2++;
            }
            else // i is fitter than j
            {
                betterCount2++;
                betterEqualCount2++;
            }
        }
    }

    iDominatedByJ = (betterEqualCount == N && betterCount > 0)?true:false;
    jDominatedByI = (betterEqualCount2 == N && betterCount2 > 0)?true:false;
}

}
