#include "fasternondominatedsetcreator.h"

using namespace errut;
using namespace std;

namespace eatk
{

FasterNonDominatedSetCreator::FasterNonDominatedSetCreator(const shared_ptr<FitnessComparison> &fitCmp, size_t numObjectives)
	: m_cmp(fitCmp), m_numObjectives(numObjectives)
{
}

FasterNonDominatedSetCreator::~FasterNonDominatedSetCreator()
{
}

template<typename Update>
void FasterNonDominatedSetCreator::buildCounts(const std::vector<std::shared_ptr<Individual>> &individuals, Update update)
{
	for (size_t i = 0 ; i < individuals.size()-1 ; i++)
	{
		for (size_t j = i+1 ; j < individuals.size() ; j++)
		{
			const Fitness &fi = individuals[i]->fitnessRef();
			const Fitness &fj = individuals[j]->fitnessRef();

			bool jDominatedByI, iDominatedByJ;
			checkDominated(fi, fj, jDominatedByI, iDominatedByJ);

			if (jDominatedByI) update(i, j);
			if (iDominatedByJ) update(j, i);
		}
	}
}

bool_t FasterNonDominatedSetCreator::calculateAllNDSets(const std::vector<std::shared_ptr<Individual>> &individuals,
                                                        size_t requestStopSize)
{
	if (individuals.size() < 1)
		return "No individuals present";

	m_dominatedCount.resize(individuals.size());
	for (auto &x : m_dominatedCount)
		x = 0;

	m_dominatesList.resize(individuals.size());
	for (auto &dl : m_dominatesList)
		dl.resize(0);

	buildCounts(individuals, [this](int i, int j)
	{
		m_dominatedCount[j]++;
		m_dominatesList[i].push_back(j);
	});

	size_t totalSize = 0;

	m_sets.clear();
	bool done = false;
	while (!done)
	{
		vector<shared_ptr<Individual>> &currentNDSet = m_tmpInd;
		vector<size_t> &selectedGenomes = m_tmpSelected;

		currentNDSet.clear();
		selectedGenomes.clear();

		for (size_t i = 0 ; i < individuals.size() ; i++)
		{
			if (m_dominatedCount[i] == 0)
			{
				m_dominatedCount[i] = -1; // Make sure we skip this next time
				currentNDSet.push_back(individuals[i]);
				selectedGenomes.push_back(i);
			}
		}

		if (selectedGenomes.size() == 0)
			done = true;
		else
		{
			for (auto i : selectedGenomes)
			{
				for (auto idx : m_dominatesList[i])
				{
					// i dominates idx
					// We're going to remove i from the list, so idx is dominated by one less
					m_dominatedCount[idx]--;
					assert(m_dominatedCount[idx] >= 0);
				}
			}

			totalSize += currentNDSet.size();
			if (totalSize >= requestStopSize)
			{
				// cerr << "Req size " << requestStopSize << " reached by totalSize " << totalSize << " for ND set " << m_sets.size() << endl;
				done = true;
			}

			m_sets.push_back({});
			swap(m_sets[m_sets.size()-1], currentNDSet);
		}
	}
	return true;
}

bool_t FasterNonDominatedSetCreator::calculateNonDomitatedSet(const std::vector<std::shared_ptr<Individual>> &individuals,
	std::vector<std::shared_ptr<Individual>> &ndSet,
	std::vector<std::shared_ptr<Individual>> &remaining)
{
	if (individuals.size() < 1)
		return "No individuals present";

	ndSet.clear();
	remaining.clear();

	m_dominatedCount.resize(individuals.size());
	for (auto &x : m_dominatedCount)
		x = 0;

	buildCounts(individuals, [this](int i, int j)
	{
		m_dominatedCount[j]++;
	});

	for (size_t i = 0 ; i < individuals.size() ; i++)
	{
		if (m_dominatedCount[i] == 0)
			ndSet.push_back(individuals[i]);
		else
			remaining.push_back(individuals[i]);
	}
	return true;
}

bool_t FasterNonDominatedSetCreator::mergeNDSets(std::vector<std::shared_ptr<Individual>> &inOut, const std::vector<std::shared_ptr<Individual>> &added)
{
	m_dominatedCount.resize(inOut.size());
	m_dominatedCount2.resize(added.size());
	for (auto &x : m_dominatedCount) x = 0;
	for (auto &x : m_dominatedCount2) x = 0;

	for (size_t i = 0 ; i < inOut.size() ; i++)
	{
		for (size_t j = 0 ; j < added.size() ; j++)
		{
			const Fitness &fi = inOut[i]->fitnessRef();
			const Fitness &fj = added[j]->fitnessRef();

			bool jDominatedByI, iDominatedByJ;
			checkDominated(fi, fj, jDominatedByI, iDominatedByJ);

			if (jDominatedByI) m_dominatedCount2[j]++;
			if (iDominatedByJ) m_dominatedCount[i]++;
		}
	}

	m_tmpInd.resize(0);
	for (size_t i = 0 ; i < inOut.size() ; i++)
		if (m_dominatedCount[i] == 0) // still not dominated
			m_tmpInd.push_back(inOut[i]);

	for (size_t j = 0 ; j < added.size() ; j++)
		if (m_dominatedCount2[j] == 0)
			m_tmpInd.push_back(added[j]);

	swap(inOut, m_tmpInd);

	return true;
}

}
