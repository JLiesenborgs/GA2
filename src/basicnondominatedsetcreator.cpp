#include "basicnondominatedsetcreator.h"

using namespace std;
using namespace errut;

namespace eatk
{

BasicNonDominatedSetCreator::BasicNonDominatedSetCreator(const shared_ptr<FitnessComparison> &fitCmp, size_t numObjectives)
	 : m_cmp(fitCmp), m_numObjectives(numObjectives)
{
}

BasicNonDominatedSetCreator::~BasicNonDominatedSetCreator()
{
}

bool_t BasicNonDominatedSetCreator::mergeNDSets(vector<shared_ptr<Individual>> &inOut, const vector<shared_ptr<Individual>> &added)
{
	// For this simple version, we'll just calculate the ND set of the joined vector
	vector<shared_ptr<Individual>> all = inOut;
	for (auto &i : added)
		all.push_back(i);
	
	inOut.clear();
	return calculateNonDomitatedSet(all, inOut, m_tmpND);
}

bool_t BasicNonDominatedSetCreator::calculateAllNDSets(const vector<shared_ptr<Individual>> &individuals)
{
	m_sets.clear();

	const vector<shared_ptr<Individual>> *pIn = &individuals;
	vector<shared_ptr<Individual>> *pND = &m_tmpND;
	vector<shared_ptr<Individual>> *pRem = &m_tmpRem[0];
	size_t rIdx = 0;
	bool_t r;

	while (true)
	{
		if (pIn->size() == 0)
			break;

		pND->clear();
		pRem->clear();

		if (!(r = calculateNonDomitatedSet(*pIn, *pND, *pRem)))
			return "Error calulating one of the sets: " + r.getErrorString();

		size_t newIdx = m_sets.size();
		m_sets.resize(newIdx+1);
		swap(m_sets[newIdx], *pND);

		pIn = pRem;
		rIdx++;
		pRem = &m_tmpRem[rIdx%2];
	}
	return true;
}

bool_t BasicNonDominatedSetCreator::calculateNonDomitatedSet(const vector<shared_ptr<Individual>> &individuals,
												   vector<shared_ptr<Individual>> &ndSet,
												   vector<shared_ptr<Individual>> &remaining)
{
	ndSet.clear();
	remaining.clear();

	const size_t N = m_numObjectives;

	auto dominatesFunction = FitnessComparison::getDominatesFunction(m_cmp, 0, N);
	auto isDominated = [&dominatesFunction](auto &i, auto &j)
	{
		return dominatesFunction(j->fitnessRef(), i->fitnessRef());
	};

	for (auto &i : individuals)
	{
		bool found = false;

		for (auto &j : individuals)
		{
			if (i.get() == j.get())
				continue;

			if (isDominated(i, j))
			{
				found = true;
				break;
			}
		}

		if (!found)
			ndSet.push_back(i);
		else
			remaining.push_back(i);
	}
	return true;
}

}
