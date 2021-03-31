#include "fitnessbasedduplicateremoval.h"

using namespace std;
using namespace errut;

namespace eatk
{

FitnessBasedDuplicateRemoval::FitnessBasedDuplicateRemoval(const shared_ptr<FitnessComparison> &fitComp, size_t numObjectives)
    : m_cmp(fitComp), m_numObjectives(numObjectives)
{
}

FitnessBasedDuplicateRemoval::~FitnessBasedDuplicateRemoval()
{

}

bool_t FitnessBasedDuplicateRemoval::check(const vector<shared_ptr<Individual>> &individuals)
{
    return true;
}

bool_t FitnessBasedDuplicateRemoval::removeDuplicates(vector<shared_ptr<Individual>> &individuals)
{
    vector<shared_ptr<Individual>> undoubled;

    for (auto &i : individuals)
    {
        bool gotDouble = false;
        for (auto &i2 : undoubled)
        {
            size_t count = 0;
            for (size_t c = 0 ; c < m_numObjectives ; c++)
            {
                if ((!m_cmp->isFitterThan(i->fitnessRef(), i2->fitnessRef(), c)) && 
                    (!m_cmp->isFitterThan(i2->fitnessRef(), i->fitnessRef(), c)) )
                    count++;
            }
            if (count == m_numObjectives)
                gotDouble = true;
        }

        if (!gotDouble)
            undoubled.push_back(i);
    }
    swap(individuals, undoubled);
    return true;
}

}
