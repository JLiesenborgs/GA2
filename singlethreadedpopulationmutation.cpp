#include "singlethreadedpopulationmutation.h"

using namespace std;
using namespace errut;

SingleThreadedPopulationMutation::SingleThreadedPopulationMutation(shared_ptr<GenomeMutation> mutation)
    : m_mutation(mutation)
{
}

SingleThreadedPopulationMutation::~SingleThreadedPopulationMutation()
{
}

bool_t SingleThreadedPopulationMutation::check(const vector<shared_ptr<Population>> &populations)
{
    for (auto &pop : populations)
    {
        for (auto &i : pop->m_individuals)
        {
            bool_t r = m_mutation->check(*i->m_genome);
            if (!r)
                return "Error in mutation check: " + r.getErrorString();
        }
    }

    return true;
}

bool_t SingleThreadedPopulationMutation::mutate(const vector<shared_ptr<Population>> &populations)
{
    for (auto &pop : populations)
    {
        for (auto &i : pop->m_individuals)
        {
            bool isChanged = false;
            bool_t r = m_mutation->mutate(*i->m_genome, isChanged);
            if (!r)
                return "Error in mutation: " + r.getErrorString();
            
            if (isChanged)
                i->m_fitness->setCalculated(false);
        }
    }
    return true;
}