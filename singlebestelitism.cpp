#include "singlebestelitism.h"
#include <string>
#include <cassert>

using namespace errut;
using namespace std;

SingleBestElitism::SingleBestElitism(bool eliteWithoutMutation, bool eliteWithMutation)
    : m_eliteWithoutMutation(eliteWithoutMutation),
      m_eliteWithMutation(eliteWithMutation)
{
}

SingleBestElitism::~SingleBestElitism()
{
}

bool_t SingleBestElitism::check(const shared_ptr<SelectionPopulation> &selPop)
{
    if (!m_eliteWithMutation && !m_eliteWithoutMutation)
        return "No elitism was enabled";
    return true;
}

bool_t SingleBestElitism::introduceElites(const shared_ptr<SelectionPopulation> &selPop,
                                shared_ptr<Population> &population,
                                size_t targetPopulationSize)
{
    if (population->size() != 0)
        return "Expecting empty population, but got size " + to_string(population->size());
    
    const vector<shared_ptr<Individual>> &best = selPop->getBestIndividuals();
    if (best.size() == 0)
    {
        population->setGenomesToSkipMutation(0);
        return true;    
    }

    if (best.size() > 1)
        return "Expecting a single best genome, but got " + to_string(best.size());
    
    assert(best.front().get());
    const Individual &i = *best.front();

    auto copyBest = [&i, &population]()
    {
        population->append(i.createCopy());
    };

    if (m_eliteWithoutMutation)
    {
        copyBest();
        population->setGenomesToSkipMutation(1);
    }
    if (m_eliteWithMutation)
        copyBest();

    return true;
}