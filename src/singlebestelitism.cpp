#include "singlebestelitism.h"
#include <string>
#include <cassert>

using namespace errut;
using namespace std;

SingleBestElitism::SingleBestElitism(bool eliteWithoutMutation, const shared_ptr<GenomeMutation> &mutation)
    : m_eliteWithoutMutation(eliteWithoutMutation),
      m_mutation(mutation)
{
}

SingleBestElitism::~SingleBestElitism()
{
}

bool_t SingleBestElitism::check(const shared_ptr<SelectionPopulation> &selPop)
{
    if (!m_mutation.get() && !m_eliteWithoutMutation)
        return "No elitism was enabled";
    // TODO: check mutation?
    return true;
}

bool_t SingleBestElitism::introduceElites(size_t generation, const shared_ptr<SelectionPopulation> &selPop,
                                shared_ptr<Population> &population,
                                size_t targetPopulationSize)
{
    if (population->size() != 0)
        return "Expecting empty population, but got size " + to_string(population->size());
    
    const vector<shared_ptr<Individual>> &best = selPop->getBestIndividuals();
    if (best.size() == 0)
        return true;    

    if (best.size() > 1)
        return "Expecting a single best genome, but got " + to_string(best.size());
    
    assert(best.front().get());
    const Individual &i = *best.front();

    auto copyBest = [&i, &population]()
    {
        auto copy = i.createCopy();
        population->append(copy);
        return copy;
    };

    if (m_eliteWithoutMutation)
        copyBest();

    if (m_mutation.get()) // add copy with mutation
    {
        auto copy = copyBest();
        bool isChanged = false;
        bool_t r = m_mutation->mutate(copy->genomeRef(), isChanged);
        if (isChanged)
        {
            copy->fitness()->setCalculated(false);
            copy->setLastMutationGeneration(generation);
        }
    }
    return true;
}