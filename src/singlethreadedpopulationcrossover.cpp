#include "singlethreadedpopulationcrossover.h"
#include <cassert>
#include <iostream>

using namespace errut;
using namespace std;

namespace mogal2
{

SingleThreadedPopulationCrossover::SingleThreadedPopulationCrossover(double cloneFraction,
                                    bool keepExistingPopulation,
                                    shared_ptr<SelectionPopulation> selectionPop,
                                    shared_ptr<ParentSelection> parentSelection,
                                    shared_ptr<GenomeCrossover> genomeCrossover,
                                    shared_ptr<GenomeMutation> genomeMutation,
                                    shared_ptr<Elitism> elitism,
                                    shared_ptr<PopulationCrossoverIteration> popIteration,
                                    shared_ptr<RandomNumberGenerator> rng)
    : m_keepExistingPopulation(keepExistingPopulation),
      m_cloneFraction(cloneFraction), m_selectionPop(selectionPop), 
      m_parentSelection(parentSelection), m_genomeCrossover(genomeCrossover),
      m_genomeMutation(genomeMutation), m_elitism(elitism), m_popIteration(popIteration),
      m_rng(rng)
{
}

SingleThreadedPopulationCrossover::~SingleThreadedPopulationCrossover()
{
}

bool_t SingleThreadedPopulationCrossover::check(const vector<shared_ptr<Population>> &populations)
{
    bool_t r;

    for (auto pop : populations)
    {
        if (!(r = m_selectionPop->check(*pop)))
            return "Error checking selection preprocessing: " + r.getErrorString();
    }
    if (m_elitism.get())
    {
        if (!(r = m_elitism->check(m_selectionPop)))
            return "Error checking compatibility of elitism with selection preprocessing: " + r.getErrorString();
    }
    if (!m_popIteration.get())
        return "No new population iteration was set";

    // TODO: check crossover, mutation, check parent selection
    return true;
}

bool_t SingleThreadedPopulationCrossover::createNewPopulation(size_t generation,
                                                              vector<shared_ptr<Population>> &populations,
                                                              size_t targetPopulationSize)
{
    bool_t r;
    vector<shared_ptr<Genome>> parentGenomes, offspring;
    vector<shared_ptr<Individual>> parentIndividuals, cloneIndividual;
    
    if (m_genomeCrossover->getNumberOfParents() < 2)
        return "Number of parents should be at least 2";

    parentGenomes.resize(m_genomeCrossover->getNumberOfParents());
    parentIndividuals.resize(m_genomeCrossover->getNumberOfParents());
    cloneIndividual.resize(1);

    for (auto &population : populations)
    {
        // Here, some pruning could take place
        if (!(r = m_selectionPop->processPopulation(population, targetPopulationSize)))
            return "Error in selection preprocessing: " + r.getErrorString();

        if (!(r = onSelectionPopulationProcessed(generation, m_selectionPop)))
            return "Error inspecting selection preprocessing: " + r.getErrorString();

        assert(population->size() > 0);
        auto refFitness = population->individual(0)->fitness();

        auto newPopulation = make_shared<Population>();
        
        // introduce elitist solutions ; should itself introduce mutations if needed
        if (m_elitism.get())
        {
            if (!(r = m_elitism->introduceElites(generation, m_selectionPop, newPopulation, targetPopulationSize)))
                return "Can't introduce elitist solutions: " + r.getErrorString();
        }

        auto appendNewIndividual = [this, generation, &refFitness, &newPopulation](auto &i) -> bool_t
        {
            bool isChanged = false;
            bool_t r = m_genomeMutation->mutate(i->genomeRef(), isChanged);
            if (!r)
                return "Error mutating genome: " + r.getErrorString();

            if (isChanged)
            {
                i->fitness()->setCalculated(false);
                i->setLastMutationGeneration(generation);
            }
            newPopulation->append(i);
            return true;
        };

        m_popIteration->startNewIteration(newPopulation, targetPopulationSize);
        while(m_popIteration->iterate(newPopulation))
        {
            double x = m_rng->getRandomDouble();
            if (x < m_cloneFraction) // TODO: can we do this more efficiently?
            {
                if (!(r = m_parentSelection->selectParents(*m_selectionPop, cloneIndividual)))
                    return "Error in clone parent selection: " + r.getErrorString();

                auto i = cloneIndividual[0]->createCopy();
                if (!(r = appendNewIndividual(i)))
                    return r;
            }
            else
            {
                if (!(r = m_parentSelection->selectParents(*m_selectionPop, parentIndividuals)))
                    return "Error in parent selection: " + r.getErrorString();

                for (size_t i = 0 ; i < parentIndividuals.size() ; i++)
                    parentGenomes[i] = parentIndividuals[i]->genome();

                if (!(r = m_genomeCrossover->generateOffspring(parentGenomes, offspring)))
                    return "Error generating offspring: " + r.getErrorString();

                for (auto &g : offspring)
                {
                    auto f = refFitness->createCopy(false);
                    auto i = parentIndividuals[0]->createNew(g, f, generation);
                    if (!(r = appendNewIndividual(i)))
                        return r;
                }                
            }
        }

        if (m_keepExistingPopulation)
        {
            for (auto &i : population->individuals())
                newPopulation->append(i);
            population->clear();
        }

        std::swap(newPopulation, population);
    }

    return true;
}

}
