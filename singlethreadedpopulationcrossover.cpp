#include "singlethreadedpopulationcrossover.h"
#include <cassert>
#include <iostream>

using namespace errut;
using namespace std;

SingleThreadedPopulationCrossover::SingleThreadedPopulationCrossover(double cloneFraction,
                                    shared_ptr<SelectionPopulation> selectionPop,
                                    shared_ptr<ParentSelection> parentSelection,
                                    shared_ptr<GenomeCrossover> genomeCrossover,
                                    shared_ptr<Elitism> elitism,
                                    shared_ptr<PopulationCrossoverIteration> popIteration,
                                    shared_ptr<RandomNumberGenerator> rng)
    : m_cloneFraction(cloneFraction), m_selectionPop(selectionPop), 
      m_parentSelection(parentSelection), m_genomeCrossover(genomeCrossover),
      m_elitism(elitism), m_popIteration(popIteration),
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

    // TODO: check crossover, check parent selection
    return true;
}

bool_t SingleThreadedPopulationCrossover::createNewPopulation(vector<shared_ptr<Population>> &populations,
                                                              size_t targetPopulationSize)
{
    bool_t r;
    vector<shared_ptr<Genome>> parents, cloneParent, offspring;
    
    parents.resize(m_genomeCrossover->getNumberOfParents());
    cloneParent.resize(1);

    for (auto &population : populations)
    {
        // Here, some pruning could take place
        if (!(r = m_selectionPop->processPopulation(population, targetPopulationSize)))
            return "Error in selection preprocessing: " + r.getErrorString();

        if (!(r = onSelectionPopulationProcessed(m_selectionPop)))
            return "Error inspecting selection preprocessing: " + r.getErrorString();

        assert(population->size() > 0);
        auto refFitness = population->individual(0)->fitness();

        auto newPopulation = make_shared<Population>();
        
        // introduce elitist solutions, this also marks how many genomes should not
        // be mutated
        if (m_elitism.get())
        {
            if (!(r = m_elitism->introduceElites(m_selectionPop, newPopulation, targetPopulationSize)))
                return "Can't introduce elitist solutions: " + r.getErrorString();
        }

        m_popIteration->startNewIteration(newPopulation, targetPopulationSize);
        while(m_popIteration->iterate())
        {
            double x = m_rng->getRandomDouble();
            if (x < m_cloneFraction) // TODO: can we do this more efficiently?
            {
                if (!(r = m_parentSelection->selectParents(*m_selectionPop, cloneParent)))
                    return "Error in clone parent selection: " + r.getErrorString();

                auto f = refFitness->createCopy(false);
                auto g = cloneParent[0]->createCopy(true);
                auto ind = make_shared<Individual>(g, f);
                newPopulation->append(ind);
            }
            else
            {
                if (!(r = m_parentSelection->selectParents(*m_selectionPop, parents)))
                    return "Error in parent selection: " + r.getErrorString();

                if (!(r = m_genomeCrossover->generateOffspring(parents, offspring)))
                    return "Error generating offspring: " + r.getErrorString();

                for (auto &g : offspring)
                {
                    auto f = refFitness->createCopy(false);
                    auto ind = make_shared<Individual>(g, f);
                    newPopulation->append(ind);
                }                
            }
        }

        std::swap(newPopulation, population);
    }

    return true;
}
