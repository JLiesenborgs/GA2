#include "singlethreadedpopulationcrossover.h"
#include <cassert>
#include <iostream>

using namespace errut;
using namespace std;

SingleThreadedPopulationCrossover::SingleThreadedPopulationCrossover(double cloneFraction,
                                    shared_ptr<SelectionPopulation> selectionPop,
                                    shared_ptr<ParentSelection> parentSelection,
                                    shared_ptr<GenomeCrossover> genomeCrossover,
                                    shared_ptr<RandomNumberGenerator> rng)
    : m_cloneFraction(cloneFraction), m_selectionPop(selectionPop), 
      m_parentSelection(parentSelection), m_genomeCrossover(genomeCrossover),
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

    // TODO: check crossover, check parent selection
    return true;
}

bool_t SingleThreadedPopulationCrossover::createNewPopulation(vector<shared_ptr<Population>> &populations,
                                                              int targetPopulationSize)
{
    bool_t r;
    vector<shared_ptr<Genome>> parents, cloneParent, offspring;
    
    parents.resize(m_genomeCrossover->getNumberOfParents());
    cloneParent.resize(1);

    for (auto &population : populations)
    {
        assert(population->m_individuals.size() > 0);
        auto &refFitness = population->m_individuals[0]->m_fitness;

        // Here, some pruning could take place
        if (!(r = m_selectionPop->processPopulation(population, targetPopulationSize)))
            return "Error in selection preprocessing: " + r.getErrorString();

        auto newPopulation = make_shared<Population>();
        
        // TODO: allow something else here? Different sizes?

        const size_t popSize = (int)population->m_individuals.size();
        for (size_t i = 0 ; i < popSize ; i++)
        {
            double x = m_rng->getRandomDouble();
            if (x < m_cloneFraction) // TODO: can we do this more efficiently?
            {
                if (!(r = m_parentSelection->selectParents(*m_selectionPop, cloneParent)))
                    return "Error in clone parent selection: " + r.getErrorString();

                auto f = refFitness->createCopy(false);
                auto g = cloneParent[0]->createCopy(true);
                auto ind = make_shared<Individual>(g, f);
                newPopulation->m_individuals.push_back(ind);
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
                    newPopulation->m_individuals.push_back(ind);
                }                
            }
        }

        std::swap(newPopulation, population);
    }

    return true;
}
