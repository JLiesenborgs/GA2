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

bool_t SingleThreadedPopulationCrossover::createNewPopulations(vector<shared_ptr<Population>> &populations)
{
    bool_t r;
    vector<shared_ptr<Genome>> parents, offspring;

    for (auto &population : populations)
    {
        assert(population->m_individuals.size() > 0);
        auto &refFitness = population->m_individuals[0]->m_fitness;

        if (!(r = m_selectionPop->processPopulation(population)))
            return "Error in selection preprocessing: " + r.getErrorString();

        // Clone some
        auto newPopulation = make_shared<Population>();
        for (auto &i : population->m_individuals)
        {
            double x = m_rng->getRandomDouble();
            if (x < m_cloneFraction) // TODO: can we do this more efficiently?
            {
                // TODO: can we just copy the individual?
                //       may not be a good idea: perhaps some other shared_ptr copy is made, and mutation would then affect both
                auto f = refFitness->createCopy(true); // TODO: can we set this to true? Yes, but need care in mutation
                auto g = i->m_genome->createCopy(true);
                auto ind = make_shared<Individual>(g, f);
                newPopulation->m_individuals.push_back(ind);
            }
        }

        // Do crossover
        // TODO: allow something else here? Different sizes?

        const size_t popSize = (int)population->m_individuals.size();
        for (size_t i = newPopulation->m_individuals.size() ; i < popSize ; i++)
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

        std::swap(newPopulation, population);
    }

    return true;
}
