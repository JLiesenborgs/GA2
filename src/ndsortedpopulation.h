#pragma once

#include "mogal2config.h"
#include "crossovermutation.h"

namespace mogal2
{

class NDSortedPopulation : public SelectionPopulation
{
public:
    NDSortedPopulation() { }
    ~NDSortedPopulation() { }

    virtual int getNumberOfSets() const = 0;
    virtual int getSetSize(int s) const = 0;
    virtual std::shared_ptr<Individual> getIndividual(int s, int i) const = 0;
};

}
