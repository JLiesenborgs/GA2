#pragma once

#include "crossovermutation.h"
#include <iostream>

class ProbeSystem
{
public:
    enum EventType
    {
        FitnessCalculated,
        AlgorithmDone,
        SelectionPreProcessed
    };

    ProbeSystem() : m_generation(0) { }
    virtual ~ProbeSystem() { }

    void setGeneration(size_t n) { m_generation = n; }
    size_t getGeneration() const { return m_generation; }

    virtual errut::bool_t inspect(EventType eventType, const std::vector<std::shared_ptr<Population>> &populations) { return true; }
    virtual errut::bool_t inspect(EventType eventType, std::shared_ptr<Population> &population) { return true; }
    virtual errut::bool_t inspect(EventType eventType, std::shared_ptr<SelectionPopulation> &selPop) { return true; }
    virtual errut::bool_t inspect(EventType eventType, const std::vector<std::shared_ptr<Individual>> &bestIndividuals) { return true; }
private:
    size_t m_generation;
};

