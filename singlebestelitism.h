#pragma once

#include "crossovermutation.h"

class SingleBestElitism : public Elitism
{
public:
    SingleBestElitism(bool eliteWithoutMutation, const std::shared_ptr<GenomeMutation> &mutation);
    ~SingleBestElitism();

    errut::bool_t check(const std::shared_ptr<SelectionPopulation> &selPop) override;
	errut::bool_t introduceElites(const std::shared_ptr<SelectionPopulation> &selPop,
								  std::shared_ptr<Population> &population,
								  size_t targetPopulationSize) override;
private:
    bool m_eliteWithoutMutation;
    std::shared_ptr<GenomeMutation> m_mutation;
};
