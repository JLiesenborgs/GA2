#pragma once

#include "crossovermutation.h"

class SingleThreadedPopulationMutation : public PopulationMutation
{
public:
    SingleThreadedPopulationMutation(std::shared_ptr<GenomeMutation> mutation);
    ~SingleThreadedPopulationMutation();

    errut::bool_t check(const std::vector<std::shared_ptr<Population>> &populations) override;
	errut::bool_t mutate(const std::vector<std::shared_ptr<Population>> &populations) override;
private:
    std::shared_ptr<GenomeMutation> m_mutation;
};