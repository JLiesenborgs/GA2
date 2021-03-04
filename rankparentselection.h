#include "crossovermutation.h"
#include "randomnumbergenerator.h"

class RankParentSelection : public ParentSelection
{
public:
    RankParentSelection(double beta, std::shared_ptr<RandomNumberGenerator> rng);
    ~RankParentSelection();

	errut::bool_t check(const SelectionPopulation &pop) override;
    errut::bool_t selectParents(const SelectionPopulation &pop, std::vector<std::shared_ptr<Genome>> &parents) override;
private:
    const double m_beta;
    std::shared_ptr<RandomNumberGenerator> m_rng;
};
