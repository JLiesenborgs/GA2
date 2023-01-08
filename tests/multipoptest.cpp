#include "evolutionaryalgorithm.h"
#include "mersennerandomnumbergenerator.h"
#include "valuefitness.h"
#include "vectorgenomefitness.h"
#include "singlepopulationcrossover.h"
#include "rankparentselection.h"
#include "simplesortedpopulation.h"
#include "singlebestelitism.h"
#include "remainingtargetpopulationsizeiteration.h"
#include "singlethreadedpopulationfitnesscalculation.h"
#include "uniformvectorgenomecrossover.h"
#include "vectorgenomeflipmutation.h"
#include "multipopulationevolver.h"
#include "fasternondominatedsetcreator.h"
#include "fitnessbasedduplicateremoval.h"
#include <iostream>
#include <limits>

using namespace errut;
using namespace std;
using namespace eatk;

class Creation : public IndividualCreation
{
public:
	Creation(const shared_ptr<RandomNumberGenerator> rng, size_t numItems)
	 : m_rng(rng), m_numItems(numItems) { }

	shared_ptr<Genome> createInitializedGenome() override
	{
		auto g = make_shared<VectorGenome<int>>(m_numItems);
		
		for (auto &x : g->getValues())
			x = (m_rng->getRandomDouble() < 0.5)?0:1;

		return g;
	}
	
	shared_ptr<Fitness> createEmptyFitness() override
	{
		return make_shared<ValueFitness<double>>();
	}
private:
	shared_ptr<RandomNumberGenerator> m_rng;
	size_t m_numItems;
};

class Item
{
public:
	Item(double val, double w) : value(val), weight(w) { }
	double value, weight;
};

class MyStop : public FixedGenerationsStopCriterion
{
public:
	MyStop(size_t n) : FixedGenerationsStopCriterion(n) { }
	bool_t analyze(const PopulationEvolver &evolver, size_t generationNumber, bool &shouldStop) override
	{
		auto &currentBest = evolver.getBestIndividuals();
		if (currentBest.size() > 0)
		{
			assert(currentBest.size() == 1);
			cout << generationNumber << "| " << currentBest[0]->toString() << endl;
		}
		return FixedGenerationsStopCriterion::analyze(evolver, generationNumber, shouldStop);
	}
};

class MyGA : public EvolutionaryAlgorithm
{
protected:
	bool_t onAlgorithmDone(size_t generation, const std::vector<std::shared_ptr<Individual>> &bestIndividuals) override
	{
		cout << "Ending after " << generation << " generations, best are: " << endl;
		for (auto &i : bestIndividuals)
			cout << i->toString() << endl;
		return true;
	}
};

class KnapsackFitnessCalculation : public GenomeFitnessCalculation
{
public:
	KnapsackFitnessCalculation(const vector<Item> &items, double weightLimit)
	 : m_items(items), m_weightLimit(weightLimit) { }
	bool_t calculate(const Genome &genome, Fitness &fitness) override
	{
		const VectorGenome<int> &vg = static_cast<const VectorGenome<int> &>(genome);
		const vector<int> &indices = vg.getValues();
		ValueFitness<double> &vf = static_cast<ValueFitness<double> &>(fitness);

		double totalValue = 0;
		double totalWeight = 0;

		assert(m_items.size() == vg.getValues().size());
		for (size_t i = 0 ; i < vg.getValues().size() ; i++)
		{
			totalValue += vg.getValues()[i]*m_items[i].value;
			totalWeight += vg.getValues()[i]*m_items[i].weight;
		}

		if (totalWeight > m_weightLimit)
		{
			// In case all genomes are bad, we do want some incentive to go below the weight
			// limit
			totalValue = -(totalWeight-m_weightLimit);
		}
		vf.setValue(totalValue);
		return true;
	}
private:
	vector<Item> m_items;
	double m_weightLimit;
};

class MyExchange : public SequentialRandomIndividualExchange
{
public:
	MyExchange(const std::shared_ptr<RandomNumberGenerator> &rng, size_t iterations) : SequentialRandomIndividualExchange(rng, iterations) { }
protected:
	void onExchange(size_t generation, size_t srcPop, size_t srcIndividualIdx, size_t dstPop, size_t dstIndividualIdx) override
	{
		cout << "Generation " << generation << ": migrating " << srcIndividualIdx << " from pop " << srcPop << " to " << dstIndividualIdx << " in pop " << dstPop << endl;
	}
};

int main(int argc, char const *argv[])
{
	size_t numItems = 28;
	vector<Item> items;
	double targetWeight = 0; // fill in later

	// Initialize the items randomly, but with fixed seeds to that different runs
	// use the same items
	{
		auto rng0 = make_shared<MersenneRandomNumberGenerator>(1337);
		double allWeight = 0;
		
		for (size_t i = 0 ; i < numItems ; i++)
		{
			double value = (int)(rng0->getRandomDouble()*100.0)+1;
			double weight = (int)(rng0->getRandomDouble()*10.0)+1;

			items.push_back({value, weight});
			allWeight += weight;
		}

		targetWeight = allWeight/4.0;
	}

	random_device rd;
	unsigned int seed = rd();
	auto rng = make_shared<MersenneRandomNumberGenerator>(seed);

	auto creation = make_shared<Creation>(rng, items.size());
	auto calcSingle = make_shared<SingleThreadedPopulationFitnessCalculation>(make_shared<KnapsackFitnessCalculation>(items, targetWeight));
	auto fitnessComp = make_shared<ValueFitnessComparison<double,false>>();

	vector<size_t> popSizes { 32, 32, 32 };
	vector<shared_ptr<PopulationEvolver>> evolvers;

	// Just make several versions with same parameters for each population
	for (size_t i = 0 ; i < popSizes.size() ; i++)
	{
		auto mutation = make_shared<VectorGenomeFlipMutation<int>>(rng, 1.0/items.size());
		auto cross = make_shared<SinglePopulationCrossover>(0.1, false,
				make_shared<SimpleSortedPopulation>(fitnessComp),
				make_shared<RankParentSelection>(2.5, rng),
				make_shared<UniformVectorGenomeCrossover<int>>(rng, false),
				mutation,
				make_shared<SingleBestElitism>(true, mutation),
				make_shared<RemainingTargetPopulationSizeIteration>(),
				rng
			);
		evolvers.push_back(cross);
	}

	//shared_ptr<BestIndividualMerger> merger = make_shared<SingleObjectiveBestIndividualMerger>(fitnessComp);

	shared_ptr<NonDominatedSetCreator> ndCreator = make_shared<FasterNonDominatedSetCreator>(fitnessComp, 1);
	shared_ptr<FitnessBasedDuplicateRemoval> dupRemoval = make_shared<FitnessBasedDuplicateRemoval>(fitnessComp, 1);
	shared_ptr<BestIndividualMerger> merger = make_shared<MultiObjectiveBestIndividualMerger>(ndCreator, dupRemoval);
	MultiPopulationEvolver multiPopEvolver(evolvers, merger);
	
	MyStop stop(1000);
	MyGA ga;

	auto migrationCheck = make_shared<UniformProbabilityMigrationCheck>(rng, 0.2f, 50);
	auto migrationExchange = make_shared<MyExchange>(rng, 1);

	BasicMigrationStrategy migration(migrationCheck, migrationExchange);
	
	auto r = ga.run(*creation, multiPopEvolver, *calcSingle, stop, migration, popSizes);
	if (!r)
	{
		cerr << "Error: " << r.getErrorString() << endl;
		return -1;
	}
	cout << "Done" << endl;
	return 0;
}

