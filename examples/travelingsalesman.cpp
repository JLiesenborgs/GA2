#define _USE_MATH_DEFINES // For M_PI
#include "evolutionaryalgorithm.h"
#include "mersennerandomnumbergenerator.h"
#include "valuefitness.h"
#include "vectorgenomefitness.h"
#include "singlepopulationcrossover.h"
#include "rankparentselection.h"
#include "tournamentparentselection.h"
#include "simplesortedpopulation.h"
#include "singlebestelitism.h"
#include "remainingtargetpopulationsizeiteration.h"
#include "singlethreadedpopulationfitnesscalculation.h"
#include "permutationordercrossover.h"
#include "permutationswapmutation.h"
#include "trackbestonlyselectionpopulation.h"
#include <iostream>

using namespace errut;
using namespace std;
using namespace eatk;

class Creation : public IndividualCreation
{
public:
	Creation(const shared_ptr<RandomNumberGenerator> rng, size_t numCities)
	 : m_rng(rng), m_numCities(numCities) { }

	shared_ptr<Genome> createInitializedGenome() override
	{
		auto g = make_shared<VectorGenome<int>>(m_numCities);
		vector<int> indices(m_numCities);
		for (size_t i = 0 ; i < m_numCities ; i++)
			indices[i] = i;
		
		for (size_t i = 0 ; i < m_numCities ; i++)
		{
			uint32_t idx = m_rng->getRandomUint32() % (uint32_t)indices.size();
			g->getValues()[i] = indices[idx];

			indices[idx] = indices[indices.size()-1];
			indices.resize(indices.size()-1);
		}
		return g;
	}
	
	shared_ptr<Fitness> createEmptyFitness() override
	{
		return make_shared<ValueFitness<double>>();
	}
private:
	shared_ptr<RandomNumberGenerator> m_rng;
	size_t m_numCities;
};

class City
{
public:
	City(double X, double Y) : x(X), y(Y) { }
	double x, y;
};

class MyStop : public FixedGenerationsStopCriterion
{
public:
	MyStop(size_t n) : FixedGenerationsStopCriterion(n) { }
	bool_t analyze(const vector<shared_ptr<Individual>> &currentBest, size_t generationNumber, bool &shouldStop) override
	{
		if (currentBest.size() > 0)
		{
			assert(currentBest.size() == 1);
			cout << generationNumber << "| " << currentBest[0]->toString() << endl;
		}
		return FixedGenerationsStopCriterion::analyze(currentBest, generationNumber, shouldStop);
	}
};

class MyGA : public EvolutionaryAlgorithm
{
public:
	double getBestFitness() const { return static_cast<ValueFitness<double>&>(m_best->fitnessRef()).getValue(); }
	const VectorGenome<int> &getBestGenome() const { return static_cast<const VectorGenome<int>&>(m_best->genomeRef()); }
protected:
	bool_t onAlgorithmDone(size_t generation, const std::vector<std::shared_ptr<Individual>> &bestIndividuals) override
	{
		cout << "Ending after " << generation << " generations, best are: " << endl;
		for (auto &i : bestIndividuals)
			cout << i->toString() << endl;

		m_best = bestIndividuals[0]->createCopy();
		return true;
	}

	shared_ptr<Individual> m_best;
};

class TSPFitnessCalculation : public GenomeFitnessCalculation
{
public:
	TSPFitnessCalculation(const vector<City> &cities) : m_cities(cities) { }
	bool_t calculate(const Genome &genome, Fitness &fitness) override
	{
		const VectorGenome<int> &vg = static_cast<const VectorGenome<int> &>(genome);
		const vector<int> &indices = vg.getValues();
		ValueFitness<double> &vf = static_cast<ValueFitness<double> &>(fitness);

		double totalDist = 0;
		assert(indices.size() == m_cities.size());
		for (size_t i = 0 ; i < m_cities.size() ; i++)
		{
			size_t j = (i+1)%m_cities.size();
			City c1 = m_cities[indices[i]];
			City c2 = m_cities[indices[j]];
			double dx = c1.x - c2.x;
			double dy = c1.y - c2.y;
			totalDist += sqrt(dx*dx + dy*dy);
		}
		vf.setValue(totalDist);
		return true;
	}
private:
	vector<City> m_cities;
};

void usage()
{
	cout << R"XYZ(
travelingsalesman circle numcities

or

travelingsalesman random numcities seed
)XYZ";
	exit(-1);
}

vector<City> getRandomCities(unsigned int seed, size_t numCities)
{
	vector<City> cities;
	MersenneRandomNumberGenerator rng(seed);

	for (size_t i = 0 ; i < numCities ; i++)
	{
		double x = rng.getRandomDouble();
		double y = rng.getRandomDouble();
		cities.push_back({x, y});
	}
	return cities;
}

int main(int argc, char const *argv[])
{
	if (argc < 3)
		usage();

	string type(argv[1]);
	size_t numCities = stoul(argv[2]);
	vector<City> cities;

	if (type == "circle")
	{
		if (argc != 3)
			usage();

		for (size_t i = 0 ; i < numCities ; i++)
		{
			double theta = (double)i/(double)numCities * 2.0*M_PI;
			cities.push_back(City { cos(theta), sin(theta) } );
		}
	}
	else if (type == "random")
	{
		if (argc != 4)
			usage();

		cities = getRandomCities((unsigned int)stoul(argv[3]), numCities);
	}

	random_device rd;
	unsigned int seed = rd();
	auto rng = make_shared<MersenneRandomNumberGenerator>(seed);
	auto creation = make_shared<Creation>(rng, cities.size());
	auto cmp = make_shared<ValueFitnessComparison<double>>();
	//auto selection = make_shared<RankParentSelection>(2.5, rng);
	auto selection = make_shared<TournamentParentSelection>(rng, 3, cmp);
	
	auto calcSingle = make_shared<SingleThreadedPopulationFitnessCalculation>(make_shared<TSPFitnessCalculation>(cities));
	
	auto mutation = make_shared<PermutationSwapMutation>(rng, 0.5/cities.size());
	//auto mutation = nullptr;
	
	// auto sortedPop = make_shared<SimpleSortedPopulation>(cmp);
	auto sortedPop = make_shared<TrackBestOnlySelectionPopulation>(cmp);

	auto cross = make_shared<SinglePopulationCrossover>(0.1, false,
			sortedPop,
			selection,
			make_shared<PermutationOrderCrossover>(rng, false),
			mutation,
			make_shared<SingleBestElitism>(true, mutation),
			make_shared<RemainingTargetPopulationSizeIteration>(),
			rng
		);
	
	MyStop stop(1000);
	MyGA ga;

	auto r = ga.run(*creation,
					*cross,
					*calcSingle, stop, 128);
	if (!r)
	{
		cerr << "Error: " << r.getErrorString() << endl;
		return -1;
	}
	cout << "Done" << endl;

	const VectorGenome<int> &cityOrder = ga.getBestGenome();
	for (int idx : cityOrder.getValues())
		cerr << cities[idx].x << " " << cities[idx].y << endl;

	int idx = cityOrder.getValues()[0];
	cerr << cities[idx].x << " " << cities[idx].y << endl;

	return 0;
}
