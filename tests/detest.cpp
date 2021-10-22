#include "vectordifferentialevolution.h"
#include "evolutionaryalgorithm.h"
#include "mersennerandomnumbergenerator.h"
#include "valuefitness.h"
#include "differentialevolutionevolver.h"
#include "singlethreadedpopulationfitnesscalculation.h"
#include <random>

using namespace std;
using namespace eatk;
using namespace errut;

template<class VT, class FT>
class VectorValueFitnessCalculation : public GenomeFitnessCalculation
{
public:
	VectorValueFitnessCalculation() { };

	errut::bool_t calculate(const Genome &genome, Fitness &fitness)
	{
		const VectorGenome<VT> &vg = static_cast<const VectorGenome<VT> &>(genome);
		ValueFitness<FT> &vf = static_cast<ValueFitness<FT> &>(fitness);

		vf.setValue(calculate(vg.getValues()));
		return true;
	}

	virtual FT calculate(const vector<VT> &x) = 0;
};

template<class VT, class FT>
class Rosenbrock : public VectorValueFitnessCalculation<VT,FT>
{
public:
	FT calculate(const vector<VT> &x) override
	{
		assert(x.size() == 2);
		VT x1 = x[0];
		VT x2 = x[1];

		return (FT)((VT)100.0 * (x1*x1 - x2) * (x1*x1 - x2) + ((VT)1 - x1)*((VT)1 - x1));
	}
};

template<class VT, class FT>
class Zimmermann : public VectorValueFitnessCalculation<VT,FT>
{
public:
	FT calculate(const vector<VT> &x) override
	{
		assert(x.size() == 2);
		VT x1 = x[0];
		VT x2 = x[1];

		auto h1 = (VT)9-x1-x2;
		auto h2 = (x1-(VT)3)*(x1-(VT)3) + (x2-(VT)2)*(x2-(VT)2) - (VT)16;
		auto h3 = x1*x2 - (VT)14;

		auto p = [](auto delta) { return (VT)100 * ((VT)1 + delta); };
		auto sgn = [] (auto x) { return (x >= 0)?(VT)1:(VT)0; };

		return std::max(h1,
		
			std::max(
				std::max(p(h2)*sgn(h2),p(h3)*sgn(h3)), 
				std::max(p(-x1)*sgn(-x1),p(-x2)*sgn(-x2))
			)
		);
	}
};

template<class T>
class ValueToReachStop : public StopCriterion
{
public:
	ValueToReachStop(T value) : m_value(value) { }
	bool_t analyze(const std::vector<std::shared_ptr<Individual>> &currentBest, size_t generationNumber, bool &shouldStop)
	{
		if (currentBest.size() != 1)
			return "Expecting current best size 1";
		
		cerr << generationNumber << " " << currentBest[0]->toString() << endl;
		const ValueFitness<T> &vf = static_cast<const ValueFitness<T> &>(currentBest[0]->fitnessRef());
		if (vf.getValue() <= m_value)
			shouldStop = true;
		return true;			
	}
private:
	T m_value;
};

class MyEA : public EvolutionaryAlgorithm
{
private:
	bool_t onAlgorithmDone(size_t generation, const std::vector<std::shared_ptr<Individual>> &bestIndividuals)
	{
		cout << "EA done after " << generation << " generations, best:" << endl;
		for (auto &ind : bestIndividuals)
			cout << ind->toString() << endl;
		return true;
	}
};

int main(int argc, char const *argv[])
{
	random_device rndDev;
	shared_ptr<RandomNumberGenerator> rng = make_shared<MersenneRandomNumberGenerator>(rndDev());
	VectorDifferentialEvolutionIndividualCreation<double,double> creation(
		// { -2.048, -2.048}, { 2.048, 2.048 },
		{ 0.0, 0.0 }, { 100.0, 100.0 },
		rng);

	double F = 0.9;
	float CR = 0.1;
	auto mut = make_shared<VectorDifferentialEvolutionMutation<double>>(F);
	auto cross = make_shared<VectorDifferentialEvolutionCrossover<double>>(CR, rng);
	auto comp = make_shared<ValueFitnessComparison<double>>();
	DifferentialEvolutionEvolver evolver(rng, mut, cross, comp);
	
	size_t popSize = 10;
	ValueToReachStop<double> stop(1e-6);
	//auto calc = make_shared<Rosenbrock<double,double>>();
	auto calc = make_shared<Zimmermann<double,double>>();
	
	SingleThreadedPopulationFitnessCalculation popCalc(calc);

	MyEA ea;

	cout << "Running EA" << endl;
	bool_t r = ea.run(creation, evolver, popCalc, stop, popSize, popSize, popSize*2);
	if (!r)
	{
		cerr << r.getErrorString() << endl;
		return -1;
	}
	cout << "Finishing..." << endl;
	return 0;
}
