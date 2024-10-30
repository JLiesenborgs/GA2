#include "goodmanweareevolver.h"
#include "evolutionaryalgorithm.h"
#include "mersennerandomnumbergenerator.h"
#include "valuefitness.h"
#include "vectordifferentialevolution.h"
#include "stopcriterion.h"
#include "multithreadedpopulationfitnesscalculation.h"
#include "singlethreadedpopulationfitnesscalculation.h"
#include <random>

using namespace std;
using namespace eatk;
using namespace errut;

inline void printOnce(const string &s)
{
	static bool shown = false;
	if (shown)
		return;
	shown = true;
	cout << s << endl;
}

class MyGW : public GoodmanWeareEvolver
{
public:
	MyGW(const std::shared_ptr<RandomNumberGenerator> &rng, ProbType probType, double a = 2.0)
		: GoodmanWeareEvolver(rng, probType, a) {} 

	void onSamples(const std::vector<std::shared_ptr<Individual>> &samples)
	{
		for (auto &i : samples)
			cerr << i->genome()->toString() << endl;
	}
};

template<class T>
class ProbCalc : public GenomeFitnessCalculation
{
public:
	ProbCalc(GoodmanWeareEvolver::ProbType probType) : m_probType(probType) { }
	bool_t calculate(const Genome &genome0, Fitness &fitness0)
	{
		assert(dynamic_cast<ValueFitness<T>*>(&fitness0));
		assert(dynamic_cast<const VectorGenome<T>*>(&genome0));
		const VectorGenome<T> &g = static_cast<const VectorGenome<T> &>(genome0);
		ValueFitness<T> &f = static_cast<ValueFitness<T> &>(fitness0);

		const vector<T> &x = g.getValues();
		vector<T> mu = { 3, 5 };
		vector<T> sigma = { 2, 1 };

		assert(x.size() == mu.size() && x.size() == sigma.size());
		T prob;
		if (m_probType == GoodmanWeareEvolver::Log)
		{
			printOnce("GoodmanWeareEvolver::Log");
			prob = 0;
			for (size_t i = 0 ; i < x.size() ; i++)
				prob += -(x[i]-mu[i])*(x[i]-mu[i])/((T)2.0*sigma[i]*sigma[i]);
		}
		else if (m_probType == GoodmanWeareEvolver::NegativeLog)
		{
			printOnce("GoodmanWeareEvolver::NegativeLog");
			prob = 0;
			for (size_t i = 0 ; i < x.size() ; i++)
				prob += (x[i]-mu[i])*(x[i]-mu[i])/((T)2.0*sigma[i]*sigma[i]);
		}
		else if (m_probType == GoodmanWeareEvolver::Regular)
		{
			printOnce("GoodmanWeareEvolver::Regular");
			prob = 1;
			for (size_t i = 0 ; i < x.size() ; i++)
				prob *= std::exp(-(x[i]-mu[i])*(x[i]-mu[i])/((T)2.0*sigma[i]*sigma[i]));
		}
		else
			throw runtime_error("Invalid prob type");
		f.setValue(prob);
		return true;
	}
private:
	GoodmanWeareEvolver::ProbType m_probType;
};

template<class T>
int templateMain(GoodmanWeareEvolver::ProbType probType)
{
	random_device rndDev;
	unsigned int seed = rndDev();
	if (getenv("SEED"))
		seed = (unsigned int)stoul(getenv("SEED"));

	shared_ptr<RandomNumberGenerator> rng = make_shared<MersenneRandomNumberGenerator>(seed);
	VectorDifferentialEvolutionIndividualCreation<T,T> creation({-10, -10}, {10, 10} , rng);

	EvolutionaryAlgorithm ea;
	MyGW gw(rng, probType);
	auto prob = make_shared<ProbCalc<T>>(probType);
	SingleThreadedPopulationFitnessCalculation popCalc(prob);
	FixedGenerationsStopCriterion stop(1000);
	size_t N = 512;

	bool_t r;
	if (!(r = ea.run(creation, gw, popCalc, stop, N, N, N*3/2)))
		throw runtime_error(r.getErrorString());

	cout << "Best is:" << endl;
	cout << gw.getBestIndividuals()[0]->genome()->toString() << endl;

	return 0;
}

int main(void)
{
	return templateMain<double>(GoodmanWeareEvolver::NegativeLog);
}
