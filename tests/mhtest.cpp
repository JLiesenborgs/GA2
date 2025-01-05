#include "metropolishastingsevolver.h"
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

class MyMH : public MetropolisHastingsEvolver
{
public:
	MyMH(const std::shared_ptr<RandomNumberGenerator> &rng,
		 const std::vector<double> &stepScales,
		 ProbType probType,
		 size_t burnin)
		: MetropolisHastingsEvolver(rng, stepScales, probType), m_burnin(burnin) {} 

	void onSamples(const std::vector<std::shared_ptr<Individual>> &samples)
	{
		m_gen++;
		if (m_gen < m_burnin)
			return;

		for (auto &i : samples)
			cerr << i->genome()->toString() << endl;
	}

	size_t m_burnin = 0;
	size_t m_gen = 0;
};

template<class T>
class ProbCalc : public GenomeFitnessCalculation
{
public:
	ProbCalc(MetropolisHastingsEvolver::ProbType probType) : m_probType(probType) { }
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
		if (m_probType == MetropolisHastingsEvolver::Log)
		{
			printOnce("MetropolisHastingsEvolver::Log");
			prob = 0;
			for (size_t i = 0 ; i < x.size() ; i++)
				prob += -(x[i]-mu[i])*(x[i]-mu[i])/((T)2.0*sigma[i]*sigma[i]);
		}
		else if (m_probType == MetropolisHastingsEvolver::NegativeLog)
		{
			printOnce("MetropolisHastingsEvolver::NegativeLog");
			prob = 0;
			for (size_t i = 0 ; i < x.size() ; i++)
				prob += (x[i]-mu[i])*(x[i]-mu[i])/((T)2.0*sigma[i]*sigma[i]);
		}
		else if (m_probType == MetropolisHastingsEvolver::Regular)
		{
			printOnce("MetropolisHastingsEvolver::Regular");
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
	MetropolisHastingsEvolver::ProbType m_probType;
};

template<class T>
int templateMain(MetropolisHastingsEvolver::ProbType probType)
{
	random_device rndDev;
	unsigned int seed = rndDev();
	if (getenv("SEED"))
		seed = (unsigned int)stoul(getenv("SEED"));

	shared_ptr<RandomNumberGenerator> rng = make_shared<MersenneRandomNumberGenerator>(seed);
	VectorDifferentialEvolutionIndividualCreation<T,T> creation({-10, -10}, {10, 10} , rng);

	EvolutionaryAlgorithm ea;
	size_t totalGen = 10000;
	size_t burnin = totalGen/2;
	MyMH mh(rng,
			{0.1, 0.1},
			probType,
			burnin);

	//mh.setAnnealingExponent(0.5);

	auto prob = make_shared<ProbCalc<T>>(probType);
	SingleThreadedPopulationFitnessCalculation popCalc(prob);
	FixedGenerationsStopCriterion stop(totalGen);
	size_t N = 512;

	bool_t r;
	if (!(r = ea.run(creation, mh, popCalc, stop, N, N, N*2)))
		throw runtime_error(r.getErrorString());

	cout << "Best is:" << endl;
	cout << mh.getBestIndividuals()[0]->genome()->toString() << endl;

	return 0;
}

int main(void)
{
	return templateMain<double>(MetropolisHastingsEvolver::NegativeLog);
}
