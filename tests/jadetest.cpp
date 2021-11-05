#include "vectordifferentialevolution.h"
#include "evolutionaryalgorithm.h"
#include "mersennerandomnumbergenerator.h"
#include "valuefitness.h"
#include "jadeevolver.h"
#include "singlethreadedpopulationfitnesscalculation.h"
#include <random>

using namespace std;
using namespace eatk;
using namespace errut;

template<class VT, class FT>
class VectorValueFitnessCalculation : public GenomeFitnessCalculation
{
public:
	VectorValueFitnessCalculation() : m_evaluations(0) { };

	errut::bool_t calculate(const Genome &genome, Fitness &fitness)
	{
		const VectorGenome<VT> &vg = static_cast<const VectorGenome<VT> &>(genome);
		ValueFitness<FT> &vf = static_cast<ValueFitness<FT> &>(fitness);

		m_evaluations++;

		vf.setValue(calculate(vg.getValues()));
		return true;
	}

	virtual FT calculate(const vector<VT> &x) = 0;
protected:
	size_t m_evaluations;
};

class BaseCalculation : public VectorValueFitnessCalculation<double,double>
{
public:
	size_t getNumberOfEvaluations() const { return m_evaluations; }
	void resetEvaluationCount() { m_evaluations = 0; }
};

class f1_Sphere : public BaseCalculation
{
public:
	f1_Sphere(size_t D) : m_D(D) { }
	double calculate(const vector<double> &x) override
	{
		assert(x.size() == m_D);
		double sumSquared = 0;
		for (auto v : x)
			sumSquared += v*v;
		return sumSquared;
	}
private:
	const size_t m_D;
};

class f3 : public BaseCalculation
{
public:
	f3(size_t D) : m_D(D) { }
	double calculate(const vector<double> &x) override
	{
		assert(x.size() == m_D);
		double sum = 0;
		
		for (size_t i = 0 ; i < m_D ; i++)
		{
			double sub = 0;
			for (size_t j = 0 ; j <= i ; j++)
				sub += x[j];
			sum += sub*sub;
		}
		return sum;
	}
private:
	const size_t m_D;
};

template<class T>
class ValueToReachStop : public StopCriterion
{
public:
	ValueToReachStop(T value, size_t maxGenerations) : m_value(value), m_maxGen(maxGenerations),m_numGen(0)
	{
		m_dump = (getenv("DUMPPROGRESS"))?true:false;
	}

	bool_t analyze(const std::vector<std::shared_ptr<Individual>> &currentBest, size_t generationNumber, bool &shouldStop)
	{
		if (currentBest.size() != 1)
			return "Expecting current best size 1";
		
		if (m_dump)
			cerr << generationNumber << " " << currentBest[0]->toString() << endl;

		m_numGen = generationNumber;
		const ValueFitness<T> &vf = static_cast<const ValueFitness<T> &>(currentBest[0]->fitnessRef());
		if (vf.getValue() <= m_value)
		{
			shouldStop = true;
			m_best = currentBest[0]->createCopy();

			if (getenv("DUMPFINAL"))
				cerr << "FOUND: " << m_best->toString() << endl;
		}
		else
		{
			if (generationNumber >= m_maxGen)
			{
				shouldStop = true;
				m_numGen = generationNumber;
			}
		}
		return true;			
	}

	size_t getNumberOfGenerations() const { return m_numGen; }
	bool converged() const { return (m_best.get())?true:false; }
	const Individual &bestSolution() const { return *m_best; }
private:
	const T m_value;
	const size_t m_maxGen;
	size_t m_numGen;
	shared_ptr<Individual> m_best;
	bool m_dump;
};

struct Test
{
	string name;
	size_t popSize;
	bool archive;
	vector<double> bottom;
	vector<double> top;
	double VTR; // Value to reach
	size_t maxGenerations;
	shared_ptr<BaseCalculation> calculator;
};

class MyJADE : public JADEEvolver
{
public:
	MyJADE(const std::shared_ptr<RandomNumberGenerator> &rng,
		const std::shared_ptr<DifferentialEvolutionMutation> &mut,
		const std::shared_ptr<DifferentialEvolutionCrossover> &cross,
		const std::shared_ptr<FitnessComparison> &fitComp,
		bool archive)
		: JADEEvolver(rng, mut, cross, fitComp, 0, 0.05, 0.1, archive)
	{
		m_logFCR =  (getenv("LOGFCR"))?true:false;
	}
protected:
	void onMutationCrossoverSettings(double muF, double muCR) const override
	{
		if (m_logFCR)
			cerr << muF << " " << muCR << endl;
	}
private:
	bool m_logFCR;
};

class MyEA : public EvolutionaryAlgorithm
{
public:
	MyEA()
	{
		m_dumpPop = (getenv("DUMPPOP"))?true:false;
	}
private:
	bool_t onFitnessCalculated(size_t generation, std::shared_ptr<Population> &population) override
	{
		if (m_dumpPop)
		{
			cerr << "Generation " << generation << ":" << endl;
			for (size_t i = 0 ; i < population->individuals().size() ; i++)
				cerr << "  [" << i << "] = " << population->individual(i)->toString() << endl;
			cerr << endl;
		}
		return true;
	}

	bool m_dumpPop;
};

int main(int argc, char const *argv[])
{
	random_device rndDev;
	unsigned int seed = rndDev();
	if (getenv("SEED"))
		seed = (unsigned int)stoul(getenv("SEED"));

	size_t numRuns = 20;
	if (getenv("NUMRUNS"))
		numRuns = stoi(getenv("NUMRUNS"));

	cout << "# Seed = " << seed << endl;
	cout << "# Runs = " << numRuns << endl;

	shared_ptr<RandomNumberGenerator> rng = make_shared<MersenneRandomNumberGenerator>(seed);

	vector<Test> tests {
		// { "f1_Sphere_30_arch", 100, true, vector<double>(30, -100), vector<double>(30, 100), 1e-8, 100000, make_shared<f1_Sphere>(30) },
		// { "f1_Sphere_30_noarch", 100, false, vector<double>(30, -100), vector<double>(30, 100), 1e-8, 100000, make_shared<f1_Sphere>(30) },
		// { "f1_Sphere_100_arch", 400, true, vector<double>(100, -100), vector<double>(100, 100), 1e-8, 100000, make_shared<f1_Sphere>(100) },
		// { "f1_Sphere_100_noarch", 400, false, vector<double>(100, -100), vector<double>(100, 100), 1e-8, 100000, make_shared<f1_Sphere>(100) },
		// { "f3_30_arch", 100, true, vector<double>(30, -100), vector<double>(30, 100), 1e-8, 100000, make_shared<f3>(30) },
		{ "f3_30_noarch", 100, false, vector<double>(30, -100), vector<double>(30, 100), 1e-8, 100000, make_shared<f3>(30) },
	};

	auto comp = make_shared<ValueFitnessComparison<double>>();

	for (const auto &test : tests)
	{
		cout << test.name << ": "; 

		size_t successCount = 0;
		size_t totalCalc = 0;

		for (size_t run = 0 ; run < numRuns ; run++)
		{
			auto mut = make_shared<VectorDifferentialEvolutionMutation<double>>();
			auto cross = make_shared<VectorDifferentialEvolutionCrossover<double>>(rng);
			
			VectorDifferentialEvolutionIndividualCreation<double,double> creation(test.bottom, test.top, rng);
			ValueToReachStop<double> stop(test.VTR, test.maxGenerations);
		
			MyJADE evolver(rng, mut, cross, comp, test.archive);
			SingleThreadedPopulationFitnessCalculation popCalc(test.calculator);

			MyEA ea;
			bool_t r = ea.run(creation, evolver, popCalc, stop, test.popSize, test.popSize, test.popSize*2);
			if (!r)
			{
				cerr << r.getErrorString() << endl;
				return -1;
			}

			if (!stop.converged())
			{
				// cout << "  Failed to converge after " << stop.getNumberOfGenerations() << " generations" << endl;
			}
			else
			{
				//cout << "  Found after " << stop.getNumberOfGenerations() << " generations, " << test.calculator->getNumberOfEvaluations()
				// 	 << " evaluations: " << stop.bestSolution().toString() << endl;

				successCount++;
				totalCalc += test.calculator->getNumberOfEvaluations();
			}

			test.calculator->resetEvaluationCount();
		}

		cout << successCount << "/" << numRuns;
		cout << ", nfe " << totalCalc/numRuns << ", NP = " << test.popSize << ((test.archive)?", archive":", noarchive");
		if (successCount != numRuns)
			cout << " !!";
		cout << endl;
	}

	return 0;
}
