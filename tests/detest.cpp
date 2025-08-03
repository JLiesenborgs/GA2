#define _USE_MATH_DEFINES
#include "vectordifferentialevolution.h"
#include "evolutionaryalgorithm.h"
#include "mersennerandomnumbergenerator.h"
#include "valuefitness.h"
#include "differentialevolutionevolver.h"
#include "singlethreadedpopulationfitnesscalculation.h"
#include "testfunctions.h"
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

template <class T>
class TemplateCalculation : public BaseCalculation
{
public:
	template<typename... Args>
	TemplateCalculation(Args&&... args) : m_f(std::forward<Args>(args)...) { }

	double calculate(const vector<double> &x) override { return m_f.calculate1(x); }
	T m_f;
};

typedef TemplateCalculation<testfunctions::Sphere> f1_Sphere;
typedef TemplateCalculation<testfunctions::Rosenbrock> f2_Rosenbrock;
typedef TemplateCalculation<testfunctions::Mod3rdDeJong> f3_Mod3rdDeJong;
typedef TemplateCalculation<testfunctions::QuarticWithNoise> f4_Mod4thDeJong_quartic;
typedef TemplateCalculation<testfunctions::Foxholes> f5_Foxholes;
typedef TemplateCalculation<testfunctions::Corana> f6_Corana;
typedef TemplateCalculation<testfunctions::Griewank> f7_Griewangk;
typedef TemplateCalculation<testfunctions::Zimmermann> f8_Zimmermann;
typedef TemplateCalculation<testfunctions::k4_Poly> f9_k4_Poly;
typedef TemplateCalculation<testfunctions::k8_Poly> f9_k8_Poly;
typedef TemplateCalculation<testfunctions::HyperEllipsoid> f11_HyperEllipsoid;
typedef TemplateCalculation<testfunctions::Katsuura> f12_Katsuura;

template<class T>
class ValueToReachStop : public StopCriterion
{
public:
	ValueToReachStop(T value, size_t maxGenerations) : m_value(value), m_maxGen(maxGenerations),m_numGen(0)
	{
		m_dump = (getenv("DUMPPROGRESS"))?true:false;
	}

	bool_t analyze(const PopulationEvolver &evolver, size_t generationNumber, bool &shouldStop) override
	{
		auto &currentBest = evolver.getBestIndividuals();
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
	double F;
	double CR;
	vector<double> bottom;
	vector<double> top;
	double VTR; // Value to reach
	size_t maxGenerations;
	shared_ptr<BaseCalculation> calculator;
};

class MyEA : public EvolutionaryAlgorithm
{
public:
	MyEA()
	{
		m_dumpPop = (getenv("DUMPPOP"))?true:false;
	}
private:
	bool_t onFitnessCalculated(size_t generation, const std::shared_ptr<Population> &population) override
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
		{ "f1_Sphere", 10, 0.5, 0.3, { -5.12, -5.12, -5.12 }, { 5.12, 5.12, 5.12 }, 1e-6, 100000, make_shared<f1_Sphere>(3, pair(-5.12, 5.12)) },
		{ "f2_Rosenbrock", 10, 0.9, 0.9, { -2.048, -2.048 }, { 2.048, 2.048 }, 1e-6, 100000, make_shared<f2_Rosenbrock>(pair(-2.048,2.048)) },
		{ "f3_Mod3rdDeJong_step", 10, 0.9, 0, { -5.12, -5.12, -5.12, -5.12, -5.12 }, { 5.12, 5.12, 5.12, 5.12, 5.12 }, 1e-6, 100000, make_shared<f3_Mod3rdDeJong>(pair(-5.12,5.12)) },
		//{ "f4_Mod4thDeJong_quartic", 10, 0.9, 0, vector<double>(30, -1.28), vector<double>(30, 1.28), 0, 100000, make_shared<f4_Mod4thDeJong_quartic>(30, rng, true, pair(-1.28,1.28)) },
		{ "f5_Foxholes", 15, 0.9, 0, { -65.536, -65.536}, { 65.536, 65.536 }, 0.998005, 100000, make_shared<f5_Foxholes>(pair(-65.536,65.536)) },
		{ "f6_Corana", 10, 0.5, 0, { -1000, -1000, -1000, -1000}, { 1000, 1000, 1000, 1000}, 1e-6, 100000, make_shared<f6_Corana>(pair(-1000.0, 1000.0)) },
		{ "f7_Griewangk", 25, 0.5, 0.2, vector<double>(10, -400), vector<double>(10, 400), 1e-6, 100000, make_shared<f7_Griewangk>(10, pair(-400.0,400.0)) },
		{ "f8_Zimmermann", 10, 0.9, 0.9, { 0.0, 0.0 }, { 100.0, 100.0 }, 1e-6, 100000, make_shared<f8_Zimmermann>(pair(0.0,100.0)) },
		{ "f9_k4_Poly", 60, 0.6, 1, vector<double>(9,-100), vector<double>(9,100), 1e-6, 100000, make_shared<f9_k4_Poly>(pair(-100.0,100.0)) },
		{ "f9_k8_Poly", 100, 0.6, 1, vector<double>(17,-1000), vector<double>(17,1000), 1e-6, 100000, make_shared<f9_k8_Poly>(pair(-1000.0,1000.0)) },
		{ "f11_HyperEllipsoid_30", 20, 0.5, 0.1, vector<double>(30, -1), vector<double>(30, 1), 1e-10, 100000, make_shared<f11_HyperEllipsoid>(30, pair(-1.0,1.0)) },
		{ "f11_HyperEllipsoid_100", 20, 0.5, 0.1, vector<double>(100, -1), vector<double>(100, 1), 1e-10, 100000, make_shared<f11_HyperEllipsoid>(100, pair(-1.0,1.0)) },
		{ "f12_Katsuura_10", 15, 0.5, 0.1, vector<double>(10, -1000), vector<double>(10, 1000), 1.05, 100000, make_shared<f12_Katsuura>(10, pair(-1000.0,1000.0)) },
		{ "f12_Katsuura_30", 15, 0.5, 0.1, vector<double>(30, -1000), vector<double>(30, 1000), 1.05, 100000, make_shared<f12_Katsuura>(30, pair(-1000.0,1000.0)) },
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
		
			DifferentialEvolutionEvolver evolver(rng, mut, test.F, cross, test.CR, comp);
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
		cout << ", nfe " << totalCalc/numRuns << ", F = " << test.F << ", CR = " << test.CR << ", NP = " << test.popSize;
		if (successCount != numRuns)
			cout << " !!";
		cout << endl;
	}

	return 0;
}
