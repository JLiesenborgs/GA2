#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <vector>
#include <random>
#include <cassert>
#include <limits>
#include <algorithm>
#include <tuple>
#include <memory>
#include <string>
#include "testfunctions.h"

using namespace std;

struct Problem
{
	virtual double evaluate(const vector<double> &trial) = 0;
	virtual vector<pair<double,double>> getIPR() = 0;
};

template <class T>
struct ProblemTemplate : public Problem
{
	template<typename... Args>
	ProblemTemplate(Args&&... args) : m_f(std::forward<Args>(args)...) { }
	double evaluate(const vector<double> &trial) override { return m_f.calculate1(trial); }

	vector<pair<double,double>> getIPR() override
	{
		auto [ lower, upper ] = m_f.getInitialParameterRange();
		if (lower.size() != upper.size())
			throw runtime_error("Internal error: sizes of lower and upper IPR vectors differ");
		
		vector<pair<double,double>> ipr;
		for (size_t i = 0 ; i < lower.size() ; i++)
			ipr.push_back({lower[i],upper[i]});

		return ipr;
	}

	T m_f;
};

typedef ProblemTemplate<eatk::testfunctions::Sphere> f1_Sphere;
typedef ProblemTemplate<eatk::testfunctions::Schwefel_2_22> f2;
typedef ProblemTemplate<eatk::testfunctions::Schwefel_1_2> f3;
typedef ProblemTemplate<eatk::testfunctions::Schwefel_2_21> f4;
typedef ProblemTemplate<eatk::testfunctions::GeneralizedRosenbrock> f5;
typedef ProblemTemplate<eatk::testfunctions::Step_2_Function> f6;
typedef ProblemTemplate<eatk::testfunctions::QuarticWithNoise> f7;
typedef ProblemTemplate<eatk::testfunctions::ModifiedSchwefel_2_26> f8;
typedef ProblemTemplate<eatk::testfunctions::Rastrigin> f9;
typedef ProblemTemplate<eatk::testfunctions::AckleyFunction1> f10;
typedef ProblemTemplate<eatk::testfunctions::Griewank> f11;
typedef ProblemTemplate<eatk::testfunctions::GeneralizedPenalizedFunction1> f12;
typedef ProblemTemplate<eatk::testfunctions::GeneralizedPenalizedFunction2> f13;
typedef ProblemTemplate<eatk::testfunctions::Branin> f14_Branin;
typedef ProblemTemplate<eatk::testfunctions::GoldsteinPrice> f15_GoldsteinPrice;
typedef ProblemTemplate<eatk::testfunctions::Hartman3D> f16_Hartman3D;
typedef ProblemTemplate<eatk::testfunctions::Hartman6D> f17_Hartman6D;
typedef ProblemTemplate<eatk::testfunctions::Shekel> f_Shekel;

inline vector<double> generateIndividual(mt19937 &rng, const vector<pair<double,double>> &IPR)
{
	vector<double> x(IPR.size());
	for (size_t i = 0 ; i < x.size(); i++)
	{
		assert(IPR[i].first < IPR[i].second);
		x[i] = (uniform_real_distribution<>(IPR[i].first, IPR[i].second))(rng);
	}
	
	return x;
}

inline double clip(double x, double x0, double x1)
{
	if (x < x0)
		return x0;
	if (x > x1)
		return x1;
	return x;
}

inline double meanA(const vector<double> &P)
{
	if (P.size() == 0)
		return 0;

	double sum = 0;
	for (auto x : P)
		sum += x;
	sum /= P.size();
	return sum;
}

inline double meanL(const vector<double> &P)
{
	if (P.size() == 0)
		return 0;

	double sum = 0, sum2 = 0;
	for (auto x : P)
	{
		sum += x;
		sum2 += x*x;
	}
	return sum2/sum;
}

tuple<bool, vector<double>, double, size_t, size_t, vector<pair<double,double>>>
jade(mt19937 &rng, Problem &problem, size_t NP, const vector<pair<double,double>> &IPR,
         double VTR, size_t maxGen, bool useArchive, bool enforceBounds = false,
		 double p = 0.05, double c = 0.1)
{
	assert(IPR.size() > 0);
	assert(NP > 4);
	assert(p >= 0 && p <= 1);
	
	bool dumpPop = false;
	if (getenv("DUMPPOP"))
		dumpPop = true;

	size_t evaluations = 0;
	double muF = 0.5, muCR = 0.5;
	vector<pair<vector<double>,double>> P(NP), Q(NP);
	vector<pair<double,double>> muLog;

	for (auto &ind : P)
	{
		ind.first = generateIndividual(rng, IPR);
		ind.second = problem.evaluate(ind.first);
		evaluations++;
	}
	sort(P.begin(), P.end(), [](auto &i, auto &j) { return i.second < j.second; });

	auto logPop = [&P, dumpPop](size_t gen)
	{
		if (!dumpPop)
			return;

		cerr << "Generation " << gen << endl;
		for (size_t i = 0 ; i < P.size() ; i++)
		{
			cerr << "(" << i << ") [";
			for (auto x : P[i].first)
				cerr << " " << x;
			cerr << " ] => " << P[i].second << endl;
		}
	};

	logPop(0);

	vector<vector<double>> A;
	vector<double> SF, SCR;

	for (size_t g = 0 ; g < maxGen ; g++)
	{
		SCR.clear();
		SF.clear();

		muLog.push_back({ muF, muCR});
		for (size_t i = 0 ; i < NP ; i++)
		{
			normal_distribution<> randn(muCR, 0.1);
			cauchy_distribution<> randc(muF, 0.1);
			double CRi = clip(randn(rng), 0, 1);
			double Fi = clip(randc(rng), 0, 1);

			size_t xBest = (size_t)((uniform_real_distribution<>(0, p))(rng) * NP);
			size_t xR1 = 0, xR2 = 0;
			do { xR1 = (uniform_int_distribution<>(0, NP-1))(rng); } while (xR1 == i);
			do { xR2 = (uniform_int_distribution<>(0, NP+A.size()-1))(rng); } while (xR2 == i || xR2 == xR1);

			vector<double> v(IPR.size());
			auto &xi = P[i].first;
			if (xR2 < NP) // choose r2 from main population
			{
				for(size_t j = 0 ; j < v.size() ; j++)
				{
					v[j] = xi[j] + Fi*(P[xBest].first[j] - xi[j]) + Fi*(P[xR1].first[j] - P[xR2].first[j]);
					assert(!isnan(v[j]));
				}
			}
			else // choose r2 from archive
			{
				assert(xR2 - NP < A.size());
				for(size_t j = 0 ; j < v.size() ; j++)
				{
					v[j] = xi[j] + Fi*(P[xBest].first[j] - xi[j]) + Fi*(P[xR1].first[j] - A[xR2-NP][j]);
					assert(!isnan(v[j]));
				}
			}

			vector<double> u(v.size());
			size_t jrand = (uniform_int_distribution<>(0, u.size()-1))(rng);
			for (size_t j = 0 ; j < u.size() ; j++)
			{
				if (j == jrand || (uniform_real_distribution<>(0,1))(rng) < CRi)
					u[j] = v[j];
				else
					u[j] = xi[j];
			}

			if (enforceBounds)
			{
				for (size_t j = 0 ; j < u.size() ; j++)
				{
					double xlow = IPR[j].first;
					double xhigh = IPR[j].second;
					if (u[j] < xlow)
						u[j] = (xlow + xi[j])/2.0;
					if (u[j] > xhigh)
						u[j] = (xhigh + xi[j])/2.0;
				}
			}

			double score = problem.evaluate(u);
			evaluations++;
			if (score < VTR)
			{
				return { true, u, score, g, evaluations, muLog };
			}

			if (P[i].second <= score)
			{
				// Nothing to do, just keep individual
				Q[i] = P[i];
			}
			else
			{
				if (useArchive)
					A.push_back(xi);

				Q[i].second = score;
				Q[i].first = u;
				SF.push_back(Fi);
				SCR.push_back(CRi);
			}
		}
		swap(P, Q);

		// Randomly remove solutions from A to trim size
		while(A.size() > NP)
		{
			size_t idx = (uniform_int_distribution<>(0, A.size()-1))(rng);
			if (idx != A.size()-1) // move last one to selected position
				A[idx] = A[A.size()-1];
			A.resize(A.size()-1);
		}

		// What to do if these sets are empty? Use 0 as mean? Don't update?
		muCR = (1.0-c)*muCR + c*meanA(SCR);
		muF = (1.0-c)*muF + c*meanL(SF);

		sort(P.begin(), P.end(), [](auto &i, auto &j) { return i.second < j.second; });
		logPop(g+1);
	}

	return { false, P[0].first, P[0].second, maxGen, evaluations, muLog };
}

struct Test
{
	string name;
	shared_ptr<Problem> problem;
	double VTR;
	size_t NP;
	size_t maxGen;
	bool enforceBounds;
};

class AvgStd
{
public:
	AvgStd() : m_sum(0), m_sum2(0), m_count(0) { }
	void process(double x)
	{
		m_sum += x;
		m_sum2 += x*x;
		m_count++;
	}
	string toString() const
	{
		double avg = m_sum/m_count;
		double std = std::sqrt(m_sum2/m_count - avg*avg);
		return to_string(avg) + " (" + to_string(std) + ")";
	}
	size_t getCount() const { return m_count; }
private:
	double m_sum;
	double m_sum2;
	size_t m_count;
};

int main(int argc, char const *argv[])
{
	random_device rd;
	unsigned int seed = rd();
	if (getenv("SEED"))
		seed = (unsigned int)stoul(getenv("SEED"));

	auto rng = make_shared<mt19937>(seed);

	string matchName;
	if (argc > 1)
		matchName = argv[1];

	bool archive = true;
	if (getenv("NOARCH"))
		archive = false;

	vector<Test> tests {
		{ "f1_Sphere_30", make_shared<f1_Sphere>(30, pair(-100.0, 100.0)), 1e-8, 100, 100000, false },
		{ "f1_Sphere_100", make_shared<f1_Sphere>(100, pair(-100.0,100.0)), 1e-8, 400, 100000, false },
		{ "f2_30", make_shared<f2>(30, pair(-10.0,10.0)), 1e-8, 100, 100000, false },
		{ "f2_100", make_shared<f2>(100, pair(-10.0,10.0)), 1e-8, 400, 100000, false },
		{ "f3_30", make_shared<f3>(30, pair(-100.0,100.0)), 1e-8, 100, 100000, false },
		{ "f3_100", make_shared<f3>(100, pair(-100.0,100.0)), 1e-8, 400, 100000, false },
		{ "f4_30", make_shared<f4>(30, pair(-100.0,100.0)), 1e-8, 100, 100000, false },
		{ "f4_100", make_shared<f4>(100, pair(-100.0,100.0)), 1e-8, 400, 100000, false },
		{ "f5_30", make_shared<f5>(30, pair(-100.0,100.0)), 1e-8, 100, 100000, false },
		{ "f5_100", make_shared<f5>(100, pair(-100.0,100.0)), 1e-8, 400, 100000, false },
		{ "f6_30", make_shared<f6>(30, pair(-100.0,100.0)), 1e-8, 100, 100000, false },
		{ "f6_100", make_shared<f6>(100, pair(-100.0,100.0)), 1e-8, 400, 100000, false },
		//{ "f7_30", make_shared<f7>(30, rng, pair(-1.28,1.28)), 1e-2, 100, 100000, false },
		//{ "f7_100", make_shared<f7>(100, rng, pair(-1.28,1.28)), 1e-2, 400, 100000, false },
		{ "f8_30", make_shared<f8>(30, pair(-500.0,500.0)), 1e-8, 100, 100000, true },
		{ "f8_100", make_shared<f8>(100, pair(-500.0,500.0)), 1e-8, 400, 100000, true },
		{ "f9_30", make_shared<f9>(30, pair(-5.12,5.12)), 1e-8, 100, 100000, false },
		{ "f9_100", make_shared<f9>(100, pair(-5.12,5.12)), 1e-8, 400, 100000, false },
		{ "f10_30", make_shared<f10>(30, pair(-32.0,32.0)), 1e-8, 100, 100000, false },
		{ "f10_100", make_shared<f10>(100, pair(-32.0,32.0)), 1e-8, 400, 100000, false },
		{ "f11_30", make_shared<f11>(30, pair(-600.0,600.0)), 1e-8, 100, 100000, false },
		{ "f11_100", make_shared<f11>(100, pair(-600.0,600.0)), 1e-8, 400, 100000, false },
		{ "f12_30", make_shared<f12>(30, pair(-50.0,50.0)), 1e-8, 100, 100000, false },
		{ "f12_100", make_shared<f12>(100, pair(-50.0,50.0)), 1e-8, 400, 100000, false },
		{ "f13_30", make_shared<f13>(30, pair(-50.0,50.0)), 1e-8, 100, 100000, false },
		{ "f13_100", make_shared<f13>(100, pair(-50.0,50.0)), 1e-8, 400, 100000, false },
		{ "f14_Branin", make_shared<f14_Branin>(vector{-5.0,0.0}, vector{10.0,15.0}), 0.3978873577 + 1e-8, 30, 100000, false },
		{ "f15_GoldsteinPrice", make_shared<f15_GoldsteinPrice>(pair(-2.0,2.0)), 3.0 + 1e-8, 30, 100000, false },
		{ "f16_Hartman3D", make_shared<f16_Hartman3D>(pair(0.0,1.0)), -3.862779787 + 1e-8, 30, 100000, false },
		{ "f17_Hartman6D", make_shared<f17_Hartman6D>(pair(0.0,1.0)), -3.322368011 + 1e-8, 30, 100000, false },
		{ "f18_Shekel5", make_shared<f_Shekel>(5, pair(0.0,10.0)), -10.153199679 + 1e-8, 30, 100000, false },
		{ "f19_Shekel7", make_shared<f_Shekel>(7, pair(0.0,10.0)), -10.402915336 + 1e-8, 30, 100000, false },
		{ "f20_Shekel10", make_shared<f_Shekel>(10, pair(0.0,10.0)), -10.536443153 + 1e-8, 30, 100000, false },
	};
	
	size_t numRuns = 20;
	if (getenv("NUMRUNS"))
		numRuns = stoi(getenv("NUMRUNS"));

	bool printBest = false;
	if (getenv("DUMPFINAL"))
		printBest = true;

	size_t testCount = 0;
	for (auto &test : tests)
	{
		if (matchName.length() > 0 && matchName != test.name)
			continue;

		testCount++;

		AvgStd nfe;
		AvgStd genAvg;

		cout << test.name << ": ";
		for (size_t r = 0 ; r < numRuns ; r++)
		{
			bool converged;
			vector<double> best;
			double score;
			size_t generations, evaluations;
			vector<pair<double,double>> muLog;

			auto IPR = test.problem->getIPR();

			tie(converged, best, score, generations, evaluations, muLog) = jade(*rng, *test.problem,
				test.NP, IPR, test.VTR, test.maxGen, archive, test.enforceBounds);

			if (converged)
			{
				nfe.process(evaluations);
				genAvg.process(generations);
				if (printBest)
				{
					cout << "FOUND: [";
					for (auto x : best)
						cout << " " << x;
					cout << " ] => " << score << endl;
				}
			}

			if (getenv("LOGFCR"))
				for (auto mu : muLog)
					cerr << mu.first << " " << mu.second << endl;
		}

		cout << " gen " << genAvg.toString();
		cout << " NFE " << nfe.toString();
		cout << " " << nfe.getCount() << "/" << numRuns;
		if (archive)
			cout << " archive";
		else
			cout << " no archive";

		cout << endl;
	}

	if (testCount == 0)
	{
		cerr << "No tests run, known are:" << endl;
		for (const auto &test : tests)
			cerr << "  " << test.name << endl;
	}
	return 0;
}
