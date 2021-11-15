#include <iostream>
#include <vector>
#include <random>
#include <cassert>
#include <limits>
#include <algorithm>
#include <tuple>
#include <memory>

using namespace std;

struct Problem
{
	virtual double evaluate(const vector<double> &trial) = 0;
};

template<size_t D>
struct f1_Sphere : public Problem
{
	double evaluate(const vector<double> &x) override
	{
		assert(x.size() == D);
		double sumSquared = 0;
		for (auto v : x)
			sumSquared += v*v;
		return sumSquared;
	}
};

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
	double sum = 0;
	for (auto x : P)
		sum += x;
	sum /= P.size();
	return sum;
}

inline double meanL(const vector<double> &P)
{
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
	
	size_t evaluations = 0;
	double muF = 0.5, muCR = 0.5;
	vector<pair<vector<double>,double>> P(NP);
	vector<pair<double,double>> muLog;

	for (auto &ind : P)
	{
		ind.first = generateIndividual(rng, IPR);
		ind.second = problem.evaluate(ind.first);
		evaluations++;
	}
	sort(P.begin(), P.end(), [](auto &i, auto &j) { return i.second < j.second; });

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
					v[j] = xi[j] + Fi*(P[xBest].first[j] - xi[j]) + Fi*(P[xR1].first[j] - P[xR2].first[j]);
			}
			else // choose r2 from archive
			{
				assert(xR2 - NP < A.size());
				for(size_t j = 0 ; j < v.size() ; j++)
					v[j] = xi[j] + Fi*(P[xBest].first[j] - xi[j]) + Fi*(P[xR1].first[j] - A[xR2-NP][j]);
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
			}
			else
			{
				if (useArchive)
					A.push_back(xi);

				P[i].second = score;
				P[i].first = u;
				SF.push_back(Fi);
				SCR.push_back(CRi);
			}
		}

		// Randomly remove solutions from A to trim size
		while(A.size() > NP)
		{
			size_t idx = (uniform_int_distribution<>(0, A.size()-1))(rng);
			if (idx != A.size()-1) // move last one to selected position
				A[idx] = A[A.size()-1];
			A.resize(A.size()-1);
		}

		muCR = (1.0-c)*muCR + c*meanA(SCR);
		muF = (1.0-c)*muF + c* meanL(SF);

		sort(P.begin(), P.end(), [](auto &i, auto &j) { return i.second < j.second; });
	}

	return { false, P[0].first, P[0].second, maxGen, evaluations, muLog };
}

struct Test
{
	string name;
	shared_ptr<Problem> problem;
	vector<pair<double,double>> IPR;
	double VTR;
	size_t NP;
	size_t maxGen;
	bool enforceBounds;
};

int main(int argc, char const *argv[])
{
	random_device rd;
	mt19937 rng(rd());

	bool archive = true;
	if (getenv("NOARCH"))
		archive = false;

	vector<Test> tests {
		{ "f1_Sphere_30", make_shared<f1_Sphere<30>>(), vector<pair<double,double>>(30,{-100,100}), 1e-8, 100, 100000, false },
		{ "f1_Sphere_100", make_shared<f1_Sphere<100>>(), vector<pair<double,double>>(100,{-100,100}), 1e-8, 400, 100000, false },
	};
	
	size_t numRuns = 20;
	if (getenv("NUMRUNS"))
		numRuns = stoi(getenv("NUMRUNS"));

	bool printBest = false;
	if (getenv("DUMPFINAL"))
		printBest = true;

	for (auto &test : tests)
	{
		size_t successCount = 0;
		size_t totalEvals = 0;

		cout << test.name << ": ";
		for (size_t r = 0 ; r < numRuns ; r++)
		{
			bool converged;
			vector<double> best;
			double score;
			size_t generations, evaluations;
			vector<pair<double,double>> muLog;

			tie(converged, best, score, generations, evaluations, muLog) = jade(rng, *test.problem,
				test.NP, test.IPR, test.VTR, test.maxGen, archive, test.enforceBounds);

			if (converged)
			{
				successCount++;
				totalEvals += evaluations;
				if (printBest)
				{
					cout << "FOUND: [";
					for (auto x : best)
						cout << " " << x;
					cout << " ] => " << score << endl;
				}
			}
		}

		double avgEvals = (double)totalEvals / (double)successCount;
		cout << " NFE " << avgEvals;
		cout << " " << successCount << "/" << numRuns;
		if (archive)
			cout << " archive";
		else
			cout << " no archive";

		cout << endl;
	}
	return 0;
}
