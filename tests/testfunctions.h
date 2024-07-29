#include "eatkconfig.h"
#include "randomnumbergenerator.h"
#include <vector>
#include <stdexcept>
#include <cassert>
#include <cmath>
#include <limits>
#include <memory>

namespace eatk
{

namespace testfunctions
{

class TestFunction
{
public:
	typedef std::pair<double,double> double2_t;
	static inline constexpr double infinity() { return std::numeric_limits<double>::infinity(); }

	TestFunction(size_t dim, size_t numObjectives = 1) : m_dim(dim), m_numObjectives(numObjectives) { }
	virtual ~TestFunction() { }

	size_t getDimensions() const { return m_dim; }
	size_t getNumObjectives() const { return m_numObjectives; }

	double calculate1(const std::vector<double> &x)
	{
		std::vector<double> result = calculate(x);
		if (result.size() != 1)
			throw std::runtime_error("Expecting single objective result, but got " + std::to_string(result.size()));
		return result[0];
	}

	std::vector<double> calculate(const std::vector<double> &x)
	{
		if (x.size() != m_dim)
			throw std::runtime_error("Invalid input dimension " + std::to_string(x.size()) + ", expecting "+ std::to_string(m_dim));

		// TODO? Always check bounds?

		std::vector<double> result = calculateInternal(x);
		if (result.size() != m_numObjectives)
			throw std::runtime_error("Got result for " + std::to_string(result.size()) + " objectives, but was expecting for " + std::to_string(m_numObjectives));
		return result;
	}

	std::pair<std::vector<double>,std::vector<double>> getBounds()
	{
		auto bounds = getBoundsInternal();
		if (bounds.first.size() != m_dim || bounds.second.size() != m_dim)
			throw std::runtime_error("Invalid bounds dimension");
		for (size_t i = 0 ; i < m_dim ; i++)
			if (bounds.first[i] >= bounds.second[i])
				throw std::runtime_error("All lower bounds must be strictly smaller than upper bounds");
		return bounds;
	}

	std::pair<std::vector<double>,std::vector<double>> getInitialParameterRange()
	{
		auto ipr = getInitialParameterRangeInternal();
		if (ipr.first.size() != m_dim || ipr.second.size() != m_dim)
			throw std::runtime_error("Invalid parameter range dimension");
		for (size_t i = 0 ; i < m_dim ; i++)
			if (ipr.first[i] >= ipr.second[i])
				throw std::runtime_error("All lower parameter values must be strictly smaller than upper ones");
		return ipr;
	}
protected:
	std::vector<double> asVector(double x) const
	{
		return std::vector<double>(m_dim, x);
	}

	virtual std::vector<double> calculateInternal(const std::vector<double> &x) = 0;
	virtual std::pair<std::vector<double>, std::vector<double>> getBoundsInternal() = 0;
	virtual std::pair<std::vector<double>, std::vector<double>> getInitialParameterRangeInternal() = 0;
protected:
	const size_t m_dim, m_numObjectives;
};

class TestFunctionSimpleRanges : public TestFunction
{
public:
	TestFunctionSimpleRanges(size_t dim, std::pair<double, double> ipr,
							 std::pair<double, double> bounds, size_t numObjectives = 1)
		: TestFunction(dim, numObjectives),
		  m_ipr(ipr), m_bounds(bounds)
	{

	}
	
	static inline constexpr double2_t unbounded() { return { -std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity() }; }
protected:
	std::pair<std::vector<double>, std::vector<double>> getBoundsInternal() override
	{
		return { asVector(m_bounds.first), asVector(m_bounds.second) };
	}

	std::pair<std::vector<double>, std::vector<double>> getInitialParameterRangeInternal() override
	{
		return { asVector(m_ipr.first), asVector(m_ipr.second) };
	}
private:
	double2_t m_ipr, m_bounds;
};

// Also f1 from JADE article
// AKA Modified 1st De Jong
class Sphere : public TestFunctionSimpleRanges
{
public:
	Sphere(size_t dim, double2_t ipr, double2_t bounds = unbounded())
		: TestFunctionSimpleRanges(dim, ipr, bounds) { }
	std::vector<double> calculateInternal(const std::vector<double> &x) override
	{
		double sumSquared = 0;
		for (auto v : x)
			sumSquared += v*v;
		return { sumSquared };
	}
};

class Rosenbrock : public TestFunctionSimpleRanges
{
public:
	Rosenbrock(double2_t ipr, double2_t bounds = unbounded())
		: TestFunctionSimpleRanges(2, ipr, bounds) { }
	std::vector<double> calculateInternal(const std::vector<double> &x) override
	{
		double x1 = x[0];
		double x2 = x[1];

		return { (100.0 * (x1*x1 - x2) * (x1*x1 - x2) + (1.0 - x1)*(1.0 - x1)) };
	}
};

// Also f5 from JADE article
class GeneralizedRosenbrock : public TestFunctionSimpleRanges
{
public:
	GeneralizedRosenbrock(size_t dim, double2_t ipr, double2_t bounds = unbounded())
		: TestFunctionSimpleRanges(dim, ipr, bounds) { }
	std::vector<double> calculateInternal(const std::vector<double> &x) override
	{
		double sum = 0;
		for (size_t i = 0 ; i < x.size()-1 ; i++)
			sum += 100.0*(x[i+1]-x[i]*x[i])*(x[i+1]-x[i]*x[i]) + (x[i] - 1.0)*(x[i] - 1.0);
		return { sum };
	}
};

class Mod3rdDeJong : public TestFunctionSimpleRanges
{
public:
	Mod3rdDeJong(double2_t ipr, double2_t bounds = unbounded())
		: TestFunctionSimpleRanges(5, ipr, bounds) { }
	std::vector<double> calculateInternal(const std::vector<double> &x) override
	{
		assert(x.size() == 5);
		double r = 30.0;
		for (auto v : x)
		{
			if (std::abs(v) < 5.12)
				r += std::floor(v);
			else if (v > 5.12)
				r += 30.0*(v - 5.12);
			else
				r += 30.0*(5.12 - v);
		}
		return { r };
	}
};

// Also f2 from JADE article
class Schwefel_2_22 : public TestFunctionSimpleRanges
{
public:
	Schwefel_2_22(size_t dim, double2_t ipr, double2_t bounds = unbounded())
		: TestFunctionSimpleRanges(dim, ipr, bounds) { }

	std::vector<double> calculateInternal(const std::vector<double> &x) override
	{
		double sum = 0, prod = 1;
		for (auto v : x)
		{
			sum += std::abs(v);
			prod *= std::abs(v);
		}
		return { sum + prod };
	}
};

// Also f3 from JADE article
// AKA Double Sum, AKA Rotated Hyper Ellipsoid
class Schwefel_1_2 : public TestFunctionSimpleRanges
{
public:
	Schwefel_1_2(size_t dim, double2_t ipr, double2_t bounds = unbounded())
		: TestFunctionSimpleRanges(dim, ipr, bounds) { }

	std::vector<double> calculateInternal(const std::vector<double> &x) override
	{		
		double sum = 0;
		
		for (size_t i = 0 ; i < x.size() ; i++)
		{
			double sub = 0;
			for (size_t j = 0 ; j <= i ; j++)
				sub += x[j];
			sum += sub*sub;
		}
		return { sum };
	}
};

// Also f4 from JADE article
// AKA Max Mod function
class Schwefel_2_21 : public TestFunctionSimpleRanges
{
public:
	Schwefel_2_21(size_t dim, double2_t ipr, double2_t bounds = unbounded())
		: TestFunctionSimpleRanges(dim, ipr, bounds) { }

	std::vector<double> calculateInternal(const std::vector<double> &x) override
	{
		double max = -std::numeric_limits<double>::max();
		for (auto v : x)
			if (std::abs(v) > max)
				max = std::abs(v);
		return { max };
	}
};

// Also f6 from JADE article
// From "A Literature Survey of Benchmark Functions For Global Optimization Problems"
class Step_2_Function : public TestFunctionSimpleRanges
{
public:
	Step_2_Function(size_t dim, double2_t ipr, double2_t bounds = unbounded())
		: TestFunctionSimpleRanges(dim, ipr, bounds) { }

	std::vector<double> calculateInternal(const std::vector<double> &x) override
	{ 
		double sum = 0;
		for (auto v : x)
			sum += std::floor(v+0.5)*std::floor(v+0.5);
		return { sum };
	}
};

// Also f7 from JADE article
// AKA 4th De Jong
class QuarticWithNoise : public TestFunctionSimpleRanges
{
public:
	QuarticWithNoise(size_t dim, const std::shared_ptr<RandomNumberGenerator> &rng,
					 double2_t ipr, double2_t bounds = unbounded())
		: TestFunctionSimpleRanges(dim, ipr, bounds), m_rng(rng) { }

	std::vector<double> calculateInternal(const std::vector<double> &x) override
	{
		double sum = 0;
		for (size_t i = 0 ; i < x.size() ; i++)
			sum += (i+1.0)*x[i]*x[i]*x[i]*x[i];
		return { sum + m_rng->getRandomDouble() };
	}
private:
	std::shared_ptr<RandomNumberGenerator> m_rng;
};

// Also f8 of JADE article
class ModifiedSchwefel_2_26 : public TestFunctionSimpleRanges
{
public:
	ModifiedSchwefel_2_26(size_t dim, double2_t ipr, double2_t bounds = unbounded())
		: TestFunctionSimpleRanges(dim, ipr, bounds) { }
	
	std::vector<double> calculateInternal(const std::vector<double> &x) override
	{
		double sum = 0;
		for (auto v : x)
			sum += -v*std::sin(std::sqrt(std::abs(v)));
		return { sum + 418.98288727243369*m_dim };
	}
};

// Also f9 from JADE article
class Rastrigin : public TestFunctionSimpleRanges
{
public:
	Rastrigin(size_t dim, double2_t ipr, double2_t bounds = unbounded())
		: TestFunctionSimpleRanges(dim, ipr, bounds) { }

	std::vector<double> calculateInternal(const std::vector<double> &x) override
	{
		double sum = 0;
		for (auto v : x)
			sum += v*v -10.0*std::cos(2.0*M_PI*v) + 10.0;
		return { sum };
	}
};

// Also f10 from JADE article
class AckleyFunction1 : public TestFunctionSimpleRanges
{
public:
	AckleyFunction1(size_t dim, double2_t ipr, double2_t bounds = unbounded())
		: TestFunctionSimpleRanges(dim, ipr, bounds) { }

	std::vector<double> calculateInternal(const std::vector<double> &x) override
	{
		double sum1 = 0, sum2 = 0;

		for (auto v : x)
		{
			sum1 += v*v;
			sum2 += std::cos(2.0*M_PI*v);
		}

		return { -20.0*std::exp(-0.2*std::sqrt(sum1/m_dim)) - std::exp(sum2/m_dim) + 20.0 + std::exp(1) };
	}
};

// Also f11 in JADE article
class Griewank : public TestFunctionSimpleRanges
{
public:
	Griewank(size_t dim, double2_t ipr, double2_t bounds = unbounded())
		: TestFunctionSimpleRanges(dim, ipr, bounds) { }
	
	std::vector<double> calculateInternal(const std::vector<double> &x) override
	{
		double sum = 0, prod = 1;
		for (size_t i = 0 ; i < x.size() ; i++)
		{
			sum += x[i]*x[i];
			prod *= std::cos(x[i]/std::sqrt(i+1.0));
		}
		return { sum/4000.0 - prod + 1.0 };
	}
};

inline double u(double z, double a, double k, double m)
{
	if (z > a)
		return k*std::pow(z-a, m);
	if (z < -a)
		return k*std::pow(-z-a, m);
	return 0;
};

// Also f12 in JADE article
class GeneralizedPenalizedFunction1 : public TestFunctionSimpleRanges
{
public:
	GeneralizedPenalizedFunction1(size_t dim, double2_t ipr, double2_t bounds = unbounded())
		: TestFunctionSimpleRanges(dim, ipr, bounds) { }

	std::vector<double> calculateInternal(const std::vector<double> &x) override
	{
		auto y = [&x](size_t i) { return 1.0+0.25*(x[i] + 1.0); };

		double sum = 10.0*std::sin(M_PI*y(0))*std::sin(M_PI*y(0));

		for (size_t i = 0 ; i < m_dim-1 ; i++)
			sum += (y(i)-1.0)*(y(i)-1.0)*(1.0+10.0*std::sin(M_PI*y(i+1))*std::sin(M_PI*y(i+1)));

		sum += (y(m_dim-1)-1.0)*(y(m_dim-1)-1.0);
		sum *= M_PI/m_dim;

		for (size_t i = 0 ; i < x.size() ; i++)
			sum += u(x[i], 10, 100, 4);
		return { sum };
	}
};

// Also f13 in JADE article
class GeneralizedPenalizedFunction2 : public TestFunctionSimpleRanges
{
public:
	GeneralizedPenalizedFunction2(size_t dim, double2_t ipr, double2_t bounds = unbounded())
		: TestFunctionSimpleRanges(dim, ipr, bounds) { }

	std::vector<double> calculateInternal(const std::vector<double> &x) override
	{
		double sum = std::sin(3*M_PI*x[0])*std::sin(3*M_PI*x[0]);

		for (size_t i = 0 ; i < m_dim-1 ; i++)
			sum += (x[i] - 1.0)*(x[i] - 1.0)*(1.0+std::sin(3*M_PI*x[i+1])*std::sin(3*M_PI*x[i+1]));

		sum += (x[m_dim-1]-1.0)*(x[m_dim-1]-1.0)*(1.0*std::sin(2*M_PI*x[m_dim-1])*std::sin(2*M_PI*x[m_dim-1]));
		sum *= 0.1;

		for (size_t i = 0 ; i < x.size() ; i++)
			sum += u(x[i], 5, 100, 4);

		return { sum };
	}
};

// Also f14 in JADE article
class Branin : public TestFunction
{
public:
	Branin(std::vector<double> lowerIpr, std::vector<double> upperIpr,
		   std::vector<double> lowerBound = { -infinity(), -infinity() },
		   std::vector<double> upperBound = { infinity(), infinity() })
		: TestFunction(2), m_iprLow(lowerIpr), m_iprHigh(upperIpr), 
		  m_boundLow(lowerBound), m_boundHigh(upperBound) { }

	std::vector<double> calculateInternal(const std::vector<double> &x) override
	{
		const double a = 1.0;
		const double b = 5.1/(4.0*M_PI*M_PI);
		const double c = 5.0/M_PI;
		const double r = 6.0;
		const double s = 10.0;
		const double t = 1.0/(8.0*M_PI);

		return { a*std::pow(x[1] - b*x[0]*x[0] + c*x[0] - r,2) + s*(1.0-t)*std::cos(x[0]) + s };
	}

	std::pair<std::vector<double>, std::vector<double>> getBoundsInternal() override
	{
		return { m_boundLow, m_boundHigh };
	}

	std::pair<std::vector<double>, std::vector<double>> getInitialParameterRangeInternal() override
	{
		return { m_iprLow, m_iprHigh };
	}
private:
	std::vector<double> m_boundLow, m_boundHigh;
	std::vector<double> m_iprLow, m_iprHigh;
};

// Also f15 in JADE article
class GoldsteinPrice : public TestFunctionSimpleRanges
{
public:
	GoldsteinPrice(double2_t ipr, double2_t bounds = unbounded())
		: TestFunctionSimpleRanges(2, ipr, bounds) { }

	std::vector<double> calculateInternal(const std::vector<double> &x) override
	{
		return { (1.0+std::pow(x[0]+x[1]+1.0,2)*(19.0-14.0*x[0]+3.0*x[0]*x[0]-14.0*x[1]+6.0*x[0]*x[1]+3.0*x[1]*x[1]))*(
			30.0+std::pow(2.0*x[0]-3.0*x[1],2)*(18.0-32.0*x[0]+12.0*x[0]*x[0]+48.0*x[1]-36.0*x[0]*x[1]+27.0*x[1]*x[1])
		) };
	}
};

class Hartman3D : public TestFunctionSimpleRanges
{
public:
	Hartman3D(double2_t ipr, double2_t bounds = unbounded())
		: TestFunctionSimpleRanges(3, ipr, bounds)
	{
		alpha = { 1.0, 1.2, 3.0, 3.2 };
		A = {
			{ 3.0, 10.0, 30.0},
			{ 0.1, 10.0, 35.0},
			{ 3.0, 10.0, 30.0},
			{ 0.1, 10.0, 35.0}
		};

		P = {
			{ 3689, 1170, 2673 },
			{ 4699, 4387, 7470 },
			{ 1091, 8732, 5547 },
			{  381, 5743, 8828 }
		};
		for (auto &v : P)
			for (auto &x : v)
				x *= 1e-4;
	}
	
	std::vector<double> calculateInternal(const std::vector<double> &x) override
	{
		double sum = 0;
		for (size_t i = 0 ; i < 4 ; i++)
		{
			double s = 0;
			for (size_t j = 0 ; j < 3 ; j++)
				s += A[i][j]*std::pow(x[j]-P[i][j],2);

			sum += alpha[i]*std::exp(-s);
		}
		return { -sum };
	}
private:
	std::vector<double> alpha;
	std::vector<std::vector<double>> A, P;
};

class Hartman6D : public TestFunctionSimpleRanges
{
public:
	Hartman6D(double2_t ipr, double2_t bounds = unbounded())
		: TestFunctionSimpleRanges(6, ipr, bounds)
	{
		alpha = { 1.0, 1.2, 3.0, 3.2 };
		A = {
			{ 10, 3, 17, 3.5, 1.7, 8 },
			{ 0.05, 10, 17, 0.1, 8, 14 },
			{ 3, 3.5, 1.7, 10, 17, 8 },
			{ 17, 8, 0.05, 10, 0.1, 14 }
		};

		P = {
			{ 1312, 1696, 5569, 124, 8283, 5886 },
			{ 2329, 4135, 8307, 3736, 1004, 9991 },
			{ 2348, 1451, 3522, 2883, 3047, 6650 },
			{ 4047, 8828, 8732, 5743, 1091, 381 }
		};
		for (auto &v : P)
			for (auto &x : v)
				x *= 1e-4;
	}

	std::vector<double> calculateInternal(const std::vector<double> &x) override
	{
		double sum = 0;
		for (size_t i = 0 ; i < 4 ; i++)
		{
			double s = 0;
			for (size_t j = 0 ; j < 6 ; j++)
				s += A[i][j]*std::pow(x[j]-P[i][j],2);

			sum += alpha[i]*std::exp(-s);
		}
		return { -sum };
	}
private:
	std::vector<double> alpha;
	std::vector<std::vector<double>> A, P;
};

class Shekel : public TestFunctionSimpleRanges
{
public:
	Shekel(size_t m, double2_t ipr, double2_t bounds = unbounded())
		: TestFunctionSimpleRanges(4, ipr, bounds), m_m(m)
	{
		beta = { 1, 2, 2, 4, 4, 6, 3, 7, 5, 5 };
		for (auto &b : beta)
			b *= 0.1;
		
		C = {
			{ 4, 1, 8, 6, 3, 2, 5, 8, 6, 7 },
			{ 4, 1, 8, 6, 7, 9, 3, 1, 2, 3.6 },
			{ 4, 1, 8, 6, 3, 2, 5, 8, 6, 7 },
			{ 4, 1, 8, 6, 7, 9, 3, 1, 2, 3.6 }
		};
	}

	std::vector<double> calculateInternal(const std::vector<double> &x) override
	{
		assert(m_m <= 10);
		double sum = 0;
		for (size_t i = 0 ; i < m_m ; i++)
		{
			double s = beta[i];
			for (size_t j = 0 ; j < 4 ; j++)
				s += std::pow(x[j] - C[j][i],2);

			sum += 1.0/s;
		}
		return { -sum };
	}
private:
	const size_t m_m;
	std::vector<double> beta;
	std::vector<std::vector<double>> C;
};

} // end namespace testfunctions

} // end namespace eatk
