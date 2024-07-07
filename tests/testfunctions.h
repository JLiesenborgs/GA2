#include "eatkconfig.h"
#include <vector>
#include <stdexcept>
#include <cassert>
#include <cmath>

namespace eatk
{

namespace testfunctions
{

class TestFunction
{
public:
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
private:
    size_t m_dim, m_numObjectives;
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
protected:
    std::pair<std::vector<double>, std::vector<double>> getBoundsInternal() override
    {
        return { asVector(m_ipr.first), asVector(m_ipr.second) };
    }

    std::pair<std::vector<double>, std::vector<double>> getInitialParameterRangeInternal() override
    {
        return { asVector(m_bounds.first), asVector(m_bounds.second) };
    }
private:
    std::pair<double,double> m_ipr, m_bounds;
};

class Sphere : public TestFunctionSimpleRanges
{
public:
    Sphere(size_t dim, std::pair<double,double> ipr, std::pair<double,double> bounds) : TestFunctionSimpleRanges(dim, ipr, bounds) { }
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
    Rosenbrock(std::pair<double,double> ipr, std::pair<double,double> bounds) : TestFunctionSimpleRanges(2, ipr, bounds) { }
    std::vector<double> calculateInternal(const std::vector<double> &x) override
    {
        double x1 = x[0];
		double x2 = x[1];

		return { (100.0 * (x1*x1 - x2) * (x1*x1 - x2) + (1.0 - x1)*(1.0 - x1)) };
    }
};

class Mod3rdDeJong : public TestFunctionSimpleRanges
{
public:
    Mod3rdDeJong(std::pair<double,double> ipr, std::pair<double,double> bounds) : TestFunctionSimpleRanges(5, ipr, bounds) { }
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


} // end namespace testfunctions

} // end namespace eatk