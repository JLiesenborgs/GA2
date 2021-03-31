#include "mersennerandomnumbergenerator.h"

namespace eatk
{

MersenneRandomNumberGenerator::MersenneRandomNumberGenerator(unsigned int seed) : m_rng(seed)
{

}

MersenneRandomNumberGenerator::~MersenneRandomNumberGenerator()
{

}

double MersenneRandomNumberGenerator::getRandomDouble()
{
    std::uniform_real_distribution<double> u(0.0, 1.0);
    return u(m_rng);
}

double MersenneRandomNumberGenerator::getRandomDouble(double min, double max)
{
    std::uniform_real_distribution<double> u(min, max);
    return u(m_rng);
}

float MersenneRandomNumberGenerator::getRandomFloat()
{
    std::uniform_real_distribution<float> u(0.0, 1.0);
    return u(m_rng);
}

float MersenneRandomNumberGenerator::getRandomFloat(float min, float max)
{
    std::uniform_real_distribution<float> u(min, max);
    return u(m_rng);
}

uint32_t MersenneRandomNumberGenerator::getRandomUint32()
{
    return m_rng();
}

}
