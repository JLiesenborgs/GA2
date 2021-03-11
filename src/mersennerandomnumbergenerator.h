#pragma once

#include "mogal2config.h"
#include "randomnumbergenerator.h"
#include <random>

namespace mogal2
{

class MersenneRandomNumberGenerator : public RandomNumberGenerator
{
public:
    MersenneRandomNumberGenerator(unsigned int seed);
    ~MersenneRandomNumberGenerator();

    double getRandomDouble() override;
    double getRandomDouble(double min, double max) override;
    float getRandomFloat() override;
    float getRandomFloat(float min, float max) override;
    uint32_t getRandomUint32() override;
private:
    std::mt19937 m_rng;
};

}
