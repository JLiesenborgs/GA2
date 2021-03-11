#pragma once

#include "mogal2config.h"
#include <stdint.h>

namespace mogal2
{

class RandomNumberGenerator
{
public:
    RandomNumberGenerator() { }
    virtual ~RandomNumberGenerator() { }

    virtual double getRandomDouble() = 0;
    virtual double getRandomDouble(double min, double max)
    {
        double x = getRandomDouble();
        double diff = max-min;
        return x*diff + min;
    }
    virtual float getRandomFloat() = 0;
    virtual float getRandomFloat(float min, float max)
    {
        float x = getRandomFloat();
        float diff = max-min;
        return x*diff + min;
    }
    virtual uint32_t getRandomUint32() = 0;
};

}