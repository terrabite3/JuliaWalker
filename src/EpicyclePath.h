#pragma once

#include "Complex.h"
#include <vector>

struct Epicycle
{
    double radius;
    double angularVelocity;
    double phase;
};

std::vector<Complex> epicyclePath(const std::vector<Epicycle>& epicycles, double t0, double t1, int steps);