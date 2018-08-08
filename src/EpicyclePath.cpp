#include "EpicyclePath.h"


std::vector<Complex> epicyclePath(const std::vector<Epicycle>& epicycles, double t0, double t1, int steps)
{
    std::vector<Complex> path;

    for (int i = 0; i < steps; ++i)
    {
        double t = t0 + (t1 - t0) * (1.0 * i / steps);

        Complex sum(0, 0);
        for (Epicycle e : epicycles)
        {
            sum += Complex(e.radius * cos(t * e.angularVelocity + e.phase), e.radius * sin(t * e.angularVelocity + e.phase));
        }
        path.emplace_back(sum);
    }

    return path;
}
