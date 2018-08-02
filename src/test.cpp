
#include <stdio.h>
#include <malloc.h>
#include <memory>
#include <random>
#include <chrono>
#include "ImageWriter.h"

const double PI = 3.1415926535897323;

struct PixelPoint
{
    int x; 
    int y;
};

template<typename T>
T clip(T val, T min, T max)
{
    if (val < min)
        return min;
    if (val > max)
        return max;
    return val;
}




PixelFloat rainbow(double val, double period)
{
    double sixthPeriod = period / 6;

    val = fmod(val, period);

    PixelFloat result{ 0, 0, 0, 1 };

    if (val >= 0 && val < sixthPeriod)
    {
        result.r = 1.0;
        result.g = val / sixthPeriod;
    }
    if (val >= sixthPeriod && val < 2 * sixthPeriod)
    {
        result.r = 1.0 - (val - sixthPeriod) / sixthPeriod;
        result.g = 1.0;
    }
    if (val >= 2 * sixthPeriod && val < 3 * sixthPeriod)
    {
        result.g = 1.0;
        result.b = (val - 2 * sixthPeriod) / sixthPeriod;
    }
    if (val >= 3 * sixthPeriod && val < 4 * sixthPeriod)
    {
        result.g = 1.0 - (val - 3 * sixthPeriod) / sixthPeriod;
        result.b = 1.0;
    }
    if (val >= 4 * sixthPeriod && val < 5 * sixthPeriod)
    {
        result.b = 1.0;
        result.r = (val - 4 * sixthPeriod) / sixthPeriod;
    }
    if (val >= 5 * sixthPeriod && val < 6 * sixthPeriod)
    {
        result.b = 1.0 - (val - 5 * sixthPeriod) / sixthPeriod;
        result.r = 1.0;
    }

    return result;
}



struct MandelbrotParameters
{
    double left;
    double right;
    double top;
    double bottom;
    int maxIt;
};

PixelPoint inverseTransform(double x, double y, MandelbrotParameters params, int width, int height)
{

    int px = (x - params.left) / (params.right - params.left) * width;
    px = clip(px, 0, width - 1);

    int py = (y - params.top) / (params.bottom - params.top) * height;
    py = clip(py, 0, height - 1);

    return { px, py };
}


int mandelbrot(double x0, double y0, int maxIt)
{
    double x = 0;
    double y = 0;
    for (int it = 0; it < maxIt; ++it)
    {
        double xTemp = x * x - y * y + x0;
        y = 2 * x*y + y0;
        x = xTemp;

        if (x*x + y * y > 2 * 2)
        {
            return it;
        }
    }
    return -1;
}


std::shared_ptr<BufferFloat> renderMandelbrot(int width, int height, MandelbrotParameters params)
{
    auto buffer = std::make_unique<BufferFloat>(width, height);

    for (int py = 0; py < height; ++py)
    {
        for (int px = 0; px < width; ++px)
        {
            double x0 = params.left + (px * (params.right - params.left) / width);
            double y0 = params.top + (py * (params.bottom - params.top) / height);

            int val = mandelbrot(x0, y0, params.maxIt);

            if (val == -1)
            {
                buffer->set(px, py, { 0, 0, 0, 1 });
            }
            else
            {
                buffer->set(px, py, { val / 100.0, val / 100.0, val / 100.0, 1 });
            }
        }
    }

    return buffer;
}

void brownianWalk(std::shared_ptr<BufferFloat> buffer, MandelbrotParameters params, int iterations)
{
    const double outwardForce = 0.005;
    const double maxOutwardForce = 0.05;

    const double brownianVariance = 0.01;
    const double maxBrownianForce = 0.1;

    const double slopeForce = 0.1;
    const double maxSlopForce = 0.2;

    const double coriolisForce = 0.01;
    const double maxCoriolisForce = 0.1;

    const double pointsOnLine = 100;


    double x = 0.0;
    double y = 0.0;

    int pxPrev = buffer->getWidth() * 2 / 3;
    int pyPrev = buffer->getHeight() / 2;


    static std::default_random_engine generator(std::chrono::system_clock::now().time_since_epoch().count());
    static std::normal_distribution<double> normDist(0, brownianVariance);

    for (int i = 0; i < iterations; ++i)
    {

        double xPrev = x;
        double yPrev = y;

        double dx = 0;
        double dy = 0;

        // brownian motion
        dx += normDist(generator);
        dy += normDist(generator);

        int escape = mandelbrot(x, y, params.maxIt);
        // away from origin
        if (escape == -1)
        {
            // main lobe
            if (x > -0.75)
            {
                if (x != 0)
                    dx += clip(outwardForce / x, -maxOutwardForce, maxOutwardForce);
                if (y != 0)
                    dy += clip(outwardForce / y, -maxOutwardForce, maxOutwardForce);
            }
            // secondary lobe
            else
            {
                if (x + 0.75 != 0)
                    dx += clip(outwardForce / (x + 0.75), -maxOutwardForce, maxOutwardForce);
                if (y != 0)
                    dy += clip(outwardForce / y, -maxOutwardForce, maxOutwardForce);
            }
        }

        // down slope
        if (escape > 0)
        {
            dx -= clip(slopeForce * (1.0 / escape) * x, -maxSlopForce, maxSlopForce);
            dy -= clip(slopeForce * (1.0 / escape) * y, -maxSlopForce, maxSlopForce);
        }

        // coriolis
        // +x -> +y
        // +y -> -x
        // -x -> -y
        // -y -> +x
        dx -= coriolisForce * y;
        dy += coriolisForce * x;


        x += dx;
        y += dy;


        // draw a line
        for (int j = 0; j < pointsOnLine; ++j)
        {
            double xInter = xPrev * (j / pointsOnLine) + x * (1 - j / pointsOnLine);
            double yInter = yPrev * (j / pointsOnLine) + y * (1 - j / pointsOnLine);

            PixelPoint pInter = inverseTransform(xInter, yInter, params, buffer->getWidth(), buffer->getHeight());

            buffer->set(pInter.x, pInter.y, rainbow(i, 500));
        }

    }
}


int main(int argc, char *argv[])
{
    // https://www.marksmath.org/visualization/julia_sets/



    // Make sure that the output filename argument has been provided
    if (argc != 2) {
        fprintf(stderr, "Please specify output file\n");
        return 1;
    }

    // Specify an output image size
    int width = 800;
    int height = 600;


    MandelbrotParameters params = { -2, 1, -1.125, 1.125, 1000 };


    auto buffer = renderMandelbrot(width, height, params);

    brownianWalk(buffer, params, 1000);

    // Save the image to a PNG file
    // The 'title' string is stored as part of the PNG file
    printf("Saving PNG\n");
    int result = writeImage(argv[1], *buffer, "This is my test image");

    // Free up the memorty used to store the image
    //free(buffer);

    return result;
}
