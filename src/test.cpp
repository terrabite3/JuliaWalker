
#include <stdio.h>
#include <malloc.h>
#include <memory>
#include <random>
#include <chrono>
#include "ImageWriter.h"
#include <string>
#include <iostream>
#include "EpicyclePath.h"
#include "Complex.h"
#include <thread>

const double PI = 3.1415926535897323;

struct PixelPoint
{
    int x; 
    int y;
};


double clip(double val, double min, double max)
{
    if (val < min)
        return min;
    if (val > max)
        return max;
    return val;
}




Pixel rainbow(double val, double period)
{
    double sixthPeriod = period / 6;

    val = fmod(val, period);

    double r = 0.0;
    double g = 0.0;
    double b = 0.0;

    if (val >= 0 && val < sixthPeriod)
    {
        r = 1.0;
        g = val / sixthPeriod;
    }
    if (val >= sixthPeriod && val < 2 * sixthPeriod)
    {
        r = 1.0 - (val - sixthPeriod) / sixthPeriod;
        g = 1.0;
    }
    if (val >= 2 * sixthPeriod && val < 3 * sixthPeriod)
    {
        g = 1.0;
        b = (val - 2 * sixthPeriod) / sixthPeriod;
    }
    if (val >= 3 * sixthPeriod && val < 4 * sixthPeriod)
    {
        g = 1.0 - (val - 3 * sixthPeriod) / sixthPeriod;
        b = 1.0;
    }
    if (val >= 4 * sixthPeriod && val < 5 * sixthPeriod)
    {
        b = 1.0;
        r = (val - 4 * sixthPeriod) / sixthPeriod;
    }
    if (val >= 5 * sixthPeriod && val < 6 * sixthPeriod)
    {
        b = 1.0 - (val - 5 * sixthPeriod) / sixthPeriod;
        r = 1.0;
    }

    return Pixel::fromNormalizedFloat(r, g, b);
}



struct MandelbrotParameters
{
    double left;
    double right;
    double top;
    double bottom;
    int maxIt;
};

PixelPoint inverseTransform(Complex z, MandelbrotParameters params, int width, int height)
{

    int px = (z.r - params.left) / (params.right - params.left) * width;
    px = clip(px, 0, width - 1);

    int py = (z.i - params.top) / (params.bottom - params.top) * height;
    py = clip(py, 0, height - 1);

    return { px, py };
}


int mandelbrot(const Complex& c, int maxIt)
{
    Complex z(0, 0);
    for (int it = 0; it < maxIt; ++it)
    {
        z = z * z + c;

        if (z.r*z.r + z.i * z.i > 2 * 2)
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

            Complex z0(x0, y0);

            int val = mandelbrot(z0, params.maxIt);

            if (val == -1)
            {
                buffer->set(px, py, Pixel(0, 0, 0));
            }
            else
            {
                uint8_t grayValue = clip(val * 255 / 100, 0, 255);
                buffer->set(px, py, Pixel(grayValue, grayValue, grayValue));
            }
        }
    }

    return buffer;
}

int julia(const Complex& z0, const Complex& c, int maxIt)
{
    Complex z = z0;
    for (int it = 0; it < maxIt; ++it)
    {
        z = z * z + c;

        if (z.r*z.r + z.i*z.i > 2*2)
        {
            return it;
        }
    }
    return -1;
}

std::shared_ptr<BufferFloat> renderJulia(int width, int height, const Complex& c, MandelbrotParameters params)
{
    auto buffer = std::make_unique<BufferFloat>(width, height);

    for (int py = 0; py < height; ++py)
    {
        for (int px = 0; px < width; ++px)
        {
            double x0 = params.left + (px * (params.right - params.left) / width);
            double y0 = params.top + (py * (params.bottom - params.top) / height);
            Complex z0(x0, y0);

            int val = julia(z0, c, params.maxIt);

            if (val == -1)
            {
                buffer->set(px, py, Pixel(0, 0, 0));
            }
            else
            {
                buffer->set(px, py, rainbow(val, 100));
            }
        }
    }

    return buffer;
}

std::vector<Complex> brownianWalk(MandelbrotParameters params, double length, int iterations)
{
    std::vector<Complex> result;
    result.emplace_back(0, 0);

    double stepCoef = length / iterations;

    const double outwardForce = 0.5 * stepCoef;
    const double maxOutwardForce = 5 * stepCoef;

    const double brownianVariance = 1 * stepCoef;
    //const double maxBrownianForce = 0.1;

    const double slopeForce = 10 * stepCoef;
    const double maxSlopForce = 20 * stepCoef;

    const double coriolisForce = 1 * stepCoef;
    const double maxCoriolisForce = 10 * stepCoef;



    double x = 0.5;
    double y = 1.0;

    //int pxPrev = buffer->getWidth() * 2 / 3;
    //int pyPrev = buffer->getHeight() / 2;


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

        int escape = mandelbrot(Complex(x, y), params.maxIt);
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


        result.emplace_back(x, y);

    }

    return result;
}

void drawPath(BufferFloat& buffer, const std::vector<Complex>& path, const MandelbrotParameters& params)
{
    const double pointsOnLine = 1;
    
    Complex dpPrev = path[0];

    int i = 0;
    for (Complex dp : path)
    {
        // draw a line
        for (int j = 0; j < pointsOnLine; ++j)
        {
            Complex interpolatePoint = dpPrev * (j / pointsOnLine) + dp * (1 - j / pointsOnLine);

            PixelPoint inter = inverseTransform(interpolatePoint, params, buffer.getWidth(), buffer.getHeight());

            buffer.set(inter.x, inter.y, rainbow(i, 100));
        }

        dpPrev = dp;
        ++i;
    }
}



void downsampleBuffer(std::shared_ptr<BufferFloat> & buffer, int supersample)
{
    if (supersample == 1) return;

    int downsampledWidth = buffer->getWidth() / supersample;
    int downsampledHeight = buffer->getHeight() / supersample;

    auto tempBuffer = std::make_shared<BufferFloat>(downsampledWidth, downsampledHeight);


    for (int y = 0; y < downsampledHeight; ++y)
    {
        for (int x = 0; x < downsampledWidth; ++x)
        {
            double r = 0, g = 0, b = 0, a = 0;
            for (int j = 0; j < supersample; ++j)
            {
                for (int i = 0; i < supersample; ++i)
                {
                    Pixel subpixel = buffer->get(x * supersample + i, y * supersample + j);
                    r += subpixel.r;
                    g += subpixel.g;
                    b += subpixel.b;
                }
            }
            r /= supersample * supersample;
            g /= supersample * supersample;
            b /= supersample * supersample;

            tempBuffer->set(x, y, Pixel{ (uint8_t)r, (uint8_t)g, (uint8_t)b });
        }
    }


    buffer = tempBuffer;
}

struct ThreadArgs
{
    std::string filename;
    int width;
    int height;
    int supersample;
    Complex c;
    MandelbrotParameters juliaParams;
};

void runJob(std::vector<ThreadArgs> argList)
{
    for (auto args : argList)
    {
        auto juliaBuffer = renderJulia(args.width * args.supersample, args.height * args.supersample, args.c, args.juliaParams);

        downsampleBuffer(juliaBuffer, args.supersample);

        writeImage(args.filename, *juliaBuffer);

        std::string doneString = "Finished frame " + args.filename + "\n";

        std::cout << doneString;
    }
}


int main(int argc, char *argv[])
{
    // https://www.marksmath.org/visualization/julia_sets/
    // https://www.desmos.com/calculator/44ka3v8igm



    // Make sure that the output filename argument has been provided
    if (argc != 2) {
        fprintf(stderr, "Please specify output file\n");
        return 1;
    }

    // Specify an output image size
    int supersample = 2;
    int width = 800;
    int height = 600;


    MandelbrotParameters params = { -2, 1, -1.125, 1.125, 1000 };
    //MandelbrotParameters params = { -1.5, 1.5, -1.125, 1.125, 1000 };


    auto buffer = renderMandelbrot(width, height, params);

    double length = 0.3;
    int steps = 3142;
    //auto path = brownianWalk(params, length, steps);
    double distance = 0.7885;
    //std::vector<Complex> path = circlePath<double>(distance, 2.26, 2.28, steps);
    std::vector<Epicycle> epicycles {
        Epicycle{ .5, 2 * PI, -PI/2 },
        Epicycle{ .25, 4 * PI, 0  }
    };
    auto path = epicyclePath(epicycles, 0.25, 0.75, steps);

    drawPath(*buffer, path, params);


    printf("Saving PNG\n");
    int result = writeImage(argv[1], *buffer);

    MandelbrotParameters juliaParams = { -1.5, 1.5, -1.125, 1.125, 1000 };


    int maxThreads = 8;
    std::vector<std::vector<ThreadArgs>> argsLists;
    for (int i = 0; i < maxThreads; ++i)
    {
        argsLists.emplace_back();
    }

    for (int i = 0; i < path.size(); ++i)
    {
        ThreadArgs args;
        args.filename = "frame" + std::to_string(i) + ".png";
        args.width = width;
        args.height = height;
        args.supersample = supersample;
        args.c = path[i];
        args.juliaParams = juliaParams;

        argsLists[i % maxThreads].emplace_back(args);
    }

    std::vector<std::thread> runningJobs;
    for (auto argList : argsLists)
    {
        runningJobs.emplace_back(runJob, argList);
    }

    for (auto& job : runningJobs)
    {
        job.join();
    }


    return result;
}
