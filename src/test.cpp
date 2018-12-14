
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
#include <SDL.h>
#undef main

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
            // Swapping to rotate 90 degrees
            double y0 = params.left + (px * (params.right - params.left) / width);
            double x0 = params.top + (py * (params.bottom - params.top) / height);
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

            buffer.set(inter.x, inter.y, rainbow(i, path.size()));
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
    //std::string filename;
    uint8_t* pixels;
    int width;
    int height;
    int supersample;
    Complex c;
    MandelbrotParameters juliaParams;
};

void runJob(ThreadArgs args)
{
    auto juliaBuffer = renderJulia(args.width * args.supersample, args.height * args.supersample, args.c, args.juliaParams);

    downsampleBuffer(juliaBuffer, args.supersample);


    int i = 0;
    for (int y = 0; y < args.height; ++y)
    {
        for (int x = 0; x < args.width; ++x)
        {
            auto px = juliaBuffer->get(x, y);
            args.pixels[i++] = px.r;
            args.pixels[i++] = px.g;
            args.pixels[i++] = px.b;
            //pixelPointer2[i++] = 255;
        }
    }

}

//
//int main(int argc, char *argv[])
//{
//    // https://www.marksmath.org/visualization/julia_sets/
//    // https://www.desmos.com/calculator/44ka3v8igm
//
//
//
//    //// Make sure that the output filename argument has been provided
//    //if (argc != 2) {
//    //    fprintf(stderr, "Please specify output file\n");
//    //    return 1;
//    //}
//
//    bool quit = false;
//    SDL_Event event;
//
//    SDL_Init(SDL_INIT_VIDEO);
//
//    SDL_Window * window = SDL_CreateWindow("SDL2 Displaying Image",
//        SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, 640, 480, 0);
//
//
//    // Specify an output image size
//    int supersample = 2;
//    int width = 800;
//    int height = 600;
//
//
//    MandelbrotParameters params = { -2, 1, -1.125, 1.125, 1000 };
//    //MandelbrotParameters params = { -1.5, 1.5, -1.125, 1.125, 1000 };
//
//
//    auto buffer = renderMandelbrot(width, height, params);
//
//    double repeatPeriodYears = 0.2;
//    double year = 1.0 / repeatPeriodYears;
//    double day = year / 365;
//    double hour = day / 24;
//    double minute = hour / 60;
//    double second = minute / 60;
//
//    double start = 60 * day * repeatPeriodYears;
//    double length = 1 * hour;
//    int steps = 600;
//
//    //auto path = brownianWalk(params, length, steps);
//    double distance = 0.7885;
//    //std::vector<Complex> path = circlePath<double>(distance, 2.26, 2.28, steps);
//    std::vector<Epicycle> epicycles {
//        Epicycle{ .5, 2 * PI, 0 },
//        Epicycle{ .25, 4 * PI, PI  },
//        Epicycle{ 0.01, 1000, 0 }
//    };
//    auto path = epicyclePath(epicycles, start, start + length, steps);
//
//    drawPath(*buffer, path, params);
//
//    printf("Saving PNG\n");
//    int result = writeImage("output.png", *buffer);
//
//
//    SDL_Renderer * renderer = SDL_CreateRenderer(window, -1, 0);
//
//    //SDL_Surface * image = SDL_CreateRGBSurface(0, width, height, 32, 0, 0, 0, 0);
//
//
//    uint8_t* pixels = new uint8_t[width * height * 4];
//
//    int i = 0;
//    for (int y = 0; y < height; ++y)
//    {
//        for (int x = 0; x < width; ++x)
//        {
//            auto px = buffer->get(x, y);
//            pixels[i++] = px.r;
//            pixels[i++] = px.g;
//            pixels[i++] = px.b;
//            pixels[i++] = 255;
//        }
//    }
//
//    SDL_Surface * image = SDL_CreateRGBSurfaceFrom(pixels, width, height, 32, 4 * width, 0, 0, 0, 0);
//
//    SDL_Texture * texture = SDL_CreateTextureFromSurface(renderer, image);
//
//
//
//
//
//
//    while (!quit)
//    {
//        SDL_WaitEvent(&event);
//
//        switch (event.type)
//        {
//        case SDL_QUIT:
//            quit = true;
//            break;
//        }
//
//
//        SDL_RenderCopy(renderer, texture, NULL, NULL);
//        SDL_RenderPresent(renderer);
//
//    }
//
//
//
//
//    MandelbrotParameters juliaParams = { -1.5, 1.5, -1.125, 1.125, 1000 };
//
//
//    int maxThreads = 8;
//    std::vector<std::vector<ThreadArgs>> argsLists;
//    for (int i = 0; i < maxThreads; ++i)
//    {
//        argsLists.emplace_back();
//    }
//
//    for (int i = 0; i < path.size(); ++i)
//    {
//        ThreadArgs args;
//        args.filename = "frame" + std::to_string(i) + ".png";
//        args.width = width;
//        args.height = height;
//        args.supersample = supersample;
//        args.c = path[i];
//        args.juliaParams = juliaParams;
//
//        argsLists[i % maxThreads].emplace_back(args);
//    }
//
//    std::vector<std::thread> runningJobs;
//    for (auto argList : argsLists)
//    {
//        runningJobs.emplace_back(runJob, argList);
//    }
//
//    for (auto& job : runningJobs)
//    {
//        job.join();
//    }
//
//
//    //return result;
//    return 0;
//}


int main()
{
    int width = 1920;
    int height = 1080;



    bool quit = false;
    SDL_Event event;

    SDL_Init(SDL_INIT_VIDEO);



    auto win = SDL_CreateWindow("Julia Set", 0, 0, width, height, SDL_WINDOW_SHOWN);
    if (win == nullptr)
    {
        std::cout << "SDL_CreateWindow Error: " << SDL_GetError() << std::endl;
        SDL_Quit();
        return {};
    }

    // Create the Renderer that will render our image into the window.
    SDL_Renderer* ren = SDL_CreateRenderer(win, -1, SDL_RENDERER_ACCELERATED | SDL_RENDERER_PRESENTVSYNC);
    if (ren == nullptr)
    {
        SDL_DestroyWindow(win);
        std::cout << "SDL_CreateRenderer Error: " << SDL_GetError() << std::endl;
        SDL_Quit();
        return {};
    }

    SDL_Texture* tex = SDL_CreateTexture(ren, SDL_PIXELFORMAT_RGB24, SDL_TEXTUREACCESS_STREAMING, width, height);
    if (tex == nullptr)
    {
        std::cout << "SDL_CreateTexture Error: " << SDL_GetError() << std::endl;
        SDL_Quit();
        return {};
    }




    int supersample = 2;

    double repeatPeriodYears = 0.2;
    double year = 1.0 / repeatPeriodYears;
    double day = year / 365;
    double hour = day / 24;
    double minute = hour / 60;
    double second = minute / 60;

    double start = 60 * day * repeatPeriodYears;
    double length = 16 * hour;
    int steps = 600;

    //auto path = brownianWalk(params, length, steps);
    double distance = 0.7885;
    //std::vector<Complex> path = circlePath<double>(distance, 2.26, 2.28, steps);
    std::vector<Epicycle> epicycles {
        Epicycle{ .5, 2 * PI, 0 },
        Epicycle{ .25, 4 * PI, PI  },
        Epicycle{ 0.01, 1000, 0 }
    };

    auto path = epicyclePath(epicycles, start, start + length, steps);




    int maxThreads = 12;
    std::vector<uint8_t*> pixelBuffers;
    int eachBufferLength = width * height * 3 / maxThreads;
    for (int i = 0; i < maxThreads; ++i)
    {
        pixelBuffers.emplace_back(new uint8_t[eachBufferLength]);
    }


    int pathIndex = 0;



    while (!quit)
    {
        if (SDL_WaitEventTimeout(&event, 1))
        {
            switch (event.type)
            {
            case SDL_QUIT:
                quit = true;
                break;
            }
        }


        double y = 1.125;
        double x = y * width / height;

        MandelbrotParameters juliaParams = { -x, x, -y, y, 100 };
        auto c = path[pathIndex++];

        std::cout << "c = " << c.r << " + " << c.i << "i\n";





        std::vector<std::thread> runningJobs;
        for (int i = 0; i < maxThreads; ++i)
        {
            //auto pixBuf = ;
            //pixelBuffers.emplace_back(pixBuf);


            ThreadArgs args;
            //args.filename = "frame" + std::to_string(i) + ".png";
            args.pixels = pixelBuffers[i];
            args.width = width;
            args.height = height / maxThreads;
            args.supersample = supersample;
            args.c = c;
            args.juliaParams = juliaParams;


            double top = args.juliaParams.top + i * (args.juliaParams.bottom - args.juliaParams.top) / maxThreads;
            double bottom = args.juliaParams.top + (i + 1) * (args.juliaParams.bottom - args.juliaParams.top) / maxThreads;


            args.juliaParams.top = top;
            args.juliaParams.bottom = bottom;


            runningJobs.emplace_back(runJob, args);
            //argsLists[i % maxThreads].emplace_back(args);
        }

        for (auto& job : runningJobs)
        {
            job.join();
        }





        // Lock the texture so we can write to it.
        void* pixels = nullptr;
        auto pitch = 0;
        if (SDL_LockTexture(tex, nullptr, &pixels, &pitch) != 0)
        {
            SDL_DestroyTexture(tex);
            SDL_DestroyRenderer(ren);
            SDL_DestroyWindow(win);
            std::cout << "SDL_LockTexture Error: " << SDL_GetError() << std::endl;
            SDL_Quit();
            return {};
        }







        // Push the pixels to the texture.
        for (int i = 0; i < maxThreads; ++i)
        {
            void* subpixels = (void*)((uint8_t*)pixels + i * eachBufferLength);
            memcpy(subpixels, pixelBuffers[i], eachBufferLength);
        }
        //memcpy(pixels, pixelPointer2, pitch * height);

        // Unlock the texture so that it updates.
        SDL_UnlockTexture(tex);

        // ReSharper disable once CppExpressionWithoutSideEffects

        //First clear the renderer
        SDL_RenderClear(ren);
        //Draw the texture
        SDL_RenderCopy(ren, tex, nullptr, nullptr);
        //Update the screen
        SDL_RenderPresent(ren);
    }

    return 0;
}