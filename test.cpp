
#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include <png.h>
#include <vector>
#include <cassert>
#include <memory>
#include <iostream>
#include <random>
#include <chrono>

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

//struct Point
//{
//    Point(double x, double y):
//        x(x),
//        y(y)
//    {}
//
//    void translate(double dx, double dy)
//    {
//        x += dx;
//        y += dy;
//    }
//
//    void scale(double scale)
//    {
//        x *= scale;
//        y *= scale;
//    }
//
//    double x;
//    double y;
//};


//std::ostream& operator<<(std::ostream& os, const Point& p)
//{
//    os << "(" << p.x << ", " << p.y << ")";
//    return os;
//}




//struct Line
//{
//    Line(Point start, Point end, double hue):
//        start(start),
//        end(end),
//        hue(hue)
//    {}
//
//
//    Point start;
//    Point end;
//    double hue;   // 0 = red, 0.333 = green, 0.667 = blue, 1 = red
//};

//std::ostream& operator<<(std::ostream& os, const Line& line)
//{
//    os << "Line " << line.start << " to " << line.end << " hue " << line.hue;
//    return os;
//}


struct PixelFloat
{
    double r = 0;
    double g = 0;
    double b = 0;
    double a = 0;
};

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

class BufferFloat
{
public:
    BufferFloat(uint32_t width, uint32_t height) :
        m_width(width),
        m_height(height)
    {
        assert(width > 0);
        assert(height > 0);
        m_data = new PixelFloat[width * height];
    }

    ~BufferFloat()
    {
        delete[] m_data;
    }

    void set(uint32_t x, uint32_t y, PixelFloat data)
    {
        assert(x < m_width);
        assert(y < m_height);
        m_data[y * m_width + x] = data;
    }

    const PixelFloat& get(uint32_t x, uint32_t y) const
    {
        assert(x < m_width);
        assert(y < m_height);
        return m_data[y * m_width + x];
    }

    uint32_t getWidth() const
    {
        return m_width;
    }

    uint32_t getHeight() const
    {
        return m_height;
    }

private:
    const uint32_t m_width;
    const uint32_t m_height;

    PixelFloat* m_data;
};


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


int writeImage(char* filename, std::shared_ptr<BufferFloat> buffer, char* title)
{
    int code = 0;
    FILE *fp = NULL;
    png_structp png_ptr = NULL;
    png_infop info_ptr = NULL;
    png_bytep row = NULL;

    auto width = buffer->getWidth();
    auto height = buffer->getHeight();

    // Open file for writing (binary mode)
    fp = fopen(filename, "wb");
    if (fp == NULL) {
        fprintf(stderr, "Could not open file %s for writing\n", filename);
        code = 1;
        goto finalise;
    }

    // Initialize write structure
    png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    if (png_ptr == NULL) {
        fprintf(stderr, "Could not allocate write struct\n");
        code = 1;
        goto finalise;
    }

    // Initialize info structure
    info_ptr = png_create_info_struct(png_ptr);
    if (info_ptr == NULL) {
        fprintf(stderr, "Could not allocate info struct\n");
        code = 1;
        goto finalise;
    }

    // Setup Exception handling
    if (setjmp(png_jmpbuf(png_ptr))) {
        fprintf(stderr, "Error during png creation\n");
        code = 1;
        goto finalise;
    }

    png_init_io(png_ptr, fp);

    // Write header (8 bit colour depth)
    png_set_IHDR(png_ptr, info_ptr, width, height,
        8, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
        PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);

    // Set title
    if (title != NULL) {
        png_text title_text;
        title_text.compression = PNG_TEXT_COMPRESSION_NONE;
        title_text.key = "Title";
        title_text.text = title;
        png_set_text(png_ptr, info_ptr, &title_text, 1);
    }

    png_write_info(png_ptr, info_ptr);

    // Allocate memory for one row (3 bytes per pixel - RGB)
    row = (png_bytep)malloc(3 * width * sizeof(png_byte));


    ////////////////////////////////////////////////////////////////

    // Draw a diagonal line
    for (int y = 0; y < height; ++y)
    {
        for (int x = 0; x < width; ++x)
        {
            auto outPixel = &(row[x * 3]);

            auto inPixel = buffer->get(x, y);

            outPixel[0] = (png_byte)(inPixel.r * 255);
            outPixel[1] = (png_byte)(inPixel.g * 255);
            outPixel[2] = (png_byte)(inPixel.b * 255);

        }

        png_write_row(png_ptr, row);
    }


    ////////////////////////////////////////////////////////////////


    // End write
    png_write_end(png_ptr, NULL);

finalise:
    if (fp != NULL) fclose(fp);
    if (info_ptr != NULL) png_free_data(png_ptr, info_ptr, PNG_FREE_ALL, -1);
    if (png_ptr != NULL) png_destroy_write_struct(&png_ptr, (png_infopp)NULL);
    if (row != NULL) free(row);

    return code;
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
    int result = writeImage(argv[1], move(buffer), "This is my test image");

    // Free up the memorty used to store the image
    //free(buffer);

    return result;
}
