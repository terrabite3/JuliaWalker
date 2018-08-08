#pragma once

#include <stdint.h>
#include <string>
#include <cassert>

struct Pixel
{
    uint8_t r = 0;
    uint8_t g = 0;
    uint8_t b = 0;

    Pixel() :
        r(0), g(0), b(0) {}
    explicit Pixel(uint8_t r, uint8_t g, uint8_t b) :
        r(r), g(g), b(b) {}

    static Pixel fromNormalizedFloat(double r, double g, double b)
    {
        assert(r >= 0.0);
        assert(r <= 1.0);
        assert(g >= 0.0);
        assert(g <= 1.0);
        assert(b >= 0.0);
        assert(b <= 1.0);
        return Pixel {
            (uint8_t)(r * 255),
            (uint8_t)(g * 255),
            (uint8_t)(b * 255),
        };
    }
};


class BufferFloat
{
public:
    BufferFloat(uint32_t width, uint32_t height);
    ~BufferFloat();

    void set(uint32_t x, uint32_t y, Pixel data);
    const Pixel& get(uint32_t x, uint32_t y) const;
    uint32_t getWidth() const;
    uint32_t getHeight() const;

private:
    const uint32_t m_width;
    const uint32_t m_height;

    Pixel* m_data;
};


int writeImage(std::string filename, const BufferFloat& buffer);