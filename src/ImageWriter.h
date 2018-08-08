#pragma once

#include <stdint.h>
#include <string>

struct PixelFloat
{
    double r = 0;
    double g = 0;
    double b = 0;
    double a = 0;
};


class BufferFloat
{
public:
    BufferFloat(uint32_t width, uint32_t height);
    ~BufferFloat();

    void set(uint32_t x, uint32_t y, PixelFloat data);
    const PixelFloat& get(uint32_t x, uint32_t y) const;
    uint32_t getWidth() const;
    uint32_t getHeight() const;

private:
    const uint32_t m_width;
    const uint32_t m_height;

    PixelFloat* m_data;
};


int writeImage(std::string filename, const BufferFloat& buffer);