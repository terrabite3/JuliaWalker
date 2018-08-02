#include "ImageWriter.h"
#include <cassert>
#include <png.h>
#include <cstdlib>

BufferFloat::BufferFloat(uint32_t width, uint32_t height) :
    m_width(width),
    m_height(height)
{
    assert(width > 0);
    assert(height > 0);
    m_data = new PixelFloat[width * height];
}

BufferFloat::~BufferFloat()
{
    delete[] m_data;
}

void BufferFloat::set(uint32_t x, uint32_t y, PixelFloat data)
{
    assert(x < m_width);
    assert(y < m_height);
    m_data[y * m_width + x] = data;
}

const PixelFloat& BufferFloat::get(uint32_t x, uint32_t y) const
{
    assert(x < m_width);
    assert(y < m_height);
    return m_data[y * m_width + x];
}

uint32_t BufferFloat::getWidth() const
{
    return m_width;
}

uint32_t BufferFloat::getHeight() const
{
    return m_height;
}




int writeImage(char* filename, const BufferFloat& buffer, char* title)
{
    int code = 0;
    FILE *fp = NULL;
    png_structp png_ptr = NULL;
    png_infop info_ptr = NULL;
    png_bytep row = NULL;

    auto width = buffer.getWidth();
    auto height = buffer.getHeight();

    // Open file for writing (binary mode)
    fp = fopen(filename, "wb");
    if (fp == NULL)
    {
        fprintf(stderr, "Could not open file %s for writing\n", filename);
        code = 1;
        goto finalise;
    }

    // Initialize write structure
    png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    if (png_ptr == NULL)
    {
        fprintf(stderr, "Could not allocate write struct\n");
        code = 1;
        goto finalise;
    }

    // Initialize info structure
    info_ptr = png_create_info_struct(png_ptr);
    if (info_ptr == NULL)
    {
        fprintf(stderr, "Could not allocate info struct\n");
        code = 1;
        goto finalise;
    }

    // Setup Exception handling
    if (setjmp(png_jmpbuf(png_ptr)))
    {
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
    if (title != NULL)
    {
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

            auto inPixel = buffer.get(x, y);

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