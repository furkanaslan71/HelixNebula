//
// Created by furkan on 12/29/25.
//

#include "image_io.h"

#define STB_IMAGE_IMPLEMENTATION
#include "io/stb_image.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "../io/stb_image_write.h"

#define TINYEXR_IMPLEMENTATION
#include "external/tinyexr-release/tinyexr.h"

unsigned char* readSDR(const std::string& file_path, int* width, int* height, int* channels)
{
    unsigned char* data = stbi_load(file_path.c_str(), width, height, channels, 0);
    if (!data)
    {
        throw std::runtime_error("Error loading SDR image!!!");
    }
    return data;
}

void readHDR(float*& out_img, const std::string& file_path, int* width, int* height, int* channels)
{
    const char* err = nullptr;
    int ret = LoadEXR(&out_img, width, height, file_path.c_str(), &err);
    if (ret != TINYEXR_SUCCESS) {
        if (err) {
            fprintf(stderr, "ERR : %s\n", err);
            FreeEXRErrorMessage(err); // release memory of error message.
        }
    }
    *channels = 4;
}

void saveImage(const std::string& outputDir,
                            const std::string& fileName,
                           std::vector<std::vector<Color>>& image, ImageType img_type)
{
    if (image.empty() || image[0].empty())
    {
        std::cerr << "Error: Image buffer is empty or has zero width/height."
          << std::endl;
        return;
    }

    int height = image.size();
    int width = image[0].size();
    int channels = 3;
    if (img_type == ImageType::HDR)
        channels = 3;

    std::filesystem::path outputPath =
      std::filesystem::path(outputDir) / fileName;
    std::string fullPath = outputPath.string();

    std::vector<unsigned char> data(width * height * channels);

    for (int y = 0; y < height; ++y)
    {
        for (int x = 0; x < width; ++x)
        {
            Color c = image[y][x];
            c.clamp();
            unsigned char r = static_cast<unsigned char>
                      (std::round(c.r));
            unsigned char g = static_cast<unsigned char>
                (std::round(c.g));
            unsigned char b = static_cast<unsigned char>
                (std::round(c.b));
            int index = (y * width + x) * channels;
            data[index + 0] = r;
            data[index + 1] = g;
            data[index + 2] = b;
        }
    }

    int success = 0;
    if (img_type == ImageType::SDR)
    {
        success = stbi_write_png(fullPath.c_str(),
      width, height, channels, data.data(), width * channels);
    }
    else
    {
        // write hdr image
        success = stbi_write_png(fullPath.c_str(),
      width, height, channels, data.data(), width * channels);
    }

    if (success)
    {
        std::cout << "Successfully saved image: " << fullPath << std::endl;
    }
    else
    {
        std::cerr << "Error: Could not save image: " << fullPath << std::endl;
    }
}

void saveEXR(const std::string& full_path, const std::vector<std::vector<Color>>& image)
{
    int height = static_cast<int>(image.size());
    int width = static_cast<int>(image[0].size());

    EXRHeader header;
    InitEXRHeader(&header);

    EXRImage exr_image;
    InitEXRImage(&exr_image);

    exr_image.num_channels = 3;
    exr_image.width = width;
    exr_image.height = height;

    // TinyEXR requires Planar format (separate R, G, B arrays)
    std::vector<float> r_channel(width * height);
    std::vector<float> g_channel(width * height);
    std::vector<float> b_channel(width * height);

    for (int y = 0; y < height; ++y)
    {
        for (int x = 0; x < width; ++x)
        {
            int idx = y * width + x;
            // Casting your Color (double) to float for EXR storage
            r_channel[idx] = static_cast<float>(image[y][x].r);
            g_channel[idx] = static_cast<float>(image[y][x].g);
            b_channel[idx] = static_cast<float>(image[y][x].b);
        }
    }

    // EXR expects BGR order for channel names often,
    // but the actual data pointers must match the header order
    float* image_ptr[3];
    image_ptr[0] = b_channel.data(); // B
    image_ptr[1] = g_channel.data(); // G
    image_ptr[2] = r_channel.data(); // R

    exr_image.images = (unsigned char**)image_ptr;

    header.num_channels = 3;
    header.channels = (EXRChannelInfo*)malloc(sizeof(EXRChannelInfo) * 3);
    strncpy(header.channels[0].name, "B", 255);
    strncpy(header.channels[1].name, "G", 255);
    strncpy(header.channels[2].name, "R", 255);

    header.pixel_types = (int*)malloc(sizeof(int) * 3);
    header.requested_pixel_types = (int*)malloc(sizeof(int) * 3);

    for (int i = 0; i < 3; i++)
    {
        header.pixel_types[i] = TINYEXR_PIXELTYPE_FLOAT;      // Input is 32-bit float
        header.requested_pixel_types[i] = TINYEXR_PIXELTYPE_HALF; // Save as 16-bit to save space
    }

    const char* err = nullptr;
    int ret = SaveEXRImageToFile(&exr_image, &header, full_path.c_str(), &err);
    if (ret != TINYEXR_SUCCESS)
    {
        std::cerr << "TinyEXR Save Error: " << err << std::endl;
        FreeEXRErrorMessage(err);
    }

    free(header.channels);
    free(header.pixel_types);
    free(header.requested_pixel_types);
}