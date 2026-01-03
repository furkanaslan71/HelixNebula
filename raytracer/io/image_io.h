//
// Created by furkan on 12/29/25.
//

#ifndef SAVE_IMAGE_H
#define SAVE_IMAGE_H

#include <vector>
#include "core/color.h"
#include <filesystem>
#include <iostream>

enum class ImageType : int {
    SDR,
    HDR
};

unsigned char* readSDR(const std::string& file_path, int* width, int* height, int* channels);

void readHDR(float*& out_img, const std::string& file_path, int* width, int* height, int* channels);

void saveImage(const std::string& outputDir, const std::string& fileName,
    std::vector<std::vector<Color>>& image, ImageType img_type);



#endif //SAVE_IMAGE_H
