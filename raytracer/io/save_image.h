//
// Created by furkan on 12/29/25.
//

#ifndef SAVE_IMAGE_H
#define SAVE_IMAGE_H

#include <vector>
#include "core/color.h"
#include <filesystem>
#include <iostream>

void saveImage(const std::string& outputDir,
                            const std::string& fileName,
                           std::vector<std::vector<Color>>& image);


#endif //SAVE_IMAGE_H
