//
// Created by furkan on 1/2/26.
//

#ifndef TONEMAPPING_H
#define TONEMAPPING_H

#include "parser/parser.hpp"
#include "core/color.h"


enum class TMO {
    Photographic,
    Filmic,
    ACES
};

struct Tonemap {
    Tonemap() = default;
    Tonemap(const Tonemap_& tm);

    TMO tmo;
    glm::vec2 TMOOptions;
    float saturation;
    float gamma;
    std::string extension;
};

//float logAvgLuminance(const std::vector<std::vector<Color>>& image);

//float scaledLuminance(float inv_log_avg, float key, const Color& color);

void tonemap(const std::vector<std::vector<Color>>& input_img,
    std::vector<std::vector<Color>>& output_img, const Tonemap& tonemap);



#endif //TONEMAPPING_H
