#pragma once

#define DEBUG_PIXEL 0
#define MULTI_THREADING 1
#define INTERSECTION_STACK_SIZE 64
#define BACKFACE_CULLING 0
#define COMMAND_LINE_INPUT 0


#if DEBUG_PIXEL
#define WIDTH 770
#define HEIGHT 68
#endif


#define FS std::filesystem

#include <string>
#include "glm_config.h"

enum class RendererParams : int {
    ImportanceSampling,
    NEE,
    MIS_BALANCE,
    MIS_POWER,
    MIS_01,
    RussianRoulette
};

struct CameraContext {
    glm::vec3 forward;
    glm::vec3 right;
    glm::vec3 up;
    float tan_half_fov_x;
    float tan_half_fov_y;

    std::unordered_map<RendererParams, bool> options = {
        {RendererParams::ImportanceSampling, false},
        {RendererParams::NEE,                false},
        {RendererParams::MIS_BALANCE,        false},
        {RendererParams::MIS_POWER,          false},
        {RendererParams::MIS_01,             false},
        {RendererParams::RussianRoulette,    false},
    };
    int splitting_factor;
};


inline std::string getFileExtension(const std::string& filename)
{
    size_t dotPosition = filename.rfind('.');

    if (dotPosition != std::string::npos && dotPosition != 0)
        return filename.substr(dotPosition + 1);

    return "";
}