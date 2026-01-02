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

struct CameraContext {
    glm::vec3 forward;
    glm::vec3 right;
    glm::vec3 up;
    float tan_half_fov_x;
    float tan_half_fov_y;
};


inline std::string getFileExtension(const std::string& filename)
{
    size_t dotPosition = filename.rfind('.');

    if (dotPosition != std::string::npos && dotPosition != 0)
        return filename.substr(dotPosition + 1);

    return "";
}