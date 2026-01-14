//
// Created by furkan on 1/14/26.
//

#ifndef RENDERCONTEXT_H
#define RENDERCONTEXT_H

#include "core/color.h"
#include "texture_mapping/texture_data.h"


enum class BackgroundType : bool {
    Color,
    Texture
};

struct RenderContext {
    union {
        Color background_color;
        Texture* background_tex;
    }background_info;
    float shadow_ray_epsilon;
    float intersection_test_epsilon;
    int max_recursion_depth;
    BackgroundType b_type;
    Image* env_map;

    RenderContext() = default;
};

#endif //RENDERCONTEXT_H
