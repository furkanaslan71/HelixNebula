//
// Created by furkan on 12/30/25.
//

#ifndef TEXTURE_LOOKUP_H
#define TEXTURE_LOOKUP_H

#include <core/hittable.h>
#include "texture_data.h"
#include "perlin.h"
#include "light/light.h"

struct CameraContext;

glm::vec3 fetch_single_sample(Image* image, glm::vec2 uv);

glm::vec3 fetch_interpolated_sample(Image* image, glm::vec2 uv, Interpolation interp);

glm::vec3 lookupImageTexture(Texture* texture, glm::vec2 uv);

glm::vec3 lookupTexture(Texture* texture, glm::vec2 uv, const glm::vec3& , bool replace_all = false);

glm::vec3 lookupNormalMap(Texture* texture, const HitRecord& rec);

glm::vec3 lookupBumpMap(Texture* texture, const HitRecord& rec);

glm::vec3 lookupBackgroundTex(Texture* texture, const glm::vec3& dir, const CameraContext& cam_context);

glm::vec3 lookupEnvMap(const EnvironmentLight& env_map, const glm::vec3& dir);


#endif //TEXTURE_LOOKUP_H
