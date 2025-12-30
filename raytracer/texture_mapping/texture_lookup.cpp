//
// Created by furkan on 12/30/25.
//

#include "texture_lookup.h"


glm::vec3 fetch_single_sample(Image* image, int x, int y)
{
    int index = (y * image->width + x) * image->channels;
    glm::vec3 result = { image->data[index + 0],
        image->data[index + 1],
        image->data[index + 2] };
    return result;
}

glm::vec3 fetch_interpolated_sample(Image* image, glm::vec2 uv, Interpolation interp)
{
    float u = uv.x - std::floor(uv.x);
    float v = uv.y - std::floor(uv.y);

    // UV -> pixel
    float i = u * (image->width - 1);
    float j = v * (image->height - 1);
    glm::vec3 color{};
    if (interp == Interpolation::nearest)
    {
        int x = round(i);
        int y = round(j);
        color = fetch_single_sample(image, x, y);
    }
    else if (interp == Interpolation::bilinear)
    {
        int p = floor(i);
        int q = floor(j);
        float dx = i - p;
        float dy = j - q;
        color = fetch_single_sample(image, p, q) * (1.f - dx) * (1.f - dy) +
            fetch_single_sample(image, p + 1, q) * (dx) * (1.f - dy) +
            fetch_single_sample(image, p, q + 1) * (1.f - dx) * (dy) +
            fetch_single_sample(image, p + 1, q + 1) * (dx) * (dy);
    }
    return color;
}

glm::vec3 lookupImageTexture(Texture* texture, glm::vec2 uv)
{
    // returns 0 - 255
    glm::vec3 color = fetch_interpolated_sample(texture->texture_data.image.image, uv, texture->interp);
    return color;
}

glm::vec3 lookupTexture(Texture* texture, glm::vec2 uv)
{
    switch (texture->type)
    {
        case(TextureType::image):
        {
            // returns 0 - 1
            glm::vec3 sample = lookupImageTexture(texture, uv);
            sample = glm::clamp(sample / texture->texture_data.image.normalizer, 0.0f, 1.0f);
            return sample;
        }
        case (TextureType::perlin):
        {
            return glm::vec3(0.0f);
            break;
        }
        default:
        {
            throw std::runtime_error("Unknown texture type");
        }
    }
}

glm::vec3 lookupNormalMap(Texture* texture, const HitRecord& rec)
{
    glm::vec3 rgb_sample = lookupTexture(texture, rec.uv) * 255.0f;
    glm::vec3 tangent_normal = rgb_sample / 127.5f - glm::vec3(1.0f);
    tangent_normal = glm::normalize(tangent_normal);

    glm::vec3 T, B;

    glm::vec3 res = rec.surface_tangents.u * tangent_normal.x +
        rec.surface_tangents.v * tangent_normal.y +
        rec.normal * tangent_normal.z;
    res = glm::normalize(res);
    return res;
}

glm::vec3 lookupBumpMap(Texture* texture, const HitRecord& rec)
{
    return glm::vec3(0.0f);
}

