//
// Created by furkan on 12/30/25.
//

#include "texture_lookup.h"


glm::vec3 fetch_single_sample(Image* image, int x, int y)
{
    //x = glm::clamp(x, 0, image->width  - 1);
    //y = glm::clamp(y, 0, image->height - 1);
    int index = (y * image->width + x) * image->channels;
    glm::vec3 result;
    if (image->type == ImageType::SDR)
    {
        result = { image->data.sdr[index + 0],
        image->data.sdr[index + 1],
        image->data.sdr[index + 2] };
    }
    else
    {
        result = { image->data.hdr[index + 0],
        image->data.hdr[index + 1],
        image->data.hdr[index + 2] };
    }
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
    glm::vec3 color = fetch_interpolated_sample(texture->texture_data.image.image, uv, texture->interp);
    if (texture->texture_data.image.image->type == ImageType::SDR)
        color = glm::clamp(color, 0.0f, 255.0f);
    return color;
}

glm::vec3 lookupTexture(Texture* texture, glm::vec2 uv, const glm::vec3& hit_point, bool replace_all)
{
    switch (texture->type)
    {
        case(TextureType::image):
        {
            // returns 0 - 1
            glm::vec3 sample = lookupImageTexture(texture, uv);
            if (texture->texture_data.image.image->type == ImageType::SDR)
            {
                sample = glm::clamp(sample / texture->texture_data.image.normalizer, 0.0f, 1.0f);
                if (replace_all)
                    sample = sample * 255.0f;
            }
            return sample;
        }
        case (TextureType::perlin):
        {
            auto sample = glm::vec3(perlinNoise(hit_point,
                texture->texture_data.perlin.noise_scale,
                texture->texture_data.perlin.noise_conversion,
                texture->texture_data.perlin.num_octaves));

            if (replace_all)
                sample *= 255.0f;
            return sample;
        }
        case (TextureType::checkerboard):
        {
            float scale = texture->texture_data.checkerboard.scale;
            float offset = texture->texture_data.checkerboard.offset;
            glm::vec3 black_color = texture->texture_data.checkerboard.black_color;
            glm::vec3 white_color = texture->texture_data.checkerboard.white_color;

            int x = static_cast<int>(std::floor((hit_point.x + offset) * scale));
            int y = static_cast<int>(std::floor((hit_point.y + offset) * scale));
            int z = static_cast<int>(std::floor((hit_point.z + offset) * scale));

            bool is_x_even = (x % 2) == 0;
            bool is_y_even = (y % 2) == 0;
            bool is_z_even = (z % 2) == 0;

            if ((is_x_even ^ is_y_even) ^ is_z_even)
                if (replace_all)
                    return black_color * 255.0f;
                else
                    return black_color;

            if (replace_all)
                return white_color * 255.0f;
            else
                return white_color;
        }
        default:
        {
            throw std::runtime_error("Unknown texture type");
        }
    }
}

glm::vec3 lookupNormalMap(Texture* texture, const HitRecord& rec)
{
    glm::vec3 rgb_sample = lookupTexture(texture, rec.uv, rec.point) * 255.0f;
    glm::vec3 tangent_normal = rgb_sample / 127.5f - glm::vec3(1.0f);
    tangent_normal = glm::normalize(tangent_normal);

    glm::vec3 T, B;

    glm::vec3 res = rec.surface_tangents.u * tangent_normal.x +
        rec.surface_tangents.v * tangent_normal.y +
        rec.normal * tangent_normal.z;
    res = glm::normalize(res);
    return res;
}

float height_function(Texture* texture, glm::vec2 uv)
{
    const glm::vec3 rgb_sample = lookupImageTexture(texture, uv);
    float grayscale = (rgb_sample.x + rgb_sample.y + rgb_sample.z) / 3.0f;
    grayscale /= texture->texture_data.bump_image.normalizer;
    return grayscale * texture->texture_data.bump_image.bump_factor;
}

glm::vec3 lookupBumpMap(Texture* texture, const HitRecord& rec)
{
    if (texture->type == TextureType::image)
    {
        float h = height_function(texture, rec.uv);
        float delta_u = 1.0f / texture->texture_data.bump_image.image->width;
        float delta_v = 1.0f / texture->texture_data.bump_image.image->height;
        float h_u = height_function(texture, rec.uv + glm::vec2(delta_u, 0.0f)) - h;
        float h_v = height_function(texture, rec.uv + glm::vec2(0.0f, delta_v)) - h;
        glm::vec3 grad = h_u * rec.surface_tangents.u + h_v * rec.surface_tangents.v;
        return glm::normalize(rec.normal - grad);
    }
    if (texture->type == TextureType::perlin)
    {
        float invEps = 1e3f;
        float eps = 1e-3f;

        float h = perlinNoise(rec.point,
            texture->texture_data.bump_perlin.noise_scale,
            texture->texture_data.bump_perlin.noise_conversion,
            texture->texture_data.bump_perlin.num_octaves);

        float h_x = perlinNoise(rec.point + glm::vec3(eps, 0.0f, 0.0f),
            texture->texture_data.bump_perlin.noise_scale,
            texture->texture_data.bump_perlin.noise_conversion,
            texture->texture_data.bump_perlin.num_octaves);

        float h_y = perlinNoise(rec.point + glm::vec3(0.0f, eps, 0.0f),
            texture->texture_data.bump_perlin.noise_scale,
            texture->texture_data.bump_perlin.noise_conversion,
            texture->texture_data.bump_perlin.num_octaves);

        float h_z = perlinNoise(rec.point + glm::vec3(0.0f, 0.0f, eps),
            texture->texture_data.bump_perlin.noise_scale,
            texture->texture_data.bump_perlin.noise_conversion,
            texture->texture_data.bump_perlin.num_octaves);

        auto g = invEps * (glm::vec3(h_x, h_y, h_z) - glm::vec3(h, h, h));

        glm::vec3 g_parallel = glm::dot(rec.normal, g) * rec.normal;
        glm::vec3 g_perp = g - g_parallel;
        return glm::normalize(rec.normal - (g_perp * texture->texture_data.bump_perlin.bump_factor));
    }
    throw std::runtime_error("Unsupported texture type for bump mapping");
}

glm::vec3 lookupBackgroundTex(Texture* texture, const glm::vec3& dir, const CameraContext& cam_context)
{
    float cos_theta = glm::dot(dir, cam_context.forward);
    if (cos_theta <= 0.0001f) return glm::vec3(0, 0, 0);

    float x = glm::dot(dir, cam_context.right) / cos_theta;
    float y = glm::dot(dir, cam_context.up) / cos_theta;

    float u = (x / (2.0f * cam_context.tan_half_fov_x)) + 0.5f;
    float v = (y / (2.0f * cam_context.tan_half_fov_y)) + 0.5f;

    v = 1.0f - v;

    u = glm::clamp(u, 0.0f, 1.0f);
    v = glm::clamp(v, 0.0f, 1.0f);

    return lookupImageTexture(texture, glm::vec2(u, v));
}

