//
// Created by furkan on 1/2/26.
//

#include "tonemapping.h"

#define LUM(v) 0.2126 * v.r + 0.7152 * v.g + 0.0722 * v.b

bool degamma = false;

Tonemap::Tonemap(const Tonemap_ &tm)
    : TMOOptions(tm.TMOOptions), saturation(tm.saturation), gamma(tm.gamma), extension(tm.extension)
{
    if (tm.TMO == "Photographic")
    {
        tmo = TMO::Photographic;
    }
    else if (tm.TMO == "Filmic")
    {
        tmo = TMO::Filmic;
    }
    else if (tm.TMO == "ACES")
    {
        tmo = TMO::ACES;
    }
    else
    {
        throw std::runtime_error("Unsupported TMO type");
    }
}

float logAvgLuminance(const std::vector<std::vector<Color>>& image)
{
    Expects(!image.empty());
    constexpr float eps = 1e-6f;
    int N = image.size() * image[0].size();
    float power = 0.0f;
    for (const auto& v : image)
    {
        for (const auto& c : v)
        {
            power += glm::log(eps + LUM(c));
        }
    }
    return glm::exp(power / N);
}

float scaledLuminance(float inv_log_avg, float key, const Color& color, float Yi)
{
    return key * inv_log_avg * Yi;
}

float gammaCorrect(float inv_g, float color_channel)
{
    if (degamma)
        return glm::clamp(color_channel, 0.0f, 1.0f) * 255.0f;
    return glm::clamp(glm::pow(color_channel, inv_g), 0.0f, 1.0f) * 255.0f;
}

void reinhard(const std::vector<std::vector<Color>>& input_img,
    std::vector<std::vector<Color>>& output_img, const Tonemap& tonemap)
{
    float log_avg_lum = logAvgLuminance(input_img);
    float inv_log_avg_lum = 1.0f / log_avg_lum;
    float inv_g = 1.0f / tonemap.gamma;

    float inv_Lwhite_2 = 0.0f;
    if (tonemap.TMOOptions.y > 0)
    {
        std::vector<float> luminances;
        for (const auto& row : input_img)
            for (const auto& col : row)
                luminances.push_back(LUM(col));

        std::sort(luminances.begin(), luminances.end());
        size_t index = std::round(luminances.size() * ((100.0f - tonemap.TMOOptions.y) / 100.0f));
        index = std::clamp(index, (size_t)0, luminances.size() - 1);

        float L_white_raw = luminances[index];
        float L_white_scaled = L_white_raw * tonemap.TMOOptions.x * inv_log_avg_lum;
        inv_Lwhite_2 = 1.0f / (L_white_scaled * L_white_scaled);
    }

    for (size_t i = 0; i < input_img.size(); i++)
    {
        for (size_t j = 0; j < input_img[i].size(); j++)
        {
            float Yi = LUM(input_img[i][j]);
            float L_scaled = scaledLuminance(inv_log_avg_lum, tonemap.TMOOptions.x, input_img[i][j], Yi);

            float Yo;
            if (tonemap.TMOOptions.y > 0)
            {
                Yo = (L_scaled * (1.0f + L_scaled * inv_Lwhite_2)) / (1.0f + L_scaled);
            }
            else
            {
                Yo = L_scaled / (1.0f + L_scaled);
            }

            // Safe division for saturation logic
            float Yi_safe = std::max(Yi, 1e-6f);
            float ratio = Yo / Yi_safe;

            // Apply saturation and gamma
            output_img[i][j].r = gammaCorrect(inv_g, glm::pow(input_img[i][j].r * ratio, tonemap.saturation));
            output_img[i][j].g = gammaCorrect(inv_g, glm::pow(input_img[i][j].g * ratio, tonemap.saturation));
            output_img[i][j].b = gammaCorrect(inv_g, glm::pow(input_img[i][j].b * ratio, tonemap.saturation));
        }
    }
}


void filmic(const std::vector<std::vector<Color>>& input_img,
    std::vector<std::vector<Color>>& output_img, const Tonemap& tonemap)
{
    float log_avg_lum = logAvgLuminance(input_img);
    float inv_log_avg_lum = 1.0f / log_avg_lum;

    // Helper for Filmic curve
    auto map_filmic = [](float L) {
        float a = 0.22f;
        float b = 0.30f;
        float c = 0.10f;
        float d = 0.20f;
        float e = 0.01f;
        float f = 0.30f;
        return ((L * (a * L + c * b) + d * e) / (L * (a * L + b) + d * f)) - (e / f);
    };

    // Determine White Point scaling
    std::vector<float> luminances;
    for (const auto& row : input_img)
        for (const auto& col : row)
            luminances.push_back(LUM(col));

    std::sort(luminances.begin(), luminances.end());
    size_t index = std::round(luminances.size() * ((100.0f - tonemap.TMOOptions.y) / 100.0f));
    float L_white_raw = luminances[std::min(index, luminances.size() - 1)];

    float L_white_scaled = L_white_raw * tonemap.TMOOptions.x * inv_log_avg_lum;
    float map_W = map_filmic(L_white_scaled);

    float inv_g = 1.0f / tonemap.gamma;

    for (size_t i = 0; i < input_img.size(); i++)
    {
        for (size_t j = 0; j < input_img[i].size(); j++)
        {
            float Yi = LUM(input_img[i][j]);
            float L_scaled = scaledLuminance(inv_log_avg_lum, tonemap.TMOOptions.x, input_img[i][j], Yi);

            float Yo = map_filmic(L_scaled) / map_W;

            Yi = std::max(Yi, 1e-6f);
            float Ro = Yo * glm::pow((input_img[i][j].r / Yi), tonemap.saturation);
            float Go = Yo * glm::pow((input_img[i][j].g / Yi), tonemap.saturation);
            float Bo = Yo * glm::pow((input_img[i][j].b / Yi), tonemap.saturation);

            output_img[i][j].r = gammaCorrect(inv_g, Ro);
            output_img[i][j].g = gammaCorrect(inv_g, Go);
            output_img[i][j].b = gammaCorrect(inv_g, Bo);
        }
    }
}

void aces(const std::vector<std::vector<Color>>& input_img,
    std::vector<std::vector<Color>>& output_img, const Tonemap& tonemap)
{
    float log_avg_lum = logAvgLuminance(input_img);
    float inv_log_avg_lum = 1.0f / log_avg_lum;

    // Helper for ACES curve
    auto map_aces = [](float L) {
        float A = 2.51f;
        float B = 0.03f;
        float C = 2.43f;
        float D = 0.59f;
        float E = 0.14f;
        return (L * (A * L + B)) / (L * (C * L + D) + E);
    };

    // Determine White Point scaling
    std::vector<float> luminances;
    for (const auto& row : input_img)
        for (const auto& col : row)
            luminances.push_back(LUM(col));

    std::sort(luminances.begin(), luminances.end());
    size_t index = std::round(luminances.size() * ((100.0f - tonemap.TMOOptions.y) / 100.0f));
    float L_white_raw = luminances[std::min(index, luminances.size() - 1)];

    // Scale white point by the key
    float L_white_scaled = L_white_raw * tonemap.TMOOptions.x * inv_log_avg_lum;
    float map_W = map_aces(L_white_scaled);

    float inv_g = 1.0f / tonemap.gamma;

    for (size_t i = 0; i < input_img.size(); i++)
    {
        for (size_t j = 0; j < input_img[i].size(); j++)
        {
            float Yi = LUM(input_img[i][j]);
            float L_scaled = scaledLuminance(inv_log_avg_lum, tonemap.TMOOptions.x, input_img[i][j], Yi);

            // Apply ACES and normalize by white point
            float Yo = map_aces(L_scaled) / map_W;

            Yi = std::max(Yi, 1e-6f);
            float Ro = Yo * glm::pow((input_img[i][j].r / Yi), tonemap.saturation);
            float Go = Yo * glm::pow((input_img[i][j].g / Yi), tonemap.saturation);
            float Bo = Yo * glm::pow((input_img[i][j].b / Yi), tonemap.saturation);

            output_img[i][j].r = gammaCorrect(inv_g, Ro);
            output_img[i][j].g = gammaCorrect(inv_g, Go);
            output_img[i][j].b = gammaCorrect(inv_g, Bo);
        }
    }
}

void tonemap(const std::vector<std::vector<Color>>& input_img,
    std::vector<std::vector<Color>>& output_img, const Tonemap& tonemap)
{
    if (tonemap.tmo == TMO::Photographic)
    {
        reinhard(input_img, output_img, tonemap);
    }
    else if (tonemap.tmo == TMO::Filmic)
    {
        filmic(input_img, output_img, tonemap);
    }
    else if (tonemap.tmo == TMO::ACES)
    {
        aces(input_img, output_img, tonemap);
    }
    else
    {
        throw std::runtime_error("Unsupported TMO");
    }
}




