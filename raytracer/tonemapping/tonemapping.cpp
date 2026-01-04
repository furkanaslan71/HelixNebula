//
// Created by furkan on 1/2/26.
//

#include "tonemapping.h"

#define LUM(v) 0.2126 * v.r + 0.7152 * v.g + 0.0722 * v.b

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
    constexpr float eps = 1e-3f;
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

float gammaCorrect(float inv_g, float color_channel, bool degamma = false)
{
    if (degamma)
        return color_channel * 255.0f;
    return glm::clamp(glm::pow(color_channel, inv_g), 0.0f, 1.0f) * 255.0f;
}

void reinhard(const std::vector<std::vector<Color>>& input_img,
    std::vector<std::vector<Color>>& output_img, const Tonemap& tonemap)
{
    float log_avg_lum = logAvgLuminance(input_img);
    float inv_log_avg_lum = 1.0f / log_avg_lum;
    if (tonemap.TMOOptions.y == 0)
    {
        for (size_t i = 0; i < input_img.size(); i++)
        {
            for (size_t j = 0; j < input_img[i].size(); j++)
            {
                float Yi = LUM(input_img[i][j]);
                float scaled_lum = scaledLuminance(inv_log_avg_lum, tonemap.TMOOptions.x, input_img[i][j], Yi);
                float Yo = scaled_lum / (scaled_lum + 1.0f);
                Yi = std::max(Yi, 1e-6f);
                float Ro = Yo * glm::pow((input_img[i][j].r / Yi), tonemap.saturation);
                float Go = Yo * glm::pow((input_img[i][j].g / Yi), tonemap.saturation);
                float Bo = Yo * glm::pow((input_img[i][j].b / Yi), tonemap.saturation);
                float inv_g = 1.0f / tonemap.gamma;
                output_img[i][j].r = gammaCorrect(inv_g, Ro);
                output_img[i][j].g = gammaCorrect(inv_g, Go);
                output_img[i][j].b = gammaCorrect(inv_g, Bo);
            }
        }
    }
    else
    {
        std::vector<float> luminances;
        for (size_t i = 0; i < input_img.size(); i++)
        {
            for (size_t j = 0; j < input_img[i].size(); j++)
            {
                luminances.push_back(LUM(input_img[i][j]));
            }
        }
        std::sort(luminances.begin(), luminances.end());

        size_t index = std::round(luminances.size() * ((100.0f - tonemap.TMOOptions.y) / 100.0f));
        Ensures(index < luminances.size());
        float L_white = luminances[index];
        float inv_Lwhite_2 = 1.0f / (L_white * L_white);

        for (size_t i = 0; i < input_img.size(); i++)
        {
            for (size_t j = 0; j < input_img[i].size(); j++)
            {
                float Yi = LUM(input_img[i][j]);
                float scaled_lum = scaledLuminance(inv_log_avg_lum, tonemap.TMOOptions.x, input_img[i][j], Yi);
                float Yo = (scaled_lum * (1 + scaled_lum * inv_Lwhite_2))/ (scaled_lum + 1.0f);
                float Ro = Yo * glm::pow((input_img[i][j].r / Yi), tonemap.saturation);
                float Go = Yo * glm::pow((input_img[i][j].g / Yi), tonemap.saturation);
                float Bo = Yo * glm::pow((input_img[i][j].b / Yi), tonemap.saturation);
                float inv_g = 1.0f / tonemap.gamma;
                output_img[i][j].r = gammaCorrect(inv_g, Ro);
                output_img[i][j].g = gammaCorrect(inv_g, Go);
                output_img[i][j].b = gammaCorrect(inv_g, Bo);
            }
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
        reinhard(input_img, output_img, tonemap);
    }
    else if (tonemap.tmo == TMO::ACES)
    {
        reinhard(input_img, output_img, tonemap);
    }
    else
    {
        throw std::runtime_error("Unsupported TMO");
    }
}




