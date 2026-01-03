#include <filesystem>

#include "texture_data.h"

Image::Image(const std::string& absolute_path, ImageType _type)
    : type(_type)
{
    if (_type == ImageType::SDR)
    {
        data.sdr = readSDR(absolute_path, &width, &height, &channels);
    }
    else if (_type == ImageType::HDR)
    {
        //data.hdr = nullptr;
        readHDR(data.hdr, absolute_path, &width, &height, &channels);
    }
}


Texture::Texture(const TextureMap_& tm)
{
    if (tm.type == "image") {
        type = TextureType::image;
    }
    else if (tm.type == "perlin") {
        type = TextureType::perlin;
    }
    else if (tm.type == "checkerboard") {
        type = TextureType::checkerboard;
    }
    else {
        throw std::runtime_error("Unknown texture type: " + tm.type);
    }

    if (tm.decal_mode == "replace_kd") {
        d_mode = DecalMode::replace_kd;
    }
    else if (tm.decal_mode == "blend_kd") {
        d_mode = DecalMode::blend_kd;
    }
    else if (tm.decal_mode == "replace_ks") {
        d_mode = DecalMode::replace_ks;
    }
    else if (tm.decal_mode == "replace_background") {
        d_mode = DecalMode::replace_background;
    }
    else if (tm.decal_mode == "replace_normal") {
        d_mode = DecalMode::replace_normal;
    }
    else if (tm.decal_mode == "bump_normal") {
        d_mode = DecalMode::bump_normal;
    }
    else if (tm.decal_mode == "replace_all") {
        d_mode = DecalMode::replace_all;
    }
    else {
        throw std::runtime_error("Unknown decal mode: " + tm.decal_mode);
    }

    if (tm.interpolation == "nearest") {
        interp = Interpolation::nearest;
    }
    else if (tm.interpolation == "bilinear") {
        interp = Interpolation::bilinear;
    }
    else if (tm.interpolation == "trilinear") {
        //trilinear is not implemented, so we treat as bilinear
        interp = Interpolation::nearest;
    }
    else {
        throw std::runtime_error("Unknown interpolation: " + tm.interpolation);
    }

}

void Texture::make_image(Image* _image, float _normalizer)
{
    new (&texture_data.image) decltype(texture_data.image){
        _image,
        _normalizer
    };
}

void Texture::make_perlin(NoiseConversion _noise_conversion,
                          float _noise_scale,
                          int _num_octaves)
{
    new (&texture_data.perlin) decltype(texture_data.perlin){
        _noise_conversion,
        _noise_scale,
        _num_octaves
    };
}

void Texture::make_bump_image(Image* _image,
                              float _normalizer,
                              float _bump_factor)
{
    new (&texture_data.bump_image) decltype(texture_data.bump_image){
        _image,
        _normalizer,
        _bump_factor
    };
}

void Texture::make_bump_perlin(NoiseConversion _noise_conversion,
                               float _bump_factor,
                               float _noise_scale,
                               int _num_octaves)
{
    new (&texture_data.bump_perlin) decltype(texture_data.bump_perlin){
        _noise_conversion,
        _bump_factor,
        _noise_scale,
        _num_octaves
    };
}

void Texture::make_checkerboard(const glm::vec3& _black_color,
                                 const glm::vec3& _white_color,
                                 float _scale,
                                 float _offset)
{
    new (&texture_data.checkerboard) decltype(texture_data.checkerboard){
        _black_color,
        _white_color,
        _scale,
        _offset
    };
}


