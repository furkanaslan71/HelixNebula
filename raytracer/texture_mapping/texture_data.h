#ifndef TEXTURE_DATA_H
#define TEXTURE_DATA_H

#include <string>
#include "parser/parser.hpp"

#define SET(x) x = tm.x


template <typename EnumT>
inline EnumT getEnum(const std::vector<std::string>& types, const std::string& target)
{
    for (size_t i = 0; i < types.size(); i++)
    {
        if (types[i] == target)
            return static_cast<EnumT>(i);
    }
    return static_cast<EnumT>(-1); // or handle error differently
}

enum class TextureType : int {
	image = 0,
	perlin = 1,
	checkerboard = 2	
};

enum class Interpolation : int {
	nearest,
	bilinear,
	trilinear
};

enum class DecalMode : int {
	replace_kd,
	blend_kd,
	replace_ks,
	replace_background,
	replace_normal,
	bump_normal,
	replace_all
};

struct Image {
	int id;
	std::string file_path;
	Image() = default;
	Image(Image_ img)
	{
		id = img.id;
		file_path = img.data;
	}
};

struct TextureMap {
	int id;
	int image_id;
	TextureType type;
	Interpolation interp;
	DecalMode d_mode;
	float bump_factor;
	float noise_scale;
	std::string noise_conversion;
	int num_octaves;
	float scale;
	float offset;
	glm::vec3 black_color;
	glm::vec3 white_color;
	float normalizer;
	TextureMap() = default;
	TextureMap(TextureMap_ tm)
	{
		id = tm.id;
		std::vector<std::string> texture_types = { "image", "perlin", "checkerboard" };
		type = getEnum<TextureType>(texture_types, tm.type);
		std::vector<std::string> interpolations = { "nearest", "bilinear", "trilinear" };
		interp = getEnum<Interpolation>(interpolations, tm.interpolation);
		std::vector<std::string> decal_modes = { "replace_kd", "blend_kd",
			"replace_ks", "replace_background", "replace_normal", "bump_normal", "replace_all" };
		d_mode = getEnum<DecalMode>(decal_modes, tm.decal_mode);
		SET(image_id);
		SET(bump_factor);
		SET(noise_scale);
		SET(noise_conversion);
		SET(num_octaves);
		SET(scale);
		SET(offset);
		SET(black_color);
		SET(white_color);
		SET(normalizer);
	}
};


struct TextureData {
	std::unordered_map<int, Image> images;
	std::unordered_map<int, TextureMap> textures;
};



#endif // !TEXTURE_DATA_H