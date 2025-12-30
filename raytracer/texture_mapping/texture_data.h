#ifndef TEXTURE_DATA_H
#define TEXTURE_DATA_H

#include <string>
#include "parser/parser.hpp"

#define SET(x) x = tm.x

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

enum class NoiseConversion : bool {
	linear,
	absval
};

struct Image {
	explicit Image(const std::string& absolute_path);
	Image() = default;

	unsigned char* data;
	int width;
	int height;
	int channels;
};

struct Texture {
	explicit Texture(const TextureMap_& tm);
	Texture() = default;
	void make_image(Image* _image, float _normalizer);
	void make_perlin(NoiseConversion _noise_conversion,
		float _noise_scale, int _num_octaves);
	void make_bump_image(Image* _image, float _normalizer, float _bump_factor);
	void make_bump_perlin(NoiseConversion _noise_conversion,
		float _bump_factor, float _noise_scale, int _num_octaves);
	void make_checkerboard(const glm::vec3& _black_color,
		const glm::vec3& _white_color, float _scale, float _offset);

	TextureType type;
	Interpolation interp;
	DecalMode d_mode;

	union {
		struct {
			Image* image;
			float normalizer;
		}image;
		struct {
			NoiseConversion noise_conversion;
			float noise_scale;
			int num_octaves;
		}perlin;
		struct {
			Image* image;
			float normalizer;
			float bump_factor;
		}bump_image;
		struct {
			NoiseConversion noise_conversion;
			float bump_factor;
			float noise_scale;
			int num_octaves;
		}bump_perlin;
		struct {
			glm::vec3 black_color;
			glm::vec3 white_color;
			float scale;
			float offset;
		}checkerboard;
	} texture_data;

	static_assert(std::is_trivially_destructible_v<decltype(texture_data)>);

};

#endif // !TEXTURE_DATA_H