#ifndef TEXTURE_FETCHER
#define TEXTURE_FETCHER

#include <vector>
#include <unordered_map>
#include <iostream>

#include "texture_data.h"
#include "glm_config.h"
#include "perlin.h"
#include "math_core/math_core.h"


enum class FetchMode : int {
	value_u_v,
	derivative_u,
	derivative_v
};

struct LookupInfo {
	DecalMode mode;
	glm::vec3 tex_val;
};

struct ImageData {
	unsigned char* data;
	int width;
	int height;
	int channels;
};

class TextureFetcher {
public:
	TextureFetcher(TextureData& _data);
	~TextureFetcher();
	LookupInfo get_lookup_info(glm::vec2 tex_coord, int texture_id, const glm::vec3& point ) const;
	float height_function(glm::vec2 tex_coord, int texture_id, FetchMode mode) const;
	glm::vec3 get_background_texture(glm::vec3 dir, int texture_id, BackgroundCameraData cam) const;

	TextureType getTexType(int texture_id) const;

	glm::vec3 getPerlinGradient(glm::vec3 point, int texture_id) const;
	TextureData& data_;

private:
	float height_helper(const ImageData& img, glm::vec2 tex_coord, Interpolation interp) const;
	unsigned char* loadTexture(const std::string& filepath, int& width, int& height, int& channels);
	glm::vec3 getTextureValue(unsigned char* textureData, int texWidth,
														int texHeight, int texChannels, float u,
														float v, Interpolation interpolation) const;



	
	std::unordered_map<int, ImageData> image_datas_;
};


#endif // !TEXTURE_FETCHER