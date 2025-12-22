#ifndef TEXTURE_FETCHER
#define TEXTURE_FETCHER

#include <vector>
#include <unordered_map>
#include <iostream>

#include "texture_data.h"
#include "glm_config.h"
#include "perlin.h"

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

private:
	unsigned char* loadTexture(const std::string& filepath, int& width, int& height, int& channels);
	glm::vec3 getTextureValue(unsigned char* textureData, int texWidth,
														int texHeight, int texChannels, float u,
														float v, Interpolation interpolation) const;

	TextureData& data_;
	std::unordered_map<int, ImageData> image_datas_;
};


#endif // !TEXTURE_FETCHER