#include "texture_fetcher.h"

#define STB_IMAGE_IMPLEMENTATION
#include "io/stb_image.h"

#define TEX_INFO textureData, texWidth, texChannels


TextureFetcher::TextureFetcher(TextureData& _data)
	: data_(_data) 
{
	for (const auto& [id, img] : data_.images)
	{
		int width, height, channels;
		unsigned char* data;

		data = loadTexture(img.file_path, width, height, channels);
		image_datas_[id] = ImageData{ data, width, height, channels };
	}
	initPerlin();
}

TextureFetcher::~TextureFetcher()
{
	for (auto& [id, img] : image_datas_)
	{
		stbi_image_free(img.data);
	}
}

unsigned char* TextureFetcher::loadTexture(const std::string& filepath, int& width, int& height, int& channels)
{
	std::string absolute_path = "D:/Furkan/GITHUB/HelixNebula/inputs2/" + filepath;
	unsigned char* data = stbi_load(absolute_path.c_str(), &width, &height, &channels, 0);
	if (!data)
	{
		std::cerr << "Error loading texture!" << std::endl;
		return nullptr;
	}
	return data;
}

inline glm::vec3 fetch(unsigned char* textureData, int texWidth, int texChannels, int x, int y)
{
	int index = (y * texWidth + x) * texChannels;
	glm::vec3 result = { textureData[index + 0],
		textureData[index + 1],
		textureData[index + 2] };
	return result;
}

glm::vec3 TextureFetcher::getTextureValue(unsigned char* textureData,
																					int texWidth,
																					int texHeight,
																					int texChannels,
																					float u,
																					float v,
																					Interpolation interpolation) const
{
	u = std::min(std::max(u, 0.0f), 1.0f);
	v = std::min(std::max(v, 0.0f), 1.0f);

	// UV -> pixel
	float i = u * (texWidth - 1.0f);
	float j = v * (texHeight - 1.0f);
	glm::vec3 color{};
	if (interpolation == Interpolation::nearest)
	{
		int x = round(i);
		int y = round(j);
		color = fetch(TEX_INFO, x, y);
	}
	else if (interpolation == Interpolation::bilinear)
	{
		int p = floor(i);
		int q = floor(j);
		float dx = i - p;
		float dy = j - q;
		color = fetch(TEX_INFO, p, q) * (1.f - dx) * (1.f - dy) +
			fetch(TEX_INFO, p + 1, q) * (dx) * (1.f - dy) +
			fetch(TEX_INFO, p, q + 1) * (1.f - dx) * (dy) +
			fetch(TEX_INFO, p + 1, q + 1) * (dx) * (dy);
	}
	return color;

	
}

LookupInfo TextureFetcher::get_lookup_info(glm::vec2 tex_coord,
																					 int texture_id,
																					 const glm::vec3& point) const
{
	const TextureMap& tex = data_.textures[texture_id];
	const DecalMode d_mode = tex.d_mode;
	glm::vec3 tex_val(0.0f);
	const TextureType type = tex.type;
	switch (type)
	{
		case TextureType::image:
		{
			const ImageData& img = image_datas_.at(tex.image_id);

			tex_val = getTextureValue(
				img.data,
				img.width,
				img.height,
				img.channels,
				tex_coord.x,
				tex_coord.y,
				tex.interp
			);
			break;
		}
		case TextureType::perlin:
			tex_val = glm::vec3(perlinNoise(point, tex.noise_scale, tex.noise_conversion, tex.num_octaves));
			return LookupInfo(d_mode, tex_val);
			break;
		case TextureType::checkerboard:
			break;
		default:
			break;
	}
	tex_val = glm::clamp(tex_val / tex.normalizer, 0.0f, 1.0f);
	return LookupInfo(d_mode, tex_val);
}