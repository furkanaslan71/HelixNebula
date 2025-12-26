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

float TextureFetcher::height_helper(const ImageData& img, glm::vec2 tex_coord, Interpolation interp) const
{
	auto tex_val = getTextureValue(
		img.data,
		img.width,
		img.height,
		img.channels,
		tex_coord.x,
		tex_coord.y,
		interp
	);
	auto res = (tex_val.x + tex_val.y + tex_val.z) / 3.0f;
	return res / 255.0f;
}

float TextureFetcher::height_function(glm::vec2 tex_coord, int texture_id, FetchMode mode) const
{
	const TextureMap& tex = data_.textures[texture_id];
	const TextureType type = tex.type;
	Interpolation interp = tex.interp;
	float res = 0;
	glm::vec3 tex_val;
	switch (type)
	{
		case TextureType::image:
		{
			const ImageData& img = image_datas_.at(tex.image_id);
			float h_value = height_helper(img, tex_coord, interp);
			if (mode == FetchMode::value_u_v)
			{
				res = h_value;
			}
			else if (mode == FetchMode::derivative_u)
			{
				float delta_u = 1.0f / img.width;
				tex_coord.x += delta_u;
				//tex_coord.x = std::min(tex_coord.x, 1.0f);
				float h_udeltau = height_helper(img, tex_coord, interp); 
				res = (h_udeltau - h_value);
			}
			else if (mode == FetchMode::derivative_v)
			{
				float delta_v = 1.0f / img.height;
				tex_coord.y += delta_v;
				//tex_coord.y = std::min(tex_coord.y, 1.0f);
				float h_vdeltav = height_helper(img, tex_coord, interp);
				res = (h_vdeltav - h_value);
			}			
			break;
		}
		default:
			break;
	}
	res *= tex.bump_factor;
	return res;
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
	u = u - std::floor(u);
	v = v - std::floor(v);

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

glm::vec3 TextureFetcher::get_background_texture(glm::vec3 dir, int texture_id, BackgroundCameraData cam) const
{
	const TextureMap& tex = data_.textures[texture_id];
	const ImageData& img = image_datas_.at(tex.image_id);

	float cos_theta = glm::dot(dir, cam.cam_forward);
	if (cos_theta <= 0.0001f) return glm::vec3(0, 0, 0);

	float x = glm::dot(dir, cam.cam_right) / cos_theta;
	float y = glm::dot(dir, cam.cam_up) / cos_theta;

	float u = (x / (2.0f * cam.tan_half_fov_x)) + 0.5f;
	float v = (y / (2.0f * cam.tan_half_fov_y)) + 0.5f;

	v = 1.0f - v;

	u = glm::clamp(u, 0.0f, 1.0f);
	v = glm::clamp(v, 0.0f, 1.0f);

	return getTextureValue(img.data, img.width, img.height, img.channels, u, v, tex.interp);
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
		{
			float scale = tex.scale;
			float offset = tex.offset;

			int x = static_cast<int>(std::floor((point.x + offset) * scale));
			int y = static_cast<int>(std::floor((point.y + offset) * scale));
			int z = static_cast<int>(std::floor((point.z + offset) * scale));

			bool is_x_even = (x % 2) == 0;
			bool is_y_even = (y % 2) == 0;
			bool is_z_even = (z % 2) == 0;

			if ((is_x_even ^ is_y_even) ^ is_z_even)
			{
				tex_val = tex.black_color;
			}
			else
			{
				tex_val = tex.white_color;
			}
			break;
		}
			
		default:
			break;
	}
	tex_val = glm::clamp(tex_val / tex.normalizer, 0.0f, 1.0f);
	return LookupInfo(d_mode, tex_val);
}

glm::vec3 TextureFetcher::getPerlinGradient(glm::vec3 point, int texture_id) const
{
	const TextureMap& tex = data_.textures[texture_id];
	float eps = 1e-3f;
	float x = point.x;
	float y = point.y;
	float z = point.z;
	glm::vec3 x_eps = { eps, 0.0f, 0.0f };
	glm::vec3 y_eps = { 0.0f, eps, 0.0f };
	glm::vec3 z_eps = { 0.0f, 0.0f , eps};
	float h = perlinNoise(point, tex.noise_scale, tex.noise_conversion, tex.num_octaves);
	float h_x = perlinNoise(point + x_eps, tex.noise_scale, tex.noise_conversion, tex.num_octaves);
	float h_y = perlinNoise(point + y_eps, tex.noise_scale, tex.noise_conversion, tex.num_octaves);
	float h_z = perlinNoise(point + z_eps, tex.noise_scale, tex.noise_conversion, tex.num_octaves);
	//eps = 1;
	return glm::vec3(
		(h_x - h) / eps,
		(h_y - h) / eps,
		(h_z - h) / eps
	);
}

TextureType TextureFetcher::getTexType(int texture_id) const
{
	const TextureMap& tex = data_.textures[texture_id];
	const TextureType type = tex.type;
	return type;
}

