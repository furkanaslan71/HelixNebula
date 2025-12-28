#include "render_manager.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "../io/stb_image_write.h"

#include <iostream>
#include <filesystem>
#include <cmath>
#include <algorithm>

RenderManager::RenderManager(const Scene& _scene,
  const MaterialManager& _material_manager,
  const RendererInfo _renderer_info,
  BaseRayTracer& _rendering_technique)
  : scene(_scene),
  material_manager(_material_manager),
  renderer_info(_renderer_info),
  technique(_rendering_technique)
{
}


void RenderManager::render() const
{
  std::string saveDir = FS::absolute(__FILE__).parent_path() / "../../outputs";
  for (const auto& cam : scene.pinhole_cameras)
  {
    std::vector<std::vector<Color>> image;
    cam.render(technique, image);
    saveImage(saveDir, cam.image_name, image);
  }
  for (const auto& cam : scene.distribution_cameras)
  {
    std::vector<std::vector<Color>> image;
    cam.render(technique, image);
    saveImage(saveDir, cam.image_name, image);
  }
}

void RenderManager::saveImage(const std::string& outputDir,
                            const std::string& fileName,
                           std::vector<std::vector<Color>>& image) const
{
  if (image.empty() || image[0].empty())
  {
    std::cerr << "Error: Image buffer is empty or has zero width/height." 
      << std::endl;
    return;
  }

  int height = image.size();
  int width = image[0].size();
  int channels = 3; // RGB

  std::filesystem::path outputPath = 
    std::filesystem::path(outputDir) / fileName;
  std::string fullPath = outputPath.string();

  std::vector<unsigned char> data(width * height * channels);

  for (int y = 0; y < height; ++y)
  {
    for (int x = 0; x < width; ++x)
    {
      Color c = image[y][x];
      c.clamp();
      unsigned char r = static_cast<unsigned char>
				(std::round(c.r));
			unsigned char g = static_cast<unsigned char>
				(std::round(c.g));
			unsigned char b = static_cast<unsigned char>
				(std::round(c.b));
      int index = (y * width + x) * channels;
      data[index + 0] = r;
      data[index + 1] = g;
      data[index + 2] = b;
    }
  }

  int success = stbi_write_png(fullPath.c_str(), 
    width, height, channels, data.data(), width * channels);

  if (success)
  {
    std::cout << "Successfully saved image: " << fullPath << std::endl;
  }
  else
  {
    std::cerr << "Error: Could not save image: " << fullPath << std::endl;
  }
}