#include <iostream>
#include <memory>
#include <chrono>
#include <filesystem>
#include <external/json.hpp>

#include "parser/parser.hpp"
#include "scene/scene.h"
#include "objects/sphere.h"
#include "objects/triangle.h"
#include "objects/plane.h"
#include "config.h"
#include "objects/tlas_box.h"
#include "objects/mesh.h"
#include "objects/geometry.h"
#include "render/raytracer.h"
#include "texture_mapping/texture_data.h"

int main(int argc, char* argv[])
{
#if COMMAND_LINE_INPUT
  // Expect exactly one argument after the executable name
  if (argc != 2)
  {
    std::cerr << "usage: " << argv[0] << " <scene_file.json>" << std::endl;
    return 1;
  }
  std::string scene_path = argv[1];
#else
  std::string scene_folder = FS::absolute(__FILE__).parent_path() / "../inputs5/akin_aydemir/teapot_roughness/";
  std::string scene_filename = "teapot_roughness";
  scene_filename += ".json";
  std::string scene_path = scene_folder + scene_filename;
#endif

  Scene_ raw_scene;
  parseScene(scene_path, raw_scene);

  std::vector<Geometry> geometries;
  std::vector<std::optional<Transformation>> transformations;
  std::vector<Material> materials;
  std::vector<Texture> textures;
  std::vector<Image> images;

  images.resize(raw_scene.images.size());
  for (const auto& img : raw_scene.images)
  {
    std::string absolute_path = scene_folder + img.data;
    std::string extension = getFileExtension(absolute_path);
    if (extension == "hdr" || extension == "exr")
      images[img.id] = Image(absolute_path, ImageType::HDR);
    else
      images[img.id] = Image(absolute_path, ImageType::SDR);
  }

  textures.resize(raw_scene.texture_maps.size());
  for (const auto& tm : raw_scene.texture_maps)
  {
    textures[tm.id] = Texture(tm);
    auto& texture = textures[tm.id];
    if (texture.type == TextureType::image)
    {
      if (texture.d_mode == DecalMode::bump_normal)
      {
        texture.make_bump_image(&images[tm.image_id], tm.normalizer, tm.bump_factor);
        continue;
      }
      texture.make_image(&images[tm.image_id], tm.normalizer);

    }
    else if (texture.type == TextureType::perlin)
    {
      NoiseConversion nc = tm.noise_conversion == "linear" ? NoiseConversion::linear : NoiseConversion::absval;
      if (texture.d_mode == DecalMode::bump_normal)
      {
        texture.make_bump_perlin(nc, tm.bump_factor, tm.noise_scale, tm.num_octaves);
        continue;
      }
      texture.make_perlin(nc, tm.noise_scale, tm.num_octaves);
    }
    else if (texture.type == TextureType::checkerboard)
    {
      texture.make_checkerboard(tm.black_color, tm.white_color, tm.scale, tm.offset);
    }
    else
    {
      throw std::runtime_error("Wrong texture parsing");
    }
  }

  Texture* bg_tex = nullptr;
  for (Texture& tex : textures)
  {
    if (tex.d_mode == DecalMode::replace_background)
    {
      bg_tex = &tex;
      break;
    }
  }

  materials.reserve(raw_scene.materials.size());
  for (auto& material : raw_scene.materials)
    materials.emplace_back(material);
  
  std::vector<TLASBox> tlas_boxes;

  int t_size = raw_scene.triangles.size();
  int s_size = raw_scene.spheres.size();
  int m_size = raw_scene.meshes.size();
  int mi_size = raw_scene.mesh_instances.size();

  transformations.reserve(t_size + s_size + mi_size);
  geometries.reserve(t_size + s_size + m_size);

  const auto& vertex_data = raw_scene.vertex_data;

  for(int i = 0; i < t_size; i++)
  {
    const auto& raw_triangle = raw_scene.triangles[i];
    
    glm::vec3 indices[3] = {
      vertex_data[raw_triangle.v0_id],
      vertex_data[raw_triangle.v1_id],
      vertex_data[raw_triangle.v2_id]
    };


    if (raw_triangle.transform_matrix.has_value())
    {
      transformations.emplace_back(
          Transformation{
              raw_triangle.transform_matrix.value(),  // safer than .value()
              glm::inverse(raw_triangle.transform_matrix.value())
          }
      );
    }
    else
    {
      transformations.emplace_back(std::nullopt);
    }

    glm::vec2 tex_coords[3];
    if (!raw_scene.tex_coord_data.empty())
    {
      tex_coords[0] = raw_scene.tex_coord_data[raw_triangle.v0_id];
      tex_coords[1] = raw_scene.tex_coord_data[raw_triangle.v1_id];
      tex_coords[2] = raw_scene.tex_coord_data[raw_triangle.v2_id];
    }
    else
    {
      tex_coords[0] = {0.f, 0.f};
      tex_coords[1] = {0.f, 0.f};
      tex_coords[2] = {0.f, 0.f};
    }

    geometries.emplace_back(std::in_place_type<Triangle>, indices, tex_coords);

    std::vector<Texture*> tex;
    for (int tex_id : raw_triangle.textures)
    {
      tex.emplace_back(&textures[tex_id]);
    }

    tlas_boxes.emplace_back(&geometries.back(), &materials[raw_triangle.material_id],
      &transformations.back(), raw_triangle.motion_blur, tex);

  }

  for (int i = 0; i < s_size; i++)
  {
    const auto& raw_sphere = raw_scene.spheres[i];
    glm::vec3 center = vertex_data[raw_sphere.center_vertex_id];

    if (raw_sphere.transform_matrix.has_value())
    {
      transformations.emplace_back(
          Transformation{
              raw_sphere.transform_matrix.value(),  // safer than .value()
              glm::inverse(raw_sphere.transform_matrix.value())
          }
      );
    }
    else
    {
      transformations.emplace_back(std::nullopt);
    }
    geometries.emplace_back(std::in_place_type<Sphere>, center, static_cast<double>(raw_sphere.radius));

    std::vector<Texture*> tex;
    for (int tex_id : raw_sphere.textures)
    {
      tex.emplace_back(&textures[tex_id]);
    }
    tlas_boxes.emplace_back(&geometries.back(), &materials[raw_sphere.material_id],
      &transformations.back(), raw_sphere.motion_blur, tex);
  }

  std::unordered_map<int, size_t> mesh_order;
  size_t index = 0;
  for (const auto& [key, val] : raw_scene.meshes)
  {
    mesh_order[key] = index++;
    geometries.emplace_back(std::in_place_type<Mesh>,val.id, val.faces, vertex_data,
      val.vertex_offset, val.texture_offset, val.smooth_shading, raw_scene.tex_coord_data);
  }
  index = 0;
  int base = t_size + s_size;
  for (const auto& [id, mi] : raw_scene.mesh_instances)
  {
    if (mi.transform_matrix.has_value())
    {
      transformations.emplace_back(
          Transformation{
              mi.transform_matrix.value(),  // safer than .value()
              glm::inverse(mi.transform_matrix.value())
          }
      );
    }
    else
    {
      transformations.emplace_back(std::nullopt);
    }

    std::vector<Texture*> tex;
    for (int tex_id : mi.textures)
    {
      tex.emplace_back(&textures[tex_id]);
    }
    tlas_boxes.emplace_back(&geometries[base + mesh_order[mi.base_mesh_id]], &materials[mi.material_id],
      &transformations.back(), mi.motion_blur, tex);
    index++;
  }
  

  std::vector<Plane> planes;
  for (const Plane_& raw_plane : raw_scene.planes)
  {
    //planes are transformed
    planes.push_back(
            Plane(raw_plane, raw_scene.vertex_data,raw_plane.motion_blur,
              &materials[raw_plane.material_id], textures)
            );
  }

  //Scene scene(raw_scene, tlas_boxes, planes);


  Image* env_map = nullptr;
  if (!raw_scene.spherical_directional_lights.empty())
    env_map = &images[raw_scene.spherical_directional_lights[0].image_id];

  initPerlin();
  RenderContext render_context;
  if (bg_tex != nullptr)
  {
    render_context.b_type = BackgroundType::Texture;
    render_context.background_info.background_tex = bg_tex;
  }
  else
  {
    render_context.b_type = BackgroundType::Color;
    render_context.background_info.background_color = raw_scene.background_color;
  }
  render_context.shadow_ray_epsilon = raw_scene.shadow_ray_epsilon;
  render_context.intersection_test_epsilon = raw_scene.intersection_test_epsilon;
  render_context.max_recursion_depth = raw_scene.max_recursion_depth;
  render_context.env_map = env_map;




  Raytracer raytracer(std::make_unique<Scene>(raw_scene, tlas_boxes, planes, env_map), render_context);


    std::cout << "Rendering started for scene file: " << scene_filename << std::endl;
  auto start = std::chrono::high_resolution_clock::now();

  raytracer.renderScene();

  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = end - start;
	std::cout << "Rendering completed in " << elapsed.count() << " seconds." << std::endl;
  return 0;
}