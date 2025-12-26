#include <iostream>
#include <memory>
#include <chrono>

#include "parser/parser.hpp"
#include "render/render_manager.h"
#include "render/base_ray_tracer.h"
#include "material/material_manager.h"
#include "scene/scene.h"
#include "objects/sphere.h"
#include "objects/triangle.h"
#include "objects/plane.h"
#include "config.h"
#include "objects/tlas_box.h"
#include "objects/mesh.h"
#include "objects/geometry.h"
#include "texture_mapping/texture_fetcher.h"


std::vector<Geometry> geometries;
std::vector<ObjectContext> object_contexes;


int main(int argc, char* argv[])
{
#if COMMAND_LINE_INPUT
  // Expect exactly one argument after the executable name
  if (argc != 2)
  {
    std::cerr << "usage: " << argv[0] << " <scene_file.json>" << std::endl;
    return 1;
  }
  std::string scene_filename = argv[1];
#else
  //std::string scene_filename = "D:/Furkan/GITHUB/HelixNebula/inputs/dragon_dynamic.json";
  //std::string scene_filename = "D:/Furkan/GITHUB/HelixNebula/inputs/cornellbox_brushed_metal.json";
  //std::string scene_filename = "D:/Furkan/GITHUB/HelixNebula/inputs/spheres_dof.json";
  //std::string scene_filename = "D:/Furkan/GITHUB/HelixNebula/inputs/metal_glass_plates.json";
  //std::string scene_filename = "D:/Furkan/GITHUB/HelixNebula/inputs/ramazan_tokay/chessboard_arealight_dof_glass_queen.json";
  //std::string scene_filename = "D:/Furkan/GITHUB/HelixNebula/inputs/tap_water/json/tap_0010.json";
  //std::string scene_filename = "D:/Furkan/GITHUB/HelixNebula/inputs/test.json";
  //std::string scene_filename = "D:/Furkan/GITHUB/HelixNebula/inputs2/galactica_dynamic.json";
  //std::string scene_filename = "D:/Furkan/GITHUB/HelixNebula/inputs2/bump_mapping_transformed.json";
  //std::string scene_filename = "D:/Furkan/GITHUB/HelixNebula/inputs2/wood_box_all.json";
  //std::string scene_filename = "D:/Furkan/GITHUB/HelixNebula/inputs2/dragon/dragon_new_ply.json";
  //std::string scene_filename = "D:/Furkan/GITHUB/HelixNebula/inputs2/veach_ajar/scene.json";
  std::string scene_filename = "D:/Furkan/GITHUB/HelixNebula/inputs2/mytap/mytap_final.json";
  //std::string scene_filename = "D:/Furkan/GITHUB/HelixNebula/inputs2/tunnel_of_doom/tunnel_of_doom_000.json";
#endif

  Scene_ raw_scene;
  
  parseScene(scene_filename, raw_scene);
  
  std::vector<TLASBox> tlas_boxes;

  int t_size = raw_scene.triangles.size();
  int s_size = raw_scene.spheres.size();
  int m_size = raw_scene.meshes.size();
  int mi_size = raw_scene.mesh_instances.size();

  object_contexes.reserve(t_size + s_size + mi_size);
  geometries.reserve(t_size + s_size + m_size);

  const auto& vertex_data = raw_scene.vertex_data;
  const auto& uv_data = raw_scene.tex_coord_data;
  TextureLookup tex 
    = uv_data.empty() ? TextureLookup::NoTexture : TextureLookup::Textured;

  for(int i = 0; i < t_size; i++)
  {
    const auto& raw_triangle = raw_scene.triangles[i];
    
    glm::vec3 indices[3] = {
      vertex_data[raw_triangle.v0_id],
      vertex_data[raw_triangle.v1_id],
      vertex_data[raw_triangle.v2_id]
    };

    std::optional<glm::mat4> inv_tr;

    if (raw_triangle.transform_matrix.has_value())
      inv_tr = glm::inverse(raw_triangle.transform_matrix.value());

    object_contexes.emplace_back(raw_triangle.transform_matrix, inv_tr, raw_triangle.motion_blur, raw_triangle.textures, raw_triangle.material_id);

    if (tex == TextureLookup::NoTexture)
    {
      geometries.emplace_back(std::in_place_type<TriangleNew<Shading::Flat, TextureLookup::NoTexture>>, indices[0], indices[1], indices[2]);
    }
    else
    {
      geometries.emplace_back(
        std::in_place_type<
        TriangleNew<Shading::Flat, TextureLookup::Textured>>,
        indices[0], indices[1], indices[2],
        TexCoords(raw_triangle, uv_data)
      );
    }
    
    tlas_boxes.emplace_back(i, &object_contexes, &geometries[i]);

  }
  for (int i = 0; i < s_size; i++)
  {
    const auto& raw_sphere = raw_scene.spheres[i];
    bool tex = !raw_sphere.textures.empty();
    glm::vec3 center = vertex_data[raw_sphere.center_vertex_id];

    std::optional<glm::mat4> inv_tr;
    if (raw_sphere.transform_matrix.has_value())
      inv_tr = glm::inverse(raw_sphere.transform_matrix.value());

    object_contexes.emplace_back(raw_sphere.transform_matrix, inv_tr, raw_sphere.motion_blur, raw_sphere.textures, raw_sphere.material_id);
    geometries.emplace_back(std::in_place_type<Sphere>, center, static_cast<double>(raw_sphere.radius), tex);
    tlas_boxes.emplace_back(i + t_size, &object_contexes, &geometries[i + t_size]);
  }
  std::unordered_map<int, size_t> mesh_order;
  size_t index = 0;
  for (const auto& [key, val] : raw_scene.meshes)
  {
    mesh_order[key] = index++;
    if (tex == TextureLookup::NoTexture)
    {
      if (val.smooth_shading)
      {
        geometries.emplace_back(
          std::in_place_type<Mesh<Shading::Smooth, TextureLookup::NoTexture>>,
          val.id, val.faces, vertex_data, uv_data, val.vertex_offset, val.texture_offset);
      }
      else
      {
        geometries.emplace_back(
          std::in_place_type<Mesh<Shading::Flat, TextureLookup::NoTexture>>,
          val.id, val.faces, vertex_data, uv_data, val.vertex_offset, val.texture_offset);
      }
      
    }
    else
    {
      if (val.smooth_shading)
      {
        geometries.emplace_back(
          std::in_place_type<Mesh<Shading::Smooth, TextureLookup::Textured>>,
          val.id, val.faces, vertex_data, uv_data, val.vertex_offset, val.texture_offset);
      }
      else
      {
        geometries.emplace_back(
          std::in_place_type<Mesh<Shading::Flat, TextureLookup::Textured>>,
          val.id, val.faces, vertex_data, uv_data, val.vertex_offset, val.texture_offset);
      }
    }
    
  }
  index = 0;
  int base = t_size + s_size;
  for (const auto& [id, mi] : raw_scene.mesh_instances)
  {
    std::optional<glm::mat4> inv_tr;
    if (mi.transform_matrix.has_value())
      inv_tr = glm::inverse(mi.transform_matrix.value());

    object_contexes.emplace_back(mi.transform_matrix, inv_tr, mi.motion_blur, mi.textures, mi.material_id);
    tlas_boxes.emplace_back(base + index, &object_contexes, &geometries[base + mesh_order[mi.base_mesh_id]]);
    index++;
  }
  

  std::vector<Plane> planes;
  for (const Plane_& raw_plane : raw_scene.planes)
  {
    //planes are transformed
    planes.push_back(
            Plane(raw_plane, raw_scene.vertex_data,raw_plane.motion_blur));
  }

	MaterialManager material_manager(raw_scene.materials);

  Scene scene(raw_scene, tlas_boxes);

  RendererInfo renderer_info(raw_scene.shadow_ray_epsilon, 
    raw_scene.intersection_test_epsilon, 
    raw_scene.max_recursion_depth);

  TextureData data;

  for (const auto& img : raw_scene.images)
  {
    data.images[img.id] = Image(img);
  }

  

  uint8_t replace_background_num = 0;
  for (const auto& tm : raw_scene.texture_maps)
  {
    if (tm.decal_mode == "replace_background")
    {
      ++replace_background_num;
      Expects(replace_background_num <= 1);
      renderer_info.background_tex_id  = tm.id;
    }
    data.textures[tm.id] = TextureMap(tm);
  }

  
  TextureFetcher tex_fetcher(data);

  BaseRayTracer ray_tracer(scene.background_color, scene.light_sources, 
    scene.world, planes, material_manager, renderer_info, tex_fetcher);

  RenderManager renderer(scene, material_manager, renderer_info, ray_tracer);

  std::cout << "Rendering started for scene file: " << scene_filename << std::endl;
  auto start = std::chrono::high_resolution_clock::now();

  renderer.render();

  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = end - start;
	std::cout << "Rendering completed in " << elapsed.count() << " seconds." << std::endl;
  return 0;
}