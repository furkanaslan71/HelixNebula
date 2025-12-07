#include <iostream>
#include <memory>
#include "../include/parser.hpp"
#include "../render/render_manager.h"
#include "../render/base_ray_tracer.h"
#include "../material/material_manager.h"
#include "../scene/scene.h"
#include "../objects/sphere.h"
#include "../objects/triangle.h"
#include "../objects/plane.h"
#include "../include/config.h"
#include "../objects/tlas_box.h"
#include "../objects/mesh.h"

#define GLOBALS local_space_objects, transform_matrices, material_ids, motion_blur

std::vector<std::shared_ptr<Hittable>> local_space_objects;
std::vector<std::optional<glm::mat4>> transform_matrices;
std::vector<int> material_ids;
std::vector<Vec3> motion_blur;

int main(int argc, char* argv[])
{
  // Expect exactly one argument after the executable name
  /*
  if (argc != 2)
  {
    std::cerr << "Usage: " << argv[0] << " <scene_file.json>" << std::endl;
    return 1;
  }

  std::string scene_filename = argv[1];
  */

  

#if DAVID_ZOOM
for(int i = 0; i < 360; i++)
{ 
  std::string scene_filename = "D:/Furkan/repos/raytracer/HelixNebula/inputs2/raven/camera_around_david/davids_camera_";
	// Append zero-padded frame number
	scene_filename += (i < 10) ? "00" : (i < 100) ? "0" : "";
	scene_filename += std::to_string(i) + ".json";
#else
	//std::string scene_filename = "D:/Furkan/repos/raytracer/HelixNebula/inputs2/raven/camera_zoom_david/davids_camera_zoom_0.json";
  //std::string scene_filename = "D:/Furkan/repos/raytracer/HelixNebula/inputs2/akif_uslu/berserker/two_berserkers.json";
  //std::string scene_filename = "D:/Furkan/repos/raytracer/HelixNebulaNew/inputs2/marching_dragons.json";
  //std::string scene_filename = "D:/Furkan/repos/raytracer/HelixNebula/inputs2/raven/dragon/dragon_new_right_ply.json";
  std::string scene_filename = "D:/Furkan/GITHUB/HelixNebula/inputs/metal_glass_plates.json";
  //std::string scene_filename = "D:/Furkan/GITHUB/HelixNebula/inputs/spheres_dof.json";
  //std::string scene_filename = "D:/Furkan/repos/raytracer/HelixNebulaNew/inputs3/ramazan_tokay/wine_glass_scene.json";
  //std::string scene_filename = "D:/Furkan/repos/raytracer/HelixNebulaNew/inputs3/tap_water/json/tap_0010.json";
#endif
  Scene_ raw_scene;
  
#if !DAVID_ZOOM
  std::cout << "Parsing scene file: " << scene_filename << std::endl;
#endif
  parseScene(scene_filename, raw_scene);

  //printSceneSummary(scene);
  //printScene(raw_scene);
  std::vector<std::shared_ptr<Hittable>> tlas_boxes;

  int t_size = raw_scene.triangles.size();
  int s_size = raw_scene.spheres.size();
  int m_size = raw_scene.meshes.size();
  int mi_size = raw_scene.mesh_instances.size();

  local_space_objects.reserve(t_size + s_size + m_size);
  transform_matrices.reserve(t_size + s_size + mi_size);
  material_ids.reserve(t_size + s_size + mi_size);
  motion_blur.reserve(t_size + s_size + mi_size);
  tlas_boxes.reserve(t_size + s_size + mi_size);

  const auto& vertex_data = raw_scene.vertex_data;

  for(int i = 0; i < t_size; i++)
  {
    const auto& raw_triangle = raw_scene.triangles[i];
    
    Vec3 indices[3] = { 
      vertex_data[raw_triangle.v0_id],
      vertex_data[raw_triangle.v1_id],
      vertex_data[raw_triangle.v2_id]
    };
    local_space_objects.emplace_back(std::make_shared<Triangle>(indices));
    transform_matrices.emplace_back(raw_triangle.transform_matrix);
    material_ids.emplace_back(raw_triangle.material_id);
    motion_blur.emplace_back(raw_triangle.motion_blur);
    tlas_boxes.emplace_back(std::make_shared<TLASBox>(i, i, GLOBALS));
  }
  for (int i = 0; i < s_size; i++)
  {
    const auto& raw_sphere = raw_scene.spheres[i];
    Vec3 center = vertex_data[raw_sphere.center_vertex_id];
    local_space_objects.emplace_back(std::make_shared<Sphere>(center, static_cast<double>(raw_sphere.radius)));
    transform_matrices.emplace_back(raw_sphere.transform_matrix);
    material_ids.emplace_back(raw_sphere.material_id);
    motion_blur.emplace_back(raw_sphere.motion_blur);
    tlas_boxes.emplace_back(std::make_shared<TLASBox>(
      i + t_size,
      i + t_size,
      GLOBALS));
  }
  std::unordered_map<int, size_t> mesh_order;
  size_t index = 0;
  for (const auto& [key, val] : raw_scene.meshes)
  {
    mesh_order[key] = index++;
    local_space_objects.emplace_back(std::make_shared<Mesh>(val.id, val.smooth_shading, val.faces, vertex_data));
  }
  index = 0;
  int base = t_size + s_size;
  for (const auto& [id, mi] : raw_scene.mesh_instances)
  {
    transform_matrices.emplace_back(mi.transform_matrix);
    material_ids.emplace_back(mi.material_id);
    motion_blur.emplace_back(mi.motion_blur);
    tlas_boxes.emplace_back(std::make_shared<TLASBox>(
      base + mesh_order[mi.base_mesh_id],
      base + index,
      GLOBALS));
    index++;
  }
  

  std::vector<Plane> planes;
  for (const Plane_& raw_plane : raw_scene.planes)
  {
    //planes are transformed
    planes.push_back(
      Plane(raw_plane, raw_scene.vertex_data,raw_plane.motion_blur)
    );
  }

	MaterialManager material_manager(raw_scene.materials);

  Scene scene(raw_scene, tlas_boxes);

  RendererInfo renderer_info(raw_scene.shadow_ray_epsilon, 
    raw_scene.intersection_test_epsilon, 
    raw_scene.max_recursion_depth);

  BaseRayTracer ray_tracer(scene.background_color, scene.light_sources, 
    scene.world, planes, material_manager, renderer_info);

  RenderManager renderer(scene, material_manager, renderer_info, ray_tracer);
#if !DAVID_ZOOM
  std::cout << "Rendering will start here in the future." << std::endl;
#endif
	//measure time with standart library
	double time_start = static_cast<double>(clock()) / CLOCKS_PER_SEC;
  renderer.render();
	double time_end = static_cast<double>(clock()) / CLOCKS_PER_SEC;

	std::cout << "Rendering completed in " << (time_end - time_start) << " seconds." << std::endl;

#if DAVID_ZOOM
	}
#endif
  return 0;
}