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

#include "../objects/geometry.h"


std::vector<Geometry> geometries;
std::vector<ObjectContext> object_contexes;
//std::vector<std::shared_ptr<Hittable>> local_space_objects;

int main(int argc, char* argv[])
{
  //Geometry sphere(std::in_place_type<Sphere>, glm::vec3{ 0,0,0 }, 1.0);
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
  //std::string scene_filename = "D:/Furkan/GITHUB/HelixNebula/inputs/dragon_dynamic.json";
  //std::string scene_filename = "D:/Furkan/GITHUB/HelixNebula/inputs/spheres_dof.json";
  //std::string scene_filename = "D:/Furkan/GITHUB/HelixNebula/inputs/metal_glass_plates.json";
  std::string scene_filename = "D:/Furkan/GITHUB/HelixNebula/inputs/focusing_dragons.json";
  //std::string scene_filename = "D:/Furkan/GITHUB/HelixNebula/inputs/ramazan_tokay/chessboard_arealight_dof_glass_queen.json";
  //std::string scene_filename = "D:/Furkan/GITHUB/HelixNebula/inputs/tap_water/json/tap_0010.json";
#endif
  Scene_ raw_scene;
  
#if !DAVID_ZOOM
  std::cout << "Parsing scene file: " << scene_filename << std::endl;
#endif
  parseScene(scene_filename, raw_scene);

  //std::vector<std::shared_ptr<Hittable>> tlas_boxes;
  std::vector<TLASBox> tlas_boxes;

  int t_size = raw_scene.triangles.size();
  int s_size = raw_scene.spheres.size();
  int m_size = raw_scene.meshes.size();
  int mi_size = raw_scene.mesh_instances.size();

  object_contexes.reserve(t_size + s_size + mi_size);
  geometries.reserve(t_size + s_size + m_size);
  //local_space_objects.reserve(t_size + s_size + m_size);
  //tlas_boxes.reserve(t_size + s_size + mi_size);

  const auto& vertex_data = raw_scene.vertex_data;

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

    object_contexes.emplace_back(raw_triangle.transform_matrix, inv_tr, raw_triangle.motion_blur, raw_triangle.material_id);
    //local_space_objects.emplace_back(std::make_shared<Triangle>(indices));
    geometries.emplace_back(std::in_place_type<Triangle>, indices);
    //tlas_boxes.emplace_back(std::make_shared<TLASBox>(i, object_contexes, local_space_objects[i]));
    tlas_boxes.emplace_back(i, &object_contexes, &geometries[i]);

  }
  for (int i = 0; i < s_size; i++)
  {
    const auto& raw_sphere = raw_scene.spheres[i];
    glm::vec3 center = vertex_data[raw_sphere.center_vertex_id];

    std::optional<glm::mat4> inv_tr;
    if (raw_sphere.transform_matrix.has_value())
      inv_tr = glm::inverse(raw_sphere.transform_matrix.value());

    object_contexes.emplace_back(raw_sphere.transform_matrix, inv_tr, raw_sphere.motion_blur, raw_sphere.material_id);
    //local_space_objects.emplace_back(std::make_shared<Sphere>(center, static_cast<double>(raw_sphere.radius)));
    geometries.emplace_back(std::in_place_type<Sphere>, center, static_cast<double>(raw_sphere.radius));
    //tlas_boxes.emplace_back(std::make_shared<TLASBox>(i + t_size, object_contexes, local_space_objects[i + t_size]));
    tlas_boxes.emplace_back(i + t_size, &object_contexes, &geometries[i + t_size]);
  }
  std::unordered_map<int, size_t> mesh_order;
  size_t index = 0;
  for (const auto& [key, val] : raw_scene.meshes)
  {
    mesh_order[key] = index++;
    //local_space_objects.emplace_back(std::make_shared<Mesh>(val.id, val.smooth_shading, val.faces, vertex_data));
    geometries.emplace_back(std::in_place_type<Mesh>, val.id, val.smooth_shading, val.faces, vertex_data);
  }
  index = 0;
  int base = t_size + s_size;
  for (const auto& [id, mi] : raw_scene.mesh_instances)
  {
    std::optional<glm::mat4> inv_tr;
    if (mi.transform_matrix.has_value())
      inv_tr = glm::inverse(mi.transform_matrix.value());

    object_contexes.emplace_back(mi.transform_matrix, inv_tr, mi.motion_blur, mi.material_id);
    //tlas_boxes.emplace_back(std::make_shared<TLASBox>(base + index, object_contexes, local_space_objects[base + mesh_order[mi.base_mesh_id]]));
    tlas_boxes.emplace_back(base + index, &object_contexes, &geometries[base + mesh_order[mi.base_mesh_id]]);
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