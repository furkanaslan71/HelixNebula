#include <fstream>
#include <iostream>
#include <sstream>

#include "parser.hpp"
#include "external/json.hpp"
#include "math_core/transformations/transformations.h"

#define M_PI 3.14159265358979323846
#define COMPOSITE_TRANSFORM calculateCompositeTransformationMatrix( \
transformations, \
scene.translations, \
scene.scalings,scene.rotations \
) \

enum class FieldRequirement {
  MUST,
  DEFAULTABLE,
  OPTIONAL
};

// For convenience
using json = nlohmann::json;

// Helper functions
glm::vec3 parseVec3(const std::string& str) {
    glm::vec3 vec;
    std::stringstream ss(str);
    ss >> vec.x >> vec.y >> vec.z;
    return vec;
}

glm::vec4 parseVec4(const std::string& str) {
    glm::vec4 vec;
    std::stringstream ss(str);
    ss >> vec.x >> vec.y >> vec.z >> vec.w;
    return vec;
}

static float readFloatLE(std::ifstream& f)
{
  float val;
  f.read(reinterpret_cast<char*>(&val), sizeof(float));
  return val;
}

static int32_t readInt32LE(std::ifstream& f)
{
  int32_t val;
  f.read(reinterpret_cast<char*>(&val), sizeof(int32_t));
  return val;
}

static uint8_t readUInt8(std::ifstream& f)
{
  uint8_t val;
  f.read(reinterpret_cast<char*>(&val), sizeof(uint8_t));
  return val;
}

bool perpendicular(const glm::vec3& a, const glm::vec3& b)
{
  return fabs(glm::dot(a, b)) < 1e-8;
}

// --- UNIVERSAL, ROBUST PLY PARSER (ASCII + basic binary) ---
namespace PlyHelpers
{
  // --- Endianness Check ---
  bool is_system_little_endian()
  {
    uint32_t i = 1;
    char* c = (char*)&i;
    return (*c);
  }

  // --- Byte Swap ---
  template <typename T>
  T byte_swap(T val)
  {
    T swapped_val;
    char* src = reinterpret_cast<char*>(&val);
    char* dst = reinterpret_cast<char*>(&swapped_val);
    for (size_t i = 0; i < sizeof(T); ++i)
    {
      dst[i] = src[sizeof(T) - 1 - i];
    }
    return swapped_val;
  }

  // --- Data Extraction from Buffer ---
  template <typename T>
  T extract_from_buffer(const char* buffer, bool file_is_little_endian)
  {
    T val;
    std::memcpy(&val, buffer, sizeof(T)); // Use memcpy to avoid alignment issues
    if (file_is_little_endian != is_system_little_endian())
    {
      val = byte_swap(val);
    }
    return val;
  }

  // --- Get Type Size ---
  size_t get_ply_type_size(const std::string& type)
  {
    if (type == "char" || type == "uchar" || type == "int8" || type == "uint8") return 1;
    if (type == "short" || type == "ushort" || type == "int16" || type == "uint16") return 2;
    if (type == "int" || type == "uint" || type == "float" || type == "int32" || type == "uint32" || type == "float32") return 4;
    if (type == "double" || type == "int64" || type == "uint64" || type == "float64") return 8;
    std::cerr << "Warning: Unknown PLY type '" << type << "', assuming 0 size." << std::endl;
    return 0;
  }

  // --- Extract Integer (any size) ---
  int extract_int_from_buffer(const char* buffer, const std::string& type, bool file_is_little_endian)
  {
    if (type == "char" || type == "int8") return static_cast<int>(extract_from_buffer<int8_t>(buffer, file_is_little_endian));
    if (type == "uchar" || type == "uint8") return static_cast<int>(extract_from_buffer<uint8_t>(buffer, file_is_little_endian));
    if (type == "short" || type == "int16") return static_cast<int>(extract_from_buffer<int16_t>(buffer, file_is_little_endian));
    if (type == "ushort" || type == "uint16") return static_cast<int>(extract_from_buffer<uint16_t>(buffer, file_is_little_endian));
    if (type == "int" || type == "int32") return static_cast<int>(extract_from_buffer<int32_t>(buffer, file_is_little_endian));
    if (type == "uint" || type == "uint32") return static_cast<int>(extract_from_buffer<uint32_t>(buffer, file_is_little_endian));
    return 0;
  }

  // --- Extract Float/Double ---
  float extract_float_from_buffer(const char* buffer, const std::string& type, bool file_is_little_endian)
  {
    if (type == "float" || type == "float32") return extract_from_buffer<float>(buffer, file_is_little_endian);
    if (type == "double" || type == "float64") return static_cast<float>(extract_from_buffer<double>(buffer, file_is_little_endian));
    return 0.0f;
  }

  // --- PLY Property Structs ---
  struct PlyProperty {
    std::string name;
    std::string type;
    size_t size = 0;
    size_t offset = 0;
  };

  struct PlyFaceProperty {
    std::string count_type = "uchar";
    std::string index_type = "int";
    size_t count_size = 1;
    size_t index_size = 4;
  };

  enum PlyFormat { UNKNOWN, ASCII, BINARY_LITTLE_ENDIAN, BINARY_BIG_ENDIAN };

} // namespace PlyHelpers


// --- UNIVERSAL, ROBUST PLY PARSER (ASCII + BINARY) ---
void parsePlyFile(const std::string& ply_filename, Mesh_& mesh, Scene_& scene)
{
  using namespace PlyHelpers;

  std::ifstream file(ply_filename, std::ios::binary);
  if (!file.is_open())
  {
    std::cerr << "Error: Could not open PLY file: " << ply_filename << std::endl;
    return;
  }

  // --- Header Parsing ---
  std::string line;
  PlyFormat format = UNKNOWN;
  long num_vertices = 0;
  long num_faces = 0;

  std::vector<PlyProperty> vertex_props;
  PlyFaceProperty face_prop;
  size_t vertex_record_size = 0;
  int x_prop_idx = -1, y_prop_idx = -1, z_prop_idx = -1;

  bool in_vertex_element = false;
  bool in_face_element = false;

  while (std::getline(file, line))
  {
    // Trim trailing CR/LF
    line.erase(line.find_last_not_of(" \t\n\r") + 1);

    std::stringstream ss(line);
    std::string token;
    ss >> token;

    if (token == "ply") continue;
    if (token == "comment") continue;

    if (token == "format")
    {
      ss >> token;
      if (token == "ascii") format = ASCII;
      else if (token == "binary_little_endian") format = BINARY_LITTLE_ENDIAN;
      else if (token == "binary_big_endian") format = BINARY_BIG_ENDIAN;
    }
    else if (token == "element")
    {
      ss >> token;
      if (token == "vertex")
      {
        ss >> num_vertices;
        in_vertex_element = true;
        in_face_element = false;
        vertex_record_size = 0;
        vertex_props.clear();
      }
      else if (token == "face")
      {
        ss >> num_faces;
        in_vertex_element = false;
        in_face_element = true;
      }
      else
      {
        in_vertex_element = false;
        in_face_element = false;
      }
    }
    else if (token == "property")
    {
      if (in_vertex_element)
      {
        PlyProperty p;
        ss >> p.type >> p.name;
        p.size = get_ply_type_size(p.type);
        p.offset = vertex_record_size;
        vertex_record_size += p.size;

        if (p.name == "x") x_prop_idx = vertex_props.size();
        else if (p.name == "y") y_prop_idx = vertex_props.size();
        else if (p.name == "z") z_prop_idx = vertex_props.size();

        vertex_props.push_back(p);
      }
      else if (in_face_element)
      {
        std::string list_token, count_type, index_type, prop_name;
        ss >> list_token >> count_type >> index_type >> prop_name;
        if (list_token == "list" && prop_name == "vertex_index")
        {
          face_prop.count_type = count_type;
          face_prop.index_type = index_type;
          face_prop.count_size = get_ply_type_size(count_type);
          face_prop.index_size = get_ply_type_size(index_type);
        }
      }
    }
    else if (token == "end_header")
    {
      break; // End of header, data follows
    }
  }

  // --- Header Validation ---
  if (format == UNKNOWN)
  {
    std::cerr << "Error: Unknown PLY format in " << ply_filename << std::endl;
    return;
  }
  if (num_vertices == 0)
  {
    std::cerr << "Error: PLY file has 0 vertices: " << ply_filename << std::endl;
    return;
  }
  if (x_prop_idx == -1 || y_prop_idx == -1 || z_prop_idx == -1)
  {
    std::cerr << "Error: PLY file missing x, y, or z vertex property: " << ply_filename << std::endl;
    return;
  }

  // --- Data Parsing ---

  // This is the critical step for combining meshes:
  // All indices from this file will be offset by the current vertex count.
  int vertex_base_index = scene.vertex_data.size();
  bool file_is_little_endian = (format == BINARY_LITTLE_ENDIAN);

#if !DAVID_ZOOM
  std::cout << "  Parsing " << num_vertices << " vertices from " << ply_filename << "..." << std::endl;
#endif
  // 1. Read Vertices
  std::vector<char> binary_buffer(vertex_record_size);
  for (long i = 0; i < num_vertices; ++i)
  {
    glm::vec3 v;

    if (format == ASCII)
    {
      if (!std::getline(file, line))
      {
        std::cerr << "Error: Unexpected end of file reading vertices from " << ply_filename << std::endl;
        return;
      }
      std::stringstream ss(line);
      float val;
      for (size_t j = 0; j < vertex_props.size(); ++j)
      {
        ss >> val; // Read value
        // Assign to correct component if it's x, y, or z
        if (j == x_prop_idx) v.x = val;
        else if (j == y_prop_idx) v.y = val;
        else if (j == z_prop_idx) v.z = val;
      }
    }
    else
    { // Binary formats
      file.read(binary_buffer.data(), vertex_record_size);
      if (!file)
      {
        std::cerr << "Error: Unexpected end of file reading binary vertices from " << ply_filename << std::endl;
        return;
      }

      const PlyProperty& p_x = vertex_props[x_prop_idx];
      v.x = extract_float_from_buffer(binary_buffer.data() + p_x.offset, p_x.type, file_is_little_endian);

      const PlyProperty& p_y = vertex_props[y_prop_idx];
      v.y = extract_float_from_buffer(binary_buffer.data() + p_y.offset, p_y.type, file_is_little_endian);

      const PlyProperty& p_z = vertex_props[z_prop_idx];
      v.z = extract_float_from_buffer(binary_buffer.data() + p_z.offset, p_z.type, file_is_little_endian);
    }

    // Add the new vertex to the GLOBAL scene vertex list
    scene.vertex_data.push_back(v);
  }
#if !DAVID_ZOOM
  std::cout << "  Parsing " << num_faces << " faces from " << ply_filename << "..." << std::endl;
#endif

  // 2. Read Faces
  std::vector<char> count_buffer(face_prop.count_size);
  std::vector<char> index_buffer(face_prop.index_size);

  for (long i = 0; i < num_faces; ++i)
  {
    std::vector<int> local_indices;
    int num_verts_in_face = 0;

    if (format == ASCII)
    {
      if (!std::getline(file, line))
      {
        std::cerr << "Error: Unexpected end of file reading faces from " << ply_filename << std::endl;
        return;
      }
      std::stringstream ss(line);
      ss >> num_verts_in_face;
      int index;
      for (int j = 0; j < num_verts_in_face; ++j)
      {
        ss >> index;
        local_indices.push_back(index);
      }
    }
    else
    { // Binary formats
      file.read(count_buffer.data(), face_prop.count_size);
      if (!file)
      {
        std::cerr << "Error: Unexpected end of file reading face count from " << ply_filename << std::endl;
        return;
      }
      num_verts_in_face = extract_int_from_buffer(count_buffer.data(), face_prop.count_type, file_is_little_endian);

      for (int j = 0; j < num_verts_in_face; ++j)
      {
        file.read(index_buffer.data(), face_prop.index_size);
        if (!file)
        {
          std::cerr << "Error: Unexpected end of file reading face index from " << ply_filename << std::endl;
          return;
        }
        local_indices.push_back(extract_int_from_buffer(index_buffer.data(), face_prop.index_type, file_is_little_endian));
      }
    }

    // --- Triangulation (same logic for all formats) ---
    if (num_verts_in_face < 3)
    {
      std::cerr << "Warning: Skipping face " << i << " (less than 3 vertices) in " << ply_filename << std::endl;
      continue;
    }

    // Triangulate the polygon (N-gon) using a triangle fan
    // Pivot vertex is v0
    int v0_global = local_indices[0] + vertex_base_index;

    for (size_t j = 1; j < local_indices.size() - 1; ++j)
    {
      int v1_global = local_indices[j] + vertex_base_index;
      int v2_global = local_indices[j + 1] + vertex_base_index;

      // Add the new face (triangle) to the MESH's face list
      // It uses the mesh's material_id and the GLOBAL vertex indices
      mesh.faces.push_back({ mesh.material_id, v0_global, v1_global, v2_global });
    }
  }
#if !DAVID_ZOOM
  std::cout << "  Finished parsing PLY file. Total vertices in scene: " << scene.vertex_data.size() << std::endl;
#endif
  file.close();
}


// Function implementation
void parseScene(const std::string& filename, Scene_& scene) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Could not open scene file: " << filename << std::endl;
        return;
    }

    json j;
    file >> j;
    const auto& scene_json = j["Scene"];

    // --- Global Scene Settings ---
    scene.background_color = parseVec3(scene_json["BackgroundColor"]);

    if (scene_json.contains("ShadowRayEpsilon"))
        scene.shadow_ray_epsilon = std::stof(scene_json["ShadowRayEpsilon"].get<std::string>());
    else
      scene.shadow_ray_epsilon = 1e-3;
			//scene.shadow_ray_epsilon = 0.0f;

    if (scene_json.contains("IntersectionTestEpsilon"))
        scene.intersection_test_epsilon = std::stof(scene_json["IntersectionTestEpsilon"].get<std::string>());
    else
        scene.intersection_test_epsilon = 1e-5;

    if (scene_json.contains("MaxRecursionDepth"))
        scene.max_recursion_depth = std::stoi(scene_json["MaxRecursionDepth"].get<std::string>());
    else
      scene.max_recursion_depth = 6;
			//scene.max_recursion_depth = 0;


    if (scene_json.contains("Transformations"))
    {
      const auto& transformations_json = scene_json["Transformations"];
      if (transformations_json.contains("Translation"))
      {
        const auto& translations_json = transformations_json["Translation"];
        if (!translations_json.is_array())
        {
          scene.translations.resize(2);
          glm::vec3 data = parseVec3(translations_json["_data"]);
          scene.translations[1] = Translation_{ data.x, data.y, data.z };
        }
        else
        {
          int size = translations_json.size();
          scene.translations.resize(size + 1);
          for (int i = 1; i < size + 1; i++)
          {
            glm::vec3 data = parseVec3(translations_json[i - 1]["_data"]);
            scene.translations[i] = Translation_{ data.x, data.y, data.z };
          }
        }
      }
      if (transformations_json.contains("Scaling"))
      {
        const auto& scalings_json = transformations_json["Scaling"];
        if (!scalings_json.is_array())
        {
          scene.scalings.resize(2);
          glm::vec3 data = parseVec3(scalings_json["_data"]);
          scene.scalings[1] = Scaling_{ data.x, data.y, data.z };
        }
        else
        {
          int size = scalings_json.size();
          scene.scalings.resize(size + 1);
          for (int i = 1; i < size + 1; i++)
          {
            glm::vec3 data = parseVec3(scalings_json[i - 1]["_data"]);
            scene.scalings[i] = Scaling_{ data.x, data.y, data.z };
          }
        }
      }
      if (transformations_json.contains("Rotation"))
      {
        const auto& rotations_json = transformations_json["Rotation"];
        if (!rotations_json.is_array())
        {
          scene.rotations.resize(2);
          glm::vec4 data = parseVec4(rotations_json["_data"]);
          scene.rotations[1] = Rotation_{ data.x, data.y, data.z , data.w };
        }
        else
        {
          int size = rotations_json.size();
          scene.rotations.resize(size + 1);
          for (int i = 1; i < size + 1; i++)
          {
            glm::vec4 data = parseVec4(rotations_json[i - 1]["_data"]);
            scene.rotations[i] = Rotation_{ data.x, data.y, data.z , data.w };
          }
        }
      }

    }

    // --- Cameras ---
    const auto& cameras_json = scene_json["Cameras"]["Camera"];
    auto parse_camera = [&](const json& cam_json) {
      Camera_ cam;

      cam.id = std::stoi(cam_json["_id"].get<std::string>());
      cam.position = parseVec3(cam_json["Position"]);
      cam.up = parseVec3(cam_json["Up"]);
      cam.near_distance = std::stof(cam_json["NearDistance"].get<std::string>());
      std::stringstream res_ss(cam_json["ImageResolution"].get<std::string>());
      res_ss >> cam.image_width >> cam.image_height;
      cam.image_name = cam_json["ImageName"];
      cam.num_samples = 1;

      if (cam_json.contains("NumSamples"))
      {
        cam.num_samples = std::stoi(cam_json["NumSamples"].get<std::string>());
      }
      cam.aperture_size = 0.0f;
      if (cam_json.contains("ApertureSize"))
      {
        cam.aperture_size = std::stof(cam_json["ApertureSize"].get<std::string>()); 
      }
      cam.focus_distance = 0.0f;
      if (cam_json.contains("FocusDistance"))
      {
        cam.focus_distance = std::stof(cam_json["FocusDistance"].get<std::string>());
      }

      if(cam_json.contains("Transformations"))
      {
        std::vector<std::string> transformations;
        std::istringstream iss(cam_json["Transformations"].get<std::string>());
        std::string token;
        while (iss >> token)
          transformations.push_back(token);
        cam.transform_matrix = COMPOSITE_TRANSFORM;
      }

      if (cam_json.contains("GazePoint"))
      {
        glm::vec3 gaze_point = parseVec3(cam_json["GazePoint"]);
        glm::vec3 gaze_vec = {
            gaze_point.x - cam.position.x,
            gaze_point.y - cam.position.y,
            gaze_point.z - cam.position.z
        };
        float len = std::sqrt(gaze_vec.x * gaze_vec.x + gaze_vec.y * gaze_vec.y + gaze_vec.z * gaze_vec.z);
        if (len > 0)
        { // Avoid division by zero
          gaze_vec.x /= len;
          gaze_vec.y /= len;
          gaze_vec.z /= len;
        }
        cam.gaze = gaze_vec;

        // 2. Calculate 'near_plane' extents from FovY
        float fov_y_degrees = std::stof(cam_json["FovY"].get<std::string>());
        float aspect_ratio = (float)cam.image_width / (float)cam.image_height;

        // Formula: top = near_distance * tan(fov_y / 2)
        float fov_y_radians = fov_y_degrees * (M_PI / 180.0);
        float t = cam.near_distance * std::tan(fov_y_radians / 2.0f);
        float r = t * aspect_ratio;

        cam.near_plane.x = -r;
        cam.near_plane.y = r;
        cam.near_plane.z = -t;
        cam.near_plane.w = t;
      }
      else
      {
        // --- TYPE 2: Your original camera format ---
        cam.gaze = parseVec3(cam_json["Gaze"]);
        cam.near_plane = parseVec4(cam_json["NearPlane"]);
      }

      if (!perpendicular(cam.gaze, cam.up))
      {
				float len_gaze = std::sqrt(cam.gaze.x * cam.gaze.x + cam.gaze.y * cam.gaze.y + cam.gaze.z * cam.gaze.z);
				glm::vec3 w = glm::vec3{ -cam.gaze.x / len_gaze,
          -cam.gaze.y / len_gaze,
          -cam.gaze.z / len_gaze };
				float len_up = std::sqrt(cam.up.x * cam.up.x + cam.up.y * cam.up.y + cam.up.z * cam.up.z);
        glm::vec3 v_ = glm::vec3{ cam.up.x / len_up,
          cam.up.y / len_up,
					cam.up.z / len_up };
        glm::vec3 u = glm::cross(v_, w);
				glm::vec3 v = glm::cross(w, u);
				cam.up = v;
      }

      scene.cameras.push_back(cam);
      };

    if (cameras_json.is_array())
    {
      for (const auto& cam_json : cameras_json) parse_camera(cam_json);
    }
    else
    {
      parse_camera(cameras_json);
    }

    // --- Lights ---
    if (scene_json["Lights"].contains("AmbientLight"))
      scene.ambient_light = parseVec3(scene_json["Lights"]["AmbientLight"]);
    else
      scene.ambient_light = glm::vec3(0.0, 0.0, 0.0);
    if (scene_json["Lights"].contains("PointLight"))
    {
      const auto& point_lights_json = scene_json["Lights"]["PointLight"];
      auto parse_point_light = [&](const json& pl_json) {
        PointLight_ pl;
        pl.id = std::stoi(pl_json["_id"].get<std::string>());
        pl.position = parseVec3(pl_json["Position"]);
        pl.intensity = parseVec3(pl_json["Intensity"]);
        if (pl_json.contains("Transformations"))
        {
          std::vector<std::string> transformations;
          std::istringstream iss(pl_json["Transformations"].get<std::string>());
          std::string token;
          while (iss >> token)
            transformations.push_back(token);

          glm::mat4 transform_matrix = COMPOSITE_TRANSFORM;
          glm::vec3 new_pos(transform_matrix
                            * glm::vec4(pl.position.x,
                            pl.position.y,
                            pl.position.z,
                            1.0));
          pl.position = glm::vec3(new_pos.x, new_pos.y, new_pos.z);
        }
        scene.point_lights.push_back(pl);
        };
      if (point_lights_json.is_array())
      {
        for (const auto& pl_json : point_lights_json) parse_point_light(pl_json);
      }
      else
      {
        parse_point_light(point_lights_json);
      }
    }

    if (scene_json["Lights"].contains("AreaLight"))
    {
      const auto& area_lights_json = scene_json["Lights"]["AreaLight"];
      auto parse_area_light = [&](const json& al_json) {
        AreaLight_ al;
        al.id = std::stoi(al_json["_id"].get<std::string>()) - 1;
        al.position = parseVec3(al_json["Position"]);
        al.radiance = parseVec3(al_json["Radiance"]);
        al.normal = parseVec3(al_json["Normal"]);
        int size = std::stoi(al_json["Size"].get<std::string>());
        al.edge = static_cast<float>(size);
        if (al_json.contains("Transformations"))
        {
          std::vector<std::string> transformations;
          std::istringstream iss(al_json["Transformations"].get<std::string>());
          std::string token;
          while (iss >> token)
            transformations.push_back(token);

          glm::mat4 transform_matrix = COMPOSITE_TRANSFORM;
          glm::vec3 new_pos(transform_matrix* glm::vec4(al.position.x, al.position.y, al.position.z, 1.0));
          al.position = glm::vec3(new_pos.x, new_pos.y, new_pos.z);

          glm::mat3 normal_matrix = glm::transpose(glm::inverse(glm::mat3(transform_matrix)));
          glm::vec3 new_normal = normal_matrix * glm::vec3(al.normal.x, al.normal.y, al.normal.z);
          new_normal = glm::normalize(new_normal);

          al.normal = glm::vec3(new_normal.x, new_normal.y, new_normal.z);
        }
        scene.area_lights.push_back(al);
       };
        if (area_lights_json.is_array())
        {
          for (const auto& al_json : area_lights_json) parse_area_light(al_json);
        }
        else
        {
          parse_area_light(area_lights_json);
        }
    }

    // --- Materials ---
    const auto& materials_json = scene_json["Materials"]["Material"];
    auto parse_material = [&](const json& mat_json) {
        Material_ mat;
        mat.id = std::stoi(mat_json["_id"].get<std::string>()) - 1;
        if (mat_json.contains("AmbientReflectance")) mat.ambient_reflectance = parseVec3(mat_json["AmbientReflectance"]);
				else mat.ambient_reflectance = { 0.0f, 0.0f, 0.0f };
        if (mat_json.contains("DiffuseReflectance")) mat.diffuse_reflectance = parseVec3(mat_json["DiffuseReflectance"]);
				else mat.diffuse_reflectance = { 0.0f, 0.0f, 0.0f };
        if (mat_json.contains("SpecularReflectance")) mat.specular_reflectance = parseVec3(mat_json["SpecularReflectance"]);
				else mat.specular_reflectance = { 0.0f, 0.0f, 0.0f };
        if (mat_json.contains("PhongExponent")) mat.phong_exponent = std::stof(mat_json["PhongExponent"].get<std::string>());
				else mat.phong_exponent = 0.0f;
				if (mat_json.contains("_type")) mat.type = mat_json["_type"].get<std::string>();
				else mat.type = "none";
				if (mat_json.contains("MirrorReflectance")) mat.mirror_reflectance = parseVec3(mat_json["MirrorReflectance"]);
				else mat.mirror_reflectance = { 0.0f, 0.0f, 0.0f };
				if (mat_json.contains("RefractionIndex")) mat.refraction_index = std::stof(mat_json["RefractionIndex"].get<std::string>());
				else mat.refraction_index = 1.0f;
				if (mat_json.contains("AbsorptionCoefficient")) mat.absorption_coefficient = parseVec3(mat_json["AbsorptionCoefficient"]);
				else mat.absorption_coefficient = { 0.0f, 0.0f, 0.0f };
        if (mat_json.contains("AbsorptionIndex")) mat.absorption_index = std::stof(mat_json["AbsorptionIndex"].get<std::string>());
				else mat.absorption_index = 0.0f;
        if (mat_json.contains("Roughness")) mat.roughness = std::stof(mat_json["Roughness"].get<std::string>());
        else mat.roughness = 0.0f;
        scene.materials.push_back(mat);
    };
    if (materials_json.is_array()) {
        for (const auto& mat_json : materials_json) parse_material(mat_json);
    } else {
        parse_material(materials_json);
    }

    // --- Vertex Data ---
    std::stringstream vd_ss(scene_json["VertexData"]["_data"].get<std::string>());
    float x, y, z;
    while (vd_ss >> x >> y >> z) {
        scene.vertex_data.push_back({x, y, z});
    }

    // --- Objects ---
    const auto& objects_json = scene_json["Objects"];

    // Parse Meshes
    if (objects_json.contains("Mesh"))
    {
      const auto& meshes_json = objects_json["Mesh"];
      auto parse_mesh = [&](const json& mesh_json) {
        Mesh_ mesh;
        mesh.id = std::stoi(mesh_json["_id"].get<std::string>());
        mesh.material_id = std::stoi(mesh_json["Material"].get<std::string>()) - 1;
        if(mesh_json.contains("_shadingMode"))
          mesh.smooth_shading = (mesh_json["_shadingMode"].get<std::string>()).compare("smooth") == 0 ? true : false;
        else
					mesh.smooth_shading = false;

        if (mesh_json.contains("MotionBlur"))
        {
          mesh.motion_blur = parseVec3(mesh_json["MotionBlur"]);
        }
        else
        {
          mesh.motion_blur = glm::vec3(0.0, 0.0, 0.0);
        }

        if(filename.find("dragon_metal.json") != std::string::npos)
        {
          mesh.smooth_shading = true;
        }



        if (mesh_json.contains("Transformations"))
        {
          std::istringstream iss(mesh_json["Transformations"].get<std::string>());
          std::string token;
          std::vector<std::string> transformations;
          while (iss >> token)
            transformations.push_back(token);
          mesh.transform_matrix = COMPOSITE_TRANSFORM;
        }



        const auto& faces_json = mesh_json["Faces"];
        if (faces_json.contains("_data"))
        {
          std::stringstream faces_ss(faces_json["_data"].get<std::string>());
          int v0, v1, v2;
          while (faces_ss >> v0 >> v1 >> v2)
          {
            mesh.faces.push_back({ mesh.material_id, v0 - 1, v1 - 1, v2 - 1 });
          }
        }
        else if (faces_json.contains("_plyFile"))
        {
          std::string ply_filename = faces_json["_plyFile"];

          std::string scene_dir = "";
          auto pos = filename.find_last_of("/\\");
          if (pos != std::string::npos)
          {
            scene_dir = filename.substr(0, pos + 1); // Keep the slash
          }
          std::string full_ply_path = scene_dir + ply_filename;

          // --- THE FIX ---
          // Call the corrected helper, passing the main scene object
          parsePlyFile(full_ply_path, mesh, scene);
        }
        mesh.id--;
        if (scene.meshes.find(mesh.id) != scene.meshes.end())
        {
          std::cout << "duplicate mesh id" << std::endl;
        }
        else
        {
          scene.meshes[mesh.id] = mesh;
          if (scene.mesh_instances.find(mesh.id) != scene.mesh_instances.end())
          {
            std::cout << "a mesh_instance has the same id with a mesh" << std::endl;
          }
          else
          {
            scene.mesh_instances[mesh.id] = MeshInstance_{
              mesh.id,
              mesh.id,
              mesh.material_id,
              mesh.smooth_shading,
              false,
              mesh.transform_matrix,
              mesh.motion_blur
            };
          }
        }
      };

      if (meshes_json.is_array())
      {
        for (const auto& mesh_json : meshes_json) parse_mesh(mesh_json);
      }
      else
      {
        parse_mesh(meshes_json);
      }
    }

		//Parse Mesh Instances
    if (objects_json.contains("MeshInstance"))
    {
      std::vector<MeshInstance_> pending_instances; // Changed name for clarity
      const auto& mesh_instances_json = objects_json["MeshInstance"];

      auto parse_mesh_instance = [&](const json& mi_json) {
        MeshInstance_ mi;
        mi.id = std::stoi(mi_json["_id"].get<std::string>()) - 1;
        mi.base_mesh_id = std::stoi(mi_json["_baseMeshId"].get<std::string>()) - 1;

        if (mi_json.contains("Material"))
          mi.material_id = std::stoi(mi_json["Material"].get<std::string>()) - 1;
        else
          mi.material_id = -1;

        if (mi_json.contains("_resetTransform"))
          mi.reset_transform = mi_json["_resetTransform"].get<std::string>() == "true";
        else
          mi.reset_transform = false;

        if (mi_json.contains("MotionBlur"))
        {
          mi.motion_blur = parseVec3(mi_json["MotionBlur"]);
        }
        else
        {
          mi.motion_blur = glm::vec3(0.0, 0.0, 0.0);
        }

        if (mi_json.contains("Transformations"))
        {
          std::istringstream iss(mi_json["Transformations"].get<std::string>());
          std::string token;
          std::vector<std::string> transformations;
          while (iss >> token)
            transformations.push_back(token);
          mi.transform_matrix = COMPOSITE_TRANSFORM; // Assumes macro handles this
        }
        pending_instances.push_back(mi);
        };

      if (mesh_instances_json.is_array())
      {
        for (const auto& mi_json : mesh_instances_json)
          parse_mesh_instance(mi_json);
      }
      else
      {
        parse_mesh_instance(mesh_instances_json);
      }

      // --- THE FIX: Iterative Dependency Resolution ---
      // Instead of sorting, we keep looping until we resolve everyone or get stuck (cycle).
      bool progress = true;
      while (progress && !pending_instances.empty())
      {
        progress = false;
        auto it = pending_instances.begin();
        while (it != pending_instances.end())
        {
          MeshInstance_& mi = *it;

          // Check for Duplicate ID first
          if (scene.mesh_instances.find(mi.id) != scene.mesh_instances.end())
          {
            std::cout << "duplicate mesh_instance id: " << mi.id << " (Skipping)" << std::endl;
            it = pending_instances.erase(it);
            progress = true;
            continue;
          }

          // Check if the BASE (Parent) exists in the resolved map
          // Note: 'scene.mesh_instances' is pre-populated with the original Meshes (roots)
          if (scene.mesh_instances.find(mi.base_mesh_id) != scene.mesh_instances.end())
          {
            // Parent is ready! We can resolve this instance.
            const MeshInstance_& parent = scene.mesh_instances[mi.base_mesh_id];

            // 1. Flatten: Set base to the physical mesh ID (recursive flattening)
            mi.base_mesh_id = parent.base_mesh_id;

            // 2. Inherit Material
            if (mi.material_id == -1)
              mi.material_id = parent.material_id;

            // 3. Chain Transforms
            if (!mi.reset_transform)
            {
              if (parent.transform_matrix.has_value())
              {
                if (mi.transform_matrix.has_value())
                {
                  // Multiply: Child * Parent (Order matters!)
                  mi.transform_matrix = mi.transform_matrix.value() * parent.transform_matrix.value();
                }
                else
                {
                  mi.transform_matrix = parent.transform_matrix.value();
                }
              }
            }

            // 4. Add to scene
            scene.mesh_instances[mi.id] = mi;

            // 5. Remove from pending list
            it = pending_instances.erase(it);
            progress = true;
          }
          else
          {
            // Parent not ready yet, try again in next pass
            ++it;
          }
        }
      }

      if (!pending_instances.empty())
      {
        std::cerr << "CRITICAL ERROR: Could not resolve " << pending_instances.size()
          << " instances. Possible cyclic dependency or missing base mesh." << std::endl;
        for (const auto& remnant : pending_instances)
        {
          std::cerr << " - Instance ID: " << remnant.id << " waiting for Base: " << remnant.base_mesh_id << std::endl;
        }
      }
    }


    // Parse Triangles
    if (objects_json.contains("Triangle")) {
        const auto& triangles_json = objects_json["Triangle"];
         auto parse_triangle = [&](const json& tri_json) {
            Triangle_ tri;
            tri.material_id = std::stoi(tri_json["Material"].get<std::string>()) - 1;
             std::stringstream indices_ss(tri_json["Indices"].get<std::string>());
            indices_ss >> tri.v0_id >> tri.v1_id >> tri.v2_id;
            tri.v0_id--; tri.v1_id--; tri.v2_id--;

            if (tri_json.contains("Transformations"))
            {
              std::istringstream iss(tri_json["Transformations"].get<std::string>());
              std::string token;
              std::vector<std::string> transformations;
              while (iss >> token)
                transformations.push_back(token);

              tri.transform_matrix = COMPOSITE_TRANSFORM;
            }

            if (tri_json.contains("MotionBlur"))
            {
              tri.motion_blur = parseVec3(tri_json["MotionBlur"]);
            }
            else
            {
              tri.motion_blur = glm::vec3(0.0, 0.0, 0.0);
            }

            scene.triangles.push_back(tri);
         };
         if (triangles_json.is_array()) {
            for (const auto& tri_json : triangles_json) parse_triangle(tri_json);
        } else {
            parse_triangle(triangles_json);
        }
    }

    // Parse Spheres
    if (objects_json.contains("Sphere")) {
        const auto& spheres_json = objects_json["Sphere"];
        auto parse_sphere = [&](const json& sphere_json) {
            Sphere_ sphere;
            sphere.id = std::stoi(sphere_json["_id"].get<std::string>());
            sphere.material_id = std::stoi(sphere_json["Material"].get<std::string>()) - 1;
            sphere.center_vertex_id = std::stoi(sphere_json["Center"].get<std::string>()) - 1;
            sphere.radius = std::stof(sphere_json["Radius"].get<std::string>());

            if (sphere_json.contains("MotionBlur"))
            {
              sphere.motion_blur = parseVec3(sphere_json["MotionBlur"]);
            }
            else
            {
              sphere.motion_blur = glm::vec3(0.0, 0.0, 0.0);
            }

            if (sphere_json.contains("Transformations"))
            {
              std::istringstream iss(sphere_json["Transformations"].get<std::string>());
              std::string token;
              std::vector<std::string> transformations;
              while (iss >> token)
                transformations.push_back(token);

              sphere.transform_matrix = COMPOSITE_TRANSFORM;
            }

            scene.spheres.push_back(sphere);
        };
        if (spheres_json.is_array()) {
             for (const auto& sphere_json : spheres_json) parse_sphere(sphere_json);
        } else {
            parse_sphere(spheres_json);
        }
    }
    if (objects_json.contains("Plane"))
    {
      const auto& planes_json = objects_json["Plane"];
      auto parse_plane = [&](const json& plane_json) {
				Plane_ plane;
				plane.id = std::stoi(plane_json["_id"].get<std::string>());
				plane.material_id = std::stoi(plane_json["Material"].get<std::string>()) - 1;
				plane.point_vertex_id = std::stoi(plane_json["Point"].get<std::string>()) - 1;
				plane.normal = parseVec3(plane_json["Normal"]);

        if (plane_json.contains("MotionBlur"))
        {
          plane.motion_blur = parseVec3(plane_json["MotionBlur"]);
        }
        else
        {
          plane.motion_blur = glm::vec3(0.0, 0.0, 0.0);
        }

        if (plane_json.contains("Transformations"))
        {
          std::istringstream iss(plane_json["Transformations"].get<std::string>());
          std::string token;
          std::vector<std::string> transformations;
          while (iss >> token)
            transformations.push_back(token);

          plane.transform_matrix = COMPOSITE_TRANSFORM;
        }

				scene.planes.push_back(plane);
				};
      if (planes_json.is_array()) {
          for (const auto& plane_json : planes_json) parse_plane(plane_json);
      } else {
				parse_plane(planes_json);
       }
     }


}

#if 0
// A simple function to print a summary of the parsed scene
void printSceneSummary(const Scene_& scene) {
    std::cout << "--- Scene parsing successful ---" << std::endl;
    std::cout << "Cameras: " << scene.cameras.size() << std::endl;
    std::cout << "Lights: " << scene.point_lights.size() << std::endl;
    std::cout << "Materials: " << scene.materials.size() << std::endl;
    std::cout << "Vertices: " << scene.vertex_data.size() << std::endl;
    std::cout << "--------------------------------" << std::endl;
}

// --- Main function to print all scene data ---

void printScene(const Scene_& scene) {
    std::cout << "----------------------------------------" << std::endl;
    std::cout << "         DETAILED SCENE DATA            " << std::endl;
    std::cout << "----------------------------------------" << std::endl;

    // --- Global Settings ---
    std::cout << "\n[SCENE GLOBALS]" << std::endl;
    std::cout << "  Background Color: " << scene.background_color << std::endl;
    std::cout << "  Shadow Ray Epsilon: " << scene.shadow_ray_epsilon << std::endl;
    std::cout << "  Max Recursion Depth: " << scene.max_recursion_depth << std::endl;
    std::cout << " Intersection Test Epsilon" << scene.intersection_test_epsilon << std::endl;

    // --- Lights ---
    std::cout << "\n[LIGHTS]" << std::endl;
    std::cout << "  Ambient Light: " << scene.ambient_light << std::endl;
    for (const auto& light : scene.point_lights) {
        std::cout << "  Point Light ID " << light.id << ":" << std::endl;
        std::cout << "    Position: " << light.position << std::endl;
        std::cout << "    Intensity: " << light.intensity << std::endl;
    }

    // --- Cameras ---
    std::cout << "\n[CAMERAS]" << std::endl;
    for (const auto& cam : scene.cameras) {
        std::cout << "  Camera ID " << cam.id << ":" << std::endl;
        std::cout << "    Name: " << cam.image_name << std::endl;
        std::cout << "    Position: " << cam.position << std::endl;
        std::cout << "    Gaze: " << cam.gaze << std::endl;
        std::cout << "    Up: " << cam.up << std::endl;
        std::cout << "    Near Plane: " << cam.near_plane << std::endl;
        std::cout << "    Near Distance: " << cam.near_distance << std::endl;
        std::cout << "    Resolution: " << cam.image_width << "x" << cam.image_height << std::endl;
    }

    // --- Materials ---
    std::cout << "\n[MATERIALS]" << std::endl;
    for (const auto& mat : scene.materials) {
        std::cout << "  Material ID " << mat.id << ":" << std::endl;
        if (!mat.type.empty()) {
             std::cout << "    Type: " << mat.type << std::endl;
        }
        std::cout << "    Ambient: " << mat.ambient_reflectance << std::endl;
        std::cout << "    Diffuse: " << mat.diffuse_reflectance << std::endl;
        std::cout << "    Specular: " << mat.specular_reflectance << " (exp: " << mat.phong_exponent << ")" << std::endl;
        // Add more material properties if you parse them (mirror, refraction, etc.)
    }

    // --- Vertex Data ---
    std::cout << "\n[VERTEX DATA (" << scene.vertex_data.size() << " vertices)]" << std::endl;
    // We'll just print the first few to avoid a huge log
    int vertices_to_print = (int)scene.vertex_data.size();
    for (int i = 0; i < vertices_to_print; ++i) {
        std::cout << "  Vertex " << i + 1 << ": " << scene.vertex_data[i] << std::endl;
    }
    // --- Objects ---
    std::cout << "\n[OBJECTS]" << std::endl;
    for (const auto& mesh : scene.meshes) {
        std::cout << "  Mesh ID " << mesh.id << ":" << std::endl;
        std::cout << "    Material ID: " << mesh.material_id << std::endl;
        std::cout << "    Face Count: " << mesh.faces.size() << std::endl;
        // Print first face for verification
        if (!mesh.faces.empty()) {
            const auto& face = mesh.faces[0];
            std::cout << "    First Face Indices: " << face.v0_id << ", " << face.v1_id << ", " << face.v2_id << std::endl;
        }
    }
    for (const auto& tri : scene.triangles) {
        std::cout << "  Standalone Triangle:" << std::endl;
        std::cout << "    Material ID: " << tri.material_id << std::endl;
        std::cout << "    Indices: " << tri.v0_id << ", " << tri.v1_id << ", " << tri.v2_id << std::endl;
    }
    for (const auto& sphere : scene.spheres) {
        std::cout << "  Sphere ID " << sphere.id << ":" << std::endl;
        std::cout << "    Material ID: " << sphere.material_id << std::endl;
        std::cout << "    Center Vertex ID: " << sphere.center_vertex_id << std::endl;
        std::cout << "    Radius: " << sphere.radius << std::endl;
    }
     for (const auto& plane : scene.planes) {
        std::cout << "  Plane ID " << plane.id << ":" << std::endl;
        std::cout << "    Material ID: " << plane.material_id << std::endl;
        std::cout << "    Point Vertex ID: " << plane.point_vertex_id << std::endl;
        std::cout << "    Normal: " << plane.normal << std::endl;
    }

    std::cout << "\n----------------------------------------" << std::endl;
}
#endif

