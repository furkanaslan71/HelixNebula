#ifndef PARSER_H
#define PARSER_H
#include <string>
#include <vector>
#include <ostream>
#include <algorithm>
#include <optional>
#include <unordered_map>

#include "config.h"
#include "glm_config.h"

struct Tonemap_ {
    std::string TMO;
    glm::vec2 TMOOptions;
    float saturation;
    float gamma;
    std::string extension;
};

struct BRDF_ {
    std::string type;
    int id;
    bool normalized = false;
    float exponent;
    bool kdfresnel = false;
};

typedef struct Camera_ {
    int id;
    glm::vec3 position;
    glm::vec3 gaze;
    glm::vec3 up;
    glm::vec4 near_plane;
    float near_distance;
    int image_width;
    int image_height;
    std::string image_name;
    //std::vector<std::string> transformations;
    std::optional<glm::mat4> transform_matrix;
    int num_samples;
    float aperture_size;
    float focus_distance;
    std::vector<Tonemap_> tone_maps;
    bool flip_x;
    std::string renderer;
    std::vector<std::string> renderer_params;
} Camera_;

typedef struct PointLight_ {
    int id;
    glm::vec3 position;
    glm::vec3 intensity;
		//std::vector<std::string> transformations;
} PointLight_;

struct AreaLight_ {
  int id;
  glm::vec3 position;
  glm::vec3 normal;
  float edge;
  glm::vec3 radiance;
  //std::vector<std::string> transformations;
};

struct DirectionalLight_ {
    int id;
    glm::vec3 direction;
    glm::vec3 radiance;
};

struct SpotLight_ {
    int id;
    glm::vec3 position;
    glm::vec3 direction;
    glm::vec3 intensity;
    float coverage_angle;
    float falloff_angle;
};

struct SphericalDirectionalLight_ {
    int id;
    std::string type;
    int image_id;
    std::string sampler;
};

typedef struct Material_ {
    int id;
    std::string type;
    glm::vec3 ambient_reflectance;
    glm::vec3 diffuse_reflectance;
    glm::vec3 specular_reflectance;
    glm::vec3 mirror_reflectance;
    float phong_exponent;
    float refraction_index;
    glm::vec3 absorption_coefficient;
    float absorption_index;
    float roughness;
    bool degamma;
    std::optional<int> brdf_id;
} Material_;

typedef struct Triangle_ {
    int material_id;
    int v0_id, v1_id, v2_id;
    std::optional<glm::mat4> transform_matrix;
    glm::vec3 motion_blur;
    std::vector<int> textures;
    std::optional<glm::vec3> radiance;
} Triangle_;

typedef struct Mesh_ {
    int id;
    int material_id;
		bool smooth_shading;
    std::vector<Triangle_> faces;
    std::optional<glm::mat4> transform_matrix;
    glm::vec3 motion_blur;
    std::vector<int> textures;
    int vertex_offset;
    int texture_offset;
    std::optional<glm::vec3> radiance;
} Mesh_;

typedef struct MeshInstance_ {
    int id;
    int base_mesh_id;
	int material_id;
	bool smooth_shading;
	bool reset_transform;
    std::optional<glm::mat4> transform_matrix;
    glm::vec3 motion_blur;
    std::vector<int> textures;
    std::optional<glm::vec3> radiance;
}MeshInstance_;

typedef struct Sphere_ {
    int id;
    int material_id;
    int center_vertex_id;
    float radius;
    std::optional<glm::mat4> transform_matrix;
    glm::vec3 motion_blur;
    std::vector<int> textures;
    std::optional<glm::vec3> radiance;
} Sphere_;

typedef struct Plane_ {
    int id;
    int material_id;
    int point_vertex_id;
    glm::vec3 normal;
    std::optional<glm::mat4> transform_matrix;
    glm::vec3 motion_blur;
    std::vector<int> textures;
    std::optional<glm::vec3> radiance;
} Plane_;

typedef struct Translation_ {
   float tx, ty, tz;
}Translation_;

typedef struct Scaling_ {
  float sx, sy, sz;
}Scaling_;

typedef struct Rotation_ {
  float angle;
  float axis_x, axis_y, axis_z;
}Rotation_;

struct Image_ {
  std::string data;
  int id;
};

struct TextureMap_ {
  int id;
  std::string type;
  int image_id;
  std::string decal_mode;
  std::string interpolation;
  float bump_factor;
  float noise_scale;
  std::string noise_conversion;
  int num_octaves;
  float scale;
  float offset;
  glm::vec3 black_color;
  glm::vec3 white_color;
  float normalizer;
};
	

typedef struct Scene_ {
    glm::vec3 background_color;
    float shadow_ray_epsilon;
    float intersection_test_epsilon;
    int max_recursion_depth;         
    std::vector<Camera_> cameras;
    glm::vec3 ambient_light;
    std::vector<PointLight_> point_lights;
    std::vector<AreaLight_> area_lights;
    std::vector<DirectionalLight_> directional_lights;
    std::vector<SphericalDirectionalLight_> spherical_directional_lights;
    std::vector<SpotLight_> spot_lights;
    std::vector<Material_> materials;
    std::vector<glm::vec3> vertex_data;
		std::vector<Translation_> translations;
		std::vector<Scaling_> scalings;
		std::vector<Rotation_> rotations;
    std::vector<glm::mat4> composite_matrices;
    std::vector<Triangle_> triangles;
    std::vector<Sphere_> spheres;
    std::vector<Plane_> planes;

    std::unordered_map<int, Mesh_> meshes;
    std::unordered_map<int, MeshInstance_> mesh_instances;

    std::vector<glm::vec2> tex_coord_data;
    std::vector<Image_> images;
    std::vector<TextureMap_> texture_maps;

    std::vector<BRDF_> brdfs;
} Scene_;

// --- Function Declaration ---

void parseScene(const std::string& filename, Scene_& scene);

inline std::ostream& operator<<(std::ostream& os, const glm::vec3& v) {
    os << "(" << v.x << ", " << v.y << ", " << v.z << ")";
    return os;
}

inline std::ostream& operator<<(std::ostream& os, const glm::vec4& v) {
    os << "(l:" << v.x << ", r:" << v.y << ", b:" << v.z << ", t:" << v.w << ")";
    return os;
}

void printSceneSummary(const Scene_& scene);
void printScene(const Scene_& scene);

#endif //PARSER_H