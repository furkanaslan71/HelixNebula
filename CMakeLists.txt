# CMakeList.txt : CMake project for HelixNebula
#
# REQUIREMENT: CMake 3.25 or higher is needed for native C++23 support
cmake_minimum_required (VERSION 3.25)



# Enable Hot Reload for MSVC compilers if supported.
if (POLICY CMP0141)
  cmake_policy(SET CMP0141 NEW)
  set(CMAKE_MSVC_DEBUG_INFORMATION_FORMAT "$<IF:$<AND:$<C_COMPILER_ID:MSVC>,$<CXX_COMPILER_ID:MSVC>>,$<$<CONFIG:Debug,RelWithDebInfo>:EditAndContinue>,$<$<CONFIG:Debug,RelWithDebInfo>:ProgramDatabase>>")
endif()

project ("HelixNebula")

# --- C++23 CONFIGURATION ---
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
# Disabling extensions ensures /permissive- on MSVC, improving IntelliSense accuracy
set(CMAKE_CXX_EXTENSIONS OFF) 

# Add source to this project's executable. 
# ONLY .cpp (source) files should be listed here, NOT .h (header) files.
add_executable (HelixNebula
    # The file that caused the 'Cannot find source file' error:
    # Assuming the file is actually in 'raytracer/src/parser.cpp' or was simply mistyped.
    # I'll assume it should be in 'src' like other .cpp files unless proven otherwise.
    "raytracer/parser/parser.cpp" 
    "raytracer/main.cpp"
    "raytracer/core/aabb.cpp"
        "raytracer/scene/scene.cpp"
    "raytracer/material/material.cpp" 
    "raytracer/objects/plane.cpp" 
    "raytracer/objects/tlas_box.cpp" 
    "raytracer/objects/mesh.cpp"
    "raytracer/camera/base_camera.cpp"
        "raytracer/core/ray.cpp"
    "raytracer/math_core/random/my_random.cpp"
        "raytracer/math_core/transformations/transformations.h"
    "raytracer/math_core/transformations/transformations.cpp" 
    "raytracer/math_core/vector_utility.h" 
    "raytracer/objects/geometry.h" 
    "raytracer/objects/geometry_concepts.h" "raytracer/glm_config.h"
        "raytracer/texture_mapping/perlin.h"
        raytracer/render/raytracer.h
        raytracer/render/raytracer.cpp
        raytracer/io/save_image.h
        raytracer/io/save_image.cpp
        raytracer/texture_mapping/texture_data.cpp
        raytracer/texture_mapping/texture_lookup.h
        raytracer/texture_mapping/texture_lookup.cpp)


# Important: Add include directories so your compiler can find your headers.
# This tells the compiler where to look for files like "vec3.h", "Ray.h", etc.
target_include_directories(HelixNebula PRIVATE
    "${CMAKE_SOURCE_DIR}/raytracer"
)