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
    "raytracer/src/parser.cpp" 
    
    # Core Source Files
    "raytracer/src/main.cpp"
    "raytracer/src/aabb.cpp" 
    "raytracer/src/bvh.cpp" 

    # Other Source Files
    "raytracer/render/render_manager.cpp"
    "raytracer/material/material_manager.cpp" 
    "raytracer/render/base_ray_tracer.cpp" 
    "raytracer/scene/scene.cpp" 
    "raytracer/material/material.cpp" 
    "raytracer/objects/plane.cpp" 
    "raytracer/objects/tlas_box.cpp" 
    "raytracer/objects/mesh.cpp" 
    "raytracer/camera/base_camera.cpp" 
    "raytracer/camera/distribution_camera.cpp" 
    "raytracer/src/ray.cpp" # Note: Usually ray.cpp would be in 'src'
    "raytracer/random/my_random.cpp"

    # Fixed typo in path: removed double slash
    "raytracer/camera/pinhole_camera.cpp" 
    
    # --- IMPORTANT ---
    # The following header files were removed from the executable list:
    # "raytracer/objects/sphere.h" 
    # "raytracer/objects/plane.h" 
    # "raytracer/include/config.h" 
    # "raytracer/objects/mesh.h" 
    # "raytracer/include/vec3.h"
    # "raytracer/camera/pinhole_camera.h" 
    # "raytracer/camera/distribution_camera.h" 
    # "raytracer/include/ray.h"
)


# Important: Add include directories so your compiler can find your headers.
# This tells the compiler where to look for files like "vec3.h", "Ray.h", etc.
target_include_directories(HelixNebula PRIVATE
    "${CMAKE_SOURCE_DIR}/raytracer/include"
    "${CMAKE_SOURCE_DIR}/raytracer/io"
    "${CMAKE_SOURCE_DIR}/raytracer/light"
    "${CMAKE_SOURCE_DIR}/raytracer/material"
    "${CMAKE_SOURCE_DIR}/raytracer/render"
    "${CMAKE_SOURCE_DIR}/raytracer/objects"
    "${CMAKE_SOURCE_DIR}/raytracer/src"
    "${CMAKE_SOURCE_DIR}/raytracer/random"
    "${CMAKE_SOURCE_DIR}/raytracer/camera"
    "${CMAKE_SOURCE_DIR}/raytracer/external/glm"
)