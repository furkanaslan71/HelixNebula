#ifndef  DISTRIBUTION_CAMERA_H
#define DISTRIBUTION_CAMERA_H

#include "base_camera.h"

class DistributionCamera : public BaseCamera {
public:
  DistributionCamera(
    const Camera_& cam,
    const std::vector<Translation_>& translations,
    const std::vector<Scaling_>& scalings,
    const std::vector<Rotation_>& rotations,
    int _recursion_depth,
    int _num_area_lights,
    std::vector<AreaLight>& _area_lights,
    std::vector<Tonemap_> _tonemaps
  ) : BaseCamera(cam,
    translations, 
    scalings, 
    rotations, 
    _recursion_depth,
    _num_area_lights,
    _area_lights,
    _tonemaps
  ),
    aperture_size(cam.aperture_size),
    focus_distance(cam.focus_distance)
  {
  }

  float aperture_size;
  float focus_distance;


private:
  void generateApertureSamples(std::vector<glm::vec3>& out_samples) const;
  glm::vec3 calculateDir(const glm::vec3& pixel_sample, const glm::vec3& a) const;

};

#endif // ! DISTRIBUTION_CAMERA_H
