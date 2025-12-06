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
    std::vector<AreaLight>& _area_lights
  ) : BaseCamera(cam,
    translations, 
    scalings, 
    rotations, 
    _recursion_depth,
    _num_area_lights,
    _area_lights
  ),
    aperture_size(cam.aperture_size),
    focus_distance(cam.focus_distance)
  {
  }

  float aperture_size;
  float focus_distance;

  void render(const BaseRayTracer& technique,
    std::vector<std::vector<Color>>& image) const override;

private:
  void generateApertureSamples(std::vector<Vec3>& out_samples) const;
  Vec3 calculateDir(const Vec3& pixel_sample, const Vec3& a) const;

};

#endif // ! DISTRIBUTION_CAMERA_H
