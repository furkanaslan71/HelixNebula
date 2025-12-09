#ifndef PINHOLE_CAMERA_H
#define PINHOLE_CAMERA_H

#include "base_camera.h"

class PinholeCamera : public BaseCamera {
public:
  PinholeCamera(
    const Camera_& cam,
    const std::vector<Translation_>& translations,
    const std::vector<Scaling_>& scalings,
    const std::vector<Rotation_>& rotations,
    int _recursion_depth,
    int _num_area_lights,
    std::vector<AreaLight>& _area_lights
  ) : BaseCamera(cam, translations, scalings, rotations, _recursion_depth,
    _num_area_lights, _area_lights)
  {}

	void render(const BaseRayTracer& rendering_technique,
		std::vector<std::vector<Color>>& image) const override;
};



#endif