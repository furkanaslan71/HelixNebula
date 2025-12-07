#include "base_camera.h"

BaseCamera::BaseCamera(
	const Camera_& cam,
	const std::vector<Translation_>& translations,
	const std::vector<Scaling_>& scalings,
	const std::vector<Rotation_>& rotations,
	int _recursion_depth,
	int _num_area_lights,
	std::vector<AreaLight>& _area_lights
) : id(cam.id),
		position(cam.position.x, cam.position.y, cam.position.z),
		gaze(cam.gaze.x, cam.gaze.y, cam.gaze.z),
		up(cam.up.x, cam.up.y, cam.up.z),
		near_distance(cam.near_distance),
		image_width(cam.image_width),
		image_height(cam.image_height),
		image_name(cam.image_name),
		num_samples(cam.num_samples),
		recursion_depth(_recursion_depth),
		num_area_lights(_num_area_lights),
		area_lights(_area_lights)
{

	gaze.normalize();

	if (cam.transform_matrix.has_value())
	{
		glm::mat4 composite_transformation_matrix = cam.transform_matrix.value();

		glm::vec4 transformed_position = composite_transformation_matrix * glm::vec4(position.x, position.y, position.z, 1.0f);
		glm::vec4 transformed_gaze = composite_transformation_matrix * glm::vec4(gaze.x, gaze.y, gaze.z, 0.0f);
		glm::vec4 transformed_up = composite_transformation_matrix * glm::vec4(up.x, up.y, up.z, 0.0f);

		position = Vec3(transformed_position.x, transformed_position.y, transformed_position.z);
		gaze = Vec3(transformed_gaze.x, transformed_gaze.y, transformed_gaze.z);
		up = Vec3(transformed_up.x, transformed_up.y, transformed_up.z);
	}

	
	double zoom_factor = gaze.length();
	if (zoom_factor < 1e-6) zoom_factor = 1.0;

	// 5. Create a new, clean orthonormal basis to prevent distortion.
	//    This corrects any skew introduced by non-uniform scaling.
	this->w = gaze.normalize() * -1.0;
	this->u = up.cross(w).normalize();
	this->v = w.cross(u); // Already normalized because w and u are orthonormal

	// 6. Proceed with calculating the view plane geometry.
	near_plane[0] = cam.near_plane.l;
	near_plane[1] = cam.near_plane.r;
	near_plane[2] = cam.near_plane.b;
	near_plane[3] = cam.near_plane.t;

	// Divide the plane dimensions by the zoom factor to narrow the field of view.
	double l = near_plane[0] / zoom_factor;
	double r = near_plane[1] / zoom_factor;
	double b = near_plane[2] / zoom_factor;
	double t = near_plane[3] / zoom_factor;

	Vec3 m = position + (gaze.normalize() * near_distance);
	q = m + u * l + v * t;
	su = u * ((r - l) / image_width);
	sv = v * ((b - t) / image_height);
}

void BaseCamera::generatePixelSamples(int i, int j, std::vector<Vec3>& out_samples) const
{
	std::vector<std::pair<float, float>> samples
		= generateJitteredSamples(num_samples);

	out_samples.clear();

	for (const auto& [x, y] : samples)
	{
		Vec3 pixel_sample = q + su * (j + x) + sv * (i + y);
		out_samples.emplace_back(pixel_sample);
	}
}