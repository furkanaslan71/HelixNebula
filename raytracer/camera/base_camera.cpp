#include "base_camera.h"

BaseCamera::BaseCamera(
	const Camera_& cam,
	const std::vector<Translation_>& translations,
	const std::vector<Scaling_>& scalings,
	const std::vector<Rotation_>& rotations,
	int _recursion_depth,
	int _num_area_lights,
	std::vector<AreaLight>& _area_lights,
	std::vector<Tonemap_>& _tonemaps
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

	gaze = glm::normalize(gaze);

	if (cam.transform_matrix.has_value())
	{
		glm::mat4 composite_transformation_matrix = cam.transform_matrix.value();

		glm::vec4 transformed_position = composite_transformation_matrix * glm::vec4(position.x, position.y, position.z, 1.0f);
		glm::vec4 transformed_gaze = composite_transformation_matrix * glm::vec4(gaze.x, gaze.y, gaze.z, 0.0f);
		glm::vec4 transformed_up = composite_transformation_matrix * glm::vec4(up.x, up.y, up.z, 0.0f);

		position = glm::vec3(transformed_position.x, transformed_position.y, transformed_position.z);
		gaze = glm::vec3(transformed_gaze.x, transformed_gaze.y, transformed_gaze.z);
		up = glm::vec3(transformed_up.x, transformed_up.y, transformed_up.z);
	}

	
	double zoom_factor = glm::length(gaze);
	if (zoom_factor < 1e-6) zoom_factor = 1.0;

	// 5. Create a new, clean orthonormal basis to prevent distortion.
	//    This corrects any skew introduced by non-uniform scaling.
	gaze = glm::normalize(gaze);
	this->w = -gaze;
	this->u = glm::normalize(glm::cross(up, w));
	this->v = glm::cross(w, u); // Already normalized because w and u are orthonormal

	// 6. Proceed with calculating the view plane geometry.
	near_plane[0] = cam.near_plane.x;
	near_plane[1] = cam.near_plane.y;
	near_plane[2] = cam.near_plane.z;
	near_plane[3] = cam.near_plane.w;

	// Divide the plane dimensions by the zoom factor to narrow the field of view.
	double l = near_plane[0] / zoom_factor;
	double r = near_plane[1] / zoom_factor;
	double b = near_plane[2] / zoom_factor;
	double t = near_plane[3] / zoom_factor;

	gaze = glm::normalize(gaze);
	glm::vec3 m = position + (gaze * (float)near_distance);
	q = m + u * (float)l + v * (float)t;
	su = u * (float)((r - l) / image_width);
	sv = v * (float)((b - t) / image_height);


	this->context.forward = -this->w;
	this->context.right = this->u;
	this->context.up = this->v;
	this->context.tan_half_fov_x = (float)((r - l) * 0.5 / near_distance);
	this->context.tan_half_fov_y = (float)((t - b) * 0.5 / near_distance);

	for (const auto& tm : _tonemaps)
	{
		tonemaps.emplace_back(tm);
	}
}

void BaseCamera::generatePixelSamples(int i, int j, std::vector<glm::vec3>& out_samples) const
{
	std::vector<std::pair<float, float>> samples
		= generateJitteredSamples(num_samples);

	out_samples.clear();

	for (const auto& [x, y] : samples)
	{
		glm::vec3 pixel_sample = q + su * (j + x) + sv * (i + y);
		out_samples.emplace_back(pixel_sample);
	}
}