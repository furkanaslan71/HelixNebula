#include "transformations.h"


Mat4 calculateCompositeTransformationMatrix(
	const std::vector<std::string> transformations,
	const std::vector<Translation_>& translations,
	const std::vector<Scaling_>& scalings,
	const std::vector<Rotation_>& rotations)

{
	auto res = glm::mat4(1.0f);
	for (const auto& transform : transformations)
	{
		auto identity = glm::mat4(1.0f);
		switch (transform[0])
		{
		case 't': // translation
		{
			int index = std::stoi(transform.substr(1));
			auto& t = translations[index];
			glm::vec3 translation = glm::vec3(t.tx, t.ty, t.tz);
			res = glm::translate(identity, translation) * res;
			break;
		}
		case 's': // scaling
		{
			int index = std::stoi(transform.substr(1));
			auto& s = scalings[index];
			glm::vec3 scaling = glm::vec3(s.sx, s.sy, s.sz);
			res = glm::scale(identity, scaling) * res;
			break;
		}
		case 'r': // rotation
		{
			int index = std::stoi(transform.substr(1));
			auto& r = rotations[index];
			glm::vec3 axis = glm::vec3(r.axis_x, r.axis_y, r.axis_z);
			res = glm::rotate(identity, glm::radians(r.angle), axis) * res;
			break;
		}
		default:
			break;
		}
	}
	return res;
}