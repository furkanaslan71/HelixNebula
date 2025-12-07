#ifndef TRANSFORMATIONS_H
#define TRANSFORMATIONS_H
#include <vector>
#include <string>

#include "../definitions.h"
#include "../include/parser.hpp"

glm::mat4 calculateCompositeTransformationMatrix(
	const std::vector<std::string> transformations,
	const std::vector<Translation_>& translations,
	const std::vector<Scaling_>& scalings,
	const std::vector<Rotation_>& rotations);

#endif //!TRANSFORMATIONS_H