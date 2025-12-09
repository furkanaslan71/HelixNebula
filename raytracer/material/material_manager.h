#ifndef MATERIAL_MANAGER_H
#define MATERIAL_MANAGER_H
#include <vector>

#include "parser/parser.hpp"
#include "material.h"

class MaterialManager {
public:
	MaterialManager() = default;
	MaterialManager(const std::vector<Material_>& materials);
	~MaterialManager() = default;
	const Material& getMaterialById(int id) const;
private:
	std::vector<Material> material_list;
};
#endif // MATERIAL_MANAGER_H