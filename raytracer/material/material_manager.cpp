#include "material_manager.h"

MaterialManager::MaterialManager(const std::vector<Material_>& materials) 
{
	material_list.resize(materials.size());
	for (const auto& mat : materials) 
	{
		material_list[mat.id] = Material(mat);
	}
}

const Material& MaterialManager::getMaterialById(int id) const 
{
	return material_list[id];
}
