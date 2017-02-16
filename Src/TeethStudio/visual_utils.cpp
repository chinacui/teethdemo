#include"visual_utils.h"
std::map<int, OpenMesh::Vec3d>CVisualUtils::rand_color_map_;
OpenMesh::Vec3d CVisualUtils::GetRandColor()
{
	static bool flag = true;
	if (flag)
	{
		srand((unsigned)time(NULL));
	}
	flag = false;
	double r = (rand() % 10000) / 11000.0;
	double g = (rand() % 10000) / 11000.0;
	double b = (rand() % 10000) / 11000.0;
	OpenMesh::Vec3d color(r, g, b);
	return color;
}

OpenMesh::Vec3d CVisualUtils::GetRandColorByTag(int id)
{
	/*if (id != -1)
		std::cerr << id << std::endl;*/
	if (rand_color_map_.find(id)==rand_color_map_.end())
	{
		
		OpenMesh::Vec3d color= GetRandColor();
		rand_color_map_[id] = color;
	}
	return rand_color_map_[id];
}