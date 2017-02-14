#ifndef CVISUAL_UTILS_H
#define CVISUAL_UTILS_H
#include"../DataColle/mesh_object.h"
#include<iostream>
#include<vector>
#include<map>

class CVisualUtils
{
protected:

	static std::map<int,OpenMesh::Vec3d>rand_color_map_;
public:
	static OpenMesh::Vec3d GetRandColor();
	static OpenMesh::Vec3d GetRandColorByTag(int id);
};
#endif