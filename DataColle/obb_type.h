#include"mesh_object.h"
#ifndef COBB_TYPE_H
#define COBB_TYPE_H
class CPqpObb;
class CObbWrapper
{
public:
	CPqpObb*obb_tree_;
	std::vector<COpenMeshT::VertexHandle>id_vh_map_;
	std::vector<COpenMeshT::FaceHandle>id_fh_map_;
	CObbWrapper()
	{
		obb_tree_ = NULL;
		id_vh_map_.clear();
		id_fh_map_.clear();
	}
	~CObbWrapper()
	{
		if (obb_tree_ != NULL)
		{
			delete obb_tree_;
		}
	}
};
#endif