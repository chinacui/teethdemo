#ifndef CTEETH_TEMPLATE_OBJECT_H
#define CTEETH_TEMPLATE_OBJECT_H
#include"mesh_object.h"
class DATACOLLE_CLASS CTeethTemplateObject :public CMeshObject
{

public:
	CTeethTemplateObject();
	CTeethTemplateObject(CTeethTemplateObject & template_object);
	void ComputeCrownMesh();
	COpenMeshT& GetCrownMesh();
	std::map<COpenMeshT::VertexHandle, COpenMeshT::VertexHandle> &GetCrown2TempVhsMap();
	OpenMesh::Vec3d & UpDir() { return updir_; }
	OpenMesh::Vec3d & FrontDir() { return front_dir_; }
	void SetType(std::string type) { type_ = type; }
	std::string GetType() { return type_; }
	std::vector<COpenMeshT::FaceHandle>&CrownFhs() { return crown_fhs_; }
	
	void SetChanged(bool flag = true);
protected:
	std::vector<COpenMeshT::FaceHandle>crown_fhs_;
	OpenMesh::Vec3d updir_, front_dir_;
	std::string type_;

	COpenMeshT crown_;
	std::map<COpenMeshT::VertexHandle, COpenMeshT::VertexHandle> crown2tempvhs_map_;

};
#endif