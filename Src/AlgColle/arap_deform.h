#ifndef CARAP_DEFORM_H
#define CARAP_DEFORM_H
#include"prereq.h"
#include <igl/arap.h>
#include"../DataColle/mesh_object.h"
#include"../DataColle/cgal_igl_converter.h"
#include<map>
//class ALGCOLLE_CLASS CARAPDeformation 
//{
//public:
//	//CARAPDeformation() : CMeshDeformation() {}
//	CARAPDeformation(COpenMeshT& mesh);
//	virtual ~CARAPDeformation();
//	//bool SetRoiVertices();//Set all of vertices on polyhedron as ROI
//	//bool SetRoiVertices(const std::vector<COpenMeshT::VertexHandle>& vertices_handle_vec);//Set specified vertices as ROI
//	bool SetDeformHandle(std::vector<COpenMeshT::VertexHandle>& deform_handles);
//	bool Deform(std::map<COpenMeshT::VertexHandle,OpenMesh::Vec3d>&handle_map);
//private:
//	COpenMeshT& mesh_;
//	Eigen::MatrixXd v_;
//	Eigen::MatrixXi f_;
//	Eigen::VectorXi handle_id_;
//	Eigen::MatrixXd handle_targt_;
//	igl::ARAPData arap_data_;
//	std::vector<COpenMeshT::VertexHandle> deform_handles_;
//	
//};
class CgalArapDeform;

class ALGCOLLE_CLASS CARAPDeformation
{
protected:
	CgalArapDeform* cgal_arap_deform_=NULL;
	Polyhedron poly_;
	COpenMeshT &mesh_;
	std::vector<COpenMeshT::VertexHandle>cgal_id_to_om_vh_map_;
	std::vector<Polyhedron::Vertex_handle>om_id_to_cgal_vh_map_;
	//std::map<COpenMeshT::VertexHandle, OpenMesh::Vec3d>deform_map_;
public:
	CARAPDeformation(COpenMeshT& mesh);
	~CARAPDeformation();
	bool SetDeformMap(std::map<COpenMeshT::VertexHandle, OpenMesh::Vec3d>& deform_handle_map);


	bool Deform();
};


#endif