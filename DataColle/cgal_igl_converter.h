#ifndef FACEMODE_DATACOLLE_CGAL_IGL_CONVERTER_H
#define FACEMODE_DATACOLLE_CGAL_IGL_CONVERTER_H
#include "prereq.h"
#include "Polyhedron_type.h"
#include <Eigen/Core>
#include"custom_openmesh_type.h"

class DATACOLLE_CLASS CConverter
{
public:
	static  bool  ConvertFromCGALToIGL(Polyhedron& mesh, Eigen::MatrixXd& v, Eigen::MatrixXi& f);
	static  bool ConvertFromCGALRoiToIGL(Polyhedron& poly, std::vector<Face_handle>&fhs, Eigen::MatrixXd& v, Eigen::MatrixXi&f);
	static  bool ConvertFromIGLToCGAL(Eigen::MatrixXd& v, Eigen::MatrixXi& f, Polyhedron &mesh, std::vector<Face_handle>&fid_fh_map, std::vector<Vertex_handle>&vid_vhs_map);
	static  bool ConvertFromOpenMeshToIGL(COpenMeshT &mesh, Eigen::MatrixXd& v, Eigen::MatrixXi& f);

	static bool ConvertFromCGALToOpenMesh(Polyhedron& mesh, COpenMeshT& openmesh, std::map<Polyhedron::Vertex_handle, COpenMeshT::VertexHandle>& vh_map, std::map<Polyhedron::Facet_handle, COpenMeshT::FaceHandle>& fh_map);
	static bool ConvertFromOpenMeshToCGAL(COpenMeshT& openmesh, Polyhedron &mesh, std::map<COpenMeshT::VertexHandle, Polyhedron::Vertex_handle>& vh_map, std::map<COpenMeshT::FaceHandle, Polyhedron::Facet_handle>& fh_map);
	
	
};
#endif