#ifndef FACEMOD_ALGCOLLE_CORRESPONDENCE_BUILDER_H
#define FACEMOD_ALGCOLLE_CORRESPONDENCE_BUILDER_H
#include "prereq.h"
#include"../DataColle/mesh_object.h"

class ALGCOLLE_CLASS CCorrespondenceBuilder
{
public:
	CCorrespondenceBuilder();

	//
	void SetParameters(double sigma_p=0.03, double sigma_n= 0.1, double sigma_c= 1)
	{
	
		sigma_p_ = sigma_p;
		sigma_n_ = sigma_n;
		sigma_c_ = sigma_c;
	};

	std::map<COpenMeshT::VertexHandle,OpenMesh::Vec3d> FindCorrespondenceMap(CMeshObject* p_polyobj, std::vector<OpenMesh::Vec3d> p_curveobj, std::vector<COpenMeshT::VertexHandle> roi_region = std::vector<COpenMeshT::VertexHandle>());
	std::map<COpenMeshT::VertexHandle, OpenMesh::Vec3d> FindCorrespondenceMapNormal(CMeshObject* p_polyobj, std::vector<OpenMesh::Vec3d>& p_curveobj, std::vector<COpenMeshT::VertexHandle> roi_region);

private:
	void ComputeEmissionProbabilities();
	void ComputeTransitionProbabilities();
	bool HMMMatching();
	CMeshObject* p_polyobj_ = NULL;
	std::vector<OpenMesh::Vec3d> p_curveobj_ ;
	std::map<COpenMeshT::VertexHandle, OpenMesh::Vec3d> roi_vertices_normal_map_;
	std::vector<std::vector<double> > ep_matrix_;
	std::vector<std::vector<double> > vdist_matrix_;
	std::vector<double> pdist_vec_;
	std::vector<OpenMesh::Vec3d> curve_normal_vec_;
	std::map<COpenMeshT::VertexHandle,OpenMesh::Vec3d> deform_handle_map_;
	double sigma_p_;
	double sigma_n_;
	double sigma_c_;
};
#endif