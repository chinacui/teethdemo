#ifndef CGEO_ALG_H
#define CGEO_ALG_H
#include"prereq.h"
#include<Eigen/Dense>
#include<iostream>
#include<vector>
#include"../DataColle/geo_primit.h"
#include"../DataColle/mesh_object.h"
class ALGCOLLE_CLASS CGeoAlg
{

protected:
	
public:
	static void FillHoles(Eigen::MatrixXd& vertexs, Eigen::MatrixXi& faces);
	static void FillHoles(COpenMeshT &mesh, bool remain_largest=false);
	static void PointSetPCA3D(std::vector<OpenMesh::Vec3d> &pts, OpenMesh::Vec3d&res_mean, std::vector<OpenMesh::Vec3d> &res_eigen_vects, std::vector<double>& res_eigen_values);//row order
	static void CutByPlane(COpenMeshT &mesh, CPlane plane, COpenMeshT &res_mesh,std::map<COpenMeshT::VertexHandle,COpenMeshT::VertexHandle>&vid_orig);//approximate return positive side of plane
	static void SeparateDisconnectedParts(COpenMeshT &mesh, std::vector<COpenMeshT*>&res_meshes, std::vector<std::map<COpenMeshT::VertexHandle, COpenMeshT::VertexHandle>>&vid_orig);
	static void GetFeatureGroupsByConnectivity(COpenMeshT&mesh, std::vector<bool>&is_feature,std::vector<std::vector<COpenMeshT::VertexHandle>>&feature_groups);
	static void SeparateMeshByVertexTag(COpenMeshT &mesh_obj, std::vector<int>&v_tag,std::vector<COpenMeshT*>&res_meshes, std::vector<std::map<COpenMeshT::VertexHandle, COpenMeshT::VertexHandle>>&vid_orig);//v_tag should be continues, start from 0
};
#endif