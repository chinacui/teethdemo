#ifndef CGEO_ALG_H
#define CGEO_ALG_H
#include"prereq.h"
#include<Eigen/Dense>
#include<iostream>
#include<vector>
#include<map>
#include"../DataColle/geo_primit.h"
#include"../DataColle/mesh_object.h"
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>


class ALGCOLLE_CLASS CGeoAlg
{

protected:
	
public:
	static bool SelfIntersectionRemoval(COpenMeshT& mesh);
	static void FillHoles(Eigen::MatrixXd& vertexs, Eigen::MatrixXi& faces);
	static void FillHoles(COpenMeshT &mesh, bool remain_largest=false);
	static void PointSetPCA3D(std::vector<OpenMesh::Vec3d> &pts, OpenMesh::Vec3d&res_mean, std::vector<OpenMesh::Vec3d> &res_eigen_vects, std::vector<double>& res_eigen_values);//row order
	static void CutByPlane(COpenMeshT &mesh, CPlane plane, COpenMeshT &res_mesh,std::map<COpenMeshT::VertexHandle,COpenMeshT::VertexHandle>&vid_orig);//approximate return positive side of plane
	static void SeparateDisconnectedParts(COpenMeshT &mesh, std::vector<COpenMeshT*>&res_meshes, std::vector<std::map<COpenMeshT::VertexHandle, COpenMeshT::VertexHandle>>&vid_orig);
	static void GetFeatureGroupsByConnectivity(COpenMeshT&mesh, std::vector<bool>&is_feature,std::vector<std::vector<COpenMeshT::VertexHandle>>&feature_groups);
	static void SeparateMeshByVertexTag(COpenMeshT &mesh_obj, std::vector<int>&v_tag,std::map<int,COpenMeshT*>&res_meshes, std::map<int,std::map<COpenMeshT::VertexHandle, COpenMeshT::VertexHandle>>&vid_orig);//v_tag should be continues, start from 0
	static bool RayMeshIntersection(OpenMesh::Vec3d  source, OpenMesh::Vec3d dir, CMeshObject &mesh_obj, COpenMeshT::FaceHandle & res_fh, OpenMesh::Vec3d &res_bary_coord);
	static bool RayMeshIntersection(OpenMesh::Vec3d  source, OpenMesh::Vec3d dir, CMeshObject &mesh_obj, COpenMeshT::VertexHandle & res_vh);
	static void ComputeClosestVertex(OpenMesh::Vec3d source, CMeshObject &mesh_obj, OpenMesh::VertexHandle &res_vh);
	static bool ComputeGeodesicPath(CMeshObject &mesh_obj,int svid,int tvid, std::vector<COpenMeshT::FaceHandle>&res_fhs, std::vector<OpenMesh::Vec3d>&res_bary_coords);
	static bool ComputeGeodesicPath(CMeshObject &mesh_obj, int svid, int tvid, std::vector<OpenMesh::Vec3d>&path);
	static bool ComputeGeodesicPath(CMeshObject &mesh_obj, COpenMeshT::VertexHandle svh, std::vector<COpenMeshT::VertexHandle>& tvhs, std::vector<OpenMesh::Vec3d>&path);
	static bool ComputeGeodesicPath(CMeshObject &mesh_obj, COpenMeshT::VertexHandle svh, std::vector<COpenMeshT::VertexHandle>& tvhs, std::vector<COpenMeshT::FaceHandle>&res_fhs, std::vector<OpenMesh::Vec3d>&res_bary_coords);
	static bool ComputeShortestPathAlongEdge(COpenMeshT &mesh, COpenMeshT::VertexHandle svh, COpenMeshT::VertexHandle tvh, std::vector<COpenMeshT::VertexHandle>&path);
	static void ComputeGeodesicDis(CMeshObject &mesh_obj, COpenMeshT::VertexHandle svh, std::vector<COpenMeshT::VertexHandle>&dst_vhs, std::vector<double>&res_dis);
	static void InitGeodesicModel(CMeshObject &mesh_obj);
	static void ComputeGradientOfScalarField(COpenMeshT &mesh, Eigen::VectorXd &scalars, Eigen::MatrixXd & res_grad);
	static void SimplifyMesh(OpenMesh::TriMesh_ArrayKernelT<COMTraits> &mesh,int edgenum);
	static void ExtractNRing(COpenMeshT &mesh, COpenMeshT::VertexHandle vh,int ringnum, std::vector<COpenMeshT::VertexHandle>&res_vhs);//except vh
	static void ComputeLocalExtremum(COpenMeshT &mesh, Eigen::VectorXd &scalars, int neighbor_num,std::vector<COpenMeshT::VertexHandle>&res_vhs);
	static void ComputeLocalExtremum(COpenMeshT &mesh, std::vector<double> &scalars, int neighbor_num, std::vector<COpenMeshT::VertexHandle>&res_vhs);
	static void LaplacianSmooth(COpenMeshT &mesh, int epoch, double step);
};
#endif