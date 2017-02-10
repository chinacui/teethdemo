#ifndef CGEO_BASE_ALG_H
#define CGEO_BASE_ALG_H
#include"prereq.h"
#include<Eigen/Dense>
#include"../DataColle/geo_primit.h"
class ALGCOLLE_CLASS CGeoBaseAlg
{

protected:

public:
	static void ComputeMeanCurvatureValue(Eigen::MatrixXd& vertexs, Eigen::MatrixXi& faces,Eigen::VectorXd &res_curvature,bool is_signed=false);
	static void ComputeMeanCurvatureValue(COpenMeshT &mesh, Eigen::VectorXd &res_curvature, bool is_signed = false);
	static void ComputeMeanCurvatureVector(COpenMeshT&mesh, Eigen::MatrixXd &res_curvature);
	static void ComputeMeanCurvatureVector(Eigen::MatrixXd& vertexs, Eigen::MatrixXi& faces, Eigen::MatrixXd &res_curvature);
	static void ComputePerVertexNormal(COpenMeshT&mesh, Eigen::MatrixXd &res_normals);
	static double ComputeAverageEdgeLength(COpenMeshT&mesh);
	static void ComputePerVertexNormal(COpenMeshT&mesh, std::map<COpenMeshT::VertexHandle, OpenMesh::Vec3d>& vh_normals);
	static void ComputePerVertexNormal(Eigen::MatrixXd& vertexs, Eigen::MatrixXi& faces, Eigen::MatrixXd &res_normals);
	static double ComputeSquaredDistance2Plane(OpenMesh::Vec3d point, double a, double b, double c, double d);
	static double ComputeSquaredDistance2Plane(OpenMesh::Vec3d point, CPlane plane);
	static bool ConvertFromOpenMeshROIToOpenMesh(COpenMeshT &mesh, std::vector<COpenMeshT::FaceHandle>&roifaces, COpenMeshT&res_mesh, std::map<COpenMeshT::FaceHandle, COpenMeshT::FaceHandle>*new_mesh_2_old_vmap = NULL);
	static void GetNeighborFaces(COpenMeshT  &mesh,std::vector<COpenMeshT::VertexHandle>&vhs,int nei_num, std::vector<COpenMeshT::FaceHandle>&res_fhs);
	static void GetNeighborVhs(COpenMeshT &mesh, std::vector<COpenMeshT::VertexHandle>&vhs, int nei_num,std::vector<COpenMeshT::VertexHandle>&res_vhs);
	static OpenMesh::Vec3d GetFacePointFromBaryCoord(COpenMeshT&mesh, COpenMeshT::FaceHandle fh, OpenMesh::Vec3d barycoord);
	static bool IsOnPositiveSide(OpenMesh::Vec3d point, CPlane plane);
	static void RemoveNonManifold(COpenMeshT &mesh);
	static double ComputeDis(OpenMesh::Vec3d a, OpenMesh::Vec3d b);
	static OpenMesh::Vec3d Transform(Eigen::Matrix4d mat, OpenMesh::Vec3d p);
	static double Min_Dist_Point2_to_Line2(OpenMesh::Vec2d p, OpenMesh::Vec2d a, OpenMesh::Vec2d b);
	static double ComputeVertexArea(COpenMeshT &mesh, COpenMeshT::VertexHandle vh);
	static double ComputeFaceArea(COpenMeshT&mesh, COpenMeshT::FaceHandle fh);
	static void NormalizeMeshSize(COpenMeshT &mesh);
	static bool IsShareCommonVertex(COpenMeshT &mesh, COpenMeshT::FaceHandle fha, COpenMeshT::FaceHandle fhb,std::vector<COpenMeshT::VertexHandle>&res_comm_vhs= std::vector<COpenMeshT::VertexHandle>());
	static bool IsShareCommonEdge(COpenMeshT &mesh, COpenMeshT::FaceHandle fha, COpenMeshT::FaceHandle fhb,COpenMeshT::EdgeHandle &res_common_edge= COpenMeshT::EdgeHandle());
	static OpenMesh::Vec3d ComputeMeshCenter(COpenMeshT &mesh);
	static bool GetCommonEdge(COpenMeshT &mesh,COpenMeshT::VertexHandle vha, COpenMeshT::VertexHandle vhb, COpenMeshT::EdgeHandle &res_h);
	static bool IsConcavity(COpenMeshT &mesh, COpenMeshT::VertexHandle vh);
	static OpenMesh::Vec3d ComputeVPosFromBaryCoord(COpenMeshT&mesh, COpenMeshT::FaceHandle fh, OpenMesh::Vec3d&bary_coord);
	static OpenMesh::Vec3d ComputeVPosFromBaryCoord(OpenMesh::Vec3d a, OpenMesh::Vec3d b, OpenMesh::Vec3d c, OpenMesh::Vec3d&bary_coord);
	static OpenMesh::Vec3d ComputeBaryCoordInTriFace(COpenMeshT &mesh,COpenMeshT::FaceHandle fh, OpenMesh::Vec3d v);
	static OpenMesh::Vec3d ComputePointFromBaryCoord(COpenMeshT &mesh, COpenMeshT::FaceHandle fh, OpenMesh::Vec3d barycoord);
	static void GetOrderedBoundary(COpenMeshT &mesh, std::vector<std::vector<COpenMeshT::VertexHandle>>&bounds);
	static void GetLargestOrderedBoundary(COpenMeshT &mesh, std::vector<COpenMeshT::VertexHandle>&res_bounds);
	static bool GetHalfedgeHandle(COpenMeshT &mesh, COpenMeshT::VertexHandle vh0, COpenMeshT::VertexHandle vh1,COpenMeshT::HalfedgeHandle &res_hh);
	static void GetEdgeVertexs(COpenMeshT&mesh, std::vector<bool>&scalars, std::vector<OpenMesh::VertexHandle>&res_vhs);
	static void GetFaceVhs(COpenMeshT&mesh, COpenMeshT::FaceHandle fh, std::vector<COpenMeshT::VertexHandle>&resvhs);
	static OpenMesh::VertexHandle GetClosestVhFromFacePoint(COpenMeshT&mesh, COpenMeshT::FaceHandle fh, OpenMesh::Vec3d barycoord);
};
#endif