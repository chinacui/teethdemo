#ifndef CGEO_BASE_ALG_H
#define CGEO_BASE_ALG_H
#include"prereq.h"
#include<Eigen/Dense>
#include"../DataColle/geo_primit.h"
class ALGCOLLE_CLASS CGeoBaseAlg
{

protected:

public:
	static void ComputeMeanCurvatureValue(Eigen::MatrixXd& vertexs, Eigen::MatrixXi& faces,Eigen::VectorXd &res_curvature);
	static void ComputeMeanCurvatureVector(Eigen::MatrixXd& vertexs, Eigen::MatrixXi& faces, Eigen::MatrixXd &res_curvature);
	static void ComputePerVertexNormal(Eigen::MatrixXd& vertexs, Eigen::MatrixXi& faces, Eigen::MatrixXd &res_normals);
	static double ComputeSquaredDistance2Plane(OpenMesh::Vec3d point, double a, double b, double c, double d);
	static double ComputeSquaredDistance2Plane(OpenMesh::Vec3d point, CPlane plane);
	static bool IsOnPositiveSide(OpenMesh::Vec3d point, CPlane plane);
	static void RemoveNonManifold(COpenMeshT &mesh);
	static double ComputeVertexArea(COpenMeshT &mesh, COpenMeshT::VertexHandle vh);
	static double ComputeFaceArea(COpenMeshT&mesh, COpenMeshT::FaceHandle fh);
	static void NormalizeMeshSize(COpenMeshT &mesh);
	static OpenMesh::Vec3d ComputeBaryCoordInTriFace(COpenMeshT &mesh,COpenMeshT::FaceHandle fh, OpenMesh::Vec3d v);
	static OpenMesh::Vec3d ComputePointFromBaryCoord(COpenMeshT &mesh, COpenMeshT::FaceHandle fh, OpenMesh::Vec3d barycoord);
};
#endif