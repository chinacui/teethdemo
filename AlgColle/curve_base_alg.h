#ifndef CCURVE_BASE_ALG_H
#define CCURVE_BASE_ALG_H
#include"../DataColle/curve_object.h"
#include"../DataColle/geo_primit.h"
#include"prereq.h"
class ALGCOLLE_CLASS CCurveBaseAlg
{
public:
	static void ProjectCurve2Plannar(std::vector<OpenMesh::Vec3d>&curve, std::vector<OpenMesh::Vec3d> proj_dir,std::vector<OpenMesh::Vec2d>&curve_2d);//proj_dir : 2d frame
	static void ComputeMeanOfCurve(std::vector<OpenMesh::Vec3d>&curve, OpenMesh::Vec3d &res_mean);
	static void ComputeBoundingBoxOf2dCurve(std::vector<OpenMesh::Vec2d>&curve,OpenMesh::Vec2d &bbox_min,OpenMesh::Vec2d &bbox_max);
	static void ComputeClosestPoint(std::vector<OpenMesh::Vec2d>&curve, OpenMesh::Vec2d p, int &res_pid);

	//static void QuadraticCurveFitting(std::vector<OpenMesh::Vec2d>&points, double&res_a, double &res_b, double &res_c, double &res_error);
};
#endif
