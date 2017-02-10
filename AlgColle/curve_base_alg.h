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
	static double ResampleCurve(std::vector<OpenMesh::Vec3d>&curve, int num, bool is_closed = false);
	//static bool CurveResampling(std::vector<OpenMesh::Vec3d> p_curveobject);
	static bool ComputeCurveNormals(std::vector<OpenMesh::Vec3d>& p_curvebobj, std::vector<OpenMesh::Vec3d> &res_normal_list);
	static bool ComputeCurveNormals(std::vector<OpenMesh::Vec2d>& p_curvebobj, std::vector<OpenMesh::Vec2d> &res_normal_list);
	static double ComputeLenOfCurve(std::vector<OpenMesh::Vec3d>&curve,bool is_closed=false);
	static double ComputeScalarOnSegment(OpenMesh::Vec3d seg0, OpenMesh::Vec3d seg1, OpenMesh::Vec3d p, double s0, double s1);
	static void ComputeClosedCurveNormal(std::vector<OpenMesh::Vec3d> &curve, OpenMesh::Vec3d &normal);
	static void ComputeArclenParamOfAnchorPoints(std::vector<OpenMesh::Vec3d>&curve, std::vector<std::pair<int, OpenMesh::Vec3d>>&anchors, std::vector<double>&res_params,bool is_closed=false);
	static bool ComputeMatchingWithAnchorsFixed(std::vector<OpenMesh::Vec3d>&src_curve, std::vector<std::pair<int, OpenMesh::Vec3d>>&src_anchor, std::vector<OpenMesh::Vec3d>&target_curve, std::vector<std::pair<int, OpenMesh::Vec3d>>&target_anchor, std::vector<std::pair<int, double>>&correspond);//assume the anchors are ordered
	static double ComputeConvexity(std::vector<OpenMesh::Vec2d>&curve, int neighbor_num, std::vector<double>&convexity, bool is_closed);
	static int ComputLocalMaximamConcavityPoints(std::vector<OpenMesh::Vec2d>&curve, int neighbor_num, int neighbor_num_for_convexity, std::vector<int>&res_points);
	static int ComputLocalMinimalConcavityPoints(std::vector<OpenMesh::Vec2d>&curve, int neighbor_num, int neighbor_num_for_convexity, std::vector<int>&res_points);
	
	static bool IsInRegion(std::vector<OpenMesh::Vec2d>&region, OpenMesh::Vec2d p);
	static double ComputeLenOfCurve(std::vector<OpenMesh::Vec2d>&curve);
	static void ComputeClosestPoint(std::vector<OpenMesh::Vec2d>&curve, OpenMesh::Vec2d p, int &res_pid);
	static void ComputeClosestPoint(std::vector<OpenMesh::Vec3d>&curve, OpenMesh::Vec3d p, int &res_pid);
	static void PolynomialFitting(std::vector<OpenMesh::Vec2d>&curve, int degree,std::vector<double>&coeffs);
	static OpenMesh::Vec2d ComputeNormalOfPolynomial(std::vector<double>&coeffs,double x);
};
#endif
