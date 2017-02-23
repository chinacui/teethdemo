#ifndef CCURVE_BASE_ALG_H
#define CCURVE_BASE_ALG_H
#include"../DataColle/curve_object.h"
#include"../DataColle/geo_primit.h"
#include"prereq.h"
class ALGCOLLE_CLASS CCurveBaseAlg
{
public:
	//proj_dir contain 2 vectors spanned a planner space
	static void ProjectCurve2Plannar(std::vector<OpenMesh::Vec3d>&curve, std::vector<OpenMesh::Vec3d> proj_dir,std::vector<OpenMesh::Vec2d>&curve_2d);//proj_dir : 2d frame
	static void ComputeMeanOfCurve(std::vector<OpenMesh::Vec3d>&curve, OpenMesh::Vec3d &res_mean);
	static void ComputeBoundingBoxOf2dCurve(std::vector<OpenMesh::Vec2d>&curve,OpenMesh::Vec2d &bbox_min,OpenMesh::Vec2d &bbox_max);
	static double ResampleCurve(std::vector<OpenMesh::Vec3d>&curve, int num, bool is_closed = false);
	//static bool CurveResampling(std::vector<OpenMesh::Vec3d> p_curveobject);
	static OpenMesh::Vec2d ComputeCurveCenter(std::vector<OpenMesh::Vec2d> &curve);
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
	static void ComputeNearestPoint(OpenMesh::Vec3d p, std::pair<OpenMesh::Vec3d,OpenMesh::Vec3d>&seg, OpenMesh::Vec3d &res_p, double &res_dis);
	static void ComputeClosestPoint(std::vector<OpenMesh::Vec3d>&curve,bool is_closed, OpenMesh::Vec3d p, OpenMesh::Vec3d &res_p, int &res_pid);
	//fit polynomial , for example y=ax^3+bx^2+cx+d
	static void PolynomialFitting(std::vector<OpenMesh::Vec2d>&curve, int degree,std::vector<double>&coeffs);
	static void SamplePointsOfPolynomial(std::vector<double>&coeff, double min_x, double max_x, double seg_len,std::vector<OpenMesh::Vec2d>&res_points);
	static OpenMesh::Vec2d ComputePointOfPolynomial(std::vector<double>&coedff, double x);
	static OpenMesh::Vec2d ComputeNormalOfPolynomial(std::vector<double>&coeffs,double x);
	static OpenMesh::Vec2d ComputeNearestPoint2Polynomial(std::vector<double>&coeff, OpenMesh::Vec2d p);
};
#endif
