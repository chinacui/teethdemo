#ifndef MORPH_SKEL_DENTAL_MESH_SEG_H
#define MORPH_SKEL_DENTAL_MESH_SEG_H
#include"prereq.h"
#include"../DataColle/mesh_object.h"
#include<vector>
#include"../DataColle/geo_primit.h"
#define DEBUG_MORPH_SKEL_DENTAL_MESH_SEG

class ALGCOLLE_CLASS CMorphSkelDentalMeshSeg
{
	enum CTag{Gingiva,Teeth,Feature,Non};
protected://variables
	CMeshObject &mesh_obj_;
	bool verbose_ ;
	//Eigen::VectorXd vertex_penalty_weight_;
protected://params
	double curvature_threshold_;
	double small_region_threshold_;
	Eigen::MatrixXd mean_curvature_vec_;
	Eigen::VectorXd mean_curvature_values_;
	CPlane cutting_plane_;
	std::vector<OpenMesh::Vec3d> edge_points_;
	std::vector<COpenMeshT::VertexHandle> edge_points_id_;
	std::vector<bool>is_edge_point_;
	std::vector<double> vertex_penalty_weight_;
	std::vector<CTag>tags_;
protected://functions
	void ComputeCuttingPlane();
	void ComputeMorphSkeleton();
	void ComputeEdgePointsFromMeanCurvatureThreshold(double thre);
	void FindOptimizePlane(CPlane ini_plane,CPlane &res_plane);
	void ComputeVertexPenaltyWeight();
	double ComputePenaltyValue(CPlane plane);
	void RemoveSmallFeatureRegions();
	void TagGingiva();
	void RemoveInnerTeethRegionInTeeth();
public:
	CMorphSkelDentalMeshSeg(CMeshObject &mesh_obj);
	~CMorphSkelDentalMeshSeg();

	void ComputeSegmentation(bool verbus=false);
	////////////////////////////////
	void TestErode();
	void TestDilate();
	void TestSkeletonize();
	void TestCurvature();
	void TestRemoveSmallFeatureRegions();
	void TestRender();
	void TestTagGingiva();
	double& TestGetCurvatureThreshold() { return curvature_threshold_; }
	double &TestGetSmallRegionThreshold() { return small_region_threshold_; }
};

#endif
