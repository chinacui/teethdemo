#include"curve_base_alg.h"
#include<opencv2\opencv.hpp>
#include"image_base_alg.h"
void CCurveBaseAlg::ProjectCurve2Plannar(std::vector<OpenMesh::Vec3d>&curve,  std::vector<OpenMesh::Vec3d> proj_dir, std::vector<OpenMesh::Vec2d>&curve_2d)
{
	std::cerr << "project" << std::endl;
	curve_2d.resize(curve.size());
	OpenMesh::Vec3d mean;
	ComputeMeanOfCurve(curve, mean);
	proj_dir[0]=proj_dir[0].normalize();
	proj_dir[1] = proj_dir[1].normalize();
	for (int i = 0; i < curve.size(); i++)
	{
		auto vdir=curve[i] - mean;
		OpenMesh::Vec2d v2d;
	
		v2d[0]=OpenMesh::dot<double,3>(vdir,proj_dir[0]);
		v2d[1] = OpenMesh::dot<double, 3>(vdir,proj_dir[1]);
		curve_2d[i] = v2d;
	}
}
void CCurveBaseAlg::ComputeClosestPoint(std::vector<OpenMesh::Vec2d>&curve, OpenMesh::Vec2d p, int &res_pid)
{
	double min_dis = std::numeric_limits<double>::max();
	int mi;
	for (int i = 0; i < curve.size(); i++)
	{
		double dis = (curve[i] - p).length();
		if (dis < min_dis)
		{
			min_dis = dis;
			res_pid = i;
		}
	}
}

void CCurveBaseAlg::ComputeBoundingBoxOf2dCurve(std::vector<OpenMesh::Vec2d>&curve, OpenMesh::Vec2d &bbox_min, OpenMesh::Vec2d &bbox_max)
{
	bbox_max = OpenMesh::Vec2d(std::numeric_limits<double>::min(), std::numeric_limits<double>::min());
	bbox_min = OpenMesh::Vec2d(std::numeric_limits<double>::max(), std::numeric_limits<double>::max());
	for (int i = 0; i < curve.size(); i++)
	{
		if (bbox_min[0] > curve[i][0])
			bbox_min[0] = curve[i][0];
		if (bbox_min[1] > curve[i][1])
			bbox_min[1] = curve[i][1];
		if (bbox_max[0] < curve[i][0])
			bbox_max[0] = curve[i][0];
		if (bbox_max[1] < curve[i][1])
			bbox_max[1] = curve[i][1];
	}
}
void CCurveBaseAlg::ComputeMeanOfCurve(std::vector<OpenMesh::Vec3d>&curve, OpenMesh::Vec3d &res_mean)
{
	res_mean = OpenMesh::Vec3d(0, 0, 0);
	for (int i = 0; i < curve.size(); i++)
	{
		res_mean = res_mean + curve[i];
	}
	res_mean = res_mean / curve.size();
}