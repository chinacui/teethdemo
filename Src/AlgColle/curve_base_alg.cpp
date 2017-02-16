#include"curve_base_alg.h"
#include<opencv2\opencv.hpp>
#include"image_base_alg.h"
#include <Eigen/Sparse>
#include"image_base_alg.h"
#include"../DataColle/Polyhedron_type.h"
#include<CGAL/squared_distance_2_1.h>
#include"geo_base_alg.h"
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
OpenMesh::Vec2d CCurveBaseAlg::ComputeNormalOfPolynomial(std::vector<double>&coeffs, double x)
{
	double px = 1;
	double py = 0;
	double dx = 1;
	for (int i = 1; i < coeffs.size(); i++)
	{
		py=i*coeffs[i] * dx;
		dx = dx*x;
	}
	OpenMesh::Vec2d dir(py, -px);
	dir.normalize();
	if (OpenMesh::dot(dir, OpenMesh::Vec2d(0, 1)) > 0)
	{
		return dir;
	}
	else
	{
		return -dir;
	}

}
void CCurveBaseAlg::PolynomialFitting(std::vector<OpenMesh::Vec2d>&points, int degree, std::vector<double>&coeffs)
{
	
	int coeff_num = degree + 1;
	Eigen::MatrixXd ox_matrix(points.size(), coeff_num);
	Eigen::VectorXd oy_matrix(points.size());
	for (int i = 0; i < points.size(); i++)
	{
		oy_matrix(i) = points[i][1];
	}
	for (size_t r = 0; r < points.size(); r++)
	{
		double nVal = 1.0f;
		for (int c = 0; c < coeff_num; c++)
		{
			ox_matrix(r, c) = nVal;
			nVal *= points[r][0];
		}
	}
	Eigen::MatrixXd ox_t_matrix = ox_matrix.transpose();
	Eigen::MatrixXd ox_t_x_matrix = ox_t_matrix*ox_matrix;
	Eigen::MatrixXd ox_t_y_matrix = ox_t_matrix*oy_matrix;

	Eigen::FullPivLU<Eigen::MatrixXd>solver;
	solver.compute(ox_t_x_matrix);
	Eigen::MatrixXd res_coeff_eigen = solver.solve(ox_t_y_matrix);
	coeffs.clear();
	std::cerr << "poly fitting coeff ";
	for (int i = 0; i < res_coeff_eigen.rows(); i++)
	{
		coeffs.push_back(res_coeff_eigen(i));
		std::cerr << res_coeff_eigen(i)<<" ";
	}
	std::cerr << std::endl;
}
void CCurveBaseAlg::ComputeClosestPoint(std::vector<OpenMesh::Vec3d>&curve, OpenMesh::Vec3d p, int &res_pid)
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
double CCurveBaseAlg::ComputeLenOfCurve(std::vector<OpenMesh::Vec2d>&curve)
{
	if (curve.size() <= 1)
		return 0;
	OpenMesh::Vec2d p = curve[0];

	double len = 0;
	for (int i = 1; i < curve.size(); i++)
	{
		len += (curve[i] - p).length();
		p = curve[i];
	}
	return len;
}
//void CCurveBaseAlg::ComputeClosedCurveNormal(std::vector<OpenMesh::Vec3d> &curve, OpenMesh::Vec3d &normal)
//{
//	double a_xy, a_zx, a_yz;
//	a_xy = a_zx = a_yz = 0;
//
//	for (int i = 0; i < curve.size(); i++)
//	{
//		int j = (i + 1) % curve.size();
//		OpenMesh::Vec3d v0, v1;
//		v0 = curve[i];
//		v1 = curve[j];
//		a_xy += v0[0] * v1[1] - v1[0] * v0[1];
//		a_zx += v0[2] * v1[0] - v1[2] * v0[0];
//		a_yz += v0[1] * v1[2] - v1[1] * v0[2];
//
//
//	}
//
//	a_xy /= 2;
//	a_zx /= 2;
//	a_yz /= 2;
//	normal = OpenMesh::Vec3d(a_yz, a_zx, a_xy);
//	normal = normal / sqrt(normal.sqrnorm());
//}
double CCurveBaseAlg::ComputeScalarOnSegment(OpenMesh::Vec3d seg0, OpenMesh::Vec3d seg1, OpenMesh::Vec3d p, double s0, double s1)
{
	OpenMesh::Vec3d diff = seg0 - seg1;
	double len = std::sqrt(diff.sqrnorm());
	OpenMesh::Vec3d diff2 = p - seg0;
	double a = std::sqrt(diff2.sqrnorm());
	a /= len;
	double b = 1 - a;
	double res = s0*b + s1 *a;
	return res;
}
void CCurveBaseAlg::ComputeArclenParamOfAnchorPoints(std::vector<OpenMesh::Vec3d>&curve, std::vector<std::pair<int, OpenMesh::Vec3d>>&anchors, std::vector<double>&res_params, bool is_closed)
{
	if (curve.size() <= 1)
		return;
	std::vector<double>lens;
	lens.push_back(0);
	for (int i = 1; i < curve.size(); i++)
	{
		double pre_len = lens.back();
		double dis=CGeoBaseAlg::ComputeDis(curve[i], curve[i - 1]);
		lens.push_back(pre_len + dis);
	}
	double tot_len = lens.back();
	if (is_closed)
	{
		tot_len += CGeoBaseAlg::ComputeDis(curve[0], curve.back());
	}
	res_params.resize(anchors.size(), 0);
	for (int i = 0; i < anchors.size(); i++)
	{
		res_params[i] = (lens[anchors[i].first] + CGeoBaseAlg::ComputeDis(anchors[i].second, curve[anchors[i].first])) / tot_len;
	}


}
bool CCurveBaseAlg::ComputeMatchingWithAnchorsFixed(std::vector<OpenMesh::Vec3d>&src_curve, std::vector<std::pair<int, OpenMesh::Vec3d>>&src_anchor, std::vector<OpenMesh::Vec3d>&target_curve, std::vector<std::pair<int, OpenMesh::Vec3d>>&target_anchor, std::vector<std::pair<int, double>>&correspond)
{
	if (src_anchor.size() != target_anchor.size())
	{
		std::cerr <<"anchors_vhs.size()!=anchor_sil.size()";
		return false;
	}

	
	
	std::vector<std::pair<int, OpenMesh::Vec3d>>tmp_target_anchor, tmp_src_anchor;

	std::vector<OpenMesh::Vec3d>tmp_target_curve, tmp_src_curve;
	//make sure src and target have same orientation
	bool is_flap = false;
	tmp_target_curve = target_curve;
	tmp_target_anchor = target_anchor;
	if (src_anchor.size() >= 3)//may need flap
	{
		std::vector<double>src_anchor_params, target_anchor_params;
		ComputeArclenParamOfAnchorPoints(src_curve, src_anchor, src_anchor_params);
		ComputeArclenParamOfAnchorPoints(target_curve, target_anchor, target_anchor_params);
		int src_anchor_orient = 0, target_anchor_orient = 0;
		if ((src_anchor_params[1] > src_anchor_params[0] && src_anchor_params[2] < src_anchor_params[1]) || (src_anchor_params[1]<src_anchor_params[0] && src_anchor_params[2]>src_anchor_params[1]) || (src_anchor_params[1] < src_anchor_params[0] && src_anchor_params[2] < src_anchor_params[1]))
		{
			src_anchor_orient = 1;
		}
		if ((target_anchor_params[1] > target_anchor_params[0] && target_anchor_params[2] < target_anchor_params[1]) || (target_anchor_params[1]<target_anchor_params[0] && target_anchor_params[2]>target_anchor_params[1]) || (target_anchor_params[1] < target_anchor_params[0] && target_anchor_params[2] < target_anchor_params[1]))
		{
			target_anchor_orient = 1;
		}
		if (src_anchor_orient != target_anchor_orient)
		{
			is_flap = true;
		}
		if (is_flap ==false)
		{
			for (int i = 0; i < src_curve.size(); i++)
			{
				tmp_src_curve.push_back(src_curve[i]);
			}
			for (int i = 0; i < target_anchor.size(); i++)
			{
				tmp_src_anchor.push_back(src_anchor[i]);
			}
		}
		else
		{
			for (int i = src_curve.size() - 1; i >= 0; i--)
			{
				tmp_src_curve.push_back(src_curve[i]);
			}
			for (int i = 0; i < src_anchor.size(); i++)
			{
				tmp_src_anchor.push_back(src_anchor[i]);
				if (src_anchor[i].second != src_curve[src_anchor[i].first])
				{
					tmp_src_anchor.back().first = tmp_src_curve.size() - 1 - (tmp_src_anchor[i].first + 1) % tmp_src_curve.size();
				}
				else
				{
					tmp_src_anchor.back().first = tmp_src_curve.size() - 1 - (tmp_src_anchor[i].first ) % tmp_src_curve.size();
				}
				
			}
		}

		if (target_anchor_orient == 1)
		{
			std::reverse(tmp_src_anchor.begin(), tmp_src_anchor.end());
			std::reverse(tmp_target_anchor.begin(), tmp_target_anchor.end());
		}

	}
	
	

	






	std::vector<std::pair<int, double>>tmp_correspond(tmp_src_curve.size());
	for (int i = 0; i < tmp_target_anchor.size(); i++)
	{

		auto t_seg0 = tmp_target_anchor[i];
		auto t_seg1 = tmp_target_anchor[(i + 1) % tmp_target_anchor.size()];
		auto s_seg0 = tmp_src_anchor[i];
		auto s_seg1 = tmp_src_anchor[(i + 1) % tmp_src_anchor.size()];
		if (s_seg0.first == s_seg1.first)
			continue;

		//compute src length
		int j = s_seg0.first;
		OpenMesh::Vec3d src_pre_p = s_seg0.second;
		double src_len = 0;
		do {
			j = (j + 1) % tmp_src_curve.size();
			OpenMesh::Vec3d diff = tmp_src_curve[j] - src_pre_p;
			src_len += std::sqrt(diff.sqrnorm());
			src_pre_p = tmp_src_curve[j];
		} while (j != s_seg1.first);

		OpenMesh::Vec3d diff = s_seg1.second - src_pre_p;
		src_len += std::sqrt(diff.sqrnorm());



		//compute target length
		int t = t_seg0.first;
		OpenMesh::Vec3d target_pre_p = t_seg0.second;
		double target_len = 0;
		if (t_seg0.first == t_seg1.first)
		{
			OpenMesh::Vec3d diff = t_seg1.second - t_seg0.second;
			target_len = std::sqrt(diff.sqrnorm());
		}
		else
		{
			do {

				t = (t + 1) % tmp_target_curve.size();
				OpenMesh::Vec3d diff = tmp_target_curve[t] - target_pre_p;
				target_len += std::sqrt(diff.sqrnorm());
				target_pre_p = tmp_target_curve[t];
			} while (t != t_seg1.first);

			OpenMesh::Vec3d diff = t_seg1.second - target_pre_p;
			target_len += std::sqrt(diff.sqrnorm());
		}
		//std::cerr << "src len: " << src_len << " target len " << target_len << std::endl;

		double src_current_length = 0, target_current_length = 0;
		double target_current_ratio = target_current_length / target_len;
		double target_pre_ratio = target_current_ratio;
		src_pre_p = s_seg0.second;
		j = s_seg0.first;
		target_pre_p = t_seg0.second;
		t = t_seg0.first;
		do {
			j = (j + 1) % tmp_src_curve.size();

			OpenMesh::Vec3d diff = tmp_src_curve[j] - src_pre_p;
			src_current_length += std::sqrt(diff.sqrnorm());
			double src_current_ratio = src_current_length / src_len;
			src_pre_p = tmp_src_curve[j];
			if (src_current_ratio == 0)
			{
				std::cerr << "error src_current_ratio == 0!!!!!!!" << std::endl;
				std::cerr << "src_pre_p " << src_pre_p << " tmp_src_curve[j] " << tmp_src_curve[j] << std::endl;
				std::cerr << "s_seg0.first " << s_seg0.first << " j " << j << std::endl;
				std::cerr << "s_seg0.second " << s_seg0.second << std::endl;
			}

			while (target_current_ratio < src_current_ratio)
			{
				target_pre_p = tmp_target_curve[t];
				target_pre_ratio = target_current_ratio;
				t = (t + 1) % tmp_target_curve.size();
				OpenMesh::Vec3d diff = tmp_target_curve[t] - target_pre_p;
				target_current_length += std::sqrt(diff.sqrnorm());
				target_current_ratio = target_current_length / target_len;

			}

			int pret = (t + tmp_target_curve.size() - 1) % tmp_target_curve.size();

			double ratio = (src_current_ratio - target_pre_ratio) / (target_current_ratio - target_pre_ratio);
			if (target_current_ratio - target_pre_ratio == 0)
			{ 
				std::cerr << "target_current_ratio - target_pre_ratio==0" << std::endl;
			}
			OpenMesh::Vec3d tmp_p = OpenMesh::Vec3d(0, 0, 0) + ((1 - ratio)*(target_pre_p - OpenMesh::Vec3d(0, 0, 0)) + ratio*(target_curve[t] - OpenMesh::Vec3d(0, 0, 0)));
			double s = ComputeScalarOnSegment(tmp_target_curve[pret], tmp_target_curve[t], tmp_p, 0, 1);
			/*if (CGeoBaseAlg::ComputeDis(tmp_target_curve[pret], tmp_target_curve[t]) == 0)
			{
				std::cerr << "error CGeoBaseAlg::ComputeDis(tmp_target_curve[pret], tmp_target_curve[t]) == 0" << std::endl;
				std::cerr << "pret " << pret << " t " << t << std::endl;
			}*/
			tmp_correspond[j] = std::make_pair(pret, s);
			if (s == 1)
			{
				tmp_correspond[j].first = (tmp_correspond[j].first+1 + tmp_target_curve.size()) % tmp_target_curve.size();
				tmp_correspond[j].second = 0;
			}
			//std::cerr << "tmp corres:" << j << " " << tmp_correspond[j].first << " " << ratio << " " << tmp_p << "t p " << target_pre_p << " t t" << target_curve[t] << " " << 1 - ratio << std::endl;

		} while (j != s_seg1.first);

	}

	correspond.clear();
	std::cerr << "tmp corres size: " << tmp_correspond.size() << std::endl;

	if (is_flap)
	{
		for (int i = 0; i < tmp_correspond.size(); i++)
		{
			correspond.push_back(tmp_correspond[tmp_correspond.size() - i - 1]);
			
		}
	}
	else
	{
		for (int i = 0; i < tmp_correspond.size(); i++)
		{
			correspond.push_back(tmp_correspond[i]);
		}
	}
	std::cerr << " corres size: " << correspond.size() << std::endl;
	std::cerr << "end" << std::endl;

}
int CCurveBaseAlg::ComputLocalMaximamConcavityPoints(std::vector<OpenMesh::Vec2d>&curve, int neighbor_num, int neighbor_num_for_convexity, std::vector<int>&res_points)
{
	res_points.clear();
	bool isccw = true;
	std::vector<OpenMesh::Vec3d> tmp_curve;
	for (int i = 0; i < curve.size(); i++)
	{
		tmp_curve.push_back(OpenMesh::Vec3d(curve[i][0], curve[i][1], 0));
	}
	OpenMesh::Vec3d normal;
	ComputeClosedCurveNormal(tmp_curve, normal);
	if (normal[2] < 0)
		isccw = false;
	std::vector < double>convexity;
	ComputeConvexity(curve, neighbor_num_for_convexity, convexity,true);
	if (!isccw)
	{
		for (int i = 0; i < convexity.size(); i++)
		{
			convexity[i] = -convexity[i];
		}
	}
	for (int i = neighbor_num; i < curve.size() - neighbor_num; i++)
	{
		if (convexity[i] <= 0)
			continue;
		double maxm = -std::numeric_limits<double>::max();
		int mi = -1;
		for (int j = i - neighbor_num; j < i + neighbor_num; j++)
		{
			if (maxm < convexity[j])
			{
				maxm = convexity[j];
				mi = j;
			}
		}
		if (mi == i)
			res_points.push_back(i);
	}
	return res_points.size();
}

int CCurveBaseAlg::ComputLocalMinimalConcavityPoints(std::vector<OpenMesh::Vec2d>&curve, int neighbor_num, int neighbor_num_for_convexity, std::vector<int>&res_points)
{
	res_points.clear();
	bool isccw = true;
	std::vector<OpenMesh::Vec3d> tmp_curve;
	for (int i = 0; i < curve.size(); i++)
	{
		tmp_curve.push_back(OpenMesh::Vec3d(curve[i][0], curve[i][1], 0));
	}
	OpenMesh::Vec3d normal;
	ComputeClosedCurveNormal(tmp_curve, normal);
	if (normal[2] < 0)
		isccw = false;
	std::vector < double>convexity;
	ComputeConvexity(curve, neighbor_num_for_convexity, convexity,true);
	if (!isccw)
	{
		for (int i = 0; i < convexity.size(); i++)
		{
			convexity[i] = -convexity[i];
		}
	}
	for (int i = neighbor_num; i < curve.size() - neighbor_num; i++)
	{
		if (convexity[i] >= 0)
			continue;
		double minimal = std::numeric_limits<double>::max();
		int mi = -1;
		for (int j = i - neighbor_num; j < i + neighbor_num; j++)
		{
			if (minimal > convexity[j])
			{
				minimal = convexity[j];
				mi = j;
			}
		}
		if (mi == i)
			res_points.push_back(i);
	}
	return res_points.size();
}
void CCurveBaseAlg::ComputeClosedCurveNormal(std::vector<OpenMesh::Vec3d> &curve, OpenMesh::Vec3d &normal)
{
	double a_xy, a_zx, a_yz;
	a_xy = a_zx = a_yz = 0;

	for (int i = 0; i < curve.size(); i++)
	{
		int j = (i + 1) % curve.size();
		OpenMesh::Vec3d v0, v1;
		v0 = curve[i];
		v1 = curve[j];
		a_xy += v0[0] * v1[1] - v1[0] * v0[1];
		a_zx += v0[2] * v1[0] - v1[2] * v0[0];
		a_yz += v0[1] * v1[2] - v1[1] * v0[2];


	}

	a_xy /= 2;
	a_zx /= 2;
	a_yz /= 2;
	normal = OpenMesh::Vec3d(a_yz, a_zx, a_xy);
	normal.normalize();

}
double CCurveBaseAlg::ComputeConvexity(std::vector<OpenMesh::Vec2d>&curve, int neighbor_num, std::vector<double>&convexity, bool is_closed)
{
	convexity.resize(curve.size(), 0);
	double sum = 0;
	if (is_closed == false)
	{
		for (int i = neighbor_num; i < curve.size() - neighbor_num; i++)
		{
			int a = i - neighbor_num;
			int b = i + neighbor_num;
			OpenMesh::Vec2d pa = curve[a];
			OpenMesh::Vec2d pb = curve[b];
			
			Line_2 line(Point_2(pa[0],pa[1]), Point_2(pb[0],pb[1]));
			
			convexity[i] = std::sqrt(CGAL::squared_distance(line, Point_2(curve[i][0],curve[i][1])));

			OpenMesh::Vec3d dir0 = OpenMesh::Vec3d(pb[0], pb[1], 0) - OpenMesh::Vec3d(pa[0], pa[1], 0);
			OpenMesh::Vec3d dir1 = OpenMesh::Vec3d(curve[i][0], curve[i][1], 0) - OpenMesh::Vec3d(pa[0], pa[1], 0);
			OpenMesh::Vec3d cross_res = OpenMesh::cross(dir1, dir0);
			if (cross_res[2] > 0)
			{
				convexity[i] = -convexity[i];
			}
			sum += std::abs(convexity[i]);
		}
	}
	else
	{
		for (int i = 0; i < curve.size(); i++)
		{
			int a = (i - neighbor_num + curve.size()) % curve.size();
			int b = (i + neighbor_num) % curve.size();
			Point_2 pa(curve[a][0],curve[a][1]);
			Point_2 pb(curve[b][0],curve[b][1]);
			Line_2 line(pa, pb);
			convexity[i] = std::sqrt(CGAL::squared_distance(line, Point_2(curve[i][0],curve[i][1])));

			Vector dir0 = Point(pb[0], pb[1], 0) - Point(pa[0], pa[1], 0);
			Vector dir1 = Point(curve[i][0], curve[i][1], 0) - Point(pa[0], pa[1], 0);
			Vector cross_res = CGAL::cross_product(dir1, dir0);
			if (cross_res[2] > 0)
			{
				convexity[i] = -convexity[i];
			}
			sum += std::abs(convexity[i]);
		}
	}

	for (int i = 0; i < curve.size(); i++)
	{
		convexity[i] /= sum;
	}
	sum = 0;
	for (int i = 0; i < curve.size(); i++)
	{
		sum += convexity[i];
	}

	return sum;
}
double CCurveBaseAlg::ComputeLenOfCurve(std::vector<OpenMesh::Vec3d>&curve, bool is_closed)
{
	if (curve.size() <= 1)
		return 0;
	OpenMesh::Vec3d p = curve[0];
	double len = 0;
	for (int i = 1; i < curve.size(); i++)
	{
		len += (curve[i] - p).length();
		p = curve[i];
	}
	if (is_closed)
	{
		len += (curve[0] - curve.back()).length();
	}
	return len;
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

double CCurveBaseAlg::ResampleCurve(std::vector<OpenMesh::Vec3d>&curve, int sample_num, bool is_closed)
{
	if (curve.size() == 0)
		return 0;
	if (is_closed)
		curve.push_back(curve[0]);
	float  sketchTotLen = 0;
	float sketchAveLen;
	int vn = curve.size();
	for (int i = 1; i < vn; i++)
	{
		sketchTotLen += std::sqrt((curve[i] - curve[i - 1]).sqrnorm());
	}
	sketchAveLen = sketchTotLen / sample_num;
	std::vector<OpenMesh::Vec3d> curve_tmp;
	curve_tmp.push_back(curve[0]);
	OpenMesh::Vec3d pre_v = curve[0];
	OpenMesh::Vec3d cur_v;
	float clen = 0;
	for (int i = 1; i < vn; i++)
	{

		clen += std::sqrt((curve[i] - curve[i - 1]).sqrnorm());
		if (clen < sketchAveLen)
		{
			if (i == vn - 1 && (!is_closed))
			{
				double len = std::sqrt((curve[i] - curve_tmp.back()).sqrnorm());
				if (len < sketchAveLen / 2)
				{
					curve_tmp[curve_tmp.size() - 1] = curve[i];
				}
				else
				{
					curve_tmp.push_back(curve[i]);
				}
			}
			continue;
		}
		else
		{
			while (clen >= sketchAveLen)
			{
				if (i == vn - 1 && clen / sketchAveLen<1.4)
				{
					if (!is_closed)
					{
						curve_tmp.push_back(curve[i]);
					}
					break;
				}

				clen -= sketchAveLen;
				cur_v = curve[i] + (curve[i - 1] - curve[i]).normalized()*clen;
				curve_tmp.push_back(cur_v);


			}
			if (i == vn - 1 && !is_closed)
			{
				curve_tmp.back() = curve.back();
			}

		}
	}


	curve.clear();
	for (int i = 0; i < curve_tmp.size(); i++)
	{
		curve.push_back(curve_tmp[i]);
	}
	if (is_closed)
	{
		if ((curve[0] - curve.back()).length() < 0.4*sketchAveLen)
		{
			curve.pop_back();
		}
	}
	vn = curve.size();
	sketchTotLen = 0;
	for (int i = 1; i < vn; i++)
	{
		sketchTotLen += (curve[i] - curve[i - 1]).length();
	}
	//curve->ComputeBoundingBox();

	if (is_closed)
	{

		sketchTotLen += (curve.front() - curve.back()).length();
		return sketchTotLen / vn;
	}
	else
		return sketchTotLen / (vn - 1);
}
//bool CCurveBaseAlg::CurveResampling(std::vector<OpenMesh::Vec3d> curve)
//{
//	//compute curve length
//	double curve_length = 0;
//	for (int j = 1; j < curve.size(); j++)
//	{
//		curve_length += sqrt((curve[j] - curve[j - 1]).squared_length());
//	}
//	double dis_threshold = curve_length / 500;
//	std::queue<std::pair<int, int> > refine_intervals;
//	std::set<int> resampled_curve;
//	refine_intervals.push(std::pair<int, int>(0, (int)curve.size() - 1));
//	resampled_curve.insert(0);
//	resampled_curve.insert((int)curve.size() - 1);
//	while (!refine_intervals.empty())
//	{
//		int interval_begin = refine_intervals.front().first;
//		int interval_end = refine_intervals.front().second;
//		refine_intervals.pop();
//		Kernel::Segment_3 seg(curve[interval_begin], curve[interval_end]);
//
//		double max_dis = 0;
//		int flag_id = -1;
//		for (int j = interval_begin + 1; j < interval_end; j++)
//		{
//			double temp_dis = sqrt(CGAL::squared_distance(curve[j], seg));
//			if (temp_dis > max_dis)
//			{
//				flag_id = j;
//				max_dis = temp_dis;
//			}
//		}
//		if (max_dis > dis_threshold)
//		{
//			refine_intervals.push(std::pair<int, int>(interval_begin, flag_id));
//			refine_intervals.push(std::pair<int, int>(flag_id, interval_end));
//			resampled_curve.insert(flag_id);
//		}
//	}
//	//curve.clear();
//	std::vector<OpenMesh::Vec3d> new_curve;
//	for (auto itr = resampled_curve.begin(); itr != resampled_curve.end(); itr++)
//	{
//		new_curve.push_back(curve[*itr]);
//	}
//	curve.clear();
//	curve = new_curve;
//	p_curveobject->SetChanged(true);
//	return true;
//}


bool CCurveBaseAlg::ComputeCurveNormals(std::vector<OpenMesh::Vec3d>&curve, std::vector<OpenMesh::Vec3d> &normal_list)
{
	OpenMesh::Vec3d plane_normal =OpenMesh::Vec3d(0,0,1);
	//for (int i = 0; i < p_curvebobj->curve_.size(); i++)
	//{
	OpenMesh::Vec3d normal_vec;
	if (curve.size() == 1)
	{
		return false;//something unpredicable will occured
	}
	for (int j = 0; j <curve.size(); j++)
	{
		OpenMesh::Vec3d temp_v;
		if (j == 0)
		{
			temp_v = curve[j + 1] - curve[j];
		}
		else if (j + 1 == curve.size())
		{
			temp_v = curve[j] - curve[j - 1];
		}
		else
		{
			temp_v = curve[j + 1] - curve[j - 1];
		}
		normal_vec=OpenMesh::cross(plane_normal, temp_v);
		normal_vec = normal_vec.normalized();
		normal_list.push_back(normal_vec);
		
	}
	//}
	return true;
}
bool CCurveBaseAlg::ComputeCurveNormals(std::vector<OpenMesh::Vec2d>& curve, std::vector<OpenMesh::Vec2d> &res_normal_list)
{
	res_normal_list.clear();
	OpenMesh::Vec3d plane_normal = OpenMesh::Vec3d(0, 0, 1);
	//for (int i = 0; i < p_curvebobj->curve_.size(); i++)
	//{
	OpenMesh::Vec2d normal_vec;
	if (curve.size() == 1)
	{
		return false;//something unpredicable will occured
	}
	for (int j = 0; j <curve.size(); j++)
	{
		OpenMesh::Vec2d temp_v;
		if (j == 0)
		{
			temp_v = curve[j + 1] - curve[j];
		}
		else if (j + 1 == curve.size())
		{
			temp_v = curve[j] - curve[j - 1];
		}
		else
		{
			temp_v = curve[j + 1] - curve[j - 1];
		}

		OpenMesh::Vec3d tmp_v3d = OpenMesh::Vec3d(temp_v[0], temp_v[1], 0);
		OpenMesh::cross(plane_normal, tmp_v3d);
		OpenMesh::Vec3d tmp_normal_vec = OpenMesh::cross(plane_normal, tmp_v3d);
		normal_vec = OpenMesh::Vec2d(tmp_normal_vec[0], tmp_normal_vec[1]);
		normal_vec = normal_vec.normalized();
		res_normal_list.push_back(normal_vec);

	}
	//}
	return true;
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
bool CCurveBaseAlg::IsInRegion(std::vector<OpenMesh::Vec2d>&region, OpenMesh::Vec2d p)
{
	int nvert = region.size();
	bool c = false;
	int i, j;
	for (i = 0, j = nvert - 1; i < nvert; j = i++) {
		if ((region[i][1]> p[1]) != (region[j][1] > p[1]) &&
			(p[0] < (region[j][0] - region[i][0]) * (p[1] - region[i][1]) / (region[j][1] - region[i][1]) + region[i][0]))
			c = !c;
	}
	return c;
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