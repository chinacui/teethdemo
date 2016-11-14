#include"dental_base_alg.h"
#include"geo_alg.h"
#include"geo_base_alg.h"
#include"numerical_base_alg.h"
#include<igl/writeDMAT.h>
#include <igl/gaussian_curvature.h>
#include"../DataColle/cgal_igl_converter.h"
#include<opencv2\opencv.hpp>
#include"image_base_alg.h"
#include"curve_base_alg.h"
#include"morphlogic_operation.h"
#include"../DataColle/data_pool.h"
void CDentalBaseAlg::ComputePCAFrameFromHighCurvaturePoints(COpenMeshT&mesh, double threshold, OpenMesh::Vec3d&mean, std::vector<OpenMesh::Vec3d>&res_frame)
{
	Eigen::MatrixXd curvatures_vec;
	CGeoBaseAlg::ComputeMeanCurvatureVector(mesh, curvatures_vec);
	Eigen::VectorXd mean_curvature_values;
	mean_curvature_values.resize(curvatures_vec.rows());
	Eigen::MatrixXd N;
	CGeoBaseAlg::ComputePerVertexNormal(mesh, N);

	for (int i = 0; i < curvatures_vec.rows(); i++)
	{
		mean_curvature_values(i) = curvatures_vec.row(i).norm();
	}
	std::vector<OpenMesh::Vec3d>feature_points;
	for (auto v_iter = mesh.vertices_begin(); v_iter != mesh.vertices_end(); v_iter++)
	{
		int vid = v_iter->idx();
		
		if (mean_curvature_values(vid) > threshold&& curvatures_vec.row(vid).dot(N.row(vid)) < 0.0 && (!mesh.is_boundary(v_iter)))
		{
			feature_points.push_back(mesh.point(v_iter));
		}
	}

	std::vector<double>eigen_values;

	CGeoAlg::PointSetPCA3D(feature_points, mean, res_frame, eigen_values);


}

double CDentalBaseAlg::ComputePenaltyValue(COpenMeshT&mesh, CPlane plane,std::vector<COpenMeshT::VertexHandle>&feature_points, Eigen::VectorXd &mean_curvature, Eigen::VectorXd& vertex_penalty_weight_)
{
	double res = 0;


	for (int i = 0; i < feature_points.size(); i++)
	{
		auto vh = feature_points[i];
		int vid = vh.idx();
		OpenMesh::Vec3d point = mesh.point(vh);
		double dis = std::sqrt(CGeoBaseAlg::ComputeSquaredDistance2Plane(point, plane.a(), plane.b(), plane.c(), plane.d()));
		res += vertex_penalty_weight_[vid] * std::pow(dis - 18, 2);
	}
	res /= feature_points.size();
	return res;
}
void CDentalBaseAlg::ComputeHoleVertMean(COpenMeshT&mesh, OpenMesh::Vec3d &res_mean)
{
	res_mean = OpenMesh::Vec3d(0, 0, 0);
	int count = 0;
	for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); viter++)
	{
		if (mesh.is_boundary(viter))
		{
			count++;
			res_mean = res_mean + mesh.point(viter);
		}
	}
	res_mean /= count;
	
}

void CDentalBaseAlg::ComputeExtremePointsOfClosedCurve(std::vector<OpenMesh::Vec2d>&curve_2d, std::vector<int>&res_extreme_pids)
{
	cv::Mat bound_cur_img;
	unsigned char back_color = 0, fore_color = 255;
	int width = 400, height = 400;
	CImageBaseAlg::Curve2dToGrayImage(curve_2d, width, height, fore_color, back_color, bound_cur_img);
		//cv::imshow("a", bound_cur_img);
	//cv::waitKey(0);
	for (int i = 0; i < bound_cur_img.rows; i++)
	{
		if (bound_cur_img.at<uchar>(i, 0) == fore_color)
			cv::floodFill(bound_cur_img, cv::Point(0, i), back_color);
		if (bound_cur_img.at<uchar>(i, width - 1) == fore_color)
			cv::floodFill(bound_cur_img, cv::Point(width - 1, i), back_color);
	}

	for (int i = 0; i < bound_cur_img.cols; i++)
	{
		if (bound_cur_img.at<uchar>(0, i) == fore_color)
			cv::floodFill(bound_cur_img, cv::Point(i, 0), back_color);
		if (bound_cur_img.at<uchar>(height - 1, i) == fore_color)
			cv::floodFill(bound_cur_img, cv::Point(i, height - 1), back_color);
	}
	CImageBaseAlg::MorphSkeleton(bound_cur_img);
	std::vector<cv::Point>extreme_points;
	CImageBaseAlg::GetExtremePoints(bound_cur_img, fore_color, extreme_points);
	if (extreme_points.size() < 2)
	{
		std::cerr << "no extreme points found" << std::endl;
		res_extreme_pids.clear();
		return;
	}
	double max_dis = -1;
	int max_disi, max_disj;
	for (int i = 0; i < extreme_points.size(); i++)
	{
		std::vector<double>dis;
		CImageBaseAlg::ShortestDis(bound_cur_img, extreme_points[i], extreme_points, dis);
		for (int j = 0; j < dis.size(); j++)
		{
			if (dis[j] > max_dis)
			{
				max_dis = dis[j];
				max_disi = i;
				max_disj = j;
			}
		}
	}
	//cv::imshow("b", bound_cur_img);
	//cv::waitKey(0);
	cv::circle(bound_cur_img, extreme_points[max_disi], 4, 200);
	cv::circle(bound_cur_img, extreme_points[max_disj], 4, 200);
	OpenMesh::Vec2d curve2d_bboxmin, curve2d_bboxmax;
	OpenMesh::Vec2d curve2d_len;
	CCurveBaseAlg::ComputeBoundingBoxOf2dCurve(curve_2d, curve2d_bboxmin, curve2d_bboxmax);
	curve2d_len = curve2d_bboxmax - curve2d_bboxmin;
	OpenMesh::Vec2d om_extreme_points[2];
	om_extreme_points[0] = OpenMesh::Vec2d(extreme_points[max_disi].x*1.0 / width*curve2d_len[0] + curve2d_bboxmin[0], extreme_points[max_disi].y*1.0 / height*curve2d_len[1] + curve2d_bboxmin[1]);
	om_extreme_points[1] = OpenMesh::Vec2d(extreme_points[max_disj].x*1.0 / width*curve2d_len[0] + curve2d_bboxmin[0], extreme_points[max_disj].y*1.0 / height*curve2d_len[1] + curve2d_bboxmin[1]);
	res_extreme_pids.resize(2);
	CCurveBaseAlg::ComputeClosestPoint(curve_2d, om_extreme_points[0], res_extreme_pids[0]);
	CCurveBaseAlg::ComputeClosestPoint(curve_2d, om_extreme_points[1], res_extreme_pids[1]);
	if (res_extreme_pids[0] > res_extreme_pids[1])
	{
		int tmp = res_extreme_pids[0];
		res_extreme_pids[0] = res_extreme_pids[1];
		res_extreme_pids[1] = tmp;
	}

}
void CDentalBaseAlg::ComputeBoundCuttingPointsOfToothMesh(CMeshObject&meshobj, std::vector<std::vector<COpenMeshT::VertexHandle>>&res_inside_vhs, std::vector<std::vector<COpenMeshT::VertexHandle>>&res_outside_vhs)
{
	COpenMeshT &mesh = meshobj.GetMesh();

	std::vector<std::vector<COpenMeshT::VertexHandle>>inside_bounds, outside_bounds;
	CDentalBaseAlg::ComputeTwoSideBoundsOfToothMesh(mesh, inside_bounds, outside_bounds);

	std::vector<std::vector<COpenMeshT::VertexHandle>>tmp_inside_vhs(inside_bounds.size()), tmp_outside_vhs(outside_bounds.size());
	std::map<COpenMeshT::VertexHandle, std::pair<int,int>> path_from;
	std::vector<bool>is_cuttingpoint(mesh.n_vertices(), false);
	int bound_num = inside_bounds.size();
	for (int t = 0; t < bound_num; t++)
	{
		tmp_inside_vhs[t].clear();
		tmp_outside_vhs[t].clear();
		std::vector<COpenMeshT::VertexHandle>& bound0 = inside_bounds[t];
		std::vector<COpenMeshT::VertexHandle>& bound1 = outside_bounds[t];

		std::vector<COpenMeshT::FaceHandle>res_fhs;

		


		std::set<COpenMeshT::VertexHandle>oppovhs_set;
		int step = 3;
		for (int i = 0; i < bound0.size(); i += step)
		{

			std::vector<OpenMesh::Vec3d> barycoords;
			CGeoAlg::ComputeGeodesicPath(meshobj, bound0[i], bound1, res_fhs, barycoords);
			/*CCurveObject *ccobj = new CCurveObject();
			std::vector<OpenMesh::Vec3d> &curve = ccobj->GetCurve();
			curve.clear();
			for (int i = 0; i < res_fhs.size(); i++)
			{
				curve.push_back(CGeoBaseAlg::ComputePointFromBaryCoord(mesh, res_fhs[i], barycoords[i]));
			}
			ccobj->SetColor(OpenMesh::Vec3d(0, 0, 1));
			ccobj->SetChanged();
			DataPool::AddCurveObject(ccobj);*/


			auto tfh = res_fhs.back();
			auto tbary = barycoords.back();
			double max_bary = -1;
			int mi;
			for (int j = 0; j < 3; j++)
			{
				if (max_bary < tbary[j])
				{
					max_bary = tbary[j];
					mi = j;
				}
			}
			std::vector<COpenMeshT::VertexHandle>fvhs;
			CGeoBaseAlg::GetFaceVhs(mesh, tfh, fvhs);
			COpenMeshT::VertexHandle tvh = fvhs[mi];
			if (tvh != bound1.back() && tvh != bound1.front() )
			{
				if (is_cuttingpoint[tvh.idx()] == false)
				{
					is_cuttingpoint[tvh.idx()] = true;
					tmp_outside_vhs[t].push_back(tvh);
				}
	
				if (path_from.find(tvh) == path_from.end())
				{
					path_from[tvh] = std::pair<int,int>();
					path_from[tvh].first = path_from[tvh].second = i;
				}
				if (path_from[tvh].first > i)
					path_from[tvh].first = i;
				if (path_from[tvh].second < i)
					path_from[tvh].second = i;
				
		
			}
			
		}
		
		for (int i = 0; i < bound1.size(); i += step)
		{

			std::vector<OpenMesh::Vec3d> barycoords;
			CGeoAlg::ComputeGeodesicPath(meshobj, bound1[i], bound0, res_fhs, barycoords);
			auto tfh = res_fhs.back();
			auto tbary = barycoords.back();
			double max_bary = -1;
			int mi;
			for (int j = 0; j < 3; j++)
			{
				if (max_bary < tbary[j])
				{
					max_bary = tbary[j];
					mi = j;
				}
			}
			std::vector<COpenMeshT::VertexHandle>fvhs;
			CGeoBaseAlg::GetFaceVhs(mesh, tfh, fvhs);
			COpenMeshT::VertexHandle tvh = fvhs[mi];
			if (tvh != bound0.back() && tvh != bound0.front())
			{
				if (is_cuttingpoint[tvh.idx()] == false)
				{
					is_cuttingpoint[tvh.idx()] = true;
					tmp_inside_vhs[t].push_back(tvh);
				}

				if (path_from.find(tvh) == path_from.end())
				{
					path_from[tvh] = std::pair<int, int>();
					path_from[tvh].first = path_from[tvh].second = i;
				}
				if (path_from[tvh].first > i)
					path_from[tvh].first = i;
				if (path_from[tvh].second < i)
					path_from[tvh].second = i;
			}
		}
	}
	for (int i = 0; i < bound_num; i++)
	{
		if (tmp_inside_vhs[i].size() != 0 || tmp_outside_vhs[i].size() != 0)
		{
			for (auto iter = tmp_inside_vhs[i].begin(); iter != tmp_inside_vhs[i].end(); iter++)
			{
				auto vh = *iter;
				auto vid_from = path_from[vh];
				std::vector<COpenMeshT::VertexHandle>vh_from;
				bool flag = true;
				for (int j = vid_from.first; j <= vid_from.second; j++)
				{
					vh_from.push_back(outside_bounds[i][j]);
					if (is_cuttingpoint[outside_bounds[i][j].idx()])
					{
						flag = false;
						break;
					}
				}
				if (flag)
				{
					std::vector<COpenMeshT::FaceHandle>res_fhs;
					std::vector<OpenMesh::Vec3d> barycoords;
					CGeoAlg::ComputeGeodesicPath(meshobj, vh, vh_from, res_fhs, barycoords);
					auto tfh = res_fhs.back();
					auto tbary = barycoords.back();
					double max_bary = -1;
					int mi;
					for (int j = 0; j < 3; j++)
					{
						if (max_bary < tbary[j])
						{
							max_bary = tbary[j];
							mi = j;
						}
					}
					std::vector<COpenMeshT::VertexHandle>fvhs;
					CGeoBaseAlg::GetFaceVhs(mesh, tfh, fvhs);
					COpenMeshT::VertexHandle tvh = fvhs[mi];
					if (tvh != outside_bounds[i].back() && tvh != outside_bounds[i].front() && is_cuttingpoint[tvh.idx()] == false)
					{
						is_cuttingpoint[tvh.idx()] = true;
						tmp_outside_vhs[i].push_back(tvh);
				
					}

				}

			}
			for (auto iter = tmp_outside_vhs[i].begin(); iter != tmp_outside_vhs[i].end(); iter++)
			{
				auto vh = *iter;
				auto vid_from = path_from[vh];
				std::vector<COpenMeshT::VertexHandle>vh_from;
				bool flag = true;
				for (int j = vid_from.first; j <= vid_from.second; j++)
				{
					vh_from.push_back(inside_bounds[i][j]);
					if (is_cuttingpoint[inside_bounds[i][j].idx()])
					{
						flag = false;
						break;
					}
				}
				if (flag)
				{
					std::vector<COpenMeshT::FaceHandle>res_fhs;
					std::vector<OpenMesh::Vec3d> barycoords;
					CGeoAlg::ComputeGeodesicPath(meshobj, vh, vh_from, res_fhs, barycoords);
					auto tfh = res_fhs.back();
					auto tbary = barycoords.back();
					double max_bary = -1;
					int mi;
					for (int j = 0; j < 3; j++)
					{
						if (max_bary < tbary[j])
						{
							max_bary = tbary[j];
							mi = j;
						}
					}
					std::vector<COpenMeshT::VertexHandle>fvhs;
					CGeoBaseAlg::GetFaceVhs(mesh, tfh, fvhs);
					COpenMeshT::VertexHandle tvh = fvhs[mi];
					if (tvh != inside_bounds[i].back() && tvh != inside_bounds[i].front() && is_cuttingpoint[tvh.idx()] == false)
					{
						is_cuttingpoint[tvh.idx()] = true;
						tmp_inside_vhs[i].push_back(tvh);
			
					}

				}
			}
		}
	}
	res_inside_vhs.clear();
	res_outside_vhs.clear();
	int ri = 0;
	for (int i = 0; i < tmp_inside_vhs.size(); i++)
	{
	
		if (tmp_inside_vhs[i].size() != 0 || tmp_outside_vhs[i].size() != 0)
		{
			res_inside_vhs.push_back(std::vector<COpenMeshT::VertexHandle>());
			res_outside_vhs.push_back(std::vector<COpenMeshT::VertexHandle>());
			for (int j = 0; j < inside_bounds[i].size(); j++)
			{
				auto vh= inside_bounds[i][j];
				if (is_cuttingpoint[vh.idx()])
				{
					res_inside_vhs[ri].push_back(vh);
				}
			}
			for (int j = 0; j < outside_bounds[i].size(); j++)
			{
				auto vh = outside_bounds[i][j];
				if (is_cuttingpoint[vh.idx()])
				{
					res_outside_vhs[ri].push_back(vh);
				}
			}
			
			ri++;
		}
	}


}
void CDentalBaseAlg::ComputeCuttingPath(CMeshObject&meshobj, std::vector<std::vector<COpenMeshT::VertexHandle>>&inside_cutting_vhs, std::vector<std::vector<COpenMeshT::VertexHandle>>&outside_cutting_vhs, std::vector<CCuttingPath>&res_cuttingpath)
{
	res_cuttingpath.clear();
	COpenMeshT &mesh = meshobj.GetMesh();
	for (int i = 0; i < inside_cutting_vhs.size(); i++)
	{
		for (int j = 0; j < outside_cutting_vhs[i].size(); j++)
		{
			mesh.set_color(outside_cutting_vhs[i][j], OpenMesh::Vec3d(1, 0, 1));
		}
		for (int j = 0; j < inside_cutting_vhs[i].size(); j++)
		{
			mesh.set_color(inside_cutting_vhs[i][j], OpenMesh::Vec3d(0, 0, 1));
		}
		std::map<COpenMeshT::VertexHandle,CCuttingPath>cuttingpath_oftvhs_map;
		for (int j = 0; j < inside_cutting_vhs[i].size(); j++)
		{
			std::vector<COpenMeshT::FaceHandle>res_fhs;
			std::vector<OpenMesh::Vec3d> res_bary_coords;
			CGeoAlg::ComputeGeodesicPath(meshobj, inside_cutting_vhs[i][j], outside_cutting_vhs[i], res_fhs, res_bary_coords);
			
			COpenMeshT::VertexHandle tvh = CGeoBaseAlg::GetClosestVhFromFacePoint(mesh,res_fhs.back(),res_bary_coords.back());

			CCuttingPath cutting_path(mesh, inside_cutting_vhs[i][j], tvh, res_fhs, res_bary_coords);

		

			//res_cuttingpath.push_back(cutting_path);
			if (cuttingpath_oftvhs_map.find(tvh) == cuttingpath_oftvhs_map.end())
			{
			/*	CCurveObject *ccobj = new CCurveObject();
				cutting_path.GetPathPoints(ccobj->GetCurve());
				ccobj->SetChanged();
				ccobj->SetColor(OpenMesh::Vec3d(1, 0, 1));
				DataPool::AddCurveObject(ccobj);*/
				cuttingpath_oftvhs_map[tvh] = cutting_path;
			}
			else
			{
		/*		CCurveObject *ccobj = new CCurveObject();
				cutting_path.GetPathPoints(ccobj->GetCurve());
				ccobj->SetChanged();
				ccobj->SetColor(OpenMesh::Vec3d(0, 0, 0));
				DataPool::AddCurveObject(ccobj);*/

				if (cuttingpath_oftvhs_map[tvh].GetStraightLength() > cutting_path.GetStraightLength())
				{
					COpenMeshT::VertexHandle presvh = cuttingpath_oftvhs_map[tvh].start_vh_;
					std::vector<COpenMeshT::VertexHandle>tmp_dst(0);
					for (int k = 0; k < outside_cutting_vhs[i].size(); k++)
					{
						if (tvh.idx()!= outside_cutting_vhs[i][k].idx())
							tmp_dst.push_back(outside_cutting_vhs[i][k]);
					}
					std::vector<COpenMeshT::FaceHandle>tmpres_fhs;
					std::vector<OpenMesh::Vec3d> tmpres_bary_coords;
					CGeoAlg::ComputeGeodesicPath(meshobj, presvh, tmp_dst, tmpres_fhs, tmpres_bary_coords);
					COpenMeshT::VertexHandle tmp_tvh = CGeoBaseAlg::GetClosestVhFromFacePoint(mesh, tmpres_fhs.back(), tmpres_bary_coords.back());
					CCuttingPath tmpcutting_path(mesh, presvh, tmp_tvh, tmpres_fhs, tmpres_bary_coords);
					if (cuttingpath_oftvhs_map.find(tmp_tvh) == cuttingpath_oftvhs_map.end()|| cuttingpath_oftvhs_map[tmp_tvh].GetStraightLength() > tmpcutting_path.GetStraightLength())
					{
						cuttingpath_oftvhs_map[tmp_tvh] = tmpcutting_path;
				
						

						
					}

					cuttingpath_oftvhs_map[tvh] = cutting_path;
				/*	CCurveObject *ccobj = new CCurveObject();
					tmpcutting_path.GetPathPoints(ccobj->GetCurve());
					ccobj->SetChanged();
					ccobj->SetColor(OpenMesh::Vec3d(1, 0, 0));
					DataPool::AddCurveObject(ccobj);*/
					
				}
				else
				{
					std::vector<COpenMeshT::VertexHandle>tmp_dst(0);
					for (int k = 0; k <outside_cutting_vhs[i].size(); k+=1)
					{
						if (tvh.idx() != outside_cutting_vhs[i][k].idx())
							tmp_dst.push_back(outside_cutting_vhs[i][k]);
					}
					for (int k = 0; k < tmp_dst.size(); k++)
					{
						std::cerr << tmp_dst[k].idx() << " ";
					}
					std::cerr <<"//////////////"<< std::endl;
					std::vector<COpenMeshT::FaceHandle>tmpres_fhs;
					std::vector<OpenMesh::Vec3d> tmpres_bary_coords;
					CGeoAlg::ComputeGeodesicPath(meshobj, inside_cutting_vhs[i][j], tmp_dst, tmpres_fhs, tmpres_bary_coords);
					COpenMeshT::VertexHandle tmp_tvh = CGeoBaseAlg::GetClosestVhFromFacePoint(mesh, tmpres_fhs.back(), tmpres_bary_coords.back());
					CCuttingPath tmpcutting_path(mesh, inside_cutting_vhs[i][j], tmp_tvh, tmpres_fhs, tmpres_bary_coords);
					if (cuttingpath_oftvhs_map.find(tmp_tvh) == cuttingpath_oftvhs_map.end())
					{
						
						cuttingpath_oftvhs_map[tmp_tvh] = tmpcutting_path;

						
					}
					/*std::cerr << "tmpcutting_path " << tmpcutting_path.path_fhs_.size() << " tvh: " << tvh.idx() << " tmp_tvh: " << tmp_tvh.idx() << std::endl;
					CCurveObject *ccobj = new CCurveObject();
					tmpcutting_path.GetPathPoints(ccobj->GetCurve());
					ccobj->SetChanged();
					ccobj->SetColor(OpenMesh::Vec3d(0,0, 1));
					DataPool::AddCurveObject(ccobj);*/
				}
			}
			
			/*CCurveObject *p_curve_obj = new CCurveObject();
			std::vector<OpenMesh::Vec3d>&res_path = p_curve_obj->GetCurve();
			res_path.clear();
			for (int k = 0; k < res_fhs.size(); k++)
			{
				res_path.push_back(CGeoBaseAlg::GetFacePointFromBaryCoord(mesh, res_fhs[k], res_bary_coords[k]));
			}
			p_curve_obj->SetChanged();
			p_curve_obj->SetColor(OpenMesh::Vec3d(0, 1, 1));
			DataPool::AddCurveObject(p_curve_obj);*/
		}

		for (auto iter = cuttingpath_oftvhs_map.begin(); iter != cuttingpath_oftvhs_map.end(); iter++)
		{
			res_cuttingpath.push_back(iter->second);
		}

	}
}
int CDentalBaseAlg::TagToothByCuttingPath(COpenMeshT&mesh, std::vector<CCuttingPath>&cutting_pathes, std::vector<int>&res_tags)
{
	res_tags.resize(mesh.n_vertices(),-1);
	std::vector<bool>is_cutting_v(mesh.n_vertices(),false);
	std::vector<COpenMeshT::VertexHandle>cutting_vhs;
	for (int pi = 0; pi < cutting_pathes.size(); pi++)
	{
		CCuttingPath &path = cutting_pathes[pi];
		int num = path.path_fhs_.size();
		for (int i = 0; i <num ; i++)
		{
			auto fh = path.path_fhs_[i];
			for (auto fviter = mesh.fv_begin(fh); fviter != mesh.fv_end(fh); fviter++)
			{
				is_cutting_v[fviter->idx()] = true;
				cutting_vhs.push_back(fviter);
			}
		}
	}
	int tagid = 0;
	for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); viter++)
	{
		int vid = viter->idx();
		if (is_cutting_v[vid] == false)
		{
			CMorphlogicOperation::FloodFill(mesh, is_cutting_v, false, viter, tagid, res_tags);
			tagid++;
		}
	}
	for (int i = 0; i < cutting_vhs.size(); i++)
	{
		for (auto vviter = mesh.vv_begin(cutting_vhs[i]); vviter != mesh.vv_end(cutting_vhs[i]); vviter++)
		{
			if (res_tags[vviter->idx()] != -1)
			{
				res_tags[cutting_vhs[i].idx()] = res_tags[vviter->idx()];
				break;
			}
		}
	}
	return tagid;
}

void CDentalBaseAlg::MergeCuttingPointsByDis(COpenMeshT&mesh, std::vector<std::vector<COpenMeshT::VertexHandle>>&inside_vhs, std::vector<std::vector<COpenMeshT::VertexHandle>>&outside_vhs, double threshold)
{
	std::cerr << "merge 1" << std::endl;
	for (int i = 0; i < inside_vhs.size(); i++)
	{
		std::vector<double>lens(inside_vhs[i].size()-1,0);
		if (lens.size() <= 1)
			continue;
		CUnionFind union_find(lens.size());
		for (int j = 1; j < inside_vhs[i].size(); j++)
		{
			lens[j - 1] = (mesh.point(inside_vhs[i][j]) - mesh.point(inside_vhs[i][j - 1])).length();
		}
		bool flag = true;
		while (flag)
		{
			int del_id=-1;
			
			CNumericalBaseAlg::GetMinValue(lens, del_id);
			std::cerr << "len:" << lens[del_id]<<std::endl;
			if(lens[del_id]>=threshold)
				flag = false;
			else
			{
				std::cerr << "del_id:" << del_id << std::endl;
				if (del_id == 0)
				{
					std::cerr << "0\n" << std::endl;
					union_find.Union(1, del_id);
				}
				else if (del_id == lens.size() - 1)
				{
					std::cerr << "1\n" << std::endl;
					union_find.Union(del_id - 1, del_id);
				}
				else
				{
					std::cerr << "2\n" << std::endl;
					int lid = union_find.Find(del_id - 1);
					int rid = union_find.Find(del_id + 1);
					if (lens[lid] < lens[rid])
					{
						union_find.Union(lid, del_id);
					}
					else
					{
						union_find.Union(rid, del_id);
					}

				}
				int rootid = union_find.Find(del_id);
				lens[rootid] += lens[del_id];
				lens[del_id] = std::numeric_limits<double>::max();
			}
			
		}
		std::cerr << "end0" << std::endl;
		std::vector<COpenMeshT::VertexHandle>tmp_vhs;
		tmp_vhs.push_back(inside_vhs[i][0]);
		for (int j = 1; j < lens.size(); j++)
		{
			if (union_find.Find(j) != union_find.Find(j - 1))
			{
				tmp_vhs.push_back(inside_vhs[i][j]);
			}
		}
		tmp_vhs.push_back(inside_vhs[i].back());
		inside_vhs[i] = tmp_vhs;
	}
	std::cerr << "merge 2" << std::endl;
	for (int i = 0; i < outside_vhs.size(); i++)
	{
		std::vector<double>lens(outside_vhs[i].size() - 1, 0);
		if (lens.size() <= 1)
			continue;
		CUnionFind union_find(lens.size());
		for (int j = 1; j < outside_vhs[i].size(); j++)
		{
			lens[j - 1] = (mesh.point(outside_vhs[i][j]) - mesh.point(outside_vhs[i][j - 1])).length();
		}
		bool flag = true;
		while (flag)
		{
			int del_id = -1;

			CNumericalBaseAlg::GetMinValue(lens, del_id);
			if (lens[del_id] >= threshold)
				flag = false;
			else
			{
				if (del_id == 0)
				{
					union_find.Union(1, del_id);
				}
				else if (del_id == lens.size() - 1)
				{
					union_find.Union(del_id - 1, del_id);
				}
				else
				{
					int lid = union_find.Find(del_id - 1);
					int rid = union_find.Find(del_id + 1);
					if (lens[lid] < lens[rid])
					{
						union_find.Union(lid, del_id);
					}
					else
					{
						union_find.Union(rid, del_id);
					}

				}
				int rootid = union_find.Find(del_id);
				lens[rootid] += lens[del_id];
				lens[del_id] = std::numeric_limits<double>::max();
			}

		}
		std::vector<COpenMeshT::VertexHandle>tmp_vhs;
		tmp_vhs.push_back(outside_vhs[i][0]);
		for (int j = 1; j < lens.size(); j++)
		{
			if (union_find.Find(j) != union_find.Find(j - 1))
			{
				tmp_vhs.push_back(outside_vhs[i][j]);
			}
		}
		tmp_vhs.push_back(outside_vhs[i].back());
		outside_vhs[i] = tmp_vhs;
	}
}
void CDentalBaseAlg::ComputeTwoSideBoundsOfToothMesh(COpenMeshT&mesh, std::vector<std::vector<COpenMeshT::VertexHandle>>&res_inside_bounds, std::vector<std::vector<COpenMeshT::VertexHandle>>&res_outside_bounds)
{
	std::vector<std::vector<COpenMeshT::VertexHandle>>bounds;
	CGeoBaseAlg::GetBoundary(mesh, bounds);
	res_inside_bounds.clear();
	res_outside_bounds.clear();
	OpenMesh::Vec3d pca_mean;
	std::vector<OpenMesh::Vec3d>pca_frame;
	CDentalBaseAlg::ComputePCAFrameFromHighCurvaturePoints(mesh, 10, pca_mean, pca_frame);
	OpenMesh::Vec3d bounds_mean(0, 0, 0);
	int bound_count = 0;
	for (int i = 0; i < bounds.size(); i++)
	{
		bound_count += bounds[i].size();
		for (int j = 0; j < bounds[i].size(); j++)
		{
			bounds_mean += mesh.point(bounds[i][j]);
		}
	}
	bounds_mean /= bound_count;
	int ri = 0;
	for (int i = 0; i < bounds.size(); i++)
	{
		//std::cerr << "bound" << i << std::endl;
		
		std::vector<COpenMeshT::VertexHandle>&bound = bounds[i];
		std::vector<OpenMesh::Vec3d>bound_curve(bound.size());
		for (int j = 0; j < bound.size(); j++)
		{
			//mesh.set_color(bound[i], OpenMesh::Vec3d(1, 0, 0));
			bound_curve[j] = mesh.point(bound[j]);
		}
		std::vector<OpenMesh::Vec3d>proj_dirs(2);
		proj_dirs[0] = pca_frame[0];
		proj_dirs[1] = pca_frame[1];
		std::vector<OpenMesh::Vec2d>curve_2d;
		CCurveBaseAlg::ProjectCurve2Plannar(bound_curve, proj_dirs, curve_2d);
		std::vector<int>extreme_pids;
		ComputeExtremePointsOfClosedCurve(curve_2d, extreme_pids);
		
		if (extreme_pids.size() < 2)
			continue;
		res_inside_bounds.push_back(std::vector<COpenMeshT::VertexHandle>());
		res_outside_bounds.push_back(std::vector<COpenMeshT::VertexHandle>());
		for (int j = extreme_pids[0]; j < extreme_pids[1]; j++)
		{
			mesh.set_color(bound[j], OpenMesh::Vec3d(1, 0, 0));
			res_inside_bounds[ri].push_back(bound[j]);
		}
		for (int j = extreme_pids[1]; j != extreme_pids[0]; j = (j +1) % bound.size())
		{
			mesh.set_color(bound[j], OpenMesh::Vec3d(0, 1, 0));
			res_outside_bounds[ri].push_back(bound[j]);
		}
		double inside2mean_dis = 0;
		for (int j = 0; j < res_inside_bounds.size(); j++)
		{
			OpenMesh::Vec3d p = mesh.point(res_inside_bounds[ri][j]);
			inside2mean_dis += (p - bounds_mean).length();
		}
		inside2mean_dis /= res_inside_bounds[ri].size();
		double outside2mean_dis = 0;
		for (int j = 0; j < res_outside_bounds.size(); j++)
		{
			OpenMesh::Vec3d p = mesh.point(res_outside_bounds[ri][j]);
			outside2mean_dis += (p - bounds_mean).length();
		}
		outside2mean_dis /= res_outside_bounds[ri].size();

		if (outside2mean_dis < inside2mean_dis)
		{
			std::vector<COpenMeshT::VertexHandle>tmp;
			tmp = res_inside_bounds[ri];
			res_inside_bounds[ri] = res_outside_bounds[ri];
			res_outside_bounds[ri] = tmp;
		}
		//std::cerr << "bound " << i <<" end"<< std::endl;
		ri++;
	}

}
void CDentalBaseAlg::ComputeGingivaVhs(COpenMeshT &mesh, std::vector<COpenMeshT::VertexHandle>&res_vhs)
{
	OpenMesh::Vec3d mean;
	std::vector<OpenMesh::Vec3d>frame;
	ComputePCAFrameFromHighCurvaturePoints(mesh, 10, mean, frame);
	OpenMesh::Vec3d dir = frame[2];
	OpenMesh::Vec3d bound_mean(0, 0, 0);
	int bcount = 0;
	for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); viter++)
	{
		if (mesh.is_boundary(viter))
		{
			OpenMesh::Vec3d p = mesh.point(viter);
			bound_mean = bound_mean + p;
			bcount++;
		}
	}
	bound_mean = bound_mean / bcount;
	OpenMesh::Vec3d bound_dir = bound_mean - mean;
	if (bound_dir[0] * dir[0] + bound_dir[1] * dir[1] + bound_dir[2] * dir[2] > 0)
	{
		dir = -dir;
	}
	dir = dir.normalize();
	mean = mean - dir*0.02;
	res_vhs.clear();
	for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); viter++)
	{
		OpenMesh::Vec3d p = mesh.point(viter);
		p = p - mean;
		if (p[0] * dir[0] + p[1] * dir[1] + p[2] * dir[2] < 0)
		{
			res_vhs.push_back(viter);
		}
	}
}
void CDentalBaseAlg::ComputeVertexPenaltyWeight(COpenMeshT&mesh, std::vector<COpenMeshT::VertexHandle>&feature_points,Eigen::VectorXd &mean_curvature,Eigen::VectorXd& res_weight)
{

	res_weight.resize(mesh.n_vertices(), 0);

	std::vector<bool>is_feature_point_(mesh.n_vertices(), false);
	for (int i = 0; i < feature_points.size(); i++)
	{
		is_feature_point_[i] = true;
	}
	for (int i = 0; i < feature_points.size(); i++)
	{

		auto vh = feature_points[i];
		int vid = vh.idx();
		std::vector<double>ncurvatures;
		double vcurvature = mean_curvature(vid);
		int ncount = 0;
		double meancurv = 0;
		int nvcount = 0;
		for (auto vviter = mesh.vv_begin(vh); vviter != mesh.vv_end(vh); vviter++)
		{
			int nvid = vviter->idx();
			if (is_feature_point_[nvid])
			{
				ncount++;
				ncurvatures.push_back(mean_curvature(nvid));
				meancurv += mean_curvature(nvid);
			}
			nvcount++;
		}

		meancurv /= nvcount;
		double devicurv = CNumericalBaseAlg::ComputeStdDeviation(ncurvatures);
		double gaussian = CNumericalBaseAlg::ComputeGaussian(meancurv, devicurv, vcurvature);
		res_weight[vid] = gaussian*ncount / feature_points.size();


	}

}
void CDentalBaseAlg::ComputeTeethFeaturePointsUsingSmoothedMesh(COpenMeshT &mesh, std::vector<COpenMeshT::VertexHandle>&res_vhs)
{
	COpenMeshT tmp_mesh=mesh;
	CGeoAlg::LaplacianSmooth(tmp_mesh, 40, 0.5);
	std::vector<COpenMeshT::VertexHandle> tmp_vhs;
	ComputeTeethFeaturePoints(tmp_mesh, tmp_vhs);
	res_vhs.clear();
	for (int i = 0; i < tmp_vhs.size(); i++)
	{
		res_vhs.push_back(mesh.vertex_handle(tmp_vhs[i].idx()));
	}
}
void CDentalBaseAlg::ComputeTeethFeaturePoints(COpenMeshT &mesh, std::vector<COpenMeshT::VertexHandle>&res_vhs)
{
	std::cerr << "compute teeth feature points" << std::endl;
	
	Eigen::MatrixXd vertexs;
	Eigen::MatrixXi faces;
	CConverter::ConvertFromOpenMeshToIGL(mesh, vertexs, faces);
	Eigen::VectorXd curvatures_value;
	CGeoBaseAlg::ComputeMeanCurvatureValue(mesh, curvatures_value,true);
	//Eigen::VectorXd gaussian_curv;
//	igl::gaussian_curvature(vertexs, faces, curvatures_value);
	/*double max_curv = -1;
	for (int i = 0; i < curvatures_value.rows(); i++)
	{
		if (max_curv < std::abs(curvatures_value(i)))
		{
			max_curv = std::abs(curvatures_value(i));
		}
	}*/
//	curvatures_value = curvatures_value / max_curv;
//	igl::writeDMAT("curvdmat0.dmat", curvatures_value);
	for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); viter++)
	{
		if (mesh.is_boundary(viter))
		{
			curvatures_value(viter->idx()) = 0;
		}
	}
	//CNumericalBaseAlg::NormalizeScalarField(curvatures_value);
	for (int i = 0; i < curvatures_value.size(); i++)
	{
		curvatures_value(i) = CNumericalBaseAlg::Sigmoid(curvatures_value(i));
	}
	//curvatures_value /= 300;
	//igl::writeDMAT("curvdmat.dmat", curvatures_value);
	OpenMesh::Vec3d pca_mean;
	std::vector<OpenMesh::Vec3d>pca_frame;
	ComputePCAFrameFromHighCurvaturePoints(mesh, 10, pca_mean, pca_frame);
	OpenMesh::Vec3d bound_mean(0, 0, 0);
	int boundcount = 0;
	for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); viter++)
	{
		if (mesh.is_boundary(viter))
		{
			boundcount++;
			bound_mean = bound_mean + mesh.point(viter);
		}	
	}
	bound_mean = bound_mean / boundcount;
	auto bound_pca_dir = (bound_mean - pca_mean);
	OpenMesh::Vec3d dir = pca_frame[2];
	if (bound_pca_dir[0]*dir[0]+bound_pca_dir[1]*dir[1]+bound_pca_dir[2]*dir[2] > 0)
	{
		dir = -dir;
	}
	dir = dir.normalize();
	std::vector<double>z_values(mesh.n_vertices());
	std::vector<double>height_f(mesh.n_vertices());
	
	
	for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); viter++)
	{
		int vid = viter->idx();
		OpenMesh::Vec3d p = mesh.point(viter);
		p = p - pca_mean;
		z_values[vid] = p[0] * dir[0] + p[1] * dir[1] + p[2] * dir[2];
		
	}
	std::vector<double>tmp_z = z_values;
	CNumericalBaseAlg::NormalizeScalarField(z_values);
	double alpha = 1;
	for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); viter++)
	{
		int vid = viter->idx();
		height_f[vid] = alpha*z_values[vid] + (1 - alpha)*(curvatures_value(vid));
	}
	CNumericalBaseAlg::NormalizeScalarField(height_f);
	//igl::writeDMAT("heightdmat.dmat", height_f);
	for (auto hiter = mesh.halfedges_begin(); hiter != mesh.halfedges_end(); hiter++)
	{
		auto vh = mesh.to_vertex_handle(hiter);
		mesh.data(*hiter).SetUV(OpenMesh::Vec2f(height_f[vh.idx()], height_f[vh.idx()]));
	}
	res_vhs.clear();
	std::vector<CVertexValue>vhs_v;
	CGeoAlg::ComputeLocalExtremum(mesh, height_f, 2, res_vhs);
	for (int i = 0; i < res_vhs.size(); i++)
	{
		vhs_v.push_back(CVertexValue(res_vhs[i], height_f[res_vhs[i].idx()]));
	}
	std::sort(vhs_v.begin(), vhs_v.end());
	res_vhs.clear();
	int maxi;
	CNumericalBaseAlg::GetMaxValue(z_values, maxi);

	std::cerr << "feature size " << vhs_v.size() << std::endl;
	for (int i=0;i<vhs_v.size();i++)
	{
		if (z_values[maxi] - z_values[vhs_v[i].vh_.idx()] < 0.35)
		{
			res_vhs.push_back(vhs_v[i].vh_);
		}
	}
	std::cerr << "res feature size " << res_vhs.size() << std::endl;
	

}