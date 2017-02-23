#include"panoramic_simulation.h"
#include"../AlgColle/geo_base_alg.h"
#include"../AlgColle/curve_base_alg.h"
#include"../DataColle/data_pool.h"
#include"../AlgColle/geo_alg.h"

void CPanoramicProjectorBase::ComputeProjParams(std::vector<OpenMesh::Vec3d>&pts, std::vector<OpenMesh::Vec2d>&res_params)
{
	res_params.resize(pts.size());
	for (int i = 0; i < res_params.size(); i++)
	{
		res_params[i] = ComputeProjParam(pts[i]);
	}
}
void CPanoramicProjectorBase::ProjectMesh(CMeshObject* mesh_obj, std::vector<OpenMesh::Vec3d> &res_pts, std::vector<COpenMeshT::VertexHandle > &res_orig_vhs)
{
	res_pts.clear();
	res_orig_vhs.clear();

	COpenMeshT &mesh = mesh_obj->GetMesh();
	for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); viter++)
	{
		OpenMesh::Vec3d p = mesh_obj->TransformPointByLocalMatrix(mesh.point(viter));
		p = Project(p);
		res_pts.push_back(p);
		res_orig_vhs.push_back(viter);
	}
}
void CPanoramicProjectorBase::ProjectMeshes(std::vector<CMeshObject*>mesh_objs, std::vector<std::vector<OpenMesh::Vec3d>>&res_pts, std::vector<std::vector<COpenMeshT::VertexHandle>>&res_orig_vhs)
{
	res_pts.resize(mesh_objs.size());
	res_orig_vhs.resize(mesh_objs.size());
	for (int i = 0; i < mesh_objs.size(); i++)
	{
		ProjectMesh(mesh_objs[i], res_pts[i], res_orig_vhs[i]);
	}
}

void CCurveSurfaceProjector::SetChanged()
{
	curve_center_ = OpenMesh::Vec3d(0, 0, 0);
	for (int i = 0; i < curve_.size(); i++)
	{
		curve_center_ += curve_[i];
	}
	curve_center_ /= curve_.size();
	proj_curve_.clear();
	for (int i = 0; i < curve_.size(); i++)
	{
		OpenMesh::Vec3d tmp_dir = curve_[i] - curve_center_;
		double zlen = OpenMesh::dot(tmp_dir, updir_);
		proj_curve_.push_back(curve_[i] - zlen*updir_);
	}
	updir_.normalize();
}
OpenMesh::Vec3d CCurveSurfaceProjector::Project(OpenMesh::Vec3d p)
{
	
	OpenMesh::Vec3d proj_p;

	OpenMesh::Vec3d tmp_dir = p - curve_center_;
	double zlen = OpenMesh::dot(tmp_dir, updir_);
	proj_p = p - zlen*updir_;
	OpenMesh::Vec3d closest_p;
	int closest_seg_id;
	CCurveBaseAlg::ComputeClosestPoint(proj_curve_, false, proj_p, closest_p, closest_seg_id);

	return closest_p + zlen*updir_;
}

OpenMesh::Vec2d CCurveSurfaceProjector::ComputeProjParam(OpenMesh::Vec3d p)
{
	




	OpenMesh::Vec3d proj_p;

	OpenMesh::Vec3d tmp_dir = p - curve_center_;
	double zlen = OpenMesh::dot(tmp_dir, updir_);
	proj_p = p - zlen*updir_;
	OpenMesh::Vec3d closest_p;
	int closest_seg_id;
	CCurveBaseAlg::ComputeClosestPoint(proj_curve_, false, proj_p, closest_p, closest_seg_id);
	double x = 0;
	for (int i = 1; i <=closest_seg_id; i++)
	{
		OpenMesh::Vec3d diff = proj_curve_[i] - proj_curve_[i - 1];
		x += diff.length();
	}
	OpenMesh::Vec3d diff = proj_curve_[closest_seg_id] - closest_p;
	x += diff.length();
	return OpenMesh::Vec2d(x, zlen);
}

void CCurveSurfaceProjector::SetParams(OpenMesh::Vec3d updir, std::vector<OpenMesh::Vec3d>&curve)
{
	updir_ = updir;
	curve_ = curve;
	updir_.normalize();
}
void CCurveSurfaceProjector::SetParams(std::vector<double>&params)
{
	double len= std::sqrt(params[0] * params[0] + params[1] * params[1]);
	if (len > 1)
	{
		params[0] / len;
		params[1] /= len;
	}
	updir_[0] = params[0];
	updir_[1] = params[1];
	updir_[2] = std::sqrt(1 - params[0] * params[0] - params[1] * params[1]);
}
void CCurveSurfaceProjector::GetParams(std::vector<double>&params)
{
	params.resize(2);
	params[0] = updir_[0];
	params[1] = updir_[1];
}



void CCircularSurfaceProjector::SetParams(std::vector<double>&params)
{

	double len = std::sqrt(params[0] * params[0] + params[1] * params[1]);
	if (len > 1)
	{
		params[0] / len;
		params[1] /= len;
	}
	updir_[0] = params[0];
	updir_[1] = params[1];
	updir_[2] = std::sqrt(1 - params[0] * params[0] - params[1] * params[1]);

	center_[0] = params[2];
	center_[1] = params[3];
	center_[2] = params[4];

	radius_ = params[5];

	

}
void CCircularSurfaceProjector::GetParams(std::vector<double>&params)
{
	params.resize(6);
	params[0] = updir_[0];
	params[1] = updir_[1];

	params[2] = center_[0];
	params[3] = center_[1];
params[4] = center_[2];

params[5] = radius_;



}

void CCircularSurfaceProjector::SetParams(OpenMesh::Vec3d updir, OpenMesh::Vec3d center, OpenMesh::Vec3d start_dir, double radius)
{
	updir_ = updir;
	center_ = center;

	radius_ = radius;
	updir_.normalize();

	OpenMesh::Vec3d p = center + start_dir;
	OpenMesh::Vec3d tmp_dir = p - center_;
	double zlen = OpenMesh::dot(tmp_dir, updir_);
	p = p - zlen*updir_;

	start_dir_ = p - center;
	start_dir_.normalize();

}
OpenMesh::Vec3d CCircularSurfaceProjector::Project(OpenMesh::Vec3d p)
{
	OpenMesh::Vec3d tmp_center;

	OpenMesh::Vec3d tmp_dir = p - center_;
	double zlen = OpenMesh::dot(tmp_dir, updir_);
	tmp_center = center_ + zlen*updir_;
	tmp_dir = p - tmp_center;
	tmp_dir.normalize();
	return tmp_center + tmp_dir*radius_;
}

OpenMesh::Vec2d CCircularSurfaceProjector::ComputeProjParam(OpenMesh::Vec3d p)
{

	OpenMesh::Vec3d proj_p;

	OpenMesh::Vec3d tmp_dir = p - center_;


	double zlen = OpenMesh::dot(tmp_dir, updir_);
	proj_p = p - zlen*updir_;

	tmp_dir = proj_p - center_;

	tmp_dir.normalize();
	double degree = std::acos(OpenMesh::dot(start_dir_, tmp_dir));
	double x_param = degree / (2.0*M_PI);

	if (OpenMesh::dot(OpenMesh::cross(tmp_dir, start_dir_), updir_) < 0)
	{
		x_param = 1.0 - x_param;
	}
	return OpenMesh::Vec2d(x_param * 2 * M_PI*radius_, zlen);
}
void CPanoramicSimulation::GeneratePanoramicBoundPointsByCircle(std::vector<CMeshObject *>crown_objs, OpenMesh::Vec3d updir, OpenMesh::Vec3d center, double radius, std::vector<std::vector<OpenMesh::Vec2d>>&res_bound_pano_params, std::vector<std::vector<COpenMeshT::VertexHandle>>&res_corres_vhs)
{
	std::vector<std::vector<OpenMesh::Vec2d>>params;
	std::vector<std::vector<COpenMeshT::VertexHandle>>vhs;
	GeneratePanoramicWholeBoundPointsByCircle(crown_objs, updir, center, radius, params, vhs);
	res_bound_pano_params.resize(crown_objs.size());
	res_corres_vhs.resize(crown_objs.size());
	for (int i = 0; i < params.size(); i++)
	{

		//std::vector<OpenMesh::Vec2d>normals;
		////std::cerr << "paramsize " << params[i].size()<<" "<<i<<" "<< params[i].size() << std::endl;
		//CCurveBaseAlg::ComputeCurveNormals(params[i], normals);
		OpenMesh::Vec2d center = CCurveBaseAlg::ComputeCurveCenter(params[i]);
	/*	double count_up = 0, count_down = 0;
		for (int j = 0; j < params[i].size(); j++)
		{
			OpenMesh::Vec2d dif = params[i][j] - center;
			dif.normalize();
			if (crown_objs[i]->GetMesh().is_boundary(vhs[i][j]) && OpenMesh::dot(dif, OpenMesh::Vec2d(1, 0)) < 0)

			{
				count_down++;

			}
			else
			{

				count_up++;


			}

		}
		bool is_upper_jaw = true;
		if (count_down > count_up)
		{
			is_upper_jaw = false;
		}*/
		std::vector<int>sub_len;
		sub_len.resize(params[i].size(), 0);
		for (int j = 0; j < params[i].size(); j++)
		{
			OpenMesh::Vec2d dif = params[i][j] - center;
			dif.normalize();
		/*	if (crown_objs[i]->GetMesh().is_boundary(vhs[i][j]))
				continue;*/
		/*	if (OpenMesh::dot(dif, OpenMesh::Vec2d(1, 0)) >0.0 && crown_objs[i]->GetMesh().is_boundary(vhs[i][j]) && is_upper_jaw)
				continue;
			if (OpenMesh::dot(dif, OpenMesh::Vec2d(-1, 0)) >0.0 && crown_objs[i]->GetMesh().is_boundary(vhs[i][j]) && (!is_upper_jaw))
				continue;*/
			if ((OpenMesh::dot(dif, OpenMesh::Vec2d(1, 0)) >0.75|| OpenMesh::dot(dif, OpenMesh::Vec2d(-1, 0)) >0.75) && crown_objs[i]->GetMesh().is_boundary(vhs[i][j]) )
				continue;
			
			int prej = (j - 1 + params[i].size()) % params[i].size();
			sub_len[j] = sub_len[prej]+1;
			
		}


		


		if (sub_len.back() != 0)
		{
			for (int j = 0; j < sub_len.size(); j++)
			{
				if (sub_len[j] != 0)
				{
					sub_len[j] += sub_len.back();
				}
				else
				{
					break;
				}
			}
		}
		int max_len = -1;
		int mi = 0;
		for (int j = 0; j < sub_len.size(); j++)
		{
			if (max_len < sub_len[j])
			{
				mi = j;
				max_len = sub_len[j];
			}
		}
	
		res_bound_pano_params[i].clear();
		res_corres_vhs[i].clear();
		int tmi=mi;
		for (int j = mi; j < mi + 5 && j < sub_len.size(); j++)
		{
			sub_len[j] = 1;
			tmi = j;
		}
		mi = tmi;
		while (sub_len[mi] != 0)
		{
			res_bound_pano_params[i].push_back(params[i][mi]);
			res_corres_vhs[i].push_back(vhs[i][mi]);
			mi = (mi - 1 + sub_len.size()) % sub_len.size();

			if (res_corres_vhs[i].size() == params[i].size())
				break;
		}
		for (int j = mi; j >= mi - 5&&j>=0; j--)
		{
			res_bound_pano_params[i].push_back(params[i][j]);
			res_corres_vhs[i].push_back(vhs[i][j]);
		}
	
	}
}

void CPanoramicSimulation::GeneratePanoramicWholeBoundPointsByCircle(std::vector<CMeshObject *>crown_objs, OpenMesh::Vec3d updir, OpenMesh::Vec3d center, double radius, std::vector<std::vector<OpenMesh::Vec2d>>&res_bound_pano_params, std::vector<std::vector<COpenMeshT::VertexHandle>>&res_corres_vhs )
{
	std::vector<std::vector<OpenMesh::Vec2d>>params;
	std::vector<std::vector<COpenMeshT::VertexHandle>>vhs;
	GeneratePanoramicPointsByCircle(crown_objs, updir, center, radius, params, vhs);
	res_bound_pano_params.resize(crown_objs.size());
	res_corres_vhs.resize(crown_objs.size());
	for (int i = 0; i < params.size(); i++)
	{
		double max_edge_len = CGeoBaseAlg::ComputeMaxEdgeLength(crown_objs[i]);
		std::vector<std::vector<int>>vids;
		CGeoAlg::AlphaShape2d(params[i], std::pow(max_edge_len, 2)*20, vids);
		if (vids.size() == 0)
		{
			std::cerr << "compute alpha shape error " << std::endl;
			continue;
		}
		res_bound_pano_params[i].clear();
		res_corres_vhs[i].clear();
		int max_size = -1;
		int mi = 0;
		for (int j=0; j < vids.size(); j++)
		{
			if (max_size < vids[j].size())
			{
				max_size = vids[j].size();
				mi = j;
			}
		}
		for (int j = 0; j < vids[mi].size(); j++)
		{
			res_bound_pano_params[i].push_back(params[i][vids[mi][j]]);
			res_corres_vhs[i].push_back(vhs[i][vids[mi][j]]);
		}
	}
	

}
CCircularSurfaceProjector CPanoramicSimulation::ConstructCircularProjector(std::vector<CMeshObject *>crown_objs, OpenMesh::Vec3d updir, OpenMesh::Vec3d center, double radius)
{
	std::vector<OpenMesh::Vec3d>center_pts;

	for (int i = 0; i < crown_objs.size(); i++)
	{
		crown_objs[i]->RestoreCurrentVPos();
		crown_objs[i]->ApplyTransform();
		crown_objs[i]->SetChanged(true);
		OpenMesh::Vec3d center_p = CGeoBaseAlg::ComputeMeshCenter(crown_objs[i]->GetMesh());
		center_pts.push_back(center_p);

	}
	OpenMesh::Vec3d start_dir;
	//using convex hull to compute front and back of dental arch
	{


		OpenMesh::Vec3d res_mean;
		std::vector<OpenMesh::Vec3d>res_frame;
		std::vector<double>res_eg_values;
		CGeoAlg::PointSetPCA3D(center_pts, res_mean, res_frame, res_eg_values);


		std::vector<OpenMesh::Vec2d>center_pts2d;
		CCurveBaseAlg::ProjectCurve2Plannar(center_pts, res_frame, center_pts2d);
		std::vector<int>convex_hull_idx;
		CGeoAlg::ComputeConvexHull(center_pts2d, convex_hull_idx);
		double max_dis = std::numeric_limits<double>::min();
		int mi = 0;
		//CCurveObject *co = new CCurveObject();
		//std::cerr << "convex hull num " << convex_hull_idx.size() << std::endl;
		for (int i = 0; i < convex_hull_idx.size(); i++)
		{
			int j = (i + 1) % convex_hull_idx.size();
			OpenMesh::Vec3d tmp_diff = center_pts[convex_hull_idx[i]] - center_pts[convex_hull_idx[j]];
			double dis = tmp_diff.length();
			if (dis > max_dis)
			{
				max_dis = dis;
				mi = i;
			}
			/*	co->GetCurve().push_back(center_pts[convex_hull_idx[i]]);
			co->SetColor(OpenMesh::Vec3d(0, 1, 0));
			DataPool::AddCurveObject(co);*/
		}


		OpenMesh::Vec3d p = center_pts[convex_hull_idx[mi]] + center_pts[convex_hull_idx[(mi + 1) % convex_hull_idx.size()]];
		start_dir = p - center;
	}



	//std::cerr << "center " << center << std::endl;
	CCircularSurfaceProjector cp;

	cp.SetParams(updir, center, start_dir, radius);
	return cp;
}
void CPanoramicSimulation::GeneratePanoramicPointsByCircle(std::vector<CMeshObject *>crown_objs, OpenMesh::Vec3d updir, OpenMesh::Vec3d center,double radius, std::vector<std::vector<OpenMesh::Vec2d>>&res_pano_params, std::vector<std::vector<COpenMeshT::VertexHandle>>&res_corres_vhs)
{

	std::vector<std::vector<OpenMesh::Vec3d>>projected_points;
	CCircularSurfaceProjector cp = ConstructCircularProjector(crown_objs, updir, center, radius);
	cp.ProjectMeshes(crown_objs, projected_points,res_corres_vhs);
	res_pano_params.resize(projected_points.size());
	for (int i = 0; i < projected_points.size(); i++)
	{

		//std::cerr << "ComputeParamOnCurveSurface" << std::endl;
		//ComputeParamOnCurveSurface(cs, projected_points[i], res_pano_params[i]);
		cp.ComputeProjParams(projected_points[i], res_pano_params[i]);
		//std::cerr << "end ComputeParamOnCurveSurface" << std::endl;
	}
	/*double  min_x = std::numeric_limits<double>::max();
	for (int i = 0; i < res_pano_params.size(); i++)
	{
		for (int j = 0; j < res_pano_params[i].size(); j++)
		{
			if (min_x > res_pano_params[i][j][0])
			{
				min_x = res_pano_params[i][j][0];
			}
		}
	}
	for (int i = 0; i < res_pano_params.size(); i++)
	{
		for (int j = 0; j < res_pano_params[i].size(); j++)
		{
			res_pano_params[i][j][0] -= min_x;
		}
	}*/
	
/*	srand((unsigned)time(NULL));
	for (int i = 0; i < res_pano_params.size(); i++)
	{
		double max_edge_len = CGeoBaseAlg::ComputeMaxEdgeLength(crown_objs[i]);
		std::vector<std::vector<int>>vids;
		CGeoAlg::AlphaShape2d(res_pano_params[i],std::pow( max_edge_len/2,2), vids);
		for (int j = 0; j < vids.size(); j++)
		{
			std::vector<OpenMesh::Vec3d>p3d;
			OpenMesh::Vec3d color;
			color[0] = rand() % 1000 / 1000.0;
			color[1] = rand() % 1000 / 1000.0;
			color[2] = rand() % 1000 / 1000.0;
			for (int t = 0; t < vids[j].size(); t++)
			{
				int vid = vids[j][t];
				p3d.push_back(OpenMesh::Vec3d(res_pano_params[i][vid][0], res_pano_params[i][vid][1], 0));
			}
			CCurveObject *co = new CCurveObject();
			co->SetCurve(p3d);
			co->RendereType() = CCurveObject::CurveType::Dots;
			co->SetColor(color);
			co->SetChanged();
			DataPool::AddCurveObject(co);
		}*/
		/*std::cerr << "max edge " << max_edge_len << std::endl;
		std::cerr <<"alpha "<< std::pow(max_edge_len / 2, 2) << std::endl;*/

		/*OpenMesh::Vec3d color;
		color[0] = rand() % 1000 / 1000.0;
		color[1] = rand() % 1000 / 1000.0;
		color[2] = rand() % 1000 / 1000.0;
		std::vector<OpenMesh::Vec3d>p3d;
	
		for (int j = 0; j < res_pano_params[i].size(); j++)
		{
			p3d.push_back(OpenMesh::Vec3d(res_pano_params[i][j][0], res_pano_params[i][j][1], 0));
			
		}
		CCurveObject *co = new CCurveObject();
		co->SetCurve(p3d);
		co->RendereType() = CCurveObject::CurveType::Dots;
		co->SetColor(color);
		co->SetChanged();
		//DataPool::AddCurveObject(co);
	}*/
	


	/*srand((unsigned)time(NULL));
	for (int i = 0; i < projected_points.size(); i++)
	{
		OpenMesh::Vec3d color;
		color[0] = rand() % 1000 / 1000.0;
		color[1] = rand() % 1000 / 1000.0;
		color[2] = rand() % 1000 / 1000.0;
		CCurveObject *co = new CCurveObject();
		co->SetCurve(projected_points[i]);
		co->RendereType() = CCurveObject::CurveType::Dots;
		co->SetColor(color);
		co->SetChanged();
		DataPool::AddCurveObject(co);
	}*/
	

}
void CPanoramicSimulation::GeneratePanoramicImageByDentalArch(std::vector<CMeshObject *>crown_objs, std::vector<OpenMesh::Vec3d>&frame,int img_width, cv::Mat &res_panoramic)
{
	//for the convenience of computation, store the current mesh pos ,set matrix of meshes to identity 
	for (int i = 0; i < crown_objs.size(); i++)
	{
		crown_objs[i]->RestoreCurrentVPos();
		crown_objs[i]->ApplyTransform();
		crown_objs[i]->SetChanged(true);
	}

	OpenMesh::Vec3d crown_center = OpenMesh::Vec3d(0, 0, 0);
	std::vector<OpenMesh::Vec3d>crown_centers;
	for (int i = 0; i < crown_objs.size(); i++)
	{
		OpenMesh::Vec3d center = CGeoBaseAlg::ComputeMeshCenter(crown_objs[i]->GetMesh());
		crown_centers.push_back(center);
		crown_center += center;
	}
	crown_center /= crown_objs.size();
	
	std::vector<OpenMesh::Vec3d>tgt_frame;
	tgt_frame.resize(3);
	
	OpenMesh::Vec3d tgt_updir(0, 1, 0);
	tgt_frame[0] = OpenMesh::Vec3d(1, 0, 0);
	tgt_frame[1] = tgt_updir;
	tgt_frame[2] = OpenMesh::Vec3d(0, 0, 1);
	OpenMesh::Vec3d tgt_center(0, 0, 0); 

	Eigen::Matrix4d  frame_trans_mat=CGeoBaseAlg::ComputeFrameTransMatrix(crown_center, frame, tgt_center, tgt_frame);
	
	for (int i = 0; i < crown_objs.size(); i++)
	{
		crown_objs[i]->GetMatrix() = frame_trans_mat*crown_objs[i]->GetMatrix();
		crown_objs[i]->ApplyTransform();
		crown_objs[i]->SetChanged(true);
		crown_centers[i]= CGeoBaseAlg::Transform(frame_trans_mat, crown_centers[i]);
	}
	crown_center = CGeoBaseAlg::Transform(frame_trans_mat, crown_center);

	
	double min_x = std::numeric_limits<double>::max();
	double max_x = std::numeric_limits<double>::min();

	for (int i = 0; i < crown_objs.size(); i++)
	{
		COpenMeshT crown_mesh = crown_objs[i]->GetMesh();
		for (auto viter = crown_mesh.vertices_begin(); viter != crown_mesh.vertices_end(); viter++)
		{
			OpenMesh::Vec3d p = crown_mesh.point(viter);
		////	p=CGeoBaseAlg::Transform(frame_trans_mat, p);
		//	crown_mesh.point(viter) = p;
			if (p[0] < min_x)
			{
				min_x = p[0];
			}
			if (p[0] > max_x)
			{
				max_x = p[0];
			}
		}
	}


	std::vector<OpenMesh::Vec2d>crown_centers_2d;

	std::vector<OpenMesh::Vec2d>crown_pts_2d;
	crown_centers_2d.resize(crown_objs.size());
	for (int i = 0; i < crown_objs.size(); i++)
	{
		//crown_centers[i]= CGeoBaseAlg::Transform(frame_trans_mat, crown_centers[i]);
		crown_centers_2d[i] = OpenMesh::Vec2d(crown_centers[i][0], crown_centers[i][2]);
		CMeshObject *crown = crown_objs[i];
		for (auto viter = crown->GetMesh().vertices_begin(); viter != crown->GetMesh().vertices_end(); viter++)
		{
			OpenMesh::Vec3d p=crown->GetMesh().point(viter);
			crown_pts_2d.push_back(OpenMesh::Vec2d(p[0], p[2]));
		}
	}
	std::vector<double>poly_coeffs;
	CCurveBaseAlg::PolynomialFitting(crown_pts_2d, 4, poly_coeffs);
	CCurveObject *cpolynomial = new CCurveObject();

	for (double x = min_x-0.5; x < max_x+0.5; x+=0.001)
	{
		OpenMesh::Vec2d p2d = CCurveBaseAlg::ComputePointOfPolynomial(poly_coeffs, x);
		OpenMesh::Vec3d p(x,crown_center[2], p2d[1]);
		cpolynomial->GetCurve().push_back(p);
	}
	DataPool::AddCurveObject(cpolynomial);
	

	/*CCurveObject *cdot = new CCurveObject();
	for (int i = 0; i < crown_centers.size(); i++)
	{
		cdot->GetCurve().push_back(crown_centers[i]);
	}
	cdot->SetChanged();
	cdot->SetColor(OpenMesh::Vec3d(0, 1, 0));
	DataPool::AddCurveObject(cdot);*/
	std::vector<OpenMesh::Vec3d>tmp_curve= cpolynomial->GetCurve();
	CCurveSurfaceProjector cs;
	cs.SetParams(OpenMesh::Vec3d(0, 1, 0), tmp_curve);
	
	cs.SetChanged();
	std::vector<std::vector<OpenMesh::Vec3d>>projected_points;
	std::vector<std::vector<OpenMesh::Vec2d>>proj_params;
	cs.ProjectMeshes(crown_objs, projected_points);
	//Project2CurveSurface(cs, crown_objs, projected_points);
	proj_params.resize(projected_points.size());
	for (int i = 0; i < projected_points.size(); i++)
	{
		
		std::cerr << "ComputeParamOnCurveSurface" << std::endl;
		//ComputeParamOnCurveSurface(cs, projected_points[i], proj_params[i]);
		cs.ComputeProjParams(projected_points[i], proj_params[i]);
		std::cerr << "end ComputeParamOnCurveSurface" << std::endl;
	}
	/*srand((unsigned)time(NULL));
	for (int i = 0; i < projected_points.size(); i++)
	{
		OpenMesh::Vec3d color;
		color[0] = rand() % 1000 / 1000.0;
		color[1] = rand() % 1000 / 1000.0;
		color[2] = rand() % 1000 / 1000.0;
		std::vector<OpenMesh::Vec3d>p3d;
		for (int j = 0; j < proj_params[i].size(); j++)
		{
			p3d.push_back(OpenMesh::Vec3d(proj_params[i][j][0], proj_params[i][j][1], 0));
		}
		CCurveObject *co = new CCurveObject();
		co->SetCurve(p3d);
		co->RendereType() = CCurveObject::CurveType::Dots;
		co->SetColor(color);
		co->SetChanged();
		DataPool::AddCurveObject(co);
	}


	srand((unsigned)time(NULL));
	for (int i = 0; i < projected_points.size(); i++)
	{
		OpenMesh::Vec3d color;
		color[0] = rand() % 1000 / 1000.0;
		color[1] = rand() % 1000 / 1000.0;
		color[2] = rand() % 1000 / 1000.0;
		CCurveObject *co = new CCurveObject();
		co->SetCurve(projected_points[i]);
		co->RendereType() = CCurveObject::CurveType::Dots;
		co->SetColor(color);
		co->SetChanged();
		DataPool::AddCurveObject(co);
	}*/
	

}
