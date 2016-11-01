#include"morph_skel_dental_mesh_seg.h"
#include"geo_base_alg.h"
#include"geo_alg.h"
#include"../DataColle/data_pool.h"
#include"../DataColle/aux_geo_utils.h"
#include"../DataColle/Polyhedron_type.h"
#include"numerical_base_alg.h"
#include"geo_base_alg.h"
#include"../DataColle/data_io.h"
#include"../DataColle/cgal_igl_converter.h"
#include"morphlogic_operation.h"
#include<queue>
CMorphSkelDentalMeshSeg::CMorphSkelDentalMeshSeg(CMeshObject &mesh_obj):mesh_obj_(mesh_obj)
{
	curvature_threshold_ =60;
	small_region_threshold_ = 10;
	small_teeth_region_percent_threshold_ = 0.01;
}
CMorphSkelDentalMeshSeg::~CMorphSkelDentalMeshSeg()
{

}
void CMorphSkelDentalMeshSeg::SetVertexTags(std::vector<COpenMeshT::VertexHandle>&vhs, std::vector<CMorphSkelDentalMeshSeg::CTag>&tags)
{

	if (tags_.size() != mesh_obj_.GetMesh().n_vertices());
		tags_.resize(mesh_obj_.GetMesh().n_vertices(),CTag::Non);
		if (is_edge_point_.size() != tags_.size())
			is_edge_point_.resize(tags.size(), false);

	for (int i = 0; i < vhs.size(); i++)
	{
		tags_[vhs[i].idx()] = tags[i];
		if(tags[i]==CTag::Feature)
		is_edge_point_[vhs[i].idx()] = true;
	}
	edge_points_.clear();
	edge_points_id_.clear();
	COpenMeshT &mesh = mesh_obj_.GetMesh();
	for (int i = 0; i < is_edge_point_.size(); i++)
	{
		if (is_edge_point_[i])
		{
			edge_points_.push_back(mesh.points()[i]);
			edge_points_id_.push_back(mesh.vertex_handle(i));
		}
	}
}
void CMorphSkelDentalMeshSeg::ComputeSegmentation(bool verbose)
{
	verbose_ = verbose;
	double a, b, c, d;
	ComputeCuttingPlane();
	mesh_obj_.SetAttrChanged();
}
void CMorphSkelDentalMeshSeg::TestRender()
{
	COpenMeshT &mesh = mesh_obj_.GetMesh();
	for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); viter++)
	{
		int vid = viter->idx();
		if (tags_[viter->idx()]==CTag::Feature)
		{
			mesh.set_color(viter, OpenMesh::Vec3d(1, 0, 0));
		}
		else if(tags_[viter->idx()] == CTag::Gingiva)
		{
			mesh.set_color(viter, OpenMesh::Vec3d(0, 1, 0));
		}
		else if (tags_[viter->idx()] == CTag::Teeth)
		{
			mesh.set_color(viter, OpenMesh::Vec3d(0, 0, 1));
		}
		else if (tags_[viter->idx()] == CTag::Base)
		{
			mesh.set_color(viter, OpenMesh::Vec3d(1, 1,0));
		}
		else
		{
			mesh.set_color(viter, OpenMesh::Vec3d(0.8, 0.8, 0.8));
		}

	}
	mesh_obj_.SetAttrChanged();
}
void CMorphSkelDentalMeshSeg::TestCurvature()
{
	ComputeEdgePointsFromMeanCurvatureThreshold(curvature_threshold_);
	for (int i = 0; i < is_edge_point_.size(); i++)
	{
		if (is_edge_point_[i])
			tags_[i] = CTag::Feature;
		else if (tags_[i] != CTag::Base)
			tags_[i] = CTag::Non;

	}
	TestRender();

}
void CMorphSkelDentalMeshSeg::TestTagGingiva()
{
	TagGingiva();
	TestRender();
	
}
void CMorphSkelDentalMeshSeg::TestRemoveSmallFeatureRegions()
{
	ResetEdgeVertexMarkFromTags();
	RemoveGingiva2GingivaFeatureEdge();
	//ResetEdgeVertexMarkFromTags();
	RemoveSmallIsolateTeethRegion();
	
	RemoveSmallFeatureRegions();
	
	TestRender();
	
}
void CMorphSkelDentalMeshSeg::TagGingiva()
{
	COpenMeshT &mesh = mesh_obj_.GetMesh();
	std::vector<bool>mmark(mesh.n_vertices(), false);
	for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); viter++)
	{
		if (mesh.is_boundary(viter)&& mesh.is_manifold(viter) &&is_edge_point_[viter->idx()] == false && (tags_[viter->idx()] == CTag::Non || tags_[viter->idx()] == CTag::Base) && mmark[viter->idx()]==false)
		{
			std::queue<COpenMeshT::VertexHandle>Q;
			Q.push(viter);
			if(tags_[viter->idx()]==CTag::Non)
			tags_[viter->idx()] = CTag::Gingiva;
			mmark[viter->idx()] == true;
			while (!Q.empty())
			{
				auto pv=Q.front();
				Q.pop();
				for (auto vviter = mesh.vv_begin(pv); vviter != mesh.vv_end(pv); vviter++)
				{
					if (is_edge_point_[vviter->idx()] == false && mmark[vviter->idx()] == false && (tags_[vviter->idx()] == CTag::Non || tags_[vviter->idx()] == CTag::Base))
					{
						Q.push(vviter);
						if (tags_[vviter->idx()] == CTag::Non)
						tags_[vviter->idx()] = CTag::Gingiva;
						mmark[vviter->idx()] = true;
					}
					
				
				}
			}
		}
		
	}
	for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); viter++)
	{
		if (tags_[viter->idx()] == CTag::Non)
			tags_[viter->idx()] = CTag::Teeth;
	}

}
void CMorphSkelDentalMeshSeg::RemoveInnerTeethRegionInTeeth()
{
	COpenMeshT&mesh = mesh_obj_.GetMesh();
	std::vector<bool>mmark(mesh.n_vertices(), 0);
	for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); viter++)
	{
		if (tags_[viter->idx()] == CTag::Teeth&&mmark[viter->idx()]==false)
		{
			std::queue<COpenMeshT::VertexHandle>Q;
			Q.push(viter);
			mmark[viter->idx()] = true;
			while (!Q.empty())
			{
				auto pv = Q.front();
				Q.pop();
				for (auto vviter = mesh.vv_begin(pv); vviter != mesh.vv_end(pv); vviter++)
				{
					if (mmark[vviter->idx()] == false)
					{
						
					}
				}
			}
		}
	}
}
void CMorphSkelDentalMeshSeg::TestSkeletonize()
{
	COpenMeshT &mesh = mesh_obj_.GetMesh();
	CMorphlogicOperation::Skeletonization(mesh, is_edge_point_);
	for (int i = 0; i < is_edge_point_.size(); i++)
	{
		if (is_edge_point_[i])
			tags_[i] = CTag::Feature;
		else if (tags_[i] != CTag::Base)
			tags_[i] = CTag::Non;
	}
	TestRender();
	
}
void CMorphSkelDentalMeshSeg::TestErode()
{
	COpenMeshT &mesh = mesh_obj_.GetMesh();
	CMorphlogicOperation::Erode(mesh, is_edge_point_);
	for (int i = 0; i < is_edge_point_.size(); i++)
	{
		if (is_edge_point_[i])
			tags_[i] = CTag::Feature;
		else if (tags_[i] != CTag::Base)
			tags_[i] = CTag::Non;
	}
	TestRender();
	
}
void CMorphSkelDentalMeshSeg::TestDilate()
{
	COpenMeshT &mesh = mesh_obj_.GetMesh();
	CMorphlogicOperation::Dilate(mesh, is_edge_point_);
	for (int i = 0; i < is_edge_point_.size(); i++)
	{
		if (is_edge_point_[i])
			tags_[i] = CTag::Feature;
		else if (tags_[i] != CTag::Base)
			tags_[i] = CTag::Non;
	}
	TestRender();

}
void CMorphSkelDentalMeshSeg::ComputeMorphSkeleton()
{
	COpenMeshT &mesh = mesh_obj_.GetMesh();
	CMorphlogicOperation::Dilate(mesh, is_edge_point_);



}
void CMorphSkelDentalMeshSeg::AdjustBaseCuttingPlane(double l)
{
	
	cutting_plane_.SetP(cutting_plane_.p() + cutting_plane_.dir()*l);
	
	std::vector<COpenMeshT::VertexHandle>v_bases;
	MarkBaseByPlane(cutting_plane_, v_bases);

	COpenMeshT&mesh = mesh_obj_.GetMesh();
	if (tags_.size() != mesh.n_vertices())
	{
		tags_.resize(mesh.n_vertices(), CTag::Non);
	}
	for (int i = 0; i < tags_.size(); i++)
	{
		if (tags_[i] == CTag::Base)
			tags_[i] = CTag::Non;
	}
	
	for (int i = 0; i < v_bases.size(); i++)
	{
		tags_[v_bases[i].idx()] = CTag::Base;
	}
	TestRender();
	
}
double CMorphSkelDentalMeshSeg::ComputePenaltyValue(CPlane plane)
{
	double res = 0;
	double u = 18;
	COpenMeshT &mesh = mesh_obj_.GetMesh();
	for (int i = 0; i < edge_points_id_.size(); i++)
	{
		auto vh = edge_points_id_[i];
		int vid = vh.idx();
		OpenMesh::Vec3d point = mesh.point(vh);
		double dis = std::sqrt(CGeoBaseAlg::ComputeSquaredDistance2Plane(point, plane.a(), plane.b(), plane.c(), plane.d()));
		res+=vertex_penalty_weight_[vid] * std::pow(dis - u, 2);
	}
	res /= edge_points_.size();
	return res;
}
void CMorphSkelDentalMeshSeg::ComputeVertexPenaltyWeight()
{
	COpenMeshT &mesh = mesh_obj_.GetMesh();
	vertex_penalty_weight_.resize(mesh.n_vertices(), 0);
	
	
	for(int i=0;i<edge_points_id_.size();i++)
	{
		
		auto vh = edge_points_id_[i];
		int vid = vh.idx();
		std::vector<double>ncurvatures;
		double vcurvature = mean_curvature_values_(vid);
		int ncount = 0;
		double meancurv = 0;
		int nvcount = 0;
		for (auto vviter = mesh.vv_begin(vh); vviter != mesh.vv_end(vh); vviter++)
		{
			int nvid = vviter->idx();
			if (is_edge_point_[nvid])
			{
				ncount++;
				ncurvatures.push_back(mean_curvature_values_(nvid));
				meancurv += mean_curvature_values_(nvid);
			}
			nvcount++;
		}
		
		meancurv/= nvcount;
		double devicurv=CNumericalBaseAlg::ComputeStdDeviation(ncurvatures);
		double gaussian=CNumericalBaseAlg::ComputeGaussian(meancurv, devicurv, vcurvature);
		vertex_penalty_weight_[vid]=gaussian*ncount / edge_points_.size();
		

	}

}
void CMorphSkelDentalMeshSeg::FindOptimizePlane(CPlane ini_plane, CPlane &res_plane)
{
	ComputeVertexPenaltyWeight();

	//find the optimal plane
	double current_penalty= ComputePenaltyValue(ini_plane) ;

	double step = 0.02;
	OpenMesh::Vec3d dir = ini_plane.dir();
	CPlane next_plane(ini_plane.p() + dir*step, dir);
	double next_penalty = ComputePenaltyValue(next_plane);
	if (next_penalty > current_penalty)
	{
		dir = -dir;
	}
	next_plane = CPlane(ini_plane.p() + dir*step, dir);
	next_penalty= ComputePenaltyValue(next_plane);
	while (next_penalty < current_penalty)
	{
		next_plane= CPlane(next_plane.p() + dir*step, dir);
		next_penalty = ComputePenaltyValue(next_plane);
	}
	next_plane.SetDir(-dir);
	res_plane = next_plane;
	
}
void CMorphSkelDentalMeshSeg::RemoveSmallFeatureRegions()
{
	COpenMeshT &mesh = mesh_obj_.GetMesh();
	std::vector<std::vector<COpenMeshT::VertexHandle>>feature_groups;
	CGeoAlg::GetFeatureGroupsByConnectivity(mesh_obj_.GetMesh(), is_edge_point_, feature_groups);
	for (int i = 0; i < feature_groups.size(); i++)
	{
		if (feature_groups[i].size() < small_region_threshold_)
		{
			for (int j = 0; j < feature_groups[i].size(); j++)
			{
				is_edge_point_[feature_groups[i][j].idx()] = false;
			}
		}
		else
		{
			edge_points_id_.clear();
			edge_points_.clear();
			
			for (int j = 0; j < feature_groups[i].size(); j++)
			{
				edge_points_id_.push_back(feature_groups[i][j]);
				edge_points_.push_back(mesh.point(feature_groups[i][j]));
			}
		}
	}
	for (auto viter=mesh.vertices_begin();viter!=mesh.vertices_end();viter++)
	{
		if (is_edge_point_[viter->idx()])
			tags_[viter->idx()] = CTag::Feature;
		else if (tags_[viter->idx()] == CTag::Feature)
		{
			auto vviter = mesh.vv_begin(viter);
			for (; vviter != mesh.vv_end(viter); vviter++)
			{
				if (tags_[vviter->idx()] != CTag::Feature)
					break;
			}
			if (vviter != mesh.vv_end(viter))
				tags_[viter->idx()] = tags_[vviter->idx()];
			else
				tags_[viter->idx()] = CTag::Non;
		}
			
	}
}
void CMorphSkelDentalMeshSeg::ComputeEdgePointsFromMeanCurvatureThreshold(double thre)
{
	Eigen::MatrixXd vertexs;
	Eigen::MatrixXi faces;
	COpenMeshT& mesh = mesh_obj_.GetMesh();
	CConverter::ConvertFromOpenMeshToIGL(mesh, vertexs, faces);
	CGeoBaseAlg::ComputeMeanCurvatureVector(vertexs, faces, mean_curvature_vec_);
	mean_curvature_values_.resize(mean_curvature_vec_.rows());
	if(tags_.size()!=mesh.n_vertices())
	tags_.resize(mesh.n_vertices(), CTag::Non);
	for (int i = 0; i < mean_curvature_vec_.rows(); i++)
	{
		mean_curvature_values_(i) = mean_curvature_vec_.row(i).norm();
	}

	Eigen::MatrixXd N;
	CGeoBaseAlg::ComputePerVertexNormal(vertexs, faces, N);
	edge_points_.clear();
	edge_points_id_.clear();
	is_edge_point_.resize(vertexs.rows(), 0);
	for (int i = 0; i < vertexs.rows(); i++)
	{
		is_edge_point_[i] = false;
		
	}
	for (auto v_iter = mesh.vertices_begin(); v_iter != mesh.vertices_end(); v_iter++)
	{
		int vid = v_iter->idx();
		if(tags_[vid]==CTag::Base)
			continue;
		if (mean_curvature_values_(vid) > thre&& mean_curvature_vec_.row(vid).dot(N.row(vid)) < 0.0&&(!mesh.is_boundary(v_iter)))
		{
			edge_points_id_.push_back(v_iter);
			edge_points_.push_back(mesh.point(v_iter));
			is_edge_point_[vid] = true;
			tags_[vid] = CTag::Feature;
		}
		
	}

}
void CMorphSkelDentalMeshSeg::AdjustSmallRegionThreshold(double thre)
{
	small_teeth_region_percent_threshold_ = thre;
}
void CMorphSkelDentalMeshSeg::RemoveSmallIsolateTeethRegion()
{
	COpenMeshT &mesh = mesh_obj_.GetMesh();
	if (tags_.size() != mesh.n_vertices())
		return;
	double totarea = mesh_obj_.GetTotSurfaceArea();
	std::vector<CTag>tmp_tags = tags_;
	for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); viter++)
	{
		int vid = viter->idx();
	
		if ((tags_[vid] == CTag::Teeth)&&(tmp_tags[vid] == CTag::Teeth ))
		{
			CTag t = CTag::Feature;
			if (t == tags_[vid])
				t = CTag::Teeth;
			double area=ComputeRegionArea(viter, &t);
			double percent = area / totarea;
			//std::cerr << "perc " << percent << std::endl;
			if (percent <= small_teeth_region_percent_threshold_)
			{
				CTag t = CTag::Feature;
				//std::cerr << "tag region" << std::endl;
				TagRegion(tags_,viter, CTag::Gingiva, &t);
				
			}
			else
			{
				CTag t = CTag::Feature;
				TagRegion(tmp_tags, viter, CTag::Gingiva, &t);
			}
		}
	
	}

}
void CMorphSkelDentalMeshSeg::RemoveGingiva2GingivaFeatureEdge()
{
	
	COpenMeshT &mesh = mesh_obj_.GetMesh();
	if (tags_.size() != mesh.n_vertices())
		return;
	std::vector<COpenMeshT::VertexHandle>tmp_edge_pid = edge_points_id_;
	edge_points_id_.clear();
	edge_points_.clear();
	for (int i = 0; i < tmp_edge_pid.size(); i++)
	{
		auto vh = tmp_edge_pid[i];
		if (!CMorphlogicOperation::IsCenterVertex(mesh, is_edge_point_, vh))
	//	if(true)
		{
			
			auto vviter = mesh.vv_begin(vh);
			for (; vviter != mesh.vv_end(vh); vviter++)
			{
				if (tags_[vviter->idx()] != CTag::Feature&&tags_[vviter->idx()] != CTag::Gingiva)
					break;
			}
			if (vviter == mesh.vv_end(vh))//Gingiva 2 Gingiva Feature Edge
			{
				//std::cerr << "sdf" << std::endl;
				is_edge_point_[tmp_edge_pid[i].idx()] = false;
				tags_[tmp_edge_pid[i].idx()] = CTag::Gingiva;
				//std::cerr << "gingiva" << std::endl;
			}
			else
			{
				edge_points_id_.push_back(tmp_edge_pid[i]);
				edge_points_.push_back(mesh.point(tmp_edge_pid[i]));
			}
		}
	
		else
		{
			edge_points_id_.push_back(tmp_edge_pid[i]);
			edge_points_.push_back(mesh.point(tmp_edge_pid[i]));
		}
	}
	 
}
void CMorphSkelDentalMeshSeg::SwitchGingivaAndTeeth(COpenMeshT::VertexHandle vh)
{
	COpenMeshT &mesh = mesh_obj_.GetMesh();
	if (tags_.size() != mesh.n_vertices()||tags_[vh.idx()]==CTag::Feature)
	{
		return;
	}

	if (tags_[vh.idx()] == CTag::Teeth)
		TagRegion(tags_,vh, CTag::Gingiva);
	else if (tags_[vh.idx()] == CTag::Gingiva)
		TagRegion(tags_,vh, CTag::Teeth);
	else
		return;
	
	RemoveGingiva2GingivaFeatureEdge();
	
	TestRender();

}
double CMorphSkelDentalMeshSeg::ComputeRegionArea(COpenMeshT::VertexHandle vh, CTag *glue_tag_)
{
	//std::cerr << "region area" << std::endl;
	COpenMeshT &mesh = mesh_obj_.GetMesh();
	if (tags_.size() != mesh.n_vertices())
		return -1;
	std::queue<COpenMeshT::VertexHandle>Q;
	Q.push(vh);
	auto origtag = tags_[vh.idx()];
	CTag gluetag = glue_tag_ != NULL ? *glue_tag_ : origtag;
	std::vector<bool>mmark(mesh.n_vertices(), false);
	mmark[vh.idx()] = true;
	double area = CGeoBaseAlg::ComputeVertexArea(mesh,vh);
	while (!Q.empty())
	{
		auto p = Q.front();
		Q.pop();
		for (auto vviter = mesh.vv_begin(p); vviter != mesh.vv_end(p); vviter++)
		{
			if ((gluetag==tags_[vviter->idx()]||origtag == tags_[vviter->idx()])&&mmark[vviter->idx()]==false)
			{
				mmark[vviter->idx()] = true;
				area += CGeoBaseAlg::ComputeVertexArea(mesh, vviter);
				Q.push(vviter);
			}
		}
	}
	return area;
}
void CMorphSkelDentalMeshSeg::ResetEdgeVertexMarkFromTags()
{
	COpenMeshT &mesh = mesh_obj_.GetMesh();
	if (tags_.size() == mesh.n_vertices())
	{
		is_edge_point_.resize(mesh.n_vertices(), false);
		edge_points_.clear();
		edge_points_id_.clear();
		for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); viter++)
		{
			int vid = viter->idx();
			if (tags_[vid] == CTag::Feature)
			{
				edge_points_.push_back(mesh.point(viter));
				edge_points_id_.push_back(viter);
				is_edge_point_[vid] = true;
			}
			else
			{
				is_edge_point_[vid] = false;
			}
		}
		
	}
}
void CMorphSkelDentalMeshSeg::TagRegion(std::vector<CTag>&tags, COpenMeshT::VertexHandle vh, CTag tag, CTag *glue_tag_)
{
	
	COpenMeshT &mesh = mesh_obj_.GetMesh();
	if (tags.size() != mesh.n_vertices())
		return;
	if (tags[vh.idx()] == tag)
		return;
	
	std::queue<COpenMeshT::VertexHandle>Q;
	Q.push(vh);
	auto origtag = tags[vh.idx()];
	CTag gludtag = glue_tag_ != NULL ? *glue_tag_ : origtag;
	tags[vh.idx()] = tag;
	while (!Q.empty())
	{
		auto p = Q.front();
		Q.pop();
		for (auto vviter = mesh.vv_begin(p); vviter != mesh.vv_end(p); vviter++)
		{
			if (origtag == tags[vviter->idx()]|| tags[vviter->idx()]== gludtag)
			{
				tags[vviter->idx()] = tag;
				Q.push(vviter);
			}
		}
	}
	ResetEdgeVertexMarkFromTags();
}
void CMorphSkelDentalMeshSeg::MarkBaseByPlane(CPlane plane, std::vector<COpenMeshT::VertexHandle>&v_bases)
{
	COpenMeshT&mesh = mesh_obj_.GetMesh();
	OpenMesh::Vec3d hole_mean(0, 0, 0);
	int count = 0;
	for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); viter++)
	{
		if (mesh.is_boundary(viter))
		{
			hole_mean += mesh.point(viter);
			count++;
		}
	}
	hole_mean /= count;
	OpenMesh::Vec3d dir = plane.p() - hole_mean;
	auto tt = dir*plane.dir();

	if ((OpenMesh::dot(dir,plane.dir())) < 0.0)
	{
		plane.SetDir(-plane.dir());
	}
	
	for (auto v_iter = mesh.vertices_begin(); v_iter != mesh.vertices_end(); v_iter++)
	{
		COpenMeshT::Point p = mesh.point(v_iter);
		if (!CGeoBaseAlg::IsOnPositiveSide(p, plane))
			v_bases.push_back(v_iter);
	}
}
void CMorphSkelDentalMeshSeg::ComputeCuttingPlane()
{
	
	

	

	ComputeEdgePointsFromMeanCurvatureThreshold(10);
	if(verbose_)
	std::cout << "num of edge points " << edge_points_.size() << std::endl;
	std::vector<OpenMesh::Vec3d>  eigen_vecs;
	std::vector<double>eigen_values;
	OpenMesh::Vec3d  mean;
	CGeoAlg::PointSetPCA3D(edge_points_, mean, eigen_vecs, eigen_values);
	OpenMesh::Vec3d dir = eigen_vecs[2];
	CPlane cplane(mean, dir);
	
	FindOptimizePlane(cplane, cutting_plane_);
	std::vector<COpenMeshT::VertexHandle>v_bases;
	MarkBaseByPlane(cutting_plane_, v_bases);
	
	COpenMeshT&mesh = mesh_obj_.GetMesh();
	if (tags_.size() != mesh.n_vertices())
	{
		tags_.resize(mesh.n_vertices(),CTag::Non);
	}
	for (int i = 0; i < v_bases.size(); i++)
	{
		tags_[v_bases[i].idx()] = CTag::Base;
	}
	//CDataIO::WriteObj("test0.obj", mesh_obj_);
	//std::map<COpenMeshT::VertexHandle, COpenMeshT::VertexHandle>vid_orig;

	//for debug, render the cutting plane
	/*CMeshObject *plane_obj = new CMeshObject();
	CAuxGeoUtils::GetPlaneMeshFromPointAndAxis(cutting_plane_.p(), eigen_vecs[0], eigen_vecs[1], eigen_vecs[2], 1.2, plane_obj->GetMesh());
	for (auto viter = plane_obj->GetMesh().vertices_begin(); viter != plane_obj->GetMesh().vertices_end(); viter++)
	{
		plane_obj->GetMesh().set_color(viter, OpenMesh::Vec3d(0, 1, 0));
	}
	plane_obj->SetChanged();
	DataPool::AddMeshObject(plane_obj);
	return;*/



	/*
	CGeoAlg::CutByPlane(mesh_obj_.GetMesh(), cutting_plane_, mesh_obj_.GetMesh(), vid_orig);
	std::vector<bool>tmp_is_edge_point = is_edge_point_;
	is_edge_point_.clear();
	
	COpenMeshT& mesh = mesh_obj_.GetMesh();
	is_edge_point_.resize(mesh.n_vertices(), 0);
	for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); viter++)
	{
		int vid = vid_orig[viter].idx();
		if (tmp_is_edge_point[vid])
		{
			is_edge_point_[viter->idx()] = true;
		}
		else
			is_edge_point_[viter->idx()] = false;
	}
	
	//CDataIO::WriteObj("test.obj", mesh_obj_);
	std::vector<COpenMeshT*>res_meshes;
	std::vector<std::map<COpenMeshT::VertexHandle, COpenMeshT::VertexHandle>>vid_origs;
	CGeoAlg::SeparateDisconnectedParts(mesh, res_meshes, vid_origs);
	int maxnum = -1;
	int mi=0;
	for (int i = 0; i < res_meshes.size(); i++)
	{
		if (maxnum < res_meshes[i]->n_vertices())
		{
		
			maxnum = res_meshes[i]->n_vertices();
			mi = i;
		}
	}

	mesh_obj_.GetMesh()=(*res_meshes[mi]);
	for (int i = 0; i < res_meshes.size(); i++)
	{
		if (i != mi)
			delete res_meshes[i];
	}
	
	COpenMeshT & tmesh = mesh_obj_.GetMesh();
	tmp_is_edge_point = is_edge_point_;
	is_edge_point_.resize(tmesh.n_vertices());
	for (auto viter = tmesh.vertices_begin(); viter != tmesh.vertices_end(); viter++)
	{
		int vid = viter->idx();
		if (tmp_is_edge_point[vid_origs[mi][viter].idx()])
		{
			tmesh.set_color(viter, OpenMesh::Vec3d(1, 0, 0));
			is_edge_point_[vid]=true;
		}
		else
		{
			tmesh.set_color(viter, OpenMesh::Vec3d(0.8, 0.8, 0.8));
			is_edge_point_[vid]=false;
		}
		

	}
	

	*/


	ComputeEdgePointsFromMeanCurvatureThreshold(curvature_threshold_);
	TestRender();

	
}