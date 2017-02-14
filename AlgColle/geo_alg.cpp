#include"geo_alg.h"

#include"../DataColle/cgal_igl_converter.h"
#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <opencv2/opencv.hpp>
#include"geo_base_alg.h"
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_vertex_base_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>

#include <CGAL/Delaunay_mesh_size_criteria_2.h>
#include<queue>
#include"../AdditionalLibs/obb/cpqp_obb_wrapper.h"
#include"../DataColle/obb_type.h"
#include "../AdditionalLibs/Remesher/UniformRemesherT.hh"
#include "../AdditionalLibs/Remesher/AdaptiveRemesherT.hh"
#include "../AdditionalLibs/Remesher/MeshSelectionT.hh"
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <igl/embree/EmbreeIntersector.h>
#include <CGAL/Aff_transformation_3.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include"../DataColle/aabb_type.h"
#include"../AdditionalLibs/XWGeodesic/xw_geodesic_wrapper.h"
#include <igl/grad.h>
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_length_cost.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Midpoint_placement.h>
#include <CGAL/boost/graph/graph_traits_TriMesh_ArrayKernelT.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include"../DataColle/curve_object.h"
#include"../DataColle/data_pool.h"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/Delaunay_triangulation_2.h>
void CGeoAlg::PointSetPCA3D(std::vector<OpenMesh::Vec3d> &pts, OpenMesh::Vec3d&res_mean, std::vector<OpenMesh::Vec3d> &res_eigen_vects, std::vector<double>& res_eigen_values)
{
	int sz = static_cast<int>(pts.size());
	cv::Mat data_pts =cv::Mat(sz, 3, CV_64FC1);
	for (int i = 0; i < data_pts.rows; ++i)
	{
		data_pts.at<double>(i, 0) = pts[i][0];
		data_pts.at<double>(i, 1) = pts[i][1];
		data_pts.at<double>(i, 2) = pts[i][2];
	}
	//Perform PCA analysis
	cv::PCA pca_analysis(data_pts, cv::Mat(), CV_PCA_DATA_AS_ROW);
	//Store the center of the object
	res_mean= OpenMesh::Vec3d(pca_analysis.mean.at<double>(0, 0),pca_analysis.mean.at<double>(0, 1), pca_analysis.mean.at<double>(0, 2));
	int dim = 3;
	res_eigen_vects.resize(dim);
	res_eigen_values.resize(dim);

	for (int i = 0; i < dim; i++)
	{
		
		for (int j = 0; j < dim; j++)
		{
			res_eigen_vects[i][j] = pca_analysis.eigenvectors.at<double>(i, j);
		}
		res_eigen_values[i]= pca_analysis.eigenvalues.at<double>(i, 0);
		//std::cerr << "eig value " << res_eigen_values[i] << std::endl;
	}

}
void CGeoAlg::ComputeClosestPoint(OpenMesh::Vec3d p, std::pair<OpenMesh::Vec3d, OpenMesh::Vec3d>seg, OpenMesh::Vec3d &res_p, double &res_dis)
{
	double x0 = seg.first[0];
	double y0 = seg.first[1];
	double z0 = seg.first[2];

	double x1 = seg.second[0];
	double y1 = seg.second[1];
	double z1 = seg.second[2];

	double x2 = p[0];
	double y2 = p[1];
	double z2 = p[2];

	double t = ((x0 - x1)*(x0 - x2) + (y0 - y1)*(y0 - y2) + (z0 - z1)*(z0 - z2))/(std::pow(x1 - x0, 2) + std::pow(y1 - y0, 2) + std::pow(z1 - z0, 2)) ;
	if (t >= 0 && t <= 1)
	{
		res_p = seg.first + (seg.second - seg.first)*t;
	}
	else if (t < 0)
	{
		res_p = seg.first;
	}
	else
	{
		res_p = seg.second;
	}
	OpenMesh::Vec3d tmp_dir = res_p - p;
	res_dis = std::sqrt(tmp_dir.sqrnorm());

}
bool CGeoAlg::ComputeClosestPointOfTwoSegments(std::pair<OpenMesh::Vec3d, OpenMesh::Vec3d>sega, std::pair<OpenMesh::Vec3d, OpenMesh::Vec3d>segb, OpenMesh::Vec3d&res_pa, OpenMesh::Vec3d &res_pb)
{
	res_pa = OpenMesh::Vec3d(0, 0, 0);
	res_pb = OpenMesh::Vec3d(0, 0, 0); 
	double x1 = sega.first[0];
	double y1 = sega.first[1];
	double z1 = sega.first[2];

	double x2 = sega.second[0];
	double y2 = sega.second[1];
	double z2 = sega.second[2];

	double x3 = segb.first[0];
	double y3 = segb.first[1];
	double z3 = segb.first[2];

	double x4 = segb.second[0];
	double y4 = segb.second[1];
	double z4 = segb.second[2];
	//ps1*s+pt1*t=pc1;
	double ps1 = (std::pow(x2 - x1, 2)+std::pow(y2-y1,2)+std::pow(z2-z1,2));
	double pt1 = -((x2-x1)*(x4-x3)+(y2-y1)*(y4-y3)+(z2-z1)*(z4-z3));
	double pc1 = (x1 - x2)*(x1 - x3) + (y1 - y2)*(y1 - y3) + (z1 - z2)*(z1 - z3);

	double ps2 = -((x2 - x1)*(x4 - x3) + (y2 - y1)*(y4 - y3) + (z2 - z1)*(z4 - z3));
	double pt2 = std::pow(x4 - x3, 2) + std::pow(y4 - y3, 2) + std::pow(z4 - z3, 2);
	double pc2 = (x1 - x3)*(x4 - x3) + (y1 - y3)*(y4 - y3) + (z1 - z3)*(z4 - z3);
	if (pt1*ps2 - ps1*pt2 == 0 || ps1*pt2 - ps2*pt1 == 0)
	{
		std::cerr << "pt1*ps2:"<< pt1*ps2  <<"  ps1*pt2: " <<  ps1*pt2 << std::endl;
		std::cerr << "ps1*pt2:"<<ps1*pt2<< " ps2*pt1 "  <<ps2*pt1 << std::endl;
		return false;
	}
		
	double t = (pc1*ps2 - ps1*pc2) / (pt1*ps2-ps1*pt2);
	double s = (pc1*pt2 - pc2*pt1) / (ps1*pt2-ps2*pt1);
	std::cerr << "s " << s << " t " << t << std::endl;
	if (s >= 0 && s <= 1 && t >= 0 && t <= 1)
	{
		res_pa = sega.first + (sega.second - sega.first)*s;
		res_pb = segb.first + (segb.second - segb.first)*t;
	
	}
	else
	{
		double dis_sega0,dis_sega1,dis_segb0,dis_segb1;
		OpenMesh::Vec3d p_sega0,p_sega1,p_segb0,p_segb1;
		ComputeClosestPoint(sega.first, segb, p_sega0, dis_sega0);
		ComputeClosestPoint(sega.second, segb, p_sega1, dis_sega1);

		ComputeClosestPoint(segb.first, sega, p_segb0, dis_segb0);
		ComputeClosestPoint(segb.second, sega, p_segb1, dis_segb1);

		if (dis_sega0 < dis_sega1&&dis_sega0 < dis_segb0&&dis_sega0 < dis_segb1)
		{
			res_pa = sega.first;
			res_pb = p_sega0;
			//std::cerr << "p_sega0" << std::endl;
		}
		else if (dis_sega1 < dis_segb0&&dis_sega1 < dis_segb1)
		{
			res_pa = sega.second;
			res_pb = p_sega1;
			//std::cerr << "p_sega1" << std::endl;
		}
		else if (dis_segb0 < dis_segb1)
		{
			res_pa = p_segb0;
			res_pb = segb.first;
			//std::cerr << "p_segb0" << std::endl;
		}
		else
		{
			res_pa = p_segb1;
			res_pb = segb.second;
			//std::cerr << "p_segb1" << std::endl;
		}
	}
	return true;
}
void CGeoAlg::FindEdgePointsPath(CMeshObject &mesh_obj, std::pair<COpenMeshT::FaceHandle, OpenMesh::Vec3d>fpa, std::pair<COpenMeshT::FaceHandle, OpenMesh::Vec3d>fpb, std::vector<std::pair<COpenMeshT::EdgeHandle, OpenMesh::Vec3d>>&res_eps, std::vector<bool>&is_vertex, int recursive_num)
{
	res_eps.clear();
	is_vertex.clear();
	
	COpenMeshT &mesh = mesh_obj.GetMesh();
	if (fpa.second == fpb.second)
	{
		std::vector<COpenMeshT::VertexHandle>comm_vhs;
		if (CGeoBaseAlg::IsShareCommonVertex(mesh, fpa.first, fpb.first, comm_vhs))
		{
			if (comm_vhs.size() == 1)
			{
				OpenMesh::Vec3d p = mesh_obj.TransformPointByLocalMatrix(mesh.point(comm_vhs[0]));
				res_eps.push_back(std::make_pair(mesh.ve_begin(comm_vhs[0]),p));
				is_vertex.push_back(true);
			}
			else if (comm_vhs.size() == 2)
			{
				OpenMesh::Vec3d p0 = mesh_obj.TransformPointByLocalMatrix(mesh.point(comm_vhs[0]));
				OpenMesh::Vec3d p1 = mesh_obj.TransformPointByLocalMatrix(mesh.point(comm_vhs[1]));
				if (CGeoBaseAlg::ComputeDis(p0, fpa.second) < CGeoBaseAlg::ComputeDis(p1, fpa.second))
				{
					res_eps.push_back(std::make_pair(mesh.ve_begin(comm_vhs[0]), p0));
				}
				else
				{
					res_eps.push_back(std::make_pair(mesh.ve_begin(comm_vhs[1]), p1));
				}
				is_vertex.push_back(true);
			}
			else if (comm_vhs.size() == 3)
			{
				OpenMesh::Vec3d p0 = mesh_obj.TransformPointByLocalMatrix(mesh.point(comm_vhs[0]));
				OpenMesh::Vec3d p1 = mesh_obj.TransformPointByLocalMatrix(mesh.point(comm_vhs[1]));
				OpenMesh::Vec3d p2 = mesh_obj.TransformPointByLocalMatrix(mesh.point(comm_vhs[2]));

				double dis0 = CGeoBaseAlg::ComputeDis(p0, fpa.second);
				double dis1 = CGeoBaseAlg::ComputeDis(p1, fpa.second);
				double dis2 = CGeoBaseAlg::ComputeDis(p2, fpa.second);
				if (dis0 == 0)
				{
					std::cerr << "!!!!!!!!!!!!!!!!!!!!!!!fpa and fpb are on the save vertex: "<<p0 << std::endl;
					
				}
				else if (dis1 == 0)
				{
					std::cerr << "!!!!!!!!!!!!!!!!!!!!!!!fpa and fpb are on the save vertex: " << p1 << std::endl;
					
				}
				else if (dis2 == 0)
				{
					std::cerr << "!!!!!!!!!!!!!!!!!!!!!!!fpa and fpb are on the save vertex: " << p2 << std::endl;
				}

			}
		}
		return;
	}

	COpenMeshT::EdgeHandle comm_edge;
	std::vector<COpenMeshT::VertexHandle>comm_vhs;
	bool flag = true;
	if (fpa.first == fpb.first)
	{
	//	std::cerr << "on the same triangle" << std::endl;
		flag = false;
		
	}
	
	else if (CGeoBaseAlg::IsShareCommonEdge(mesh,fpa.first, fpb.first, comm_edge))
	{
		//std::cerr << "share common edge" << std::endl;
		COpenMeshT::HalfedgeHandle he0 = mesh.halfedge_handle(comm_edge, 0), he1 = mesh.halfedge_handle(comm_edge, 1);
		OpenMesh::Vec3d p0 = mesh_obj.TransformPointByLocalMatrix(mesh.point(mesh.to_vertex_handle(he0)));
		OpenMesh::Vec3d p1 = mesh_obj.TransformPointByLocalMatrix(mesh.point(mesh.to_vertex_handle(he1)));
	//	std::cerr << "p0 " << p0 << " p1 " << p1 << std::endl;
	//	std::cerr << "fpa " << fpa.second << " fpb " << fpb.second << std::endl;
		OpenMesh::Vec3d res_pa, res_pb;
	//	std::cerr << "ComputeClosestPointOfTwoSegments 0" << std::endl;

		ComputeClosestPointOfTwoSegments(std::make_pair(fpa.second, fpb.second), std::make_pair(p0, p1), res_pa, res_pb);
		//std::cerr << "res pb " << res_pb << std::endl;
		//std::cerr << "ComputeClosestPointOfTwoSegments 1" << std::endl;
	
		if (CGeoBaseAlg::ComputeDis(p0, res_pb) ==0)
		{
			//res_pb = p0;
		/*	std::cerr << "is vertex ***********" << std::endl;
			std::cerr << "p0 " << p0 << " p1 " << p1 << std::endl;
			std::cerr << "fpa " << fpa.second << " fpb " << fpb.second << std::endl;*/
			//is_vertex.push_back(true);
		} 
		else if (CGeoBaseAlg::ComputeDis(p1, res_pb) ==0)
		{
			//res_pb = p1;
		/*	std::cerr << "is vertex ***********" << std::endl;
			std::cerr << "p0 " << p0 << " p1 " << p1 << std::endl;
			std::cerr << "fpa " << fpa.second << " fpb " << fpb.second << std::endl;*/
			//is_vertex.push_back(true);
		}
		else
		{
			flag = false;
			is_vertex.push_back(false);
			res_eps.push_back(std::make_pair(comm_edge, res_pb));
		}
		
	
	}
	

	if(flag)
	{
	
		if (recursive_num == 0)
		{
		//	std::cerr << "share nothing" << std::endl;
			
			if (CGeoBaseAlg::IsShareCommonVertex(mesh, fpa.first, fpb.first, comm_vhs))
			{
			//std::cerr << "share common vertex*********************** vh size "<<comm_vhs.size()<<" vid " << comm_vhs[0] .idx()<< std::endl;
				OpenMesh::Vec3d p = mesh_obj.TransformPointByLocalMatrix(mesh.point(comm_vhs[0]));
				res_eps.push_back(std::make_pair(mesh.ve_begin(comm_vhs[0]), p));
				is_vertex.push_back(true);
			//	std::cerr << "fpa face " << fpa.first << " fpb face " << fpb.first << std::endl;
			//	std::cerr << fpa.second << " " << fpb.second << std::endl;
			}
			else
			{
				std::cerr << "error recursive_num==0!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
			/*	CCurveObject *co = new CCurveObject();
				co->GetCurve().push_back(fpa.second);
				co->GetCurve().push_back(fpb.second);
				co->GetCurve().push_back(OpenMesh::Vec3d(0, 0, 0));
				std::cerr << fpa.second << std::endl;
				std::cerr << fpb.second << std::endl;
				co->SetColor(OpenMesh::Vec3d(1, 0, 0));

				DataPool::AddCurveObject(co);*/
			}
			return;
		}
		OpenMesh::Vec3d centerp = (fpa.second + fpb.second) / 2;
		COpenMeshT::FaceHandle cfh;
		OpenMesh::Vec3d cp;
		CGeoAlg::ComputeClosestPoint(centerp, mesh_obj, cfh, cp);
		std::vector<std::pair<COpenMeshT::EdgeHandle, OpenMesh::Vec3d>>tmp_eps;
		std::vector<bool>tmp_is_v;
		CGeoAlg::FindEdgePointsPath(mesh_obj, fpa, std::make_pair(cfh, cp), tmp_eps, tmp_is_v, recursive_num-1);
		for (int i = 0; i < tmp_eps.size(); i++)
		{
			res_eps.push_back(tmp_eps[i]);
			is_vertex.push_back(tmp_is_v[i]);
		}
		CGeoAlg::FindEdgePointsPath(mesh_obj, std::make_pair(cfh, cp),fpb, tmp_eps, tmp_is_v, recursive_num-1);
		for (int i = 0; i < tmp_eps.size(); i++)
		{
			res_eps.push_back(tmp_eps[i]);
			is_vertex.push_back(tmp_is_v[i]);
		}
	}
}
void CGeoAlg::FindPointsPath(CMeshObject &mesh_obj, std::vector<std::pair<COpenMeshT::FaceHandle, OpenMesh::Vec3d>>face_points_path, std::vector<OpenMesh::Vec3d>&res_points, std::map<int,COpenMeshT::FaceHandle>&res_fhs_map, std::map<int,COpenMeshT::EdgeHandle>&res_ehs_map, std::map<int, COpenMeshT::VertexHandle>&res_vhs_map, std::vector<int>&orig_res_map,bool is_closed)
{
	res_points.clear();
	res_ehs_map.clear();
	res_fhs_map.clear();
	int recursive_num = 100;
	COpenMeshT & mesh = mesh_obj.GetMesh();
	orig_res_map.resize(face_points_path.size());
	for (int i = 0; i < face_points_path.size() - 1; i++)
	{
		OpenMesh::Vec3d bary_coord = CGeoBaseAlg::ComputeBaryCoordInTriFace(mesh, face_points_path[i].first, face_points_path[i].second);
		double sum_bary = 0;
		for (int k = 0; k < 3; k++)
		{
			if (bary_coord[k] < 1e-3)
			{
				bary_coord[k] = 1e-3;
			}
			sum_bary += bary_coord[k];
		}
		for (int k = 0; k < 3; k++)
		{
			bary_coord[k] /= sum_bary;
		}
		face_points_path[i].second = CGeoBaseAlg::ComputeVPosFromBaryCoord(mesh, face_points_path[i].first, bary_coord);
	}

	for (int i = 0; i < face_points_path.size()-1; i++)
	{
		int j = i + 1;
		
		
		//std::cerr << "res points size " << res_points.size() << std::endl;
		/*if (res_points.size() == 83 || res_points.size() == 168)
		{
			std::cerr << "p " << res_points.back() << std::endl;
		}*/
		res_points.push_back(face_points_path[i].second);
		orig_res_map[i] = res_points.size() - 1;
		res_fhs_map[res_points.size() - 1] = face_points_path[i].first;
		std::vector<std::pair<COpenMeshT::EdgeHandle, OpenMesh::Vec3d>>res_eps;
		std::vector<bool>is_vertex;
	//	std::cerr << "start find edge point" << std::endl;
		//std::cerr << face_points_path[i].second << " " << face_points_path[j].second << std::endl;
		/*if (face_points_path[i].first == face_points_path[j].first)
		{
			std::cerr << "fhi " << face_points_path[i].first.idx() << " fhj " << face_points_path[j].first.idx() << std::endl;
		}*/
		FindEdgePointsPath(mesh_obj, face_points_path[i], face_points_path[j], res_eps, is_vertex, recursive_num);
		//std::cerr << "res_eps size " << res_eps.size() << std::endl;
		for (int t = 0; t < res_eps.size(); t++)
		{
			//std::cerr << "res points size " << res_points.size() << std::endl;
			if (res_points.size() == 83 || res_points.size() == 168)
			{
				std::cerr << "p " << res_points.back() << std::endl;
			}
			res_points.push_back(res_eps[t].second);
			//std::cerr << "res point " << res_eps[t].second << std::endl;
			if (is_vertex[t])
			{
				COpenMeshT::VertexHandle vh0 = mesh.to_vertex_handle(mesh.halfedge_handle(res_eps[t].first, 0));
				COpenMeshT::VertexHandle vh1 = mesh.to_vertex_handle(mesh.halfedge_handle(res_eps[t].first, 1));
				OpenMesh::Vec3d p0 = mesh.point(vh0);
				OpenMesh::Vec3d p1 = mesh.point(vh1);
				OpenMesh::Vec3d dif0 = res_eps[t].second - p0;
				OpenMesh::Vec3d dif1 = res_eps[t].second - p1;
				if (dif0.sqrnorm() > dif1.sqrnorm())
				{
					res_vhs_map[res_points.size() - 1] = vh1;
				}
				else
				{
					res_vhs_map[res_points.size() - 1] = vh0;
				}
			
			}
			else
			{
				res_ehs_map[res_points.size() - 1] = res_eps[t].first;
			}
			
			
		}
		//std::cerr << "end find edge point" << std::endl;

	}



	//std::cerr << "res points size " << res_points.size() << std::endl;
	/*if (res_points.size() == 83 || res_points.size() == 168)
	{
		std::cerr << "p " << res_points.back() << std::endl;
	}*/
	res_points.push_back(face_points_path.back().second);
	orig_res_map[face_points_path.size()-1] = res_points.size() - 1;
	res_fhs_map[res_points.size() - 1] = face_points_path.back().first;



	if (is_closed)
	{
		int i = face_points_path.size() - 1;
		int j = 0;
		std::vector<std::pair<COpenMeshT::EdgeHandle, OpenMesh::Vec3d>>res_eps;
		std::vector<bool>is_vertex;
		FindEdgePointsPath(mesh_obj, face_points_path[i], face_points_path[j], res_eps, is_vertex, recursive_num);
		for (int t = 0; t < res_eps.size(); t++)
		{
			//std::cerr << "res points size " << res_points.size() << std::endl;
			/*if (res_points.size() == 83 || res_points.size() == 168)
			{
				std::cerr << "p " << res_points.back() << std::endl;
			}*/
			res_points.push_back(res_eps[t].second);

			if (is_vertex[t])
			{
				COpenMeshT::VertexHandle vh0 = mesh.to_vertex_handle(mesh.halfedge_handle(res_eps[t].first, 0));
				COpenMeshT::VertexHandle vh1 = mesh.to_vertex_handle(mesh.halfedge_handle(res_eps[t].first, 1));
				OpenMesh::Vec3d p0 = mesh.point(vh0);
				OpenMesh::Vec3d p1 = mesh.point(vh1);
				OpenMesh::Vec3d dif0 = res_eps[t].second - p0;
				OpenMesh::Vec3d dif1 = res_eps[t].second - p1;
				if (dif0.sqrnorm() > dif1.sqrnorm())
				{
					res_vhs_map[res_points.size() - 1] = vh1;
				}
				else
				{
					res_vhs_map[res_points.size() - 1] = vh0;
				}

			}
			else
			{
				res_ehs_map[res_points.size() - 1] = res_eps[t].first;
			}
			
		}
	}
}
void CGeoAlg::SeparateMeshByVhsLoop(COpenMeshT &mesh, std::vector<COpenMeshT::VertexHandle>&vhs, std::vector<COpenMeshT*>&res_meshes, std::vector<std::map<COpenMeshT::VertexHandle, COpenMeshT::VertexHandle>>&new2orig_vhsmap)
{
	
	std::set<std::pair<COpenMeshT::VertexHandle, COpenMeshT::VertexHandle>>loop_edge_set;
	std::vector<bool>vmark;
	vmark.resize(mesh.n_vertices(), 0);
	for (int i = 0; i < vhs.size(); i++)
	{
		vmark[vhs[i].idx()] = true;
		loop_edge_set.insert(std::make_pair(vhs[i],vhs[(i+1)%vhs.size()]));
		loop_edge_set.insert(std::make_pair( vhs[(i + 1) % vhs.size()], vhs[i]));
	}
	std::queue<COpenMeshT::FaceHandle>Q;
	Q.push(mesh.faces_begin());
	std::vector<int>ftags;
	ftags.resize(mesh.n_faces(), 0);
	ftags[mesh.faces_begin()->idx()] = 1;
	while (!Q.empty())
	{
		COpenMeshT::FaceHandle fh = Q.front();
		Q.pop();
		for (auto ffiter = mesh.ff_begin(fh); ffiter != mesh.ff_end(fh); ffiter++)
		{
			if (ftags[ffiter->idx()] == 1)
				continue;
			COpenMeshT::EdgeHandle comm_edge;
			bool have_comm_edge = CGeoBaseAlg::IsShareCommonEdge(mesh, fh, ffiter, comm_edge);
			if (have_comm_edge)
			{
				COpenMeshT::HalfedgeHandle he0 = mesh.halfedge_handle(comm_edge, 0), he1 = mesh.halfedge_handle(comm_edge, 1);
				COpenMeshT::VertexHandle vh0 = mesh.to_vertex_handle(he0);
				COpenMeshT::VertexHandle vh1 = mesh.to_vertex_handle(he1);
				
				/*int bound_count = 0;
				for (auto fviter = mesh.fv_begin(ffiter); fviter != mesh.fv_end(ffiter); fviter++)
				{
					if (vmark[fviter->idx()] == true)
					{
						bound_count++;
					}

				}*/
				if (vmark[vh0.idx()] == false || vmark[vh1.idx()] == false )
				{
					ftags[ffiter->idx()] = 1;
					Q.push(ffiter);
				}
				else 
				{
					if (loop_edge_set.find(std::make_pair(vh0, vh1)) == loop_edge_set.end())
					{
						ftags[ffiter->idx()] = 1;
						Q.push(ffiter);
					}
				}

			}
			
			
		}
		/*for (auto fviter = mesh.fv_begin(fh); fviter != mesh.fv_end(fh); fviter++)
		{
			if (vmark[fviter->idx()] == false)
			{
				for (auto vfiter = mesh.vf_begin(fviter); vfiter != mesh.vf_end(fviter); vfiter++)
				{
					if (ftags[vfiter->idx()] == 0)
					{
						ftags[vfiter->idx()] = 1;
						Q.push(vfiter);
					}
					
				}
			}
		}*/
	}
	std::map<int, COpenMeshT*>sep_meshes;
	std::map<int, std::map<COpenMeshT::VertexHandle, COpenMeshT::VertexHandle>>tag_new_orig_vhsmap,tag_orig_new_vhsmap;
	SeparateMeshByFaceTag(mesh, ftags, sep_meshes, tag_new_orig_vhsmap);
	for (auto titer = tag_new_orig_vhsmap.begin(); titer != tag_new_orig_vhsmap.end(); titer++)
	{
		tag_orig_new_vhsmap[titer->first] = std::map<COpenMeshT::VertexHandle, COpenMeshT::VertexHandle>();
		for (auto iter = titer->second.begin(); iter != titer->second.end(); iter++)
		{
			tag_orig_new_vhsmap[titer->first][iter->second] = iter->first;

		}
	}
	for (auto titer = tag_orig_new_vhsmap.begin(); titer != tag_orig_new_vhsmap.end(); titer++)
	{
		auto mmap = titer->second;
		for (int i = 0; i < vhs.size(); i++)
		{
			if (mmap.find(vhs[i]) == mmap.end())
			{
				std::cerr <<titer->first<< " mmap.find(vhs[i]) == mmap.end() vh id: "<<vhs[i].idx() << std::endl;
			}
		}

	}
	res_meshes.clear();
	new2orig_vhsmap.clear();
	for (auto iter = sep_meshes.begin(); iter != sep_meshes.end(); iter++)
	{
		res_meshes.push_back(iter->second);
		new2orig_vhsmap.push_back(tag_new_orig_vhsmap[iter->first]);

	}

}
void CGeoAlg::GetSilhouetteEdges(COpenMeshT&mesh, OpenMesh::Vec3d view_dir,std::vector<std::pair<COpenMeshT::VertexHandle, COpenMeshT::VertexHandle>>&res_edges)
{
	mesh.request_vertex_normals();
	for (auto eiter = mesh.edges_begin(); eiter != mesh.edges_end(); eiter++)
	{
		COpenMeshT::VertexHandle vh0=mesh.to_vertex_handle(mesh.halfedge_handle(eiter, 0));
		COpenMeshT::VertexHandle vh1 = mesh.to_vertex_handle(mesh.halfedge_handle(eiter, 1));
		OpenMesh::Vec3d normal0=mesh.normal(vh0);
		OpenMesh::Vec3d normal1 = mesh.normal(vh1);
		if (OpenMesh::dot(normal0, view_dir)*OpenMesh::dot(normal1, view_dir) < 0)
		{
			res_edges.push_back(std::make_pair(vh0,vh1));
		}
	}
}

void CGeoAlg::DelaunayTriangulation3d(std::vector<Eigen::Vector3d>&points, Eigen::MatrixXd &res_v, Eigen::MatrixXi &res_f, std::vector<int>&res_orig_vidmap)
{
	Eigen::Vector3d axis_01 = points[1]- points[0];
	Eigen::Vector3d axis_02 = points[2]- points[0];

	Eigen::Vector3d axis_x = axis_01.normalized();
	Eigen::Vector3d axis_z = axis_01.cross(axis_02).normalized();
	Eigen::Vector3d axis_y = axis_z.cross(axis_x).normalized();

	std::vector<Eigen::Vector2d>points_2d;
	for (int i = 0; i < points.size(); i++)
	{
		Eigen::Vector3d dir = points[i] - points[0];
		double x = dir.dot(axis_x);
		double y = dir.dot(axis_y);
		points_2d.push_back(Eigen::Vector2d(x, y));
	}

	
	Eigen::MatrixXd tri_v2d;
	Eigen::MatrixXi tri_f;
	std::vector<int>vid_map;
	DelaunayTriangulation(points_2d, tri_v2d, tri_f, vid_map);
	std::vector<int>new2orig_vid_map;
	new2orig_vid_map.resize(points_2d.size());
	for (int i = 0; i < vid_map.size(); i++)
	{
		new2orig_vid_map[vid_map[i]] = i;
	}
	res_v.resize(tri_v2d.rows(), 3);
	res_f.resize(tri_f.rows(), 3);
	res_orig_vidmap = vid_map;
	for (int i = 0; i < points.size(); i++)
	{
		res_v(i, 0) = points[i].x();
		res_v(i, 1) = points[i].y();
		res_v(i, 2) = points[i].z();
	}
	for (int i = 0; i < res_f.rows(); i++)
	{
		res_f(i, 0) = new2orig_vid_map[tri_f(i, 0)];
		res_f(i, 1) = new2orig_vid_map[tri_f(i, 1)];
		res_f(i, 2) = new2orig_vid_map[tri_f(i, 2)];
	}
}

void CGeoAlg::ConstrainedTriangulation3d(std::vector<Eigen::Vector3d>&points, std::vector<std::pair<int, int>>&cons,Eigen::MatrixXd &res_v, Eigen::MatrixXi &res_f, std::vector<int>&res_orig2new_vidmap)
{
	Eigen::Vector3d axis_01 = points[1] - points[0];
	Eigen::Vector3d axis_02 = points[2] - points[0];

	Eigen::Vector3d axis_x = axis_01.normalized();
	Eigen::Vector3d axis_z = axis_01.cross(axis_02).normalized();
	Eigen::Vector3d axis_y = axis_z.cross(axis_x).normalized();

	std::vector<Eigen::Vector2d>points_2d;
	for (int i = 0; i < points.size(); i++)
	{
		Eigen::Vector3d dir = points[i] - points[0];
		double x = dir.dot(axis_x);
		double y = dir.dot(axis_y);
		points_2d.push_back(Eigen::Vector2d(x, y));
	}


	Eigen::MatrixXd tri_v2d;
	Eigen::MatrixXi tri_f;
	std::vector<int>vid_map;
	ConstrainedTriangulation(points_2d, cons, tri_v2d, tri_f, vid_map);
	std::vector<int>new2orig_vid_map;
	new2orig_vid_map.resize(tri_v2d.rows(),-1);
	for (int i = 0; i < vid_map.size(); i++)
	{
		new2orig_vid_map[vid_map[i]] = i;
	}
	res_v.resize(tri_v2d.rows(), 3);
	res_f.resize(tri_f.rows(), 3);
	res_orig2new_vidmap.resize(points.size());

	for (int i = 0; i < points.size(); i++)
	{
		res_v(i, 0) = points[i].x();
		res_v(i, 1) = points[i].y();
		res_v(i, 2) = points[i].z();
		res_orig2new_vidmap[i] = i;
	}


	for (int i = 0; i < res_f.rows(); i++)
	{
		res_f(i, 0) = new2orig_vid_map[tri_f(i, 0)];
		res_f(i, 1) = new2orig_vid_map[tri_f(i, 1)];
		res_f(i, 2) = new2orig_vid_map[tri_f(i, 2)];
	}
}

void CGeoAlg::ConstrainedTriangulation(std::vector<Eigen::Vector2d>&points, std::vector<std::pair<int, int>>&cons, Eigen::MatrixXd &res_v, Eigen::MatrixXi &res_f, std::vector<int>&orig_vidmap)
{
	
	typedef CGAL::Exact_predicates_inexact_constructions_kernel         K;
	typedef CGAL::Exact_predicates_tag                               Itag;
	typedef CGAL::Constrained_Delaunay_triangulation_2<K, CGAL::Default, Itag> CDT;
	typedef CDT::Point          TriPoint;
	
	std::vector< TriPoint > tri_points;
	std::map<TriPoint, int>id_map;
	orig_vidmap.clear();
	for (int i = 0; i < points.size(); i++)
	{
		TriPoint tp=TriPoint(points[i][0], points[i][1]);
		tri_points.push_back(tp);
		orig_vidmap.push_back(i); 
		id_map[tp] = i;
	}
	CDT dt;

	dt.insert(tri_points.begin(), tri_points.end());
	for (int i = 0; i < cons.size(); i++)
	{
		dt.insert_constraint(tri_points[cons[i].first], tri_points[cons[i].second]);
	}
	res_v.resize(dt.number_of_vertices(), 2);
	res_f.resize(dt.number_of_faces(), 3);

	for (auto viter = dt.finite_vertices_begin(); viter != dt.finite_vertices_end(); viter++)
	{
		res_v(id_map[viter->point()], 0) = viter->point().x();
		res_v(id_map[viter->point()], 1) = viter->point().y();

	}
	int fid = 0;
	for (auto fiter = dt.finite_faces_begin(); fiter != dt.finite_faces_end(); fiter++, fid++)
	{

		for (int i = 0; i < 3; i++)
		{
			res_f(fid, i) = id_map[fiter->vertex(i)->point()];

		}

	}

}
void CGeoAlg::DelaunayTriangulation(std::vector<Eigen::Vector2d>&points, Eigen::MatrixXd &res_v, Eigen::MatrixXi &res_f, std::vector<int>&orig_vidmap )
{

	typedef CGAL::Exact_predicates_inexact_constructions_kernel         K;
	typedef CGAL::Triangulation_vertex_base_with_info_2<unsigned, K>    Vb;
	typedef CGAL::Triangulation_data_structure_2<Vb>                    Tds;
	typedef CGAL::Delaunay_triangulation_2<K, Tds>                      Delaunay;
	typedef Delaunay::Point                                             TriPoint;
	std::vector< std::pair<TriPoint, unsigned> > tri_points;
	orig_vidmap.clear();
	for (int i = 0; i < points.size(); i++)
	{
		tri_points.push_back(std::make_pair(TriPoint(points[i][0], points[i][1]),i));
		orig_vidmap.push_back(i);
	}
	Delaunay dt;
	dt.insert(tri_points.begin(), tri_points.end());
	res_v.resize(tri_points.size(), 2);
	res_f.resize(dt.number_of_faces(), 3);

	
	for (auto viter = dt.finite_vertices_begin(); viter != dt.finite_vertices_end(); viter++)
	{
		res_v(viter->info(), 0) = viter->point().x();
		res_v(viter->info(), 1) = viter->point().y();

	}
	int fid = 0;
	for (auto fiter = dt.finite_faces_begin(); fiter != dt.finite_faces_end(); fiter++,fid++)
	{
		
		for (int i = 0; i < 3; i++)
		{
			res_f(fid, i) = fiter->vertex(i)->info();
		
		}
		
	}
	
	
		
	
}
void CGeoAlg::Triangulate(std::vector<Eigen::Vector2d>&boundv, std::vector<Eigen::Vector2d>&inner_pts, Eigen::MatrixXd &res_v, Eigen::MatrixXi &res_f, std::vector<int>&bound_vids, std::vector<int>&inner_vids, double max_edge_len)
{
	typedef CGAL::Delaunay_mesh_vertex_base_2<Kernel>Vbase;

	typedef CGAL::Triangulation_vertex_base_with_info_2<int, Kernel, Vbase> Vb;
	typedef CGAL::Delaunay_mesh_face_base_2<Kernel> Fb;
	typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
	typedef CGAL::Constrained_Delaunay_triangulation_2<Kernel, Tds> CDT;
	typedef CGAL::Delaunay_mesh_size_criteria_2<CDT> Criteria;
	typedef CGAL::Delaunay_mesher_2<CDT, Criteria> Mesher;
	typedef CDT::Vertex_handle CDT_Vertex_handle;
	typedef CDT::Point CDT_Point;
	std::cerr << "size of boundv " << boundv.size() << std::endl;

	std::vector<CDT_Vertex_handle> cdt_boundv_handles(boundv.size());

	CDT cdt;
	for (int i = 0; i < boundv.size(); i++)
	{
		std::cerr << boundv[i][0]<<" "<<boundv[i][1] << std::endl;
		cdt_boundv_handles[i] = cdt.insert(CDT_Point(boundv[i].x(), boundv[i].y()));
	}
	for (int i = 0; i <  boundv.size(); i++)
	{
		cdt.insert_constraint(cdt_boundv_handles[i], cdt_boundv_handles[(i + 1) % boundv.size()]);
	}
	std::vector<CDT_Vertex_handle>cdt_inner_v_handles(inner_pts.size());

	for (int i = 0; i < inner_pts.size(); i++)
	{
		std::cerr << inner_pts[i][0] << " " << inner_pts[i][1] << std::endl;
		cdt_inner_v_handles[i]= cdt.insert(CDT_Point(inner_pts[i].x(), inner_pts[i].y()));
	}
	if (max_edge_len != -1)
	{
		CGAL::refine_Delaunay_mesh_2(cdt, Criteria(125, max_edge_len));
	}
	else 
	{
		CGAL::refine_Delaunay_mesh_2(cdt, Criteria(std::numeric_limits<double>::max(), std::numeric_limits<double>::max()));
	}


	for (auto viter = cdt.finite_vertices_begin(); viter != cdt.finite_vertices_end(); viter++)
	{
		viter->info() = -1;//mark all vertices to -1
	}
	
	int n_faces = 0;
	int n_vertices = 0;
	for(auto fiter=cdt.finite_faces_begin();fiter!=cdt.finite_faces_end();fiter++)
	{
		if (fiter->is_in_domain())
		{
			for (int i = 0; i < 3; i++)
			{
				if (fiter->vertex(i)->info() == -1)
					n_vertices++;
				fiter->vertex(i)->info() = -2;//mark all vertices in domain to -2
			}
			n_faces++;
		}

		
	}
	
	res_v.resize(n_vertices, 2);
	res_f.resize(n_faces, 3);
	int vi = 0;
	for (auto viter = cdt.finite_vertices_begin(); viter != cdt.finite_vertices_end(); viter++)
	{
		if (viter->info() == -2)
		{
			res_v(vi, 0) = viter->point().x();
			res_v(vi, 1) = viter->point().y();
			viter->info() = vi;
			vi++;
		}
	}
	int fi = 0;
	for (auto fiter = cdt.finite_faces_begin(); fiter != cdt.finite_faces_end(); fiter++)
	{
		if (fiter->is_in_domain())
		{
			for (int i = 0; i < 3; i++)
			{
				int vid = fiter->vertex(i)->info();
				res_f(fi, i) = vid;
			}
			fi++;
		}
	}
	bound_vids.resize(cdt_boundv_handles.size());
	for (int i = 0; i < cdt_boundv_handles.size(); i++)
	{
		bound_vids[i] = cdt_boundv_handles[i]->info();
	}
	inner_vids.resize(cdt_inner_v_handles.size());
	for (int i = 0; i < cdt_inner_v_handles.size(); i++)
	{
		inner_vids[i] = cdt_inner_v_handles[i]->info();
	}

}
void CGeoAlg::ComputeBaryCoord(OpenMesh::Vec2d a, OpenMesh::Vec2d b, OpenMesh::Vec2d c, OpenMesh::Vec2d p, double &res_w0, double &res_w1, double &res_w2)
{
	std::vector<OpenMesh::Vec3d>ps;
	ps.resize(3);
	ps[0] = OpenMesh::Vec3d(a[0],a[1],0);
	ps[1] = OpenMesh::Vec3d(b[0], b[1], 0);
	ps[2] = OpenMesh::Vec3d(c[0], c[1], 0);
	OpenMesh::Vec3d p3d = OpenMesh::Vec3d(p[0], p[1], p[2]);
	OpenMesh::Vec3d  p0 = p3d - ps[0];
	OpenMesh::Vec3d  p1 = p3d - ps[1];
	OpenMesh::Vec3d  p2 = p3d - ps[2];
	OpenMesh::Vec3d  v01 = ps[1] - ps[0];
	OpenMesh::Vec3d  v02 = ps[2] - ps[0];
	double area = OpenMesh::cross(v01, v02).length() / 2;
	double area0 = OpenMesh::cross(p1, p2).length() / 2;
	double area1 = OpenMesh::cross(p0, p2).length() / 2;

	res_w0 = area0 / area;
	res_w1 = area1 / area;
	res_w2 = 1 - res_w0 - res_w1;
}
void CGeoAlg::ComputeBaryCoord(OpenMesh::Vec3d a, OpenMesh::Vec3d b, OpenMesh::Vec3d c, OpenMesh::Vec3d p, double &res_w0,  double &res_w1, double &res_w2)
{
	std::vector<OpenMesh::Vec3d>ps;
	ps.resize(3);
	ps[0] = a;
	ps[1] = b;
	ps[2] = c;

	OpenMesh::Vec3d  p0 = p - ps[0];
	OpenMesh::Vec3d  p1 = p - ps[1];
	OpenMesh::Vec3d  p2 = p - ps[2];
	OpenMesh::Vec3d  v01 = ps[1] - ps[0];
	OpenMesh::Vec3d  v02 = ps[2] - ps[0];
	double area = OpenMesh::cross(v01, v02).length() / 2;
	double area0 = OpenMesh::cross(p1, p2).length() / 2;
	double area1 = OpenMesh::cross(p0, p2).length() / 2;

	res_w0 = area0 / area;
	res_w1 = area1 / area;
	res_w2 = 1 - res_w0 - res_w1;
}
void CGeoAlg::ComputeBaryCoord(COpenMeshT& mesh, COpenMeshT::FaceHandle  fh, OpenMesh::Vec3d  p, double &w0, double &w1, double &w2)
{
	std::vector<OpenMesh::Vec3d>ps;
	int count = 0;
	
	for (auto viter = mesh.fv_begin(fh); viter != mesh.fv_end(fh); viter++)
	{
		ps.push_back(mesh.point(viter));
	}
	
	assert(count == 3);

	OpenMesh::Vec3d  p0 = p - ps[0];
	OpenMesh::Vec3d  p1 = p - ps[1];
	OpenMesh::Vec3d  p2 = p - ps[2];
	OpenMesh::Vec3d  v01 = ps[1] - ps[0];
	OpenMesh::Vec3d  v02 = ps[2] - ps[0];
	double area = OpenMesh::cross(v01, v02).length() / 2;
	double area0 = OpenMesh::cross(p1, p2).length() / 2;
	double area1 = OpenMesh::cross(p0, p2).length() / 2;

	w0 = area0 / area;
	w1 = area1 / area;
	w2 = 1 - w0 - w1;
}
void CGeoAlg::SplitFaceByPoints(COpenMeshT &mesh, COpenMeshT::FaceHandle fh, std::vector<OpenMesh::Vec3d> &pts, std::vector<OpenMesh::VertexHandle>&new_vhs)
{
	std::vector<OpenMesh::Vec3d>fs;
	std::vector<OpenMesh::VertexHandle>face_vhs;
	for (auto viter = mesh.fv_begin(fh); viter != mesh.fv_end(fh); viter++)
	{
		OpenMesh::Vec3d p = mesh.point(viter);
		fs.push_back(p);
		face_vhs.push_back(viter);
	}
	OpenMesh::Vec3d axis_01 = fs[1] - fs[0];
	OpenMesh::Vec3d axis_02 = fs[2] - fs[0];

	OpenMesh::Vec3d axis_x = axis_01.normalized();
	OpenMesh::Vec3d axis_z = OpenMesh::cross(axis_01,axis_02).normalized();
	OpenMesh::Vec3d axis_y = OpenMesh::cross(axis_z, axis_x).normalized();
	
	std::vector<Eigen::Vector2d>boundvs_2d,inner_vs2d;
	for (int i = 0; i < fs.size(); i++)
	{
		OpenMesh::Vec3d dir = fs[i] - fs[0];
		double x = OpenMesh::dot( dir,axis_x);
		double y = OpenMesh::dot(dir,axis_y);
		boundvs_2d.push_back(Eigen::Vector2d(x, y));
	}

	for (int i = 0; i < pts.size(); i++)
	{
		OpenMesh::Vec3d dir = pts[i] - fs[0];
		double x = OpenMesh::dot(dir, axis_x);
		double y = OpenMesh::dot(dir, axis_y);
		inner_vs2d.push_back(Eigen::Vector2d(x, y));
	}
	Eigen::MatrixXd tri_v2d;
	Eigen::MatrixXi tri_f;
	std::vector<int>bound_vids, inner_vids;
	Triangulate(boundvs_2d, inner_vs2d, tri_v2d, tri_f, bound_vids, inner_vids, -1);
	std::vector<bool>is_boundv;
	is_boundv.resize(tri_v2d.rows(), 0);

	std::map<int, COpenMeshT::VertexHandle>ps_handle_map;
	for (int i = 0; i < bound_vids.size(); i++)
	{
		is_boundv[bound_vids[i]] = true;
		ps_handle_map[bound_vids[i]] = face_vhs[i];
	}
	
	new_vhs.clear();
	for (int i = 0; i < is_boundv.size(); i++)
	{
		if (!is_boundv[i])
		{
			Eigen::Vector2d p2d = tri_v2d.row(i);
			OpenMesh::Vec3d p = p2d.x()*axis_x + p2d.y()*axis_y;
			ps_handle_map[i] = mesh.add_vertex(p);
			new_vhs.push_back(ps_handle_map[i]);
		}
	}

	mesh.delete_face(fh);
	for (int i = 0; i < tri_f.rows(); i++)
	{
		mesh.add_face(ps_handle_map[tri_f(i, 0)], ps_handle_map[tri_f(i, 1)], ps_handle_map[tri_f(i, 2)]);
	}
	
}
double CGeoAlg::ComputeAverageEdgeLength(COpenMeshT &mesh)
{
	int count = 0;
	double avelen = 0;
	for (auto hiter = mesh.halfedges_begin(); hiter != mesh.halfedges_begin(); hiter++)
	{
		OpenMesh::Vec3d p0=mesh.point(mesh.to_vertex_handle(hiter));
		OpenMesh::Vec3d p1 = mesh.point(mesh.to_vertex_handle(mesh.opposite_halfedge_handle(hiter)));
		OpenMesh::Vec3d dif = (p0 - p1);
		avelen+=dif.length();
		count++;
	}
	avelen /= count;
	return avelen;
}
void CGeoAlg::Remeshing(COpenMeshT &mesh, std::vector<COpenMeshT::VertexHandle>&roi)
{
	mesh.request_face_normals();
	mesh.request_vertex_normals();
	mesh.update_face_normals();
	mesh.update_vertex_normals();
	//Remeshing::UniformRemesherT<Surface_Mesh> uni_remesher(openmesh);
	std::vector<int>selected_indices;
	for (int i = 0; i <roi.size(); i++)
	{
		selected_indices.push_back(roi[i].idx());
	}
	MeshSelection::selectVertices(&mesh, selected_indices);
	Remeshing::AdaptiveRemesherT<COpenMeshT> adapt_remesher(mesh);
	//Remeshing::UniformRemesherT<COpenMeshT> uni_remesher(mesh);
	double ave_edge_len = ComputeAverageEdgeLength(mesh);
	adapt_remesher.remesh(ave_edge_len, ave_edge_len * 0.1, ave_edge_len, 5);
}
void CGeoAlg::Remeshing(COpenMeshT &mesh)
{
	mesh.request_face_normals();
	mesh.request_vertex_normals();
	mesh.update_face_normals();
	mesh.update_vertex_normals();
	//Remeshing::UniformRemesherT<Surface_Mesh> uni_remesher(openmesh);
	std::vector<int>selected_indices;
	for (int i = 0; i < mesh.n_vertices(); i++)
	{
		selected_indices.push_back(i);
	}
	MeshSelection::selectVertices(&mesh, selected_indices);
	Remeshing::AdaptiveRemesherT<COpenMeshT> adapt_remesher(mesh);
	//Remeshing::UniformRemesherT<COpenMeshT> uni_remesher(mesh);
	double ave_edge_len = ComputeAverageEdgeLength(mesh);
	adapt_remesher.remesh(ave_edge_len, ave_edge_len * 0.1, ave_edge_len , 5);
}
void CGeoAlg::RefineMeshBoundary(COpenMeshT &mesh)
{
	std::vector<std::vector<COpenMeshT::VertexHandle>> bounds;
	std::vector<COpenMeshT::VertexHandle>to_be_deleted;
	for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); viter++)
	{
		if (mesh.is_boundary(viter) && mesh.is_manifold(viter) == false)
		{
			to_be_deleted.push_back(viter);
		}
	}
	for (int i = 0; i < to_be_deleted.size(); i++)
	{
		mesh.delete_vertex(to_be_deleted[i]);
	}
	mesh.garbage_collection();
	std::cerr << "RefineMeshBoundary " << std::endl;
	to_be_deleted.clear();
	std::vector<int>mmark(mesh.n_vertices(), -1);
	std::vector<std::vector<OpenMesh::VertexHandle>>vhs_of_t;
	int t = 0;
	for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); viter++)
	{
		if (mesh.is_valid_handle(viter) == false)
			continue;
		int vid = viter->idx();
		if (mmark[vid] != -1)
			continue;
		std::queue<COpenMeshT::VertexHandle>Q;
		Q.push(viter);
		mmark[vid] = t;
	
		vhs_of_t.push_back(std::vector<OpenMesh::VertexHandle>());
		vhs_of_t[t].push_back(viter);
		while (!Q.empty())
		{
			auto pviter = Q.front();

			Q.pop();
			for (auto vviter = mesh.vv_begin(pviter); vviter != mesh.vv_end(pviter); vviter++)
			{
				if (mmark[vviter->idx()] == -1)
				{
					Q.push(vviter);
					mmark[vviter->idx()] = t;
					vhs_of_t[t].push_back(vviter);
				}

			}
		}
		
		t++;
	}

	std::cerr << "RefineMeshBoundary " << std::endl;
	for(int i=0;i<vhs_of_t.size();i++)
	{
		if (vhs_of_t[i].size() < 20)
		{
			for (int j = 0; j < vhs_of_t[i].size(); j++)
			{
				to_be_deleted.push_back(vhs_of_t[i][j]);
			}
			
		}
	}
	std::cerr << "RefineMeshBoundary " << std::endl;

	for (int i = 0; i < to_be_deleted.size(); i++)
	{
		mesh.delete_vertex(to_be_deleted[i]);
	}
	mesh.garbage_collection();
}
int CGeoAlg::PickMesh(OpenMesh::Vec3d  source, OpenMesh::Vec3d dir, std::map<int, std::shared_ptr<CMeshObject>>&data_pool,bool is_visiable)
{
	double min_dis = std::numeric_limits<double>::max();
	int mi=-1;
	for (auto iter = data_pool.begin(); iter != data_pool.end(); iter++)
	{
		/*auto camera = viewer_->GetCamera();
		OpenMesh::Vec3d orig, dir;
		camera.ConvertClickToLine(e->pos(), orig, dir);*/

		CMeshObject* meshobj = iter->second.get();
		if (meshobj->IsVisiable() != is_visiable||meshobj->IsPickAble()==false)
			continue;
		COpenMeshT&mesh = meshobj->GetMesh();
		COpenMeshT::FaceHandle fh;
		OpenMesh::Vec3d barycoord;
	
		//std::cerr << "start ray" << std::endl;
	//	std::cerr << source << " " << dir << std::endl;
		if (CGeoAlg::RayMeshIntersection(source, dir, *meshobj,fh, barycoord))
		{
			//std::cerr << "inter" << std::endl;
			OpenMesh::Vec3d p = CGeoBaseAlg::ComputeVPosFromBaryCoord(mesh,fh, barycoord);
			
			double dis=(p - source).length();
			if (min_dis > dis)
			{
				min_dis = dis;
				mi = iter->first;
			}
		}
	}
	return mi;

}
void CGeoAlg::ComputeGradientOfScalarField(COpenMeshT &mesh, Eigen::VectorXd &scalars, Eigen::MatrixXd & res_grad)
{
	Eigen::MatrixXd vertexs;
	Eigen::MatrixXi faces;

	CConverter::ConvertFromOpenMeshToIGL(mesh, vertexs, faces);
	Eigen::SparseMatrix<double> G;
	igl::grad(vertexs, faces, G);
	// Compute gradient of U
	res_grad = Eigen::Map<const Eigen::MatrixXd>((G*scalars).eval().data(), faces.rows(), 3);
	
}
void CGeoAlg::InitGeodesicModel(CMeshObject &mesh_obj)
{
	auto geodesic_model = mesh_obj.GetGeodesicModel();
	COpenMeshT &mesh = mesh_obj.GetMesh();

	if (geodesic_model == NULL)
	//if (true)
	{
		geodesic_model = (new CXWGeodesic());
		std::vector<CPoint3D>vertexs(mesh.n_vertices());
		std::vector<CBaseModel::CFace>faces(mesh.n_faces());
		
		for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); viter++)
		{
			auto p = mesh.point(viter);
			if (viter->idx() >= vertexs.size())
				std::cerr << "error" << std::endl;
			vertexs[viter->idx()] = CPoint3D(p[0], p[1], p[2]);
		}
		for (auto fiter = mesh.faces_begin(); fiter != mesh.faces_end(); fiter++)
		{
			CBaseModel::CFace f;
			int i = 0;

			for (auto fviter = mesh.fv_begin(fiter); fviter!=mesh.fv_end(fiter),i < 3; fviter++, i++)
			{
				f.verts[i] = fviter->idx();
			}
			faces[fiter->idx()] = f;
		}
		geodesic_model->SetModel(vertexs, faces);
		mesh_obj.SetGeodesicModel(geodesic_model);
		std::cerr << "init geo model finished" << std::endl;
	}
}
void CGeoAlg::ComputeGeodesicDis(CMeshObject &mesh_obj, COpenMeshT::VertexHandle svh, std::vector<COpenMeshT::VertexHandle>&dst_vhs, std::vector<double>&res_dis)
{
	InitGeodesicModel(mesh_obj);
	auto geodesic_model = mesh_obj.GetGeodesicModel();
	COpenMeshT &mesh = mesh_obj.GetMesh();
	
	std::vector<CGeoFacePoint>tpath;
	std::vector<int>tvids(0);
	for (int i = 0; i < dst_vhs.size(); i++)
	{
		tvids.push_back(dst_vhs[i].idx());
	}

	geodesic_model->GeodesicDis(svh.idx(), tvids, res_dis);
	

	
}
bool CGeoAlg::ComputeGeodesicPath(CMeshObject &mesh_obj, int svid, int tvid, std::vector<OpenMesh::Vec3d>&path)
{
	InitGeodesicModel(mesh_obj);
	auto geodesic_model = mesh_obj.GetGeodesicModel();
	COpenMeshT &mesh = mesh_obj.GetMesh();
	std::vector<CGeoFacePoint>tpath;
	geodesic_model->GeodesicPath(svid, tvid, tpath);
	path.resize(tpath.size());
	for (int i = 0; i < tpath.size(); i++)
	{

		
		path[i] = OpenMesh::Vec3d(tpath[i].pos_[0], tpath[i].pos_[1], tpath[i].pos_[2]);
	}

	return true;
}
void CGeoAlg::ComputeGeodesicDis(CMeshObject &mesh_obj, std::vector<COpenMeshT::VertexHandle>&svhs, std::vector<double>&res_dis)
{
	InitGeodesicModel(mesh_obj);
	auto geodesic_model = mesh_obj.GetGeodesicModel();
	COpenMeshT &mesh = mesh_obj.GetMesh();
	std::vector<std::vector<CGeoFacePoint>>tpaths;
	
	std::vector<int>svids(svhs.size());
	for (int i = 0; i < svhs.size(); i++)
	{
		svids[i] = svhs[i].idx();
	}
	geodesic_model->GeodesicDis(svids, res_dis);
}
bool CGeoAlg::ComputeGeodesicPath(CMeshObject &mesh_obj, int svid, int tvid, std::vector<COpenMeshT::FaceHandle>&res_fhs, std::vector<OpenMesh::Vec3d>&res_bary_coords)
{
	InitGeodesicModel(mesh_obj);
	auto geodesic_model = mesh_obj.GetGeodesicModel();
	COpenMeshT &mesh = mesh_obj.GetMesh();
	std::vector<CGeoFacePoint>tpath;
	geodesic_model->GeodesicPath(svid, tvid, tpath);
	for (int i = 0; i < tpath.size(); i++)
	{
		std::cerr << tpath[i].vids_[0] << std::endl;
	}
	res_fhs.resize(tpath.size());
	res_bary_coords.resize(tpath.size());

	for (int i = 0; i < tpath.size(); i++)
	{
		
		res_fhs[i] = mesh.face_handle(tpath[i].fid_);
		//std::cerr << res_fhs[i].idx() << " " << tpath[i].fid_ << std::endl;
		if (res_fhs[i].idx() != tpath[i].fid_)
		{
			std::cerr << "ddddddddddddddddddddddddddddddd" << std::endl;
		}
		std::vector<int>vids;
		for (auto viter = mesh.fv_begin(res_fhs[i]); viter!= mesh.fv_end(res_fhs[i]);viter++)
		{
			vids.push_back(viter->idx());
		}

			if (vids[0] != tpath[i].vids_[0])
			{
				std::cerr << "error vids[i]!=tpath[i].vids_[i] " << vids[i] << " " << tpath[0].vids_[0] << std::endl;
				break;
			}
				
	
		res_bary_coords[i] = OpenMesh::Vec3d(tpath[i].ls_[0], tpath[i].ls_[1], tpath[i].ls_[2]);
	}

	return true;
}
bool CGeoAlg::ComputeGeodesicPath(CMeshObject &mesh_obj, std::vector<COpenMeshT::VertexHandle>&svhs, std::vector<COpenMeshT::VertexHandle>&tvhs, std::vector<std::vector<COpenMeshT::FaceHandle>>&res_fhs, std::vector<std::vector<OpenMesh::Vec3d>>&res_bary_coords)
{
	InitGeodesicModel(mesh_obj);
	auto geodesic_model = mesh_obj.GetGeodesicModel();
	COpenMeshT &mesh = mesh_obj.GetMesh();
	std::vector<std::vector<CGeoFacePoint>>tpaths;
	std::vector<int>tvids(tvhs.size());
	for (int i = 0; i < tvhs.size(); i++)
	{
		tvids[i] = tvhs[i].idx();
	}
	std::vector<int>svids(svhs.size());
	for (int i = 0; i < svhs.size(); i++)
	{
		svids[i] = svhs[i].idx();
	}
	geodesic_model->GeodesicPath(svids, tvids, tpaths);
	res_fhs.resize(tpaths.size());
	res_bary_coords.resize(tpaths.size());
	for (int i = 0; i < tpaths.size(); i++)
	{
		res_fhs[i].resize(tpaths[i].size());
		res_bary_coords[i].resize(tpaths[i].size());
		for (int j = 0; j < tpaths[i].size(); j++)
		{
			res_fhs[i][j] = mesh.face_handle(tpaths[i][j].fid_);
			res_bary_coords[i][j] = OpenMesh::Vec3d(tpaths[i][j].ls_[0], tpaths[i][j].ls_[1], tpaths[i][j].ls_[2]);
		}
	
		//std::cerr << res_bary_coords[i] << std::endl;
	}

	return true;
}
bool CGeoAlg::ComputeGeodesicPath(CMeshObject &mesh_obj, COpenMeshT::VertexHandle svh, std::vector<COpenMeshT::VertexHandle>& tvhs, std::vector<COpenMeshT::FaceHandle>&res_fhs, std::vector<OpenMesh::Vec3d>&res_bary_coords)
{
	InitGeodesicModel(mesh_obj);
	auto geodesic_model = mesh_obj.GetGeodesicModel();
	COpenMeshT &mesh = mesh_obj.GetMesh();
	std::vector<CGeoFacePoint>tpath;
	std::vector<int>tvids(tvhs.size());
	for (int i = 0; i < tvhs.size(); i++)
	{
		tvids[i] = tvhs[i].idx();
	}
	geodesic_model->GeodesicPath(svh.idx(), tvids, tpath);
	res_fhs.resize(tpath.size());
	res_bary_coords.resize(tpath.size());
	for (int i = 0; i < tpath.size(); i++)
	{
		res_fhs[i] = mesh.face_handle(tpath[i].fid_);
		res_bary_coords[i] = OpenMesh::Vec3d(tpath[i].ls_[0], tpath[i].ls_[1], tpath[i].ls_[2]);
		//std::cerr << res_bary_coords[i] << std::endl;
	}

	return true;
}
bool CGeoAlg::ComputeGeodesicPath(CMeshObject &mesh_obj, COpenMeshT::VertexHandle svh, std::vector<COpenMeshT::VertexHandle>& tvhs, std::vector<OpenMesh::Vec3d>&path)
{
	InitGeodesicModel(mesh_obj);
	auto geodesic_model = mesh_obj.GetGeodesicModel();
	COpenMeshT &mesh = mesh_obj.GetMesh();
	std::vector<CGeoFacePoint>tpath;
	std::vector<int>tvids(tvhs.size());
	for (int i = 0; i < tvhs.size(); i++)
	{
		tvids[i] = tvhs[i].idx();
	}
	geodesic_model->GeodesicPath(svh.idx(), tvids, tpath);
	path.resize(tpath.size());
	
	for (int i = 0; i < tpath.size(); i++)
	{
		path[i] = OpenMesh::Vec3d(tpath[i].pos_[0], tpath[i].pos_[1], tpath[i].pos_[2]);
	}

	return true;
}
bool CGeoAlg::ComputeShortestPathAlongEdge(COpenMeshT &mesh, COpenMeshT::VertexHandle svh, COpenMeshT::VertexHandle tvh, std::vector<COpenMeshT::VertexHandle>&path)
{
	std::queue<COpenMeshT::VertexHandle>Q;
	std::vector<bool>mmark(mesh.n_vertices(), false);
	Q.push(tvh);
	mmark[tvh.idx()] = true;
	std::vector<double>dis(mesh.n_vertices(), -1);
	dis[tvh.idx()] = 0;
	std::vector<int>pre(mesh.n_vertices());
	pre[tvh.idx()] = -1;
	double prune_threshold = -1;
	while (!Q.empty())
	{
		auto pv = Q.front();
		auto p = mesh.point(pv);
		mmark[pv.idx()] = false;
		Q.pop();
		for (auto vviter = mesh.vv_begin(pv); vviter != mesh.vv_end(pv); vviter++)
		{
			auto np = mesh.point(vviter);
			double d=std::sqrt(std::pow(p[0] - np[0], 2) + std::pow(p[1] - np[1], 2) + std::pow(p[2] - np[2], 2));
			if ((dis[vviter->idx()] == -1&&(prune_threshold==-1||dis[pv.idx()]+d<prune_threshold)) || dis[vviter->idx()] > dis[pv.idx()] + d)
			{
				pre[vviter->idx()] = pv.idx();
				dis[vviter->idx()] = dis[pv.idx()] + d;
				if (vviter->idx() == svh.idx())
				{
					prune_threshold = dis[vviter->idx()];
				}
				if (mmark[vviter->idx()] == false)
				{
					Q.push(vviter);
					mmark[vviter->idx()] = true;
				}
			}

		}
	}
	if (dis[svh.idx()] == -1)
	{
		return false;
	}
	else
	{
		path.clear();
		auto pvh = svh;
		path.push_back(pvh);
		while (pre[pvh.idx()] != -1)
		{
			pvh = mesh.vertex_handle(pre[pvh.idx()]);
			path.push_back(pvh);
		}
		return true;
	}
}
void CGeoAlg::FillHoles(Eigen::MatrixXd & vertexs, Eigen::MatrixXi & faces)
{
	Polyhedron poly;
	std::vector<Face_handle>fhmap;
	std::vector<Vertex_handle>vhmap;
	CConverter::ConvertFromIGLToCGAL(vertexs, faces, poly, fhmap, vhmap);
	std::cout << "start filling holes" << std::endl;
	int hole_count = 0;
	BOOST_FOREACH(Polyhedron::Halfedge_handle h, halfedges(poly))
	{
		if (h->is_border())
		{
			hole_count++;
			std::vector<Polyhedron::Facet_handle>  patch_facets;
			std::vector<Polyhedron::Vertex_handle> patch_vertices;
			CGAL::Polygon_mesh_processing::triangulate_and_refine_hole(
				poly,
				h,
				std::back_inserter(patch_facets),
				std::back_inserter(patch_vertices),
				CGAL::Polygon_mesh_processing::parameters::vertex_point_map(get(CGAL::vertex_point, poly)).
				geom_traits(Kernel()));
		}
	}
	std::cout << hole_count << " holes filled" << std::endl;
	CConverter::ConvertFromCGALToIGL(poly, vertexs, faces);
}

//for mesh simplification
typedef boost::graph_traits<OpenMesh::TriMesh_ArrayKernelT<COMTraits>>::edge_descriptor edge_descriptor;
typedef boost::graph_traits<OpenMesh::TriMesh_ArrayKernelT<COMTraits>>::edge_iterator edge_iterator;
class Constrained_edge_map
{
public:
	typedef boost::read_write_property_map_tag    category;
	typedef bool                                  value_type;
	typedef bool                                  reference;
	typedef edge_descriptor                       key_type;
	Constrained_edge_map(OpenMesh::TriMesh_ArrayKernelT<COMTraits>& sm)
		: sm_(sm)
	{
		sm_.add_property(constraint);
	}
	inline friend reference get(const Constrained_edge_map& em, key_type e)
	{
		bool b = em.sm_.property(em.constraint, em.sm_.edge_handle(e.idx()));
		return b;
	}

	inline friend void put(const Constrained_edge_map& em, key_type e, value_type b)
	{
		em.sm_.property(em.constraint, em.sm_.edge_handle(e.idx())) = b;
	}
private:
	OpenMesh::TriMesh_ArrayKernelT<COMTraits>& sm_;
	OpenMesh::EPropHandleT<bool> constraint;
};
////////////////////////
void CGeoAlg::ComputeLocalExtremum(COpenMeshT &mesh, Eigen::VectorXd &scalars, int neighbor_num, std::vector<COpenMeshT::VertexHandle>&res_vhs)
{
	res_vhs.clear();
	for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); viter++)
	{
		double vs = scalars(viter->idx());
		std::vector<COpenMeshT::VertexHandle>neivhs;
		ExtractNRing(mesh, viter, neighbor_num, neivhs);
		bool flag = true;
		for (int i = 0; i < neivhs.size(); i++)
		{
			if (scalars(neivhs[i].idx()) > vs)
			{
				flag = false;
				break;
			}
		}
		if (flag)
		{
			res_vhs.push_back(viter);
		}
	}
}
void CGeoAlg::GetOrderedRegionBound(COpenMeshT&mesh, std::vector<int>&region_tags, COpenMeshT::VertexHandle start_vh, std::vector<COpenMeshT::VertexHandle>&res_bound_vhs)
{

}

void CGeoAlg::LaplacianSmooth(COpenMeshT &mesh, int epoch, double step)
{
	COpenMeshT tmp_mesh[2];
	tmp_mesh[0]=tmp_mesh[1]= mesh;
	int p = 0, q = 1;
	while (epoch--)
	{
		p = 1 - p;
		q = 1 - q;
		for (auto viter = tmp_mesh[p].vertices_begin(); viter != tmp_mesh[p].vertices_end(); viter++)
		{
			OpenMesh::Vec3d mean_v(0, 0, 0);
			int count = 0;
			for (auto vviter = tmp_mesh[p].vv_begin(viter); vviter != tmp_mesh[p].vv_end(viter); vviter++)
			{
				count++;
				mean_v += tmp_mesh[p].point(vviter);
			}
			mean_v = mean_v / count;
			auto qviter = tmp_mesh[q].vertex_handle(viter->idx());
			OpenMesh::Vec3d pv = tmp_mesh[q].point(qviter);
			pv = pv + (mean_v - pv)*step;
			tmp_mesh[q].set_point(qviter,pv);
			
		}

	}
	for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); viter++)
	{
		int vid = viter->idx();
		auto pv=tmp_mesh[p].point(tmp_mesh[p].vertex_handle(vid));
		mesh.set_point(viter, pv);
	}
}

void CGeoAlg::ComputeLocalExtremum(COpenMeshT &mesh, std::vector<double> &scalars, int neighbor_num, std::vector<COpenMeshT::VertexHandle>&res_vhs)
{
	res_vhs.clear();
	for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); viter++)
	{
		double vs = scalars[viter->idx()];
		std::vector<COpenMeshT::VertexHandle>neivhs;
		ExtractNRing(mesh, viter, neighbor_num, neivhs);
		bool flag = true;
		for (int i = 0; i < neivhs.size(); i++)
		{
			if (scalars[neivhs[i].idx()] > vs)
			{
				flag = false;
				break;
			}
		}
		if (flag)
		{
			res_vhs.push_back(viter);
		}
	}
}

void CGeoAlg::ExtractNRing(COpenMeshT &mesh, COpenMeshT::VertexHandle vh, int ringnum, std::vector<COpenMeshT::VertexHandle>&res_vhs)
{
	std::queue<std::pair<COpenMeshT::VertexHandle,int>>Q;
	res_vhs.clear();
	Q.push(std::make_pair(vh,0));
	while (!Q.empty())
	{
		auto p = Q.front();
		auto vh = p.first;
		Q.pop();
		if (p.second < ringnum)
		{
			for (auto vviter = mesh.vv_begin(vh); vviter != mesh.vv_end(vh); vviter++)
			{
				res_vhs.push_back(vviter);
				Q.push(std::make_pair(vviter, p.second + 1));
			}
		}
		
	}
}
void CGeoAlg::SimplifyMesh(OpenMesh::TriMesh_ArrayKernelT<COMTraits> &mesh, int edgenum)
{
	std::cerr << "simplify" << std::endl;
	namespace SMS = CGAL::Surface_mesh_simplification;
	Constrained_edge_map constraints_map(mesh);
	
	// For the pupose of the example we mark 10 edges as constrained edges
	//edge_iterator b, e;
	//int count = 0;
	//for (boost::tie(b, e) = edges(surface_mesh); b != e; ++b) {
	//	put(constraints_map, *b, (count++ < 100));
	//}
	// This is a stop predicate (defines when the algorithm terminates).
	// In this example, the simplification stops when the number of undirected edges
	// left in the surface mesh drops below the specified number (1000)
	SMS::Count_stop_predicate<OpenMesh::TriMesh_ArrayKernelT<COMTraits>> stop(edgenum);

	// This the actual call to the simplification algorithm.
	// The surface mesh and stop conditions are mandatory arguments.
	int r = SMS::edge_collapse
	(mesh
		, stop
		, CGAL::parameters::halfedge_index_map(get(CGAL::halfedge_index, mesh))
		.vertex_point_map(get(boost::vertex_point, mesh))
		.edge_is_constrained_map(constraints_map)
	);

	mesh.garbage_collection();
	std::cerr << "simplify end" << std::endl;
}



void CGeoAlg::FillHoles(COpenMeshT &mesh, bool remain_largest)
{
	Polyhedron poly;
	
	CConverter::ConvertFromOpenMeshToCGAL(mesh, poly, std::map<COpenMeshT::VertexHandle, Polyhedron::Vertex_handle>(), std::map<COpenMeshT::FaceHandle, Polyhedron::Facet_handle>());
	std::cout << "start filling holes" << std::endl;
	int hole_count = 0;
	std::vector<Polyhedron::Halfedge_handle>border_handles;
	std::vector<int>border_size;
	std::vector<bool>mmark(poly.size_of_halfedges(),0);
	BOOST_FOREACH(Polyhedron::Halfedge_handle h, halfedges(poly))
	{
		if (h->is_border() && mmark[h->id()] == false)
		{
			mmark[h->id()] = true;
			border_handles.push_back(h);
			int bsize = 1;
			auto nexth = h->next();
			while (mmark[nexth->id()] == false)
			{
				mmark[nexth->id()] = true;
				nexth = nexth->next();
				bsize++;
			}
			border_size.push_back(bsize);
		}
	}
	
	int max_size = -1;
	int mi;
	for (int i = 0; i < border_size.size(); i++)
	{
		if (max_size < border_size[i])
		{
			mi = i;
			max_size = border_size[i];
		}
	}
	
	for (int i = 0; i<border_handles.size(); i++)
	{
		if (remain_largest&&i == mi)
			continue;
	
			hole_count++;
			std::vector<Polyhedron::Facet_handle>  patch_facets;
			std::vector<Polyhedron::Vertex_handle> patch_vertices;
			CGAL::Polygon_mesh_processing::triangulate_and_refine_hole(
				poly,
				border_handles[i],
				std::back_inserter(patch_facets),
				std::back_inserter(patch_vertices),
				CGAL::Polygon_mesh_processing::parameters::vertex_point_map(get(CGAL::vertex_point, poly)).
				geom_traits(Kernel()));
		
	}
	std::cout << hole_count << " holes filled" << std::endl;
	CConverter::ConvertFromCGALToOpenMesh(poly, mesh, std::map<Polyhedron::Vertex_handle, COpenMeshT::VertexHandle>(), std::map<Polyhedron::Facet_handle, COpenMeshT::FaceHandle>());

}
void CGeoAlg::SplitFaceByVhs(COpenMeshT &mesh, COpenMeshT::FaceHandle fh, std::vector<COpenMeshT::VertexHandle>&vhs, std::vector<std::pair<COpenMeshT::VertexHandle, COpenMeshT::VertexHandle>>&cons, std::map<COpenMeshT::VertexHandle, COpenMeshT::EdgeHandle>&vh_eh_map)
{
	
	std::vector<Eigen::Vector3d>points;
	std::vector<COpenMeshT::VertexHandle>vh_map;
	std::map<COpenMeshT::VertexHandle, int>vh_id_map;

	for (auto fviter = mesh.fv_ccwbegin(fh); fviter != mesh.fv_ccwend(fh); fviter++)
	{
		auto vh = mesh.vertex_handle(fviter->idx());
		OpenMesh::Vec3d p = mesh.point(vh);
		points.push_back(Eigen::Vector3d(p[0], p[1], p[2]));
		vh_map.push_back(vh);
		vh_id_map[vh] = vh_map.size() - 1;
	}
	
	for (int i = 0; i < vhs.size(); i++)
	{
		OpenMesh::Vec3d p = mesh.point(vhs[i]);
		points.push_back(Eigen::Vector3d(p[0], p[1], p[2]));
		vh_map.push_back(vhs[i]);
		vh_id_map[vhs[i]] = vh_map.size()-1;
	}
	Eigen::MatrixXd res_v;
	Eigen::MatrixXi res_f;
	std::vector<int>orig_vidmap;
	//CGeoAlg::DelaunayTriangulation3d(points, res_v, res_f, orig_vidmap);
	
	std::vector<std::pair<int, int>>tmp_cons;
	for (int i = 0; i < cons.size(); i++)
	{
		tmp_cons.push_back(std::make_pair(vh_id_map[cons[i].first], vh_id_map[cons[i].second]));
	}
	CGeoAlg::ConstrainedTriangulation3d(points, tmp_cons, res_v, res_f, orig_vidmap);

	std::vector<int>new2orig_vidmap;
	new2orig_vidmap.resize(res_v.rows());
	for (int i = 0; i <orig_vidmap.size(); i++)
	{
		new2orig_vidmap[orig_vidmap[i]] = i;
	}
	std::set<int>ext_fs;
	for (int i = 0; i < res_f.rows(); i++)
	{
		std::vector<COpenMeshT::VertexHandle>edge_ps, orig_ps;
		for (int j = 0; j < 3; j++)
		{
			auto vh = vh_map[new2orig_vidmap[res_f(i, j)]];
			if (new2orig_vidmap[res_f(i, j)] < 3)
			{
				orig_ps.push_back(vh);
			}
			else if (vh_eh_map.find(vh) != vh_eh_map.end())
			{
				edge_ps.push_back(vh);
			}
		}
		if (edge_ps.size() + orig_ps.size() == 3)
		{
			if (edge_ps.size() == 1)
			{
				COpenMeshT::EdgeHandle common_edge;
				if (CGeoBaseAlg::GetCommonEdge(mesh, orig_ps[0], orig_ps[1], common_edge))
				{
					if (common_edge == vh_eh_map[edge_ps[0]])
					{
						ext_fs.insert(i);
					}
				}
			}
			else if (edge_ps.size() == 2)
			{
				if (vh_eh_map[edge_ps[0]] == vh_eh_map[edge_ps[1]])
				{
					COpenMeshT::EdgeHandle eh = vh_eh_map[edge_ps[0]];
					if (orig_ps[0] == mesh.to_vertex_handle(mesh.halfedge_handle(eh, 0)) || orig_ps[0] == mesh.to_vertex_handle(mesh.halfedge_handle(eh, 1)))
					{
						ext_fs.insert(i);
					}
				}
			}
			else if (edge_ps.size() == 3)
			{
				if (vh_eh_map[edge_ps[0]] == vh_eh_map[edge_ps[1]] && vh_eh_map[edge_ps[0]] == vh_eh_map[edge_ps[2]])
				{
					ext_fs.insert(i);
				}
			}
		}
	}
	Eigen::MatrixXi tmp_fs;
	int fi = 0;
	tmp_fs.resize(res_f.rows() - ext_fs.size(), 3);
	for (int i = 0; i < res_f.rows(); i++)
	{
		if (ext_fs.find(i) == ext_fs.end())
		{
			tmp_fs(fi, 0) = res_f(i, 0);
			tmp_fs(fi, 1) = res_f(i, 1);
			tmp_fs(fi, 2) = res_f(i, 2);
			fi++;
		}
	}
	res_f = tmp_fs;
	mesh.delete_face(fh, false);
	//std::cerr << "res fi size " << res_f.rows() << std::endl;
	for (int i = 0; i < res_f.rows(); i++)
	{
		COpenMeshT::VertexHandle vh0 = vh_map[new2orig_vidmap[res_f(i, 0)]];
		COpenMeshT::VertexHandle vh1 = vh_map[new2orig_vidmap[res_f(i, 1)]];
		COpenMeshT::VertexHandle vh2 = vh_map[new2orig_vidmap[res_f(i, 2)]];
		if (vh0.idx() < 0 || vh0.idx() >= mesh.n_vertices() || vh1.idx() < 0 || vh1.idx() >= mesh.n_vertices() || vh2.idx() < 0 || vh2.idx() >= mesh.n_vertices())
		{
			std::cerr << "vhs id are wrong" << std::endl;
		}


		mesh.add_face(vh0, vh1, vh2);
	}

}
void CGeoAlg::SmoothMesh(COpenMeshT &mesh, std::vector<COpenMeshT::VertexHandle>&roi_vhs, int num)
{
	std::vector<OpenMesh::Vec3d>new_pos;
	new_pos.resize(roi_vhs.size());
	while (num--)
	{
		for (int i = 0; i < roi_vhs.size(); i++)
		{
			OpenMesh::Vec3d vpos(0, 0, 0);
			int count = 0;
			for (auto vviter = mesh.vv_begin(roi_vhs[i]); vviter != mesh.vv_end(roi_vhs[i]); vviter++)
			{
				vpos += mesh.point(vviter);
				count++;
			}
			vpos /= count;
			new_pos[i] = vpos;
		}
		for (int i = 0; i < roi_vhs.size(); i++)
		{
			mesh.point(roi_vhs[i]) = new_pos[i];
		}
	}
	
}
void CGeoAlg::SeparateMeshByFaceTag(COpenMeshT &mesh, std::vector<int>&f_tag, std::map<int, COpenMeshT*>&res_meshes, std::map<int, std::map<COpenMeshT::VertexHandle, COpenMeshT::VertexHandle>>&vnew_orig)
{
	res_meshes.clear();
	vnew_orig.clear();
	std::map<int, std::set<COpenMeshT::VertexHandle>>t_vset;
	std::map<int, std::map<COpenMeshT::VertexHandle, COpenMeshT::VertexHandle>>orig_new_vhs_map;
	for (auto fiter = mesh.faces_begin(); fiter != mesh.faces_end(); fiter++)
	{
		int t = f_tag[fiter->idx()];
		if (res_meshes.find(t) == res_meshes.end())
		{
			res_meshes.insert(std::make_pair(t, new COpenMeshT()));
			res_meshes[t]->clear();
			t_vset[t]=std::set<COpenMeshT::VertexHandle>();
			orig_new_vhs_map[t] = std::map<COpenMeshT::VertexHandle, COpenMeshT::VertexHandle>();

			vnew_orig[t]=std::map<COpenMeshT::VertexHandle, COpenMeshT::VertexHandle>();
		}
		std::vector<COpenMeshT::VertexHandle>fvhs;
		for (auto fviter = mesh.fv_begin(fiter); fviter != mesh.fv_end(fiter); fviter++)
		{
		
			if (t_vset[t].find(fviter) == t_vset[t].end())
			{
				t_vset[t].insert(fviter);
				COpenMeshT::VertexHandle new_vh=res_meshes[t]->add_vertex(mesh.point(fviter));
				vnew_orig[t][new_vh] = fviter;
				orig_new_vhs_map[t][fviter] = new_vh;
			}
			fvhs.push_back(orig_new_vhs_map[t][fviter]);
		}
		res_meshes[t]->add_face(fvhs);


	}
	for (auto iter = res_meshes.begin(); iter != res_meshes.end(); iter++)
	{
		iter->second->garbage_collection();
	}
}
void CGeoAlg::SeparateMeshByVertexTag(COpenMeshT &mesh, std::vector<int>&v_tag, std::map<int,COpenMeshT*>&res_meshes, std::map<int,std::map<COpenMeshT::VertexHandle, COpenMeshT::VertexHandle>>&v_orig)
{
	
	res_meshes.clear();
	v_orig.clear();

	
	std::vector<COpenMeshT::VertexHandle>tmp_v_new(mesh.n_vertices());
	for (auto v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it)
	{
		int t = v_tag[v_it->idx()];
		if (res_meshes.find(t) == res_meshes.end())
		{
			res_meshes.insert(std::make_pair(t,new COpenMeshT()));
			res_meshes[t]->clear();
		}
			
		auto new_vh=res_meshes[t]->add_vertex(mesh.point(v_it));
		res_meshes[t]->set_color(new_vh, mesh.color(v_it));
		if (v_orig.find(t) == v_orig.end())
			v_orig[t] = std::map<COpenMeshT::VertexHandle, COpenMeshT::VertexHandle>();
		v_orig[t][new_vh] = *v_it;
		tmp_v_new[v_it->idx()]=new_vh;
	}
	for (auto f_it = mesh.faces_begin(); f_it != mesh.faces_end(); ++f_it)
	{
		std::vector<COpenMeshT::VertexHandle> vhandles;
		int t = -std::numeric_limits<int>::max();
		auto fh = mesh.fh_begin(f_it);
		for (; fh != mesh.fh_end(f_it); fh++)
		{
			int vid = mesh.to_vertex_handle(fh).idx();
			int pt=v_tag[vid];
			if (t != -std::numeric_limits<int>::max() && pt != t)
			{
				break;
			}
			t = pt;
			vhandles.push_back(tmp_v_new[vid]);
		}
		if (fh == mesh.fh_end(f_it))
		{
			res_meshes[t]->add_face(vhandles);
		}
		
	}
	for (auto iter=res_meshes.begin();iter!=res_meshes.end();iter++)
	{
		iter->second->garbage_collection();
	}
	
}
void CGeoAlg::CutByPlane(COpenMeshT &in_mesh, CPlane plane, COpenMeshT &res_mesh, std::map<COpenMeshT::VertexHandle, COpenMeshT::VertexHandle>&vid_orig)
{
	COpenMeshT mesh = in_mesh;
	
	res_mesh.clear();
	int t = 1;
	std::vector<int>v_tags(mesh.n_vertices(),0);
	for (auto v_iter = mesh.vertices_begin(); v_iter != mesh.vertices_end(); v_iter++)
	{
		COpenMeshT::Point p = mesh.point(v_iter);
		if(CGeoBaseAlg::IsOnPositiveSide(p, plane))
		v_tags[v_iter->idx()] = t;
	}
	std::map<int,COpenMeshT*>tmp_res_meshes;
	std::map<int,std::map<COpenMeshT::VertexHandle, COpenMeshT::VertexHandle>>tmpv_orig;
	SeparateMeshByVertexTag(mesh, v_tags, tmp_res_meshes, tmpv_orig);
	res_mesh = *tmp_res_meshes[t];
	
	delete tmp_res_meshes[0];
	vid_orig = tmpv_orig[t];
}
bool CGeoAlg::SelfIntersectionRemoval(COpenMeshT& mesh)
{
	Polyhedron poly;
	CConverter::ConvertFromOpenMeshToCGAL(mesh, poly, std::map<COpenMeshT::VertexHandle, Polyhedron::Vertex_handle>(), std::map<COpenMeshT::FaceHandle, Polyhedron::Facet_handle>());
	Polyhedron*p_poly = &poly;
	namespace PMP = CGAL::Polygon_mesh_processing;
	while (PMP::does_self_intersect(*p_poly,
		PMP::parameters::vertex_point_map(CGAL::get(CGAL::vertex_point, *p_poly))))
	{
		printf("Self-intersection detected!\n");
		std::vector<std::pair<Polyhedron::Facet_handle, Polyhedron::Facet_handle> > intersected_tris;
		PMP::self_intersections(*p_poly,
			std::back_inserter(intersected_tris),
			PMP::parameters::vertex_point_map(CGAL::get(CGAL::vertex_point, *p_poly)));
		//put the facet pair to a set (expand one ring)
		std::set<Polyhedron::Facet_handle> inter_facet_set;
		for (int i = 0; i < intersected_tris.size(); i++)
		{
			Polyhedron::Halfedge_around_vertex_circulator heiter = intersected_tris[i].first->facet_begin()->vertex_begin();
			Polyhedron::Halfedge_around_vertex_circulator heend = heiter;
			CGAL_For_all(heiter, heend)
			{
				inter_facet_set.insert(heiter->facet());
			}
			heiter = intersected_tris[i].second->facet_begin()->vertex_begin();
			heend = heiter;
			CGAL_For_all(heiter, heend)
			{
				inter_facet_set.insert(heiter->facet());
			}
		}
		int solution = 0;
		switch (solution)
		{
		case 0:
		{
			//Delete intersection faces
			if (inter_facet_set.find(NULL) != inter_facet_set.end())
			{
				inter_facet_set.erase(inter_facet_set.find(NULL));
			}
			for (auto p = inter_facet_set.begin(); p != inter_facet_set.end(); p++)
			{
				p_poly->erase_facet((*p)->facet_begin());
			}
			printf("Facets Deleted\n");
			p_poly->keep_largest_connected_components(1);

		}
		break;
		case 1:
		{
			std::set<Polyhedron::Vertex_handle> inter_ver_set;
			for (auto p = inter_facet_set.begin(); p != inter_facet_set.end(); p++)
			{
				Polyhedron::Halfedge_around_facet_circulator heiter = (*p)->facet_begin();
				Polyhedron::Halfedge_around_facet_circulator heend = heiter;
				CGAL_For_all(heiter, heend)
				{
					inter_ver_set.insert(heiter->vertex());
				}
			}
			std::vector<Polyhedron::Vertex_handle> inter_ver_vec;
			for (auto p = inter_ver_set.begin(); p != inter_ver_set.end(); p++)
			{
				inter_ver_vec.push_back(*p);
			}
			typedef CGAL::Eigen_solver_traits<Eigen::SparseLU<
				CGAL::Eigen_sparse_matrix<double>::EigenType, Eigen::COLAMDOrdering<int> >  >
				Default_solver;
			typedef CGAL::internal::Cotangent_weight_with_voronoi_area_fairing<Polyhedron>
				Default_Weight_calculator;
			if (!PMP::internal::fair(*p_poly, inter_ver_vec,
				Default_solver(),
				Default_Weight_calculator(*p_poly), 1,
				CGAL::get(CGAL::vertex_point, *p_poly)))
			{
				printf("Fairing failed!\n");
				break;
			}
		}
		break;
		default:
			break;
		}
	}

	CConverter::ConvertFromCGALToOpenMesh(poly, mesh, std::map<Polyhedron::Vertex_handle, COpenMeshT::VertexHandle>(), std::map<Polyhedron::Facet_handle, COpenMeshT::FaceHandle>());
	return true;
}
double CGeoAlg::ComputeClosestPoint(OpenMesh::Vec3d source, CMeshObject &mesh_obj, OpenMesh::FaceHandle &res_fh, OpenMesh::Vec3d& res_p)
{
	Eigen::Matrix4d invmat = mesh_obj.GetMatrix().inverse();
	CObbWrapper  *obb = mesh_obj.GetPqpObb();

	if (obb == NULL)
	{
		obb = new CObbWrapper();
		Eigen::MatrixXd V;
		Eigen::MatrixXi F;
		
		CConverter::ConvertFromOpenMeshToIGL(mesh_obj.GetMesh(), V, F, &obb->id_fh_map_,&obb->id_vh_map_);
		obb->obb_tree_ = new CPqpObb(V,F);
		mesh_obj.SetPqpObb(obb);
	}
	Eigen::Vector4d ep(source[0], source[1], source[2],1);
	ep = invmat*ep;
	CPqpObb::QueryResult res=obb->obb_tree_->QueryClosestPoint(Eigen::Vector3d(ep.x(),ep.y(),ep.z()));
	res_fh = obb->id_fh_map_[res.fid_];
	Eigen::Vector4d eres_p= mesh_obj.GetMatrix()*Eigen::Vector4d(res.closest_pnt_.x(), res.closest_pnt_.y(), res.closest_pnt_.z(),1);
	res_p = OpenMesh::Vec3d(eres_p.x(), eres_p.y(), eres_p.z());
	return res.distance;
	
}
void CGeoAlg::SeparateVhsByConnectivity(std::vector<COpenMeshT::VertexHandle>&vhs, COpenMeshT &mesh, std::vector<std::vector<COpenMeshT::VertexHandle>>&res_vhs)
{
	std::set<COpenMeshT::VertexHandle>vhs_set;
	for (int i = 0; i < vhs.size(); i++)
	{
		vhs_set.insert(vhs[i]);
	}
	res_vhs.clear();
	for (int i = 0; i < vhs.size(); i++)
	{
		if (vhs_set.find(vhs[i]) == vhs_set.end())
			continue;
		std::queue<COpenMeshT::VertexHandle>Q;
		Q.push(vhs[i]);
		res_vhs.push_back(std::vector<COpenMeshT::VertexHandle>());
		res_vhs.back().push_back(vhs[i]);
		while (!Q.empty())
		{
			COpenMeshT::VertexHandle p = Q.front();
			Q.pop();
			for (auto vviter = mesh.vv_begin(p); vviter != mesh.vv_end(p); vviter++)
			{
				if (vhs_set.find(vviter) == vhs_set.end())
				{
					continue;
				}
				else
				{
					Q.push(vviter);
					res_vhs.back().push_back(vviter);
					vhs_set.erase(vviter);
				}
			}
		}
	}
}
void CGeoAlg::ComputeClosestVertex(OpenMesh::Vec3d source, CMeshObject &mesh_obj, OpenMesh::VertexHandle &res_vh)
{
	COpenMeshT &mesh = mesh_obj.GetMesh();
	OpenMesh::VertexHandle min_vhs;
	double min_dis = std::numeric_limits<double>::max();
	for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); viter++)
	{
		OpenMesh::Vec3d p = mesh.point(viter);
		double dis = (p - source).length();
		if (dis < min_dis)
		{
			min_dis = dis;
			min_vhs = viter;
		}
	}
	res_vh = min_vhs;
}
bool CGeoAlg::RayMeshIntersection(OpenMesh::Vec3d  source, OpenMesh::Vec3d dir, CMeshObject &mesh_obj, COpenMeshT::VertexHandle & res_vh)
{
	COpenMeshT::FaceHandle  res_fh;
	OpenMesh::Vec3d res_bary_coord;
	if (!RayMeshIntersection(source, dir, mesh_obj, res_fh, res_bary_coord))
		return false;
	auto p = CGeoBaseAlg::ComputePointFromBaryCoord(mesh_obj.GetMesh(), res_fh, res_bary_coord);
	COpenMeshT &mesh = mesh_obj.GetMesh();

	double mind = std::numeric_limits<double>::max();
	for (auto fviter = mesh.fv_begin(res_fh); fviter != mesh.fv_end(res_fh); fviter++)
	{
		auto fp = mesh.point(fviter);
		double d = (fp - p).length();
		if (d < mind)
		{
			mind = d;
			res_vh = fviter;
		}
	}
	return true;
}
bool CGeoAlg::RayMeshIntersection(OpenMesh::Vec3d  source, OpenMesh::Vec3d dir, CMeshObject &mesh_obj, COpenMeshT::FaceHandle & res_fh, OpenMesh::Vec3d &res_bary_coord)
{
	
	igl::Hit hit;
	Eigen::Matrix4d invmat = mesh_obj.GetMatrix().inverse();
	Eigen::Vector4d lsource= Eigen::Vector4d(source[0], source[1], source[2], 1);
	Eigen::Vector4d ldir = Eigen::Vector4d(dir[0], dir[1], dir[2], 1);
	COpenMeshT &mesh = mesh_obj.GetMesh();
	
	ldir = invmat*ldir;
	lsource = invmat*lsource;
	
	
	
	auto aabb_tree = mesh_obj.GetAABBTree();

	if (aabb_tree == NULL)
	{
		
		std::map<COpenMeshT::VertexHandle, Polyhedron::Vertex_handle> vh_map;
		std::map<COpenMeshT::FaceHandle, Polyhedron::Facet_handle>fh_map;
		
		aabb_tree = new CAABBTree();
		CConverter::ConvertFromOpenMeshToCGAL(mesh, aabb_tree->poly, vh_map, fh_map);
		aabb_tree->tree_=new AABBTree(faces(aabb_tree->poly).first, faces(aabb_tree->poly).second, aabb_tree->poly);
		aabb_tree->id_fh_map_.resize(mesh.n_faces());
		aabb_tree->id_vh_map_.resize(mesh.n_vertices());
		for (auto iter = fh_map.begin(); iter != fh_map.end(); iter++)
		{
			aabb_tree->id_fh_map_[iter->second->id()] = iter->first;
		}
		for (auto iter = vh_map.begin(); iter != vh_map.end(); iter++)
		{
			aabb_tree->id_vh_map_[iter->second->id()] = iter->first;
		}
		mesh_obj.SetAABBTree(aabb_tree);
	}

	
	Point to(lsource[0]+ldir[0],lsource[1]+ldir[1],lsource[2]+ldir[2]);
	Point orig(lsource[0], lsource[1], lsource[2]);
	
	Ray ray_query(orig, to);

	std::list<Ray_intersection> intersections;
	
	try
	{
		aabb_tree->tree_->all_intersections(ray_query, std::back_inserter(intersections));
	}
	catch (void* e)
	{
		std::cerr<<"intersection error"<<std::endl;
	}



	std::list<Polyhedron::Face_handle> face_handles;
	aabb_tree->tree_->all_intersected_primitives(ray_query, std::back_inserter(face_handles));
	assert(intersections.size() == face_handles.size());

	if (intersections.size() == 0)
		return false;
	auto p_iter = intersections.begin();
	auto fh_iter = face_handles.begin();
	double min_dis = 9999999999999999999;
	OpenMesh::Vec3d res_p;
	while (p_iter != intersections.end())
	{
		if (fh_iter == face_handles.end())
		{
			printf("error in intersection : fhiter==face_handles\n");
		}
		assert(fh_iter != face_handles.end());
		Point* cur_p = boost::get<Point>(&(p_iter->get().first));
		Point tmp_p = *cur_p;
		double tmp_dis = CGAL::squared_distance(orig, tmp_p);
		if (tmp_dis < min_dis)
		{
			min_dis = tmp_dis;
			res_fh = aabb_tree->id_fh_map_[(*fh_iter)->id()];
			res_p = OpenMesh::Vec3d(tmp_p.x(), tmp_p.y(), tmp_p.z());
			res_bary_coord=CGeoBaseAlg::ComputeBaryCoordInTriFace(mesh, res_fh, res_p);
		}

		fh_iter++;
		p_iter++;
	}
	//std::cerr << "intersected" << std::endl;
	return true;

}
void CGeoAlg::ComputeMeshPCA(COpenMeshT&mesh, OpenMesh::Vec3d&mean, std::vector<OpenMesh::Vec3d>&res_frame)
{
	std::vector<OpenMesh::Vec3d>points;
	for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); viter++)
	{
		points.push_back(mesh.point(viter));
	}
	std::vector<double>res_eigenvalues;
	CGeoAlg::PointSetPCA3D(points, mean, res_frame, res_eigenvalues);
	/*for (int i = 0; i < res_frame.size(); i++)
	{
		res_frame[i]=res_frame[i].normalize()*res_eigenvalues[i];
	}*/

}
void CGeoAlg::GetFeatureGroupsByConnectivity(COpenMeshT&mesh, std::vector<bool>&is_feature, std::vector<std::vector<COpenMeshT::VertexHandle>>&feature_groups)
{
	std::vector<bool>mmark(mesh.n_vertices(), 0);
	int gid = 0;
	feature_groups.clear();
	for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); viter++)
	{
		if (is_feature[viter->idx()]&&mmark[viter->idx()]==false)
		{
			feature_groups.push_back(std::vector<COpenMeshT::VertexHandle>());
			
			mmark[viter->idx()] = true;
			std::queue<COpenMeshT::VertexHandle>Q;
			Q.push(viter);
			feature_groups[gid].push_back(viter);
			while (!Q.empty())
			{
				auto pviter = Q.front();
				Q.pop();
				for (auto nviter = mesh.vv_begin(pviter); nviter != mesh.vv_end(pviter); nviter++)
				{
					if (is_feature[nviter->idx()] && mmark[nviter->idx()] == false)
					{
						Q.push(nviter);
						mmark[nviter->idx()] = true;
						feature_groups[gid].push_back(nviter);
					}
				}
			}
			gid++;
		}
	}
}
void CGeoAlg::SeparateDisconnectedParts(COpenMeshT &mesh, std::vector<COpenMeshT*>&res_meshes, std::vector<std::map<COpenMeshT::VertexHandle, COpenMeshT::VertexHandle>>&vid_orig)
{
	std::vector<int>mmark(mesh.n_vertices(),-1);
	int t = 0;
	for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); viter++)
	{
		
		int vid = viter->idx();
		if (mmark[vid] != -1)
			continue;
		std::queue<COpenMeshT::VertexHandle>Q;
		Q.push(viter);
		mmark[vid] = t;
		while (!Q.empty())
		{
			auto pviter = Q.front();

			Q.pop();
			for (auto vviter = mesh.vv_begin(pviter); vviter != mesh.vv_end(pviter); vviter++)
			{
				if (mmark[vviter->idx()] == -1)
				{
					Q.push(vviter);
					mmark[vviter->idx()] = t;
				}
				
			}
		}
		t++;
	}

	std::map<int, COpenMeshT*>tmp_meshes;
	std::map<int, std::map<COpenMeshT::VertexHandle, COpenMeshT::VertexHandle>>tmp_vid_orig;
	SeparateMeshByVertexTag(mesh, mmark, tmp_meshes, tmp_vid_orig);
	res_meshes.clear();
	for (auto iter = tmp_meshes.begin(); iter != tmp_meshes.end(); iter++)
	{
		res_meshes.push_back(iter->second);
	}
	for (auto iter = tmp_vid_orig.begin(); iter != tmp_vid_orig.end(); iter++)
	{
		vid_orig.push_back(iter->second);
	}
	


}
