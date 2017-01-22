#include"geo_alg.h"

#include"../DataColle/cgal_igl_converter.h"
#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <opencv2/opencv.hpp>
#include"geo_base_alg.h"
#include<queue>
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

void CGeoAlg::SeparateMeshByVertexTag(COpenMeshT &mesh, std::vector<int>&v_tag, std::map<int,COpenMeshT*>&res_meshes, std::map<int,std::map<COpenMeshT::VertexHandle, COpenMeshT::VertexHandle>>&v_orig)
{
	int max_tag = -1;

	for (int i = 0; i < v_tag.size(); i++)
	{
		if (max_tag < v_tag[i])
		{
			max_tag = v_tag[i];
			
		}
	}
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
