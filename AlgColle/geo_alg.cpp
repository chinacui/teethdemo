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
		
	}

}
bool CGeoAlg::ComputeGeodesicPath(CMeshObject &mesh_obj, int svid, int tvid, std::vector<COpenMeshT::FaceHandle>&res_fhs, std::vector<OpenMesh::Vec3d>&res_bary_coords)
{
	auto geodesic_model = mesh_obj.GetGeodesicModel();
	COpenMeshT &mesh = mesh_obj.GetMesh();
	if (geodesic_model == NULL)
	{
		geodesic_model = new CXWGeodesic();
		std::vector<CPoint3D>vertexs(mesh.n_vertices());
		std::vector<CBaseModel::CFace>faces(mesh.n_faces());
		for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); viter++)
		{
			auto p = mesh.point(viter);
			vertexs[viter->idx()] = CPoint3D(p[0], p[1], p[2]);
		}
		for (auto fiter = mesh.faces_begin(); fiter != mesh.faces_end(); fiter++)
		{
			CBaseModel::CFace f;
			int i = 0;
		
			for (auto fviter = mesh.fv_begin(fiter); i < 3; fviter++, i++)
			{
				f.verts[i] = fviter->idx();
			}
			faces[fiter->idx()] = f;
		}
		geodesic_model->SetModel(vertexs, faces);
		std::cerr << "init model finished" << std::endl;
	}
	std::vector<CGeoFacePoint>tpath;
	geodesic_model->GeodesicPath(svid, tvid, tpath);
	res_fhs.resize(tpath.size());
	res_bary_coords.resize(tpath.size());
	for (int i = 0; i < tpath.size(); i++)
	{
		res_fhs[i] = mesh.face_handle(tpath[i].fid_);
		res_bary_coords[i] = OpenMesh::Vec3d(tpath[i].ls_[0], tpath[i].ls_[1], tpath[i].ls_[2]);
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
		if (i == mi)
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
void CGeoAlg::SeparateMeshByVertexTag(COpenMeshT &mesh, std::vector<int>&v_tag, std::vector<COpenMeshT*>&res_meshes, std::vector<std::map<COpenMeshT::VertexHandle, COpenMeshT::VertexHandle>>&v_orig)
{
	int max_tag = -1;

	for (int i = 0; i < v_tag.size(); i++)
	{
		if (max_tag < v_tag[i])
		{
			max_tag = v_tag[i];
			
		}
	}
	res_meshes.resize(max_tag + 1);
	v_orig.resize(max_tag + 1);
	for (int i = 0; i <=max_tag; i++)
	{
		res_meshes[i] = new COpenMeshT();
		res_meshes[i]->clear();
		v_orig[i].clear();
	}
	std::vector<COpenMeshT::VertexHandle>tmp_v_new(mesh.n_vertices());
	for (auto v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it)
	{
		int t = v_tag[v_it->idx()];
		auto new_vh=res_meshes[t]->add_vertex(mesh.point(v_it));
		res_meshes[t]->set_color(new_vh, mesh.color(v_it));
		v_orig[t][new_vh] = *v_it;
		tmp_v_new[v_it->idx()]=new_vh;
	}
	for (auto f_it = mesh.faces_begin(); f_it != mesh.faces_end(); ++f_it)
	{
		std::vector<COpenMeshT::VertexHandle> vhandles;
		int t = -1;
		auto fh = mesh.fh_begin(f_it);
		for (; fh != mesh.fh_end(f_it); fh++)
		{
			int vid = mesh.to_vertex_handle(fh).idx();
			int pt=v_tag[vid];
			if (t != -1 && pt != t)
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
	std::vector<COpenMeshT*>tmp_res_meshes;
	std::vector<std::map<COpenMeshT::VertexHandle, COpenMeshT::VertexHandle>>tmpv_orig;
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
	//std::cerr << "intersection start" << std::endl;
	igl::Hit hit;
	auto invmat = mesh_obj.GetMatrix().inverse();
	auto lsource= invmat*Eigen::Vector4d(source[0], source[1], source[2], 1);
	auto ldir = invmat*Eigen::Vector4d(dir[0], dir[1], dir[2], 1);
	COpenMeshT &mesh = mesh_obj.GetMesh();
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

	
	SeparateMeshByVertexTag(mesh, mmark, res_meshes, vid_orig);
	


}
