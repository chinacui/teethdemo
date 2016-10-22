#include"geo_alg.h"

#include"../DataColle/cgal_igl_converter.h"
#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <opencv2/opencv.hpp>
#include"geo_base_alg.h"
#include<queue>
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
