#include"geo_base_alg.h"
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/invert_diag.h>
#include<igl/per_vertex_normals.h>
#include"../DataColle/Polyhedron_type.h"
#include"../DataColle/cgal_igl_converter.h"
#include<queue>
#include<igl/writeOBJ.h>
#include<igl/writeDMAT.h>
#include"geo_alg.h"
void CGeoBaseAlg::ComputeMeanCurvatureValue(Eigen::MatrixXd& vertexs, Eigen::MatrixXi& faces, Eigen::VectorXd &res_curvature, bool is_signed)
{
	
	Eigen::MatrixXd HN;
	Eigen::SparseMatrix<double> L, M, Minv;
	igl::cotmatrix(vertexs, faces, L);
	igl::massmatrix(vertexs, faces, igl::MASSMATRIX_TYPE_VORONOI, M);
	igl::invert_diag(M, Minv);
	HN = -Minv*(L*vertexs);
	res_curvature = HN.rowwise().norm(); //up to sign
	igl::writeDMAT("curvature.dmat", res_curvature);
	if (is_signed)
	{
		Eigen::MatrixXd v_normal;
		ComputePerVertexNormal(vertexs, faces, v_normal);
		for (int i = 0; i < vertexs.rows(); i++)
		{
			if (v_normal(i, 0)*HN(i, 0) + v_normal(i, 1)*HN(i, 1) + v_normal(i, 2)*HN(i, 2)<0)
			{
				res_curvature(i) = -res_curvature(i);
				//std::cerr << res_curvature(i) << std::endl;
			}
		}
	}

}
OpenMesh::Vec3d CGeoBaseAlg::GetFacePointFromBaryCoord(COpenMeshT&mesh, COpenMeshT::FaceHandle fh, OpenMesh::Vec3d barycoord)
{
	OpenMesh::Vec3d resv(0, 0, 0);
	int count = 0;
	for (auto viter = mesh.fv_begin(fh); viter != mesh.fv_end(fh),count<3; viter++,count++)
	{

		resv = resv + mesh.point(viter)*barycoord[count];
	}
	//std::cerr << barycoord << std::endl;
	//std::cerr << std::endl;
	return resv ;
}
bool CGeoBaseAlg::ConvertFromOpenMeshROIToOpenMesh(COpenMeshT &mesh, std::vector<COpenMeshT::FaceHandle>&roifaces, COpenMeshT&res_mesh, std::map<COpenMeshT::FaceHandle, COpenMeshT::FaceHandle>*new_mesh_2_old_vmap)
{
	res_mesh.clear();
	if (new_mesh_2_old_vmap != NULL)
	{
		new_mesh_2_old_vmap->clear();
	}
	std::vector<bool>old_vmask(mesh.n_vertices(), false);
	std::map<COpenMeshT::VertexHandle, COpenMeshT::VertexHandle>old2new_vmap;
	for (int i = 0; i < roifaces.size(); i++)
	{
		for (auto viter = mesh.fv_begin(roifaces[i]); viter != mesh.fv_end(roifaces[i]); viter++)
		{
			if (old_vmask[viter->idx()] == false)
			{
				old2new_vmap[viter] = res_mesh.add_vertex(mesh.point(viter));
				old_vmask[viter->idx()] = true;
			}

		}
	}
	for (int i = 0; i < roifaces.size(); i++)
	{
		std::vector<COpenMeshT::VertexHandle>fvhs;
		for (auto viter = mesh.fv_begin(roifaces[i]); viter != mesh.fv_end(roifaces[i]); viter++)
		{
			fvhs.push_back(old2new_vmap[viter]);
		}
		COpenMeshT::FaceHandle fh = res_mesh.add_face(fvhs);
		if (new_mesh_2_old_vmap != NULL)
		{
			(*new_mesh_2_old_vmap)[fh] = roifaces[i];
		}
	}
	return true;
}
void CGeoBaseAlg::ComputeMeanCurvatureValue(COpenMeshT &mesh, Eigen::VectorXd &res_curvature, bool is_signed)
{
	Eigen::MatrixXd vertexs;
	Eigen::MatrixXi faces;
	CConverter::ConvertFromOpenMeshToIGL(mesh, vertexs, faces);
	//igl::writeOBJ("test.obj", vertexs, faces);
	ComputeMeanCurvatureValue(vertexs, faces, res_curvature, is_signed);
}
OpenMesh::Vec3d CGeoBaseAlg::ComputeBaryCoordInTriFace(COpenMeshT &mesh, COpenMeshT::FaceHandle fh, OpenMesh::Vec3d p)
{
	std::vector<OpenMesh::Vec3d>fps;
	for (auto fviter = mesh.fv_begin(fh); fviter != mesh.fv_end(fh); fviter++)
	{
		fps.push_back(mesh.point(fviter));
	}
	if (fps.size() != 3)
	{
		std::cerr << "error : face v num must == 3" << std::endl;
	}
	
	OpenMesh::Vec3d p0 = p - fps[0];
	OpenMesh::Vec3d p1 = p - fps[1];
	OpenMesh::Vec3d p2 = p - fps[2];
	OpenMesh::Vec3d v01 = fps[1] - fps[0];
	OpenMesh::Vec3d v02 = fps[2] - fps[0];

	double area = OpenMesh::cross(v01, v02).length() / 2;
	double area0 = OpenMesh::cross(p1, p2).length() / 2;
	double area1 = OpenMesh::cross(p0, p2).length() / 2;

	OpenMesh::Vec3d res;
	res[0] = area0 / area;
	res[1] = area1 / area;
	res[2] = 1 - res[0] - res[1];
	return res;
}
double CGeoBaseAlg::ComputeVertexArea(COpenMeshT &mesh, COpenMeshT::VertexHandle vh)
{
	double area = 0;
	for (auto vfiter = mesh.vf_begin(vh); vfiter != mesh.vf_end(vh); vfiter++)
	{
		area+=(ComputeFaceArea(mesh, vfiter)/3.0);
	}
	return area;
}
double CGeoBaseAlg::ComputeFaceArea(COpenMeshT&mesh, COpenMeshT::FaceHandle fh)
{
	std::vector<OpenMesh::Vec3d>fvs;
	for (auto fviter = mesh.fv_begin(fh); fviter != mesh.fv_end(fh); fviter++)
	{
		fvs.push_back(mesh.point(fviter));
	}
	OpenMesh::Vec3d va = fvs[1] - fvs[0];
	OpenMesh::Vec3d vb = fvs[2] - fvs[0];
	return OpenMesh::cross(va, vb).length() / 2;
}
void CGeoBaseAlg::GetFaceVhs(COpenMeshT&mesh, COpenMeshT::FaceHandle fh, std::vector<COpenMeshT::VertexHandle>&resvhs)
{
	resvhs.clear();
	for (auto viter = mesh.fv_begin(fh); viter != mesh.fv_end(fh); viter++)
	{
		resvhs.push_back(viter);
	}
}

OpenMesh::VertexHandle CGeoBaseAlg::GetClosestVhFromFacePoint(COpenMeshT&mesh, COpenMeshT::FaceHandle fh, OpenMesh::Vec3d barycoord)
{
	int maxi;
	if (barycoord[0] >= barycoord[1] && barycoord[0] >= barycoord[2])
		maxi = 0;
	else if (barycoord[1] >= barycoord[2])
		maxi = 1;
	else
		maxi = 2;
	std::vector<COpenMeshT::VertexHandle>fvhs;
	CGeoBaseAlg::GetFaceVhs(mesh, fh, fvhs);
	COpenMeshT::VertexHandle tvh = fvhs[maxi];
	return tvh;
}

void CGeoBaseAlg::GetOrderedBoundary(COpenMeshT &mesh, std::vector<std::vector<COpenMeshT::VertexHandle>>&bounds)
{
	std::vector<bool>mmark(mesh.n_vertices(), 0);
	bounds.clear();
	int bid = -1;
	for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); viter++)
	{
	
		if (mesh.is_boundary(viter)&&mmark[viter->idx()]==false)
		{
			
			bounds.push_back(std::vector<COpenMeshT::VertexHandle>());
			bid = bounds.size() - 1;
			bounds[bid].clear();
			bounds[bid].push_back(viter);
			mmark[viter->idx()] = true;
			COpenMeshT::VertexHandle vh = viter;
			COpenMeshT::VertexHandle pre_vh = vh;
			int first_vid = viter->idx();
			do {
				bool flag = false;
				for (auto nviter = mesh.vv_begin(vh); nviter != mesh.vv_end(vh); nviter++)
				{
					if (mesh.is_boundary(nviter)&&(mmark[nviter->idx()]==false||(mesh.is_manifold(nviter)==false&&pre_vh!=nviter)))
					{
						//if (nviter->idx() < 0 || nviter->idx() >= mmark.size())
						//	std::cerr<<"nviter->idx() " << nviter->idx() << std::endl;
						mmark[nviter->idx()] = true;
						//if (bid<0 || bid>=bounds.size())
						//	std::cerr << "bid " << bid << std::endl;
						bounds[bid].push_back(nviter);

						pre_vh = vh;
						vh = nviter;
						flag = true;
						break;
					}
				}
				if (!flag)
				{
					break;
				}
				
			} while (vh.idx() != first_vid);
		}
	}
	//std::cerr << "get bound end" << std::endl;
}
OpenMesh::Vec3d CGeoBaseAlg::ComputePointFromBaryCoord(COpenMeshT &mesh, COpenMeshT::FaceHandle fh, OpenMesh::Vec3d barycoord)
{
	OpenMesh::Vec3d res(0, 0, 0);
	int i = 0;
	for (auto fviter=mesh.fv_begin(fh);fviter!=mesh.fv_end(fh);fviter++,i++)
	{
		if (i >= 3)
		{
			std::cerr << "error face v num must be 3" << std::endl;
		}
		else
		{
			res = res + barycoord[i] * mesh.point(fviter);
		}
		
	}
	return res;
}
OpenMesh::Vec3d CGeoBaseAlg::ComputeMeshCenter(COpenMeshT &mesh)
{
	OpenMesh::Vec3d center(0, 0, 0);
	for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); viter++)
	{
		OpenMesh::Vec3d p=mesh.point(viter);
		center = center + p;
	}
	center = center / mesh.n_vertices();
	return center;
}
void CGeoBaseAlg::NormalizeMeshSize(COpenMeshT &mesh)
{
	OpenMesh::Vec3d mmin(std::numeric_limits<double>::max(), std::numeric_limits<double>::max(), std::numeric_limits<double>::max()), mmax(std::numeric_limits<double>::min(), std::numeric_limits<double>::min(), std::numeric_limits<double>::min());
	for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); viter++)
	{
		COpenMeshT::Point p = mesh.point(viter);
		for (int i = 0; i < 3; i++)
		{
			if (mmin[i] > p[i])
				mmin[i] = p[i];
			if (mmax[i] < p[i])
				mmax[i] = p[i];
		}
			
	}
	OpenMesh::Vec3d lens;
	lens = mmax - mmin;
	int mi = 0;
	for (int i = 1; i < 3; i++)
	{
		if (lens[i] > lens[mi])
		{
			mi = i;
		}
	}
	for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); viter++)
	{
		COpenMeshT::Point p = mesh.point(viter);
		p = (p - mmin-(lens/2.0)) /( lens[mi]/2.0)*0.6;
		//p = p - COpenMeshT::Point(lens[0]/lens[mi] / 4.0, lens[1]/lens[mi] / 4.0, lens[2]/lens[mi]/2.0);
		mesh.set_point(viter, p);
	}
	//std::cerr << "mi " << mi <<" "<< lens[1] / lens[mi]<< std::endl;
}
bool CGeoBaseAlg::IsConcavity(COpenMeshT &mesh, COpenMeshT::VertexHandle vh)
{
	OpenMesh::Vec3d vnorm = mesh.normal(vh);
	OpenMesh::Vec3d neicenter(0, 0, 0);
	int c = 0;
	for(auto vviter=mesh.vv_begin(vh);vviter!=mesh.vv_end(vh);vviter++)
	{
		c++;
		neicenter = neicenter + mesh.point(vviter);
	}
	neicenter /= c;
	auto dir = neicenter - mesh.point(vh);
	if (dir[0] * vnorm[0] + dir[1] * vnorm[1] + dir[2] * vnorm[2] < 0)
		return false;
	else
		return true;

}
void CGeoBaseAlg::RemoveNonManifold(COpenMeshT &mesh)
{
	std::vector<COpenMeshT::FaceHandle>to_bedeleted;
	std::vector<COpenMeshT::VertexHandle>vto_deleted;
	for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); viter++)
	{
		if (mesh.is_manifold(viter) == false)
		{
			vto_deleted.push_back(viter);
			for (auto vfiter = mesh.vf_begin(viter); vfiter != mesh.vf_end(viter); vfiter++)
			{
				to_bedeleted.push_back(vfiter);
				
			}
		}
	}
	
	std::cerr << "mesh v num " << mesh.n_vertices() << std::endl;
	std::cerr << "mesh f num " << mesh.n_faces() << std::endl;
	std::cerr <<"non manifold faces:"<< to_bedeleted.size() << std::endl;
	for (int i = 0; i < to_bedeleted.size(); i++)
	{
		std::cerr << "delete face" << std::endl;
		mesh.delete_face(to_bedeleted[i],false);
	}
	
	std::cerr << "non manifold vertexs :" << vto_deleted.size() << std::endl;
	for (int i = 0; i < vto_deleted.size(); i++)
		mesh.delete_vertex(vto_deleted[i], false);
	mesh.delete_isolated_vertices();
	mesh.garbage_collection();
	std::cerr << "mesh v num " << mesh.n_vertices() << std::endl;
	std::cerr << "mesh f num " << mesh.n_faces() << std::endl;

}
double CGeoBaseAlg::Min_Dist_Point2_to_Line2(OpenMesh::Vec2d p, OpenMesh::Vec2d a, OpenMesh::Vec2d b)
{
	Point_2 cgal_p(p[0], p[1]);
	Point_2 cgal_a(a[0], a[1]);
	Point_2 cgal_b(b[0], b[1]);
	Line_2 l(cgal_a, cgal_b);
	Point_2 p_proj = l.projection(cgal_p);

	double d_ppa = std::sqrt((p_proj- cgal_a).squared_length());
	double d_ppb = std::sqrt((p_proj- cgal_b).squared_length());
	double d_ab = std::sqrt((cgal_b- cgal_a).squared_length());

	double min_dist = -1;

	if (d_ppa < d_ab && d_ppb < d_ab)
	{
		// to p_proj
		min_dist = std::sqrt((cgal_p- p_proj).squared_length());
	}
	else if (d_ppa < d_ppb)
	{
		// to a
		min_dist = std::sqrt((cgal_p - cgal_a).squared_length());
	}
	//else if (d_ppb < d_ppa)
	//{
	//	// to b
	//	min_dist = CCurveBaseAlg::TwoPointDistance(p, b);
	//}
	else
	{
		// to b
		min_dist = std::sqrt((cgal_p- cgal_b).squared_length());
		//std::cerr << "Error @ CCurveBaseAlg::Min_Dist_Point2_to_Line2" << std::endl;
	}

	return min_dist;

}
double CGeoBaseAlg::ComputeAverageEdgeLength(COpenMeshT&mesh)
{

	double len = 0;
	int count = 0;
	for (auto eiter = mesh.edges_begin(); eiter != mesh.edges_end(); eiter++)
	{
		COpenMeshT::VertexHandle vh0 = mesh.to_vertex_handle(mesh.halfedge_handle(eiter, 0));
		COpenMeshT::VertexHandle vh1 = mesh.to_vertex_handle(mesh.halfedge_handle(eiter, 1));
		OpenMesh::Vec3d p0 = mesh.point(vh0);
		OpenMesh::Vec3d p1 = mesh.point(vh1);
		len+=(p0 - p1).length();
		count++;
	}
	len /= count;
	return len;
}
void CGeoBaseAlg::ComputePerVertexNormal(COpenMeshT&mesh, std::map<COpenMeshT::VertexHandle, OpenMesh::Vec3d>& vh_normals)
{
	mesh.request_vertex_normals();
	vh_normals.clear();
	for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); viter++)
	{
		vh_normals[viter] = mesh.normal(viter).normalized();
		
	}
}
void CGeoBaseAlg::ComputePerVertexNormal(COpenMeshT&mesh, Eigen::MatrixXd &res_normals)
{
	Eigen::MatrixXd vertexs;
	Eigen::MatrixXi faces;
	CConverter::ConvertFromOpenMeshToIGL(mesh, vertexs, faces);
	ComputePerVertexNormal(vertexs, faces, res_normals);
}
void CGeoBaseAlg::ComputeMeanCurvatureVector(COpenMeshT&mesh, Eigen::MatrixXd &res_curvature)
{
	Eigen::MatrixXd vertexs;
	Eigen::MatrixXi faces;
	CConverter::ConvertFromOpenMeshToIGL(mesh, vertexs, faces);
	ComputeMeanCurvatureVector(vertexs, faces, res_curvature);

}
void CGeoBaseAlg::GetEdgeVertexs(COpenMeshT&mesh, std::vector<bool>&scalars,std::vector<OpenMesh::VertexHandle>&res_vhs)
{
	res_vhs.clear();
	for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); viter++)
	{
		int vid = viter->idx();
		if (scalars[vid] == true)
		{
			for (auto vviter = mesh.vv_begin(viter); vviter != mesh.vv_end(viter); vviter++)
			{
				if (scalars[vviter->idx()]==false)
				{
					res_vhs.push_back(viter);
					break;
				}
			}
		}
	}
}
OpenMesh::Vec3d CGeoBaseAlg::ComputeVPosFromBaryCoord(COpenMeshT&mesh, COpenMeshT::FaceHandle fh, OpenMesh::Vec3d&bary_coord)
{
	OpenMesh::Vec3d p(0, 0, 0);
	int count = 0;
	for (auto viter = mesh.fv_begin(fh); viter != mesh.fv_end(fh); viter++)
	{
		OpenMesh::Vec3d v=mesh.point(viter);
		return v;
		p=p+v*bary_coord[count];
		count++;
	}
	return p;
}
void CGeoBaseAlg::ComputeMeanCurvatureVector(Eigen::MatrixXd& vertexs, Eigen::MatrixXi& faces, Eigen::MatrixXd &res_curvature)
{
	Eigen::SparseMatrix<double> L, M, Minv;
	igl::cotmatrix(vertexs, faces, L);
	igl::massmatrix(vertexs, faces, igl::MASSMATRIX_TYPE_VORONOI, M);
	igl::invert_diag(M, Minv);
	res_curvature = -Minv*(L*vertexs);
	
}

void CGeoBaseAlg::ComputePerVertexNormal(Eigen::MatrixXd& vertexs, Eigen::MatrixXi& faces, Eigen::MatrixXd &res_normals)
{
	igl::per_vertex_normals(vertexs, faces, res_normals);
}
double CGeoBaseAlg::ComputeSquaredDistance2Plane(OpenMesh::Vec3d point, CPlane plane)
{
	return ComputeSquaredDistance2Plane(point, plane.a(), plane.b(), plane.c(), plane.d());
}
double CGeoBaseAlg::ComputeSquaredDistance2Plane(OpenMesh::Vec3d point, double a, double b, double c, double d)
{
	Plane plane(a, b, c, d);
	Point p(point[0], point[1], point[2]);
	return CGAL::squared_distance(plane, p);
}
bool CGeoBaseAlg::IsOnPositiveSide(OpenMesh::Vec3d point, CPlane plane)
{
	OpenMesh::Vec3d dir = point - plane.p();
	if (OpenMesh::dot(dir,plane.dir()) > 0)
		return true;
	else
		return false;
}