#include"geo_base_alg.h"
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/invert_diag.h>
#include<igl/per_vertex_normals.h>
#include"../DataColle/Polyhedron_type.h"
void CGeoBaseAlg::ComputeMeanCurvatureValue(Eigen::MatrixXd& vertexs, Eigen::MatrixXi& faces, Eigen::VectorXd &res_curvature)
{
	Eigen::MatrixXd HN;
	Eigen::SparseMatrix<double> L, M, Minv;
	igl::cotmatrix(vertexs, faces, L);
	igl::massmatrix(vertexs, faces, igl::MASSMATRIX_TYPE_VORONOI, M);
	igl::invert_diag(M, Minv);
	HN = -Minv*(L*vertexs);
	res_curvature = HN.rowwise().norm(); //up to sign
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
		p = (p - mmin) / lens[mi];
		p = p - COpenMeshT::Point(lens[0]/lens[mi] / 4.0, lens[1]/lens[mi] / 4.0, lens[2]/lens[mi]/2.0);
		mesh.set_point(viter, p);
	}
	//std::cerr << "mi " << mi <<" "<< lens[1] / lens[mi]<< std::endl;
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