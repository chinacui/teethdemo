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
void CGeoBaseAlg::RemoveNonManifold(COpenMeshT &mesh)
{
	std::vector<COpenMeshT::FaceHandle>to_bedeleted;
	std::vector<COpenMeshT::VertexHandle>vto_deleted;
	for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); viter++)
	{
		if (mesh.is_manifold(viter) == false||mesh.is_boundary(viter))
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