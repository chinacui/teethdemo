#include"rigid_icp.h"
#include"../AdditionalLibs/libicp/libicp_wrapper.h"
#include"../DataColle/cgal_igl_converter.h"
CRigidIcp::CRigidIcp(COpenMeshT* src_mesh)
{
	src_mesh_ = src_mesh;
}
void CRigidIcp::FitTarget(COpenMeshT *tgt_mesh)
{
	Eigen::MatrixXd src_vertexs, tgt_vertexs;
	Eigen::MatrixXi src_faces, tgt_faces;
	CConverter::ConvertFromOpenMeshToIGL(*src_mesh_, src_vertexs, src_faces);
	CConverter::ConvertFromOpenMeshToIGL(*tgt_mesh, tgt_vertexs, tgt_faces);
	//perform icp
	CIcpPoint2Point icp(tgt_vertexs);
	Eigen::Matrix3d rot;
	Eigen::Vector3d trans;
	Eigen::MatrixXd X;//transform matrix
	icp.Fit(src_vertexs, rot, trans);
	X.resize(src_vertexs.rows() * 4, 3);
	X.setZero();
	rot.transposeInPlace();
	Eigen::SparseMatrix<double>D;
	D.resize(src_vertexs.rows(), src_vertexs.rows() * 4);
	std::vector < Eigen::Triplet<double>>trip_d(0);

	for (int i = 0; i < src_vertexs.rows(); i++)
	{
		for (int j = 0; j < 3; j++)
		{
			trip_d.push_back(Eigen::Triplet<double>(i, i * 4 + j, src_vertexs(i, j)));
		}
		trip_d.push_back(Eigen::Triplet<double>(i, i * 4 + 3, 1));
	}
	for (int i = 0; i < src_vertexs.rows(); i++)
	{
		for (int j = 0; j < 3; j++)
		{
			X(i * 4 + j, 0) = rot(j, 0);
			X(i * 4 + j, 1) = rot(j, 1);
			X(i * 4 + j, 2) = rot(j, 2);
		}
		X(i * 4 + 3, 0) = trans(0);
		X(i * 4 + 3, 1) = trans(1);
		X(i * 4 + 3, 2) = trans(2);
	}
	Eigen::MatrixXd transformed;
	transformed = (D*X).eval();
	for (auto viter = src_mesh_->vertices_begin(); viter != src_mesh_->vertices_end(); viter++)
	{
		int id = viter->idx();
		OpenMesh::Vec3d v(transformed(id, 0), transformed(id, 1), transformed(id, 2));
		src_mesh_->set_point(viter, v);
	}
}