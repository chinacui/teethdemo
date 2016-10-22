// Copyright 2016_9 by ChenNenglun
#include"data_io.h"
#include<vector>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include<igl/writeOBJ.h>
#include"cgal_igl_converter.h"
bool CDataIO::ReadMesh(std::string fname, CMeshObject & res_mesh_obj)
{
	if (!OpenMesh::IO::read_mesh(res_mesh_obj.GetMesh(), fname))
	{
		std::cerr << "read error\n";
		return false;
	}
	COpenMeshT &mesh = res_mesh_obj.GetMesh();
	for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); viter++)
	{
		mesh.set_color(viter, OpenMesh::Vec3d(0.8, 0.8, 0.8));
	}
	res_mesh_obj.SetChanged();
	return true;
}

bool CDataIO::WriteMesh(std::string fname, CMeshObject & res_mesh_obj)
{
	Eigen::MatrixXd V;
	Eigen::MatrixXi F;
	CConverter::ConvertFromOpenMeshToIGL(res_mesh_obj.GetMesh(), V, F);
	igl::writeOBJ(fname, V, F);
	/*if (!OpenMesh::IO::write_mesh(res_mesh_obj.GetMesh(), fname))
	{
		std::cerr << "write error\n";
		return false;
	}*/
	return true;
}


