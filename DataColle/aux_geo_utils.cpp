#include"aux_geo_utils.h"
#include<iostream>
void CAuxGeoUtils::GetPlainMeshFromPointAndAxis(OpenMesh::Vec3d p, OpenMesh::Vec3d axis_a, OpenMesh::Vec3d axis_b, OpenMesh::Vec3d axis_dir, double scale, COpenMeshT &mesh)
{
	
	//axis_a.normalize();
	double half_scale =scale/2;
	axis_a = axis_a*half_scale;
	//axis_b.normalize();
	axis_b = axis_b*half_scale;
	OpenMesh::Vec3d vs[4];
	vs[0] = p+axis_a + axis_b;
	vs[1] = p + axis_a - axis_b;
	vs[2] = p - axis_a - axis_b;
	vs[3] = p - axis_a + axis_b;
	COpenMeshT::VertexHandle vhs[4];
	for (int i = 0; i < 4; i++)
	{
		vhs[i] = mesh.add_vertex(vs[i]);
	}
	
	if (OpenMesh::dot(OpenMesh::cross(axis_b,axis_a),(axis_dir)) > 0)
	{
		mesh.add_face(vhs[0], vhs[1], vhs[2]);
		mesh.add_face(vhs[2], vhs[3], vhs[0]);
	}
	else
	{
		mesh.add_face(vhs[0], vhs[2], vhs[1]);
		mesh.add_face(vhs[2], vhs[0], vhs[3]);
	}
	
}
void CAuxGeoUtils::GenVolumeRenderAuxCubeAndCuttingPlane(double x_len, double y_len, double z_len, OpenMesh::Vec3d view_point, OpenMesh::Vec3d view_dir, OpenMesh::Vec3d view_right_dir, OpenMesh::Vec3d view_up_dir, COpenMeshT& res_mesh, COpenMeshT& res_plane)
{
	res_mesh = COpenMeshT();
	res_plane = COpenMeshT();
	res_mesh.request_halfedge_texcoords3D();
	res_mesh.request_vertex_texcoords3D();
	res_plane.request_halfedge_texcoords3D();
	res_plane.request_vertex_texcoords3D();
	OpenMesh::Vec3d vs[8];
	vs[0] = OpenMesh::Vec3d(-x_len/2,-y_len/2, z_len / 2);
	vs[1] = OpenMesh::Vec3d(x_len / 2, -y_len / 2, z_len/2);
	vs[2] = OpenMesh::Vec3d(x_len / 2, -y_len / 2, -z_len / 2);
	vs[3] = OpenMesh::Vec3d(-x_len / 2, -y_len / 2, -z_len / 2);

	vs[4] = OpenMesh::Vec3d(-x_len / 2, y_len / 2, z_len / 2);
	vs[5] = OpenMesh::Vec3d(x_len / 2, y_len / 2, z_len / 2);
	vs[6] = OpenMesh::Vec3d(x_len / 2, y_len / 2, -z_len / 2);
	vs[7] = OpenMesh::Vec3d(-x_len / 2, y_len / 2, -z_len / 2);
	
	COpenMeshT::VertexHandle vhs[8];
	for (int i = 0; i < 8; i++)
	{
		vhs[i] = res_mesh.add_vertex(vs[i]);
	}

	res_mesh.add_face(vhs[0], vhs[2], vhs[1]);
	res_mesh.add_face(vhs[0], vhs[3], vhs[2]);

	res_mesh.add_face(vhs[0], vhs[1], vhs[5]);
	res_mesh.add_face(vhs[0], vhs[5], vhs[4]);

	res_mesh.add_face(vhs[1], vhs[2], vhs[6]);
	res_mesh.add_face(vhs[1], vhs[6], vhs[5]);

	res_mesh.add_face(vhs[4], vhs[6], vhs[7]);
	res_mesh.add_face(vhs[4], vhs[5], vhs[6]);

	res_mesh.add_face(vhs[0], vhs[4], vhs[3]);
	res_mesh.add_face(vhs[4], vhs[7], vhs[3]);

	res_mesh.add_face(vhs[3], vhs[7], vhs[6]);
	res_mesh.add_face(vhs[3], vhs[6], vhs[2]);

	for (auto hiter = res_mesh.halfedges_begin(); hiter != res_mesh.halfedges_end(); hiter++)
	{
		auto vh = res_mesh.to_vertex_handle(hiter);
		auto point = res_mesh.point(vh);
		double u = point[0] >0 ? 1 : 0;
		double v = point[1] > 0 ? 1 : 0;
		double w = point[2] > 0 ? 0 :1;
		res_mesh.data(*hiter).SetUV3D(COpenMeshT::TexCoord3D(u, v, w));
	}
	OpenMesh::Vec3d uv_coord[3];
	OpenMesh::Vec3d uv_center(-x_len/2, -y_len/2, z_len/2);
	uv_coord[0] = OpenMesh::Vec3d(x_len, 0, 0);
	uv_coord[1] = OpenMesh::Vec3d(0, y_len, 0);
	uv_coord[2] = OpenMesh::Vec3d(0, 0, -z_len);


	OpenMesh::Vec3d planevs[4];
	planevs[0] = view_point - view_up_dir - view_right_dir;
	planevs[1] = view_point - view_up_dir + view_right_dir;
	planevs[2] = view_point + view_up_dir + view_right_dir;
	planevs[3] = view_point + view_up_dir - view_right_dir;
	COpenMeshT::VertexHandle planevhs[4];
	for (int i = 0; i < 4; i++)
	{
		planevhs[i] = res_plane.add_vertex(planevs[i]);
	}
	res_plane.add_face(planevhs[0], planevhs[1], planevhs[2]);
	res_plane.add_face(planevhs[0], planevhs[2], planevhs[3]);

	for (int i = 0; i < 4; i++)
	{
		OpenMesh::Vec3d dir = planevs[i] - uv_center;

		double uvx=  OpenMesh::dot(dir,uv_coord[0])/ uv_coord[0].sqrnorm();
		double uvy = OpenMesh::dot(dir,uv_coord[1]) / uv_coord[1].sqrnorm();
		double uvz = OpenMesh::dot(dir,uv_coord[2]) / uv_coord[2].sqrnorm();
		for (auto vhiter = res_plane.vih_begin(planevhs[i]); vhiter != res_plane.vih_end(planevhs[i]); vhiter++)
		{
			res_plane.data(*vhiter).SetUV3D(COpenMeshT::TexCoord3D(uvx, uvy, uvz));
		}
	}


}