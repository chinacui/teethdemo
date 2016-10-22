#include"aux_geo_utils.h"
#include<iostream>
void CAuxGeoUtils::GetPlaneMeshFromPointAndAxis(OpenMesh::Vec3d p, OpenMesh::Vec3d axis_a, OpenMesh::Vec3d axis_b, OpenMesh::Vec3d axis_dir, double scale, COpenMeshT &mesh)
{
	
	axis_a.normalize();
	double half_scale =scale/2;
	axis_a = axis_a*half_scale;
	axis_b.normalize();
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