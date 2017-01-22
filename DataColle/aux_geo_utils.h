#ifndef AUX_GEO_H
#define AUX_GEO_H
#include"prereq.h"
#include<Eigen/Dense>
#include"custom_openmesh_type.h"
class DATACOLLE_CLASS CAuxGeoUtils
{
public:
	static void GetPlainMeshFromPointAndAxis(OpenMesh::Vec3d p, OpenMesh::Vec3d axis_a, OpenMesh::Vec3d axis_b, OpenMesh::Vec3d axis_dir, double scale, COpenMeshT &mesh);
	static void GenVolumeRenderAuxCubeAndCuttingPlane(double x_len, double y_len, double z_len,OpenMesh::Vec3d view_point,OpenMesh::Vec3d view_dir, OpenMesh::Vec3d view_up_dir, OpenMesh::Vec3d view_right_dir, COpenMeshT& res_cube,COpenMeshT& res_plane);
	
};
#endif