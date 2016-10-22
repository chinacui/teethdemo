#ifndef AUX_GEO_H
#define AUX_GEO_H
#include"prereq.h"
#include<Eigen/Dense>
#include"custom_openmesh_type.h"
class DATACOLLE_CLASS CAuxGeoUtils
{
public:
	static void GetPlaneMeshFromPointAndAxis(OpenMesh::Vec3d p, OpenMesh::Vec3d axis_a, OpenMesh::Vec3d axis_b, OpenMesh::Vec3d axis_dir, double scale, COpenMeshT &mesh);
};
#endif