#ifndef CSHAPE_DEFORMATION_H
#define CSHAPE_DEFORMATION_H
#include"prereq.h"
#include <igl/arap.h>
#include <Eigen/Geometry>
#include <Eigen/StdVector>
#include"../DataColle/custom_openmesh_type.h"
#include"../DataColle/mesh_object.h"
#include"../DataColle/cgal_igl_converter.h"
class CArapShapeDeformation
{
protected:
public:
	CArapShapeDeformation(CMeshObject *meshobj);
};
#endif