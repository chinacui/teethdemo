#ifndef GEO_PRIMIT_H
#define GEO_PRIMIT_H
#include"prereq.h"
#include"base_object.h"
#include<Eigen\Core>
#include"custom_openmesh_type.h"
class DATACOLLE_CLASS CPlane:public CBaseObject
{
protected:
	double a_;
	double b_;
	double c_;
	double d_;
	OpenMesh::Vec3d p_;
	OpenMesh::Vec3d dir_;
protected:
	void ComputePlaneParamFromMeanAndDir();
public:
	
	CPlane(OpenMesh::Vec3d p, OpenMesh::Vec3d dir);
	CPlane() {};
	double a() { return a_; }
	double b() { return b_; }
	double c() { return c_; }
	double d() { return d_; }
	OpenMesh::Vec3d p() { return p_; }
	OpenMesh::Vec3d dir() { return dir_; }
	void SetP(OpenMesh::Vec3d p);
	void SetDir(OpenMesh::Vec3d dir);
};
#endif