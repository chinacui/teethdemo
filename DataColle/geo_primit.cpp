#include"geo_primit.h"

CPlane::CPlane(OpenMesh::Vec3d p, OpenMesh::Vec3d dir)
{
	p_ = p;
	dir_ = dir;
	dir_.normalize();
	ComputePlaneParamFromMeanAndDir();
}
void CPlane::ComputePlaneParamFromMeanAndDir()
{
	a_ = p_[0];
	b_ = p_[1];
	c_ = p_[2];
	d_ = -dir_[0]*p_[0] - dir_[1]*p_[1] - dir_[2]*p_[2];

}
void CPlane::SetP(OpenMesh::Vec3d p)
{
	p_ = p;
	ComputePlaneParamFromMeanAndDir();
}
void CPlane::SetDir(OpenMesh::Vec3d dir)
{
	dir_ = dir;
	dir_=dir_.normalize();
	ComputePlaneParamFromMeanAndDir();
}