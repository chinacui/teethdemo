// Copyright 2016_9 by ChenNenglun
#ifndef CBASE_OBJECT_H
#define CBASE_OBJECT_H
#include"prereq.h"
#include<Eigen/Dense>
#include"custom_openmesh_type.h"
class DATACOLLE_CLASS CBaseObject
{
protected:
	Eigen::Matrix4d mat_;//transform and rotate of the local coordinate
public:
	CBaseObject();
	~CBaseObject();
	Eigen::Matrix4d & GetMatrix();
	void SetMatrix(Eigen::Matrix4d &mat);
	void Transform(OpenMesh::Vec3d trans);
	void Rotate(OpenMesh::Vec3d axis,double angle,OpenMesh::Vec3d center=OpenMesh::Vec3d(0,0,0));
};
#endif