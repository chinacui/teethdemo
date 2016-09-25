// Copyright 2016_9 by ChenNenglun
#ifndef CBASE_OBJECT_H
#define CBASE_OBJECT_H
#include"prereq.h"
#include<Eigen/Dense>
class DATACOLLE_CLASS CBaseObject
{
protected:
	Eigen::Matrix4d mat_;//transform and rotate of the local coordinate
public:
	CBaseObject();
	~CBaseObject();
	Eigen::Matrix4d & GetMatrix();
	void SetMatrix(Eigen::Matrix4d &mat);
};
#endif