// Copyright 2016_9 by ChenNenglun
#include"base_object.h"
Eigen::Matrix4d & CBaseObject::GetMatrix()
{
	return mat_;
}
void CBaseObject::SetMatrix(Eigen::Matrix4d &mat)
{
	mat_ = mat;
}
CBaseObject::CBaseObject()
{
	mat_.setIdentity();
}
CBaseObject::~CBaseObject()
{

}