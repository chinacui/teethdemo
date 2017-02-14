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
void CBaseObject::Rotate(OpenMesh::Vec3d axis, double angle, OpenMesh::Vec3d center )
{

	Eigen::Matrix3d m;
	m = Eigen::AngleAxisd(angle, Eigen::Vector3d(axis[0], axis[1], axis[2]));
	Eigen::Quaternion<double>rot(m);
	Eigen::Matrix3d rot_mat = rot.toRotationMatrix();
	Eigen::Matrix4d rot_mat4;
	rot_mat4.setIdentity();

	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			rot_mat4(i, j) = rot_mat(i, j);
		}
	}
	for (int i = 0; i < 3; i++)
	{
		mat_(i, 3) -= center[i];
	}
	mat_ = rot_mat4*mat_;
	for (int i = 0; i < 3; i++)
	{
		mat_(i, 3) += center[i];
	}

}
void CBaseObject::Transform(OpenMesh::Vec3d trans)
{
	Eigen::Matrix4d trans_mat;
	trans_mat.setIdentity();
	for (int i = 0; i < 3; i++)
	{
		trans_mat(i, 3)= trans[i];
	}
	mat_ = trans_mat*mat_;
	
}