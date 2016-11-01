// Copyright 2016_9 by ChenNenglun
#include"camera.h"
#include"QGLViewer/quaternion.h"
#include"qquaternion.h"
#include"qvector3d.h"
CCamera::CCamera() :Camera()
{

}


CCamera::CCamera(const CCamera& c) :Camera(c)
{

	
}

CCamera & CCamera::operator=(const CCamera &b)
{

	*((Camera*)this) = (Camera)b;

	return *this;
}
void CCamera::ConvertClickToLine(QPoint p, OpenMesh::Vec3d &orig, OpenMesh::Vec3d &dir)
{
	qglviewer::Vec torig, tdir;
	convertClickToLine(p, torig, tdir);
	orig = OpenMesh::Vec3d(torig.x, torig.y, torig.z);
	dir = OpenMesh::Vec3d(tdir.x, tdir.y, tdir.z);
}
