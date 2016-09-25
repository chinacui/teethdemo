// Copyright 2016_9 by ChenNenglun
#ifndef CCAMERA_H
#define CCAMERA_H
#include"QGLViewer/camera.h"
#include<Eigen/Dense>
class CCamera :public qglviewer::Camera
{
private:
	CCamera(const qglviewer::Camera &b) { qglviewer::Camera::Camera(b); };
public:
	CCamera();
	CCamera(const CCamera& c);
	CCamera & CCamera::operator=(const CCamera &b);
	

};

#endif