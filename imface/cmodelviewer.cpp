// Copyright 2016_9_25 by ChenNenglun
#include "cmodelviewer.h"

#include <cmath>
#include <float.h>

using namespace std;


void CModelViewer::draw()
{
	glClearColor(background_color_(0), background_color_(1), background_color_(2),1);
	if(scene_!=NULL)
	scene_->Render( camera_);

}
void CModelViewer::SetBackgroundColor(Eigen::Vector3d background_color)
{
	background_color_ = background_color;
}
void CModelViewer::initializeGL()
{


	QGLViewer::initializeGL();

	initializeOpenGLFunctions();



}
CModelViewer::CModelViewer(QWidget *parent)
{
	background_color_ = Eigen::Vector3d(1, 1, 1);
	scene_ = NULL;
	
	
}
void CModelViewer::InitCamera()
{
	
	camera_.setAspectRatio(1);
	camera_.setPosition(qglviewer::Vec(0, 0, 1));
	camera_.lookAt(qglviewer::Vec(0, 0, 0));

}
void CModelViewer::SetScene(CScene* scene)
{
	scene_ = scene;
}
CModelViewer::~CModelViewer()
{

}
void CModelViewer::init()
{
	InitCamera();
	setCamera(&camera_);
	setMouseBinding(Qt::ControlModifier, Qt::LeftButton, QGLViewer::CAMERA, QGLViewer::ROTATE);
	setMouseBinding(Qt::ControlModifier, Qt::RightButton, QGLViewer::CAMERA, QGLViewer::TRANSLATE);
	setMouseBinding(Qt::ControlModifier, Qt::MidButton, QGLViewer::CAMERA, QGLViewer::ZOOM);
	setWheelBinding(Qt::ControlModifier, QGLViewer::CAMERA, QGLViewer::ZOOM);
}
