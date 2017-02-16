// Copyright 2016_9_25 by ChenNenglun
#include "cmodelviewer.h"

#include <QKeyEvent>

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
	init();
	
}
void CModelViewer::InitCamera()
{
	camera_.setType(qglviewer::Camera::Type::ORTHOGRAPHIC);
	//camera_.setAspectRatio(1);
	camera_.setPosition(qglviewer::Vec(0, 0, 10));
	camera_.lookAt(qglviewer::Vec(0, 0, 0));
	camera_.setUpVector(qglviewer::Vec(0, 1, 0));
	camera_.setZClippingCoefficient(1000);
	std::cerr << "zclip " << camera_.zClippingCoefficient() << std::endl;
}
void CModelViewer::keyReleaseEvent(QKeyEvent*e)
{
	action_manager_->KeyReleaseEvent(e);
	this->update();
}
void CModelViewer::keyPressEvent(QKeyEvent *e)
{
	action_manager_->KeyPressEvent(e);
	this->update();
	
}
void CModelViewer::SetScene(CScene* scene)
{
	scene_ = scene;
}
CModelViewer::~CModelViewer()
{

}
void CModelViewer::mousePressEvent(QMouseEvent *e)
{
	if (e->modifiers() == Qt::Modifier::CTRL)
		QGLViewer::mousePressEvent(e);

	action_manager_->MousePressEvent(e);
	this->update();
}
void CModelViewer::mouseMoveEvent(QMouseEvent *e)
{
	if (e->modifiers() == Qt::Modifier::CTRL)
		QGLViewer::mouseMoveEvent(e);
	
	action_manager_->MouseMoveEvent(e);
	this->update();
}
void CModelViewer::mouseReleaseEvent(QMouseEvent *e)
{
	if (e->modifiers() == Qt::Modifier::CTRL)
		QGLViewer::mouseMoveEvent(e);
	
	action_manager_->MouseReleaseEvent(e);
	this->update();
}
void CModelViewer::mouseDoubleClickEvent(QMouseEvent *e)
{
	action_manager_->MouseDoubleClickEvent(e);
	this->update();
}
void CModelViewer::wheelEvent(QWheelEvent *e)
{
	if (e->modifiers() == Qt::Modifier::CTRL)
		QGLViewer::wheelEvent(e);
	else
	action_manager_->WheelEvent(e);
	this->update();
}
void CModelViewer::init()
{
	InitCamera();
	setCamera(&camera_);
	setMouseBinding(Qt::ControlModifier, Qt::LeftButton, QGLViewer::CAMERA, QGLViewer::ROTATE);
	setMouseBinding(Qt::ControlModifier, Qt::RightButton, QGLViewer::CAMERA, QGLViewer::TRANSLATE);
	//setMouseBinding(Qt::ControlModifier, Qt::MidButton, QGLViewer::CAMERA, QGLViewer::ZOOM);
	setWheelBinding(Qt::ControlModifier, QGLViewer::CAMERA, QGLViewer::ZOOM);

	action_manager_ = new CActionManager();
	action_manager_->Init(this);


}
