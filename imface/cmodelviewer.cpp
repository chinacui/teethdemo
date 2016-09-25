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


	setBackgroundColor(::Qt::white);

	std::cerr << "initializeGL" << std::endl;
}
CModelViewer::CModelViewer(QWidget *parent)
{
	background_color_ = Eigen::Vector3d(1, 1, 1);
	scene_ = NULL;
	
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

}
