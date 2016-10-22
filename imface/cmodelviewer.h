// Copyright 2016_9 by ChenNenglun
#ifndef MODEL_VIEWER_H
#define MODEL_VIEWER_H

#include <QGLViewer/qglviewer.h>
#include <vector>
#include"../DataColle/mesh_object.h"
#include <QOpenGLFunctions_2_1>
#include"camera.h"
#include"scene.h"

class CModelViewer : public QGLViewer, public QOpenGLFunctions_2_1
{
protected:
	CCamera camera_;
	CScene *scene_;
	Eigen::Vector3d background_color_;
protected:
	virtual void draw();
	virtual void init();//called automatically when the viewer is initialized
	virtual void initializeGL();//called automatically when the opengl context is initialized

	void InitCamera();// init parameters of the camera_
	void keyPressEvent(QKeyEvent *e);
public:
	CModelViewer(QWidget *parent);
	~CModelViewer();
	void SetScene(CScene* scene);
	void SetBackgroundColor(Eigen::Vector3d background_color);


};

#endif // MODEL_VIEWER_H
