// Copyright 2016_9 by ChenNenglun
#ifndef SCENE_H
#define SCENE_H
#include<iostream>
#include<map>
#include"scene_mesh_object.h"
#include"scene_curve_object.h"
#include"scene_volume_object.h"
#include"../DataColle/data_pool.h"
#include <QGLViewer/qglviewer.h>
#include <QOpenGLVertexArrayObject>
#include <QOpenGLBuffer>
#include <QOpenGLShaderProgram>
#include"camera.h"
class CScene
{
protected:
	std::map<int, CSceneMeshObject*> scene_mesh_;//mesh object in scene
	std::map<int, CSceneCurveObject*>scene_curve_;
	std::map<int, CSceneVolumeObject*>scene_volumedata_;
	void UpdateScene();//called each frame, check and update the scene when datapool is updated
public:
	CScene();
	~CScene();


	
	void Render( CCamera camera);//render the scene
};
#endif
