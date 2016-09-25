#ifndef SCENE_H
#define SCENE_H
#include<iostream>
#include<map>
#include"scene_mesh_object.h"
#include"../DataColle/data_pool.h"
#include <QGLViewer/qglviewer.h>
#include <QOpenGLVertexArrayObject>
#include <QOpenGLBuffer>
#include <QOpenGLShaderProgram>
#include"camera.h"
class CScene
{
protected:
	std::map<int, CSceneMeshObject*> scene_mesh_;

	
public:
	CScene();
	~CScene();


	void UpdateScene();
	void Render( CCamera camera);
};
#endif
