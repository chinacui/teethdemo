#ifndef SCENE_MESH_OBJECT_H
#define SCENE_MESH_OBJECT_H
#include"../DataColle/mesh_object.h"
#include<iostream>
#include<vector>
#include <QGLViewer/qglviewer.h>
#include <QOpenGLVertexArrayObject>
#include <QOpenGLBuffer>
#include <QOpenGLShaderProgram>
#include<memory>
#include <QGLViewer/qglviewer.h>
#include <QOpenGLFunctions_2_1>
#include"camera.h"
class CSceneMeshObject :public QOpenGLFunctions_2_1
{
protected:
	std::weak_ptr<CMeshObject> mesh_;
	std::vector<float> vertexs_pos_;
	std::vector<float> vertex_normals_;
	std::vector<float>vertex_colors_;

	mutable QOpenGLBuffer vertex_pos_buffer_,vertex_normals_buffer_,vertex_colors_buffer_;
	mutable QOpenGLVertexArrayObject vao_;
	mutable QOpenGLShaderProgram rendering_program_;

protected:
	void ComputeRenderingElements();
	void SetGLBufferDataFromElements();
	void InitShader();
public:
	CSceneMeshObject(std::weak_ptr<CMeshObject>  mesh);
	~CSceneMeshObject();
	void UpdateRenderInfo();

	void Render(CCamera camera);
};
#endif