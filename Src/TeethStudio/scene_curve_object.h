// Copyright 2016_9 by ChenNenglun
#ifndef SCENE_CURVE_OBJECT_H
#define SCENE_CURVE_OBJECT_H
#include"../DataColle/curve_object.h"
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
#include"qopengltexture.h"
//wraper of CMeshObject, CSceneMeshObject can be drawed in opengl context
class CSceneCurveObject :public QOpenGLFunctions_2_1
{
protected:
	std::weak_ptr<CCurveObject> curve_;
	std::vector<float> vertexs_pos_;//data of opengl buffer,size:n x 3,n is the vertex number in curve_
	std::vector<float>vertex_colors_;//data of opengl buffer,size:n x 3,n is the vertex number in curve_


	mutable QOpenGLBuffer vertex_pos_buffer_,  vertex_colors_buffer_;
	mutable QOpenGLVertexArrayObject vao_;
	mutable QOpenGLShaderProgram rendering_program_;

protected:

	void ComputeRenderingElements();
	void SetGLBufferDataFromElements();
	void InitShader();
public:
	CSceneCurveObject(std::weak_ptr<CCurveObject>  curve);
	~CSceneCurveObject();
	void UpdateRenderInfo();// if the is_changed_ of mesh is seted, the opengl buffer will be updated

	void Render(CCamera camera);
};
#endif