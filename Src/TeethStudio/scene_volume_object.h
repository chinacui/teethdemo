#ifndef CSCENE_VOLUME_OBJECT_H
#define CSCENE_VOLUME_OBJECT_H
#include<iostream>
#include<vector>
#include <QGLViewer/qglviewer.h>
#include <QOpenGLVertexArrayObject>
#include <QOpenGLBuffer>
#include <QOpenGLShaderProgram>
#include<qopengltexture.h>
#include<memory>
#include <QGLViewer/qglviewer.h>
#include <QOpenGLFunctions_2_1>
#include"camera.h"
#include"qopengltexture.h"
#include"../DataColle/aux_geo_utils.h"
#include"../DataColle/custom_openmesh_type.h"
#include"../DataColle/volume_data_object.h"
#include <QOpenGLFramebufferObjectFormat>
class CSceneVolumeObject :public QOpenGLFunctions_2_1
{
protected:
	
	std::weak_ptr<CVolumeDataObject> volume_data_;
	int number_of_slices_;
	int slices_over_x_;//width
	int slices_over_y_;//height
	GLuint volume_data_texture_id_;

	COpenMeshT aux_cube_,aux_plane_;
	std::vector<float> aux_cube_vertexs_pos_, aux_plane_vertexs_pos_;//data of opengl buffer,size:n x 3,n is the vertex number in mesh_
	std::vector<float>aux_cube_tex_coords_, aux_plane_tex_coords_;

	QSize offscreen_frame_size_;

	mutable QOpenGLBuffer back_texcoord_aux_cube_vertex_pos_buffer_, back_texcoord_aux_cube_tex_coords_buffer_;
	mutable QOpenGLBuffer front_texcoord_aux_cube_vertex_pos_buffer_, front_texcoord_aux_cube_tex_coords_buffer_, front_texcoord_aux_plane_vertex_pos_buffer_, front_texcoord_aux_plane_tex_coords_buffer_;
	mutable QOpenGLBuffer depth_aux_cube_vertex_pos_buffer_, depth_aux_cube_tex_coords_buffer_;
	mutable QOpenGLBuffer volume_aux_cube_vertex_pos_buffer_;
	mutable QOpenGLVertexArrayObject back_texcoord_vao_, front_cube_texcoord_vao_, front_plane_texcoord_vao_, depth_vao_, volume_vao_;
	mutable QOpenGLShaderProgram back_texcoord_rendering_program_, front_texcoord_rendering_program_,depth_rendering_program_,volume_rendering_program_;


	QOpenGLFramebufferObject* back_texcoord_fbo_;
	QOpenGLFramebufferObject* front_texcoord_fbo_;

	void SetVolumeDataTextureFromVolumeData();
	void SetBackTexcoordBufferFromElements();
	void SetFrontTexcoordBufferFromElements();
	void SetVolumeRenderBufferFromElements();
	void RenderBackTexcoord(CCamera& camera);
	void RenderFrontTexcoord(CCamera& camera);
	void RenderVolumeData(CCamera&camera);
	void SetCommonRenderAttributes(QOpenGLShaderProgram &program, CCamera& camera);

	void GenBackTexcoordFrameBuffer();
	void GenFrontTexcoordFrameBuffer();
public:
	void ComputeRenderingElements(CCamera camera);
	void SetGLBufferDataFromElements();
	void InitShader();
	CSceneVolumeObject(std::weak_ptr<CVolumeDataObject>  volumedata);
	~CSceneVolumeObject();
	void UpdateRenderInfo();// if the is_changed_ of mesh is seted, the opengl buffer will be updated

	void Render(CCamera camera);
};
#endif