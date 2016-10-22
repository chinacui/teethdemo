// Copyright 2016_9 by ChenNenglun
#include"scene_mesh_object.h"
#include<fstream>
void CSceneMeshObject::SetGLBufferDataFromElements()
{
	rendering_program_.bind();
	
	if (vao_.isCreated())
		vao_.destroy();
	vao_.create();
	vao_.bind();

	if (vertex_pos_buffer_.isCreated())
		vertex_pos_buffer_.destroy();
	vertex_pos_buffer_.create();
	vertex_pos_buffer_.bind();
	vertex_pos_buffer_.allocate(vertexs_pos_.data(), static_cast<int>(vertexs_pos_.size() * sizeof(float)));
	auto v_pos_loc = rendering_program_.attributeLocation("a_pos");
	rendering_program_.enableAttributeArray(v_pos_loc);
	rendering_program_.setAttributeBuffer(v_pos_loc, GL_FLOAT, 0, 3);
	vertex_pos_buffer_.release();

	if (vertex_normals_buffer_.isCreated())
		vertex_normals_buffer_.destroy();
	vertex_normals_buffer_.create();
	vertex_normals_buffer_.bind();
	vertex_pos_buffer_.allocate(vertex_normals_.data(), static_cast<int>(vertex_normals_.size() * sizeof(float)));
	auto v_normal_loc = rendering_program_.attributeLocation("a_normal");
	rendering_program_.enableAttributeArray(v_normal_loc);
	rendering_program_.setAttributeBuffer(v_normal_loc, GL_FLOAT, 0, 3);
	vertex_normals_buffer_.release();
	

	if (vertex_colors_buffer_.isCreated())
		vertex_colors_buffer_.destroy();
	vertex_colors_buffer_.create();
	vertex_colors_buffer_.bind();
	vertex_colors_buffer_.allocate(vertex_colors_.data(), static_cast<int>(vertex_colors_.size() * sizeof(float)));
	auto v_color_loc = rendering_program_.attributeLocation("a_color");
	rendering_program_.enableAttributeArray(v_color_loc);
	rendering_program_.setAttributeBuffer(v_color_loc, GL_FLOAT, 0, 3);
	vertex_colors_buffer_.release();
	vao_.release();
	rendering_program_.release();
}
void CSceneMeshObject::UpdateRenderInfo()
{
	if (mesh_.lock().get()->IsChanged())
	{
		ComputeRenderingElements();
		SetGLBufferDataFromElements();
		mesh_.lock().get()->SetChanged(false);
	}

}
void CSceneMeshObject::Render( CCamera camera)
{

	UpdateRenderInfo();
	vao_.bind();
	rendering_program_.bind();
	
	//compute attributes
	QMatrix4x4 mvpMatrix;
	QMatrix4x4 mvMatrix, frame_matrix, scale_matrix;
	double mat[16];
	camera.getModelViewProjectionMatrix(mat);

	auto local_mat= mesh_.lock().get()->GetMatrix();
	QMatrix4x4 q_local_mat;
	for (int i = 0; i < local_mat.rows(); i++)
	{
		for (int j = 0; j < local_mat.cols(); j++)
		{
			q_local_mat(i, j) = local_mat(i, j);
		}
	}
	for (int i = 0; i < 16; i++)
	{
		mvpMatrix.data()[i] = (float)mat[i];
	}
	mvpMatrix = mvpMatrix*q_local_mat;
	camera.getModelViewMatrix(mat);
	for (int i = 0; i < 16; i++)
	{
		mvMatrix.data()[i] = (float)mat[i];
	}
	mvMatrix = mvMatrix*q_local_mat;
	

	QVector4D	ambient(0.1,0.1,0.1,1.0);
	QVector4D	diffuse(1, 1, 1, 1);
	QVector4D	specular(0,0,0,1);
	QVector4D	light_pos(0.0f, 0.0f, 100.0f, 1.0f);

	int mvp_mat_loc = rendering_program_.uniformLocation("mvp_matrix");
	int mv_mat_loc  = rendering_program_.uniformLocation("mv_matrix");
	int light_pos_loc = rendering_program_.uniformLocation("u_light_pos");
	int light_diff_loc = rendering_program_.uniformLocation("u_light_diff");
	int light_spec_loc = rendering_program_.uniformLocation("u_light_spec");
	int light_amb_loc = rendering_program_.uniformLocation("u_light_amb");
	int spec_power_loc = rendering_program_.uniformLocation("u_spec_power");


		
	rendering_program_.setUniformValue(light_pos_loc, light_pos);
	rendering_program_.setUniformValue(mvp_mat_loc, mvpMatrix);
	rendering_program_.setUniformValue(mv_mat_loc, mvMatrix);
	rendering_program_.setUniformValue(light_diff_loc, diffuse);
	rendering_program_.setUniformValue(light_spec_loc, specular);
	rendering_program_.setUniformValue(light_amb_loc, ambient);


	
	glDrawArrays(GL_TRIANGLES, 0, static_cast<GLsizei>(vertexs_pos_.size() / 3));
	
	rendering_program_.release();
	vao_.release();
}
void CSceneMeshObject::InitShader()
{
	std::ifstream v_ifstream("mesh_render.vert");
	std::string vertex_shader_source((std::istreambuf_iterator<char>(v_ifstream)), std::istreambuf_iterator<char>());
	v_ifstream.close();

	std::cerr << "Init shaders" << std::endl;
	QOpenGLShader *vertex_shader = new QOpenGLShader(QOpenGLShader::Vertex);
	if (!vertex_shader->compileSourceCode(vertex_shader_source.c_str()))
	{
		std::cerr << "Compiling vertex source FAILED" << std::endl;
	}

	std::ifstream f_ifstream("mesh_render.frag");
	std::string frag_shader_source((std::istreambuf_iterator<char>(f_ifstream)), std::istreambuf_iterator<char>());
	f_ifstream.close();
	QOpenGLShader *fragment_shader = new QOpenGLShader(QOpenGLShader::Fragment);
	if (!fragment_shader->compileSourceCode(frag_shader_source.c_str()))
	{
		std::cerr << "Compiling fragmentsource FAILED" << std::endl;
	}

	if (!rendering_program_.addShader(vertex_shader))
	{
		std::cerr << "adding vertex shader FAILED" << std::endl;
	}
	if (!rendering_program_.addShader(fragment_shader))
	{
		std::cerr << "adding fragment shader FAILED" << std::endl;
	}
	if (!rendering_program_.link())
	{
		std::cerr << "linking Program FAILED" << std::endl;
	}
}
void CSceneMeshObject::ComputeRenderingElements()
{
	vertexs_pos_.resize(0);
	vertex_normals_.resize(0);
	vertex_colors_.resize(0);
		
	CMeshObject* mesh = mesh_.lock().get();
	auto openmesh = mesh->GetMesh();
	

	COpenMeshT::FaceIter fiter = openmesh.faces_begin();

	for (int i = 0; i < openmesh.n_faces(); i++, fiter++)
	{
		COpenMeshT::FaceHalfedgeIter fhi = openmesh.fh_begin(*fiter), fhend = openmesh.fh_end(*fiter);
		std::vector<COMTraits::Point> vs;

		while (fhi != fhend)
		{
			int vid=openmesh.to_vertex_handle(*fhi).idx();
			auto vpos =openmesh.points()[vid];
			
			auto vcolor = openmesh.vertex_colors()[openmesh.to_vertex_handle(fhi).idx()];

			for (int j = 0; j < 3; j++)
			{
				vertexs_pos_.push_back(vpos[j]);
				vertex_colors_.push_back(vcolor[j]);
			}
			vs.push_back(vpos);
			fhi++;
		}
	


		auto  ab = vs[1] - vs[0];
		auto ac = vs[2] - vs[0];
		auto fnormal = ab%ac;
		for (int j = 0; j < 3; j++)
		{
			for (int k = 0; k < 3; k++)
			{
				vertex_normals_.push_back(fnormal[k]);
			}
		}
	}

	

	

}



CSceneMeshObject::CSceneMeshObject(std::weak_ptr<CMeshObject>  mesh)
{
	InitShader();
	this->initializeOpenGLFunctions();
	mesh_ = mesh;
	mesh_.lock().get()->SetChanged(true);
	UpdateRenderInfo();
	mesh_.lock().get()->SetChanged(false);

}
CSceneMeshObject::~CSceneMeshObject()
{

}