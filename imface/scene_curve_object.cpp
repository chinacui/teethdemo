#include"scene_curve_object.h"
#include<fstream>
#include"ui_context.h"
#include"texture_wrapper.h"
void CSceneCurveObject::SetGLBufferDataFromElements()
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
void CSceneCurveObject::UpdateRenderInfo()
{

	
	if (curve_.lock().get()->IsChanged())
	{
		ComputeRenderingElements();
		SetGLBufferDataFromElements();
		curve_.lock().get()->SetChanged(false);
	}

}
void CSceneCurveObject::Render(CCamera camera)
{

	if (curve_.lock()->GetCurve().size() < 2&& curve_.lock()->RendereType()==CCurveObject::Default)
		return;
	UpdateRenderInfo();
	vao_.bind();
	rendering_program_.bind();

	//compute attributes
	QMatrix4x4 mvpMatrix;
	QMatrix4x4 mvMatrix, frame_matrix, scale_matrix;
	double mat[16];
	camera.getModelViewProjectionMatrix(mat);

	auto local_mat = curve_.lock().get()->GetMatrix();
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


	QVector4D	ambient(0.1, 0.1, 0.1, 1.0);
	QVector4D	diffuse(1, 1, 1, 1);
	QVector4D	specular(0, 0, 0, 1);
	QVector4D	light_pos(0.0f, 0.0f, 100.0f, 1.0f);

	int mvp_mat_loc = rendering_program_.uniformLocation("mvp_matrix");
	int mv_mat_loc = rendering_program_.uniformLocation("mv_matrix");

	


	rendering_program_.setUniformValue(mvp_mat_loc, mvpMatrix);
	rendering_program_.setUniformValue(mv_mat_loc, mvMatrix);



	glLineWidth(6);
	glPointSize(6);

	if (curve_.lock().get()->RendereType() == CCurveObject::CurveType::Default)
	{
		glDrawArrays(GL_LINE_STRIP, 0, static_cast<GLsizei>(vertexs_pos_.size() / 3));
	}
	else
	{
		glPointSize(10);
		glDrawArrays(GL_POINTS, 0, static_cast<GLsizei>(vertexs_pos_.size() / 3));
	}
	
	rendering_program_.release();
	vao_.release();
}
void CSceneCurveObject::InitShader()
{
	std::ifstream v_ifstream("curve_render.vert");
	std::string vertex_shader_source((std::istreambuf_iterator<char>(v_ifstream)), std::istreambuf_iterator<char>());
	v_ifstream.close();

	//std::cerr << "Init shaders" << std::endl;
	QOpenGLShader *vertex_shader = new QOpenGLShader(QOpenGLShader::Vertex);
	if (!vertex_shader->compileSourceCode(vertex_shader_source.c_str()))
	{
		std::cerr << "Compiling vertex source FAILED" << std::endl;
	}

	std::ifstream f_ifstream("curve_render.frag");
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
void CSceneCurveObject::ComputeRenderingElements()
{
	vertexs_pos_.resize(0);

	vertex_colors_.resize(0);


	CCurveObject* pcurve = curve_.lock().get();
	auto curve = pcurve->GetCurve();
	auto color = pcurve->GetColor();


	for (int i = 0; i < curve.size(); i++)
	{
		auto vpos = curve[i];
		for (int j = 0; j < 3; j++)
		{
			vertexs_pos_.push_back(vpos[j]);
			vertex_colors_.push_back(color[j]);
		}
	}









}



CSceneCurveObject::CSceneCurveObject(std::weak_ptr<CCurveObject>  curve)
{
	InitShader();
	this->initializeOpenGLFunctions();
	curve_ = curve;
	curve_.lock().get()->SetChanged(true);
	UpdateRenderInfo();
	curve_.lock().get()->SetChanged(false);

}
CSceneCurveObject::~CSceneCurveObject()
{

}