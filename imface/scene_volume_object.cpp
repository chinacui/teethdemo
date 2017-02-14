#include"scene_volume_object.h"
#include"../DataColle/aux_geo_utils.h"

void CSceneVolumeObject::ComputeRenderingElements(CCamera camera)
{
	qreal pos[3] = { 0,0,0 };
	qreal viewupdir[3] = { 0,1,0 };
	qreal viewrightdir[3] = { 1,0,0 };
	qreal viewdir[3] = { 0,0,-1 };

	qreal wpos[3], wviewupdir[3], wviewrightdir[3], wviewdir[3];

	camera.getWorldCoordinatesOf(pos, wpos);
	camera.getWorldCoordinatesOf(viewupdir, wviewupdir);
	camera.getWorldCoordinatesOf(viewrightdir, wviewrightdir);
	camera.getWorldCoordinatesOf(viewdir, wviewdir);

	for (int i = 0; i < 3; i++)
	{
		wviewupdir[i] -= wpos[i];
		wviewrightdir[i] -= wpos[i];
		wviewdir[i] -= wpos[i];
	}

	/*std::cerr << "cviewdir " << wviewdir[0] << " " << wviewdir[1] << " " << wviewdir[2] << std::endl;
	std::cerr << "cpos " << wpos[0] << " " << wpos[1] << " " << wpos[2] << std::endl;
	std::cerr << "cupvec " << wviewupdir[0] << " " <<wviewupdir[1] << " " << wviewupdir[2] << std::endl;
*/
	auto vviewdir = camera.viewDirection();
	vviewdir.normalize();
	CAuxGeoUtils::GenVolumeRenderAuxCubeAndCuttingPlane(1,1,1, OpenMesh::Vec3d(wpos[0]+ vviewdir[0]*0.01,wpos[1] + vviewdir[1] * 0.01,wpos[2] + vviewdir[2] * 0.01), OpenMesh::Vec3d(wviewdir[0], wviewdir[1], wviewdir[2] ),OpenMesh::Vec3d(wviewupdir[0] ,wviewupdir[1] ,wviewupdir[2]), OpenMesh::Vec3d(wviewrightdir[0], wviewrightdir[1] , wviewrightdir[2] ),aux_cube_, aux_plane_);

	aux_cube_vertexs_pos_.resize(0);
	aux_cube_tex_coords_.resize(0);
	aux_plane_vertexs_pos_.resize(0);
	aux_plane_tex_coords_.resize(0);

	COpenMeshT::FaceIter fiter = aux_cube_.faces_begin();

	for (int i = 0; i < aux_cube_.n_faces(); i++, fiter++)
	{
		COpenMeshT::FaceHalfedgeIter fhi = aux_cube_.fh_begin(*fiter), fhend = aux_cube_.fh_end(*fiter);
		std::vector<COMTraits::Point> vs;

		while (fhi != fhend)
		{
			int vid = aux_cube_.to_vertex_handle(*fhi).idx();
			auto vpos = aux_cube_.points()[vid];

		

			for (int j = 0; j < 3; j++)
			{
				aux_cube_vertexs_pos_.push_back(vpos[j]);
				aux_cube_tex_coords_.push_back(aux_cube_.data(*fhi).GetUV3D()[j]);

			}
			
			vs.push_back(vpos);
			fhi++;
		}
	
	}

	 fiter = aux_plane_.faces_begin();

	for (int i = 0; i < aux_plane_.n_faces(); i++, fiter++)
	{
		COpenMeshT::FaceHalfedgeIter fhi = aux_plane_.fh_begin(*fiter), fhend = aux_plane_.fh_end(*fiter);
		std::vector<COMTraits::Point> vs;

		while (fhi != fhend)
		{
			int vid = aux_plane_.to_vertex_handle(*fhi).idx();
			auto vpos = aux_plane_.points()[vid];



			for (int j = 0; j < 3; j++)
			{
				aux_plane_vertexs_pos_.push_back(vpos[j]);
				aux_plane_tex_coords_.push_back(aux_plane_.data(*fhi).GetUV3D()[j]);

			}

			vs.push_back(vpos);
			fhi++;
		}

	}

	

	


}
void CSceneVolumeObject::SetVolumeDataTextureFromVolumeData()
{
	
	CVolumeDataObject * p_vdata_obj = volume_data_.lock().get();

	ItkVolumeDataType::Pointer pvdata = p_vdata_obj->GetVolumeData();

	auto vdata_size = p_vdata_obj->GetDataSize();
	
	int d = vdata_size[2];
	int w = vdata_size[0];
	int h = vdata_size[1];
	//int d = p_vdata_obj->data_depth_;
	//int w = p_vdata_obj->data_width_;
	//int h = p_vdata_obj->data_height_;
	//int d = 500;
	//int w = 500;
	//int h = 400;
	number_of_slices_ = d;
	slices_over_x_ = std::ceil(std::sqrt(number_of_slices_));
	slices_over_y_ = slices_over_x_;
	int channel_size = 2;
	int data_size = slices_over_x_*w*slices_over_y_*h;
	float*data = new float[data_size * channel_size];
	for (int i = 0; i < data_size * channel_size; i++)
	{
		data[i] = 0;
	}

	int silce_len = w*h; 
	ItkVolumeDataType::PixelType max_v = std::numeric_limits<ItkVolumeDataType::PixelType>::min();
	ItkVolumeDataType::PixelType min_v = std::numeric_limits<ItkVolumeDataType::PixelType>::max();
	for (int i = 0; i < d; i++)
	{
		for (int j = 0; j < w; j++)
		{
			for (int k = 0; k < h; k++)
			{
				ItkVolumeDataType::IndexType pixel_idx;
				pixel_idx[0] = j;
				pixel_idx[1] = k;
				pixel_idx[2] = i;
				ItkVolumeDataType::PixelType pixelvalue = pvdata->GetPixel(pixel_idx);
				if (pixelvalue < 0)
					continue;
				//pixelvalue = std::abs(pixelvalue);
				if (max_v < pixelvalue)
					max_v = pixelvalue;
				if (min_v > pixelvalue)
					min_v = pixelvalue;
			}
		}
	}
	std::cerr << max_v << " " << min_v << std::endl;
	for (int i = 0; i < d; i++)
	{
		for (int j = 0; j < w; j++)
		{
			for (int k = 0; k < h; k++)
			{
				//OpenMesh::Vec4d p = p_vdata_obj->vdata_[i*silce_len + j*h + k];
				ItkVolumeDataType::IndexType pixel_idx;
				pixel_idx[0] = j;
				pixel_idx[1] = k;
				pixel_idx[2] = i;
				ItkVolumeDataType::PixelType pixelvalue=pvdata->GetPixel(pixel_idx);
				//pixelvalue = std::abs(pixelvalue);
				//auto max_pixel_value=std::numeric_limits<ItkVolumeDataType::PixelType>::max();
				//double max_pixel_value = 255.0;
				//max_pixel_value = 1500;
				float fpixelvalue = (pixelvalue - min_v)*1.0 / (max_v - min_v);
				//OpenMesh::Vec4d p = OpenMesh::Vec4d(fpixelvalue, fpixelvalue, fpixelvalue , fpixelvalue );
				if (fpixelvalue < 0)
					fpixelvalue = 0;
				OpenMesh::Vec2d p = OpenMesh::Vec2d(fpixelvalue, fpixelvalue);

				int rowid = i / slices_over_x_*h + k;
				int colid = i%slices_over_x_*w + j;
				data[(rowid*slices_over_x_*w + colid) * channel_size] = p[0];
				data[(rowid*slices_over_x_*w + colid) * channel_size + 1] = p[1];
				//data[(rowid*slices_over_x_*w + colid) * channel_size + 2] = p[2];
				//data[(rowid*slices_over_x_*w + colid) * channel_size + 3] = p[3];
				
			}
		}
	}

	

	glGenTextures(1, &volume_data_texture_id_);

	glBindTexture(GL_TEXTURE_2D, volume_data_texture_id_);


	glTexImage2D(GL_TEXTURE_2D, 0, GL_LUMINANCE_ALPHA, w*slices_over_x_, h*slices_over_y_, 0, GL_LUMINANCE_ALPHA, GL_FLOAT, data);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	delete data;
	


}
void CSceneVolumeObject::SetVolumeRenderBufferFromElements()
{
	volume_rendering_program_.bind();
	if (volume_vao_.isCreated())
		volume_vao_.destroy();
	volume_vao_.create();
	volume_vao_.bind();
	if (volume_aux_cube_vertex_pos_buffer_.isCreated())
		volume_aux_cube_vertex_pos_buffer_.destroy();
	volume_aux_cube_vertex_pos_buffer_.create();
	volume_aux_cube_vertex_pos_buffer_.bind();
	volume_aux_cube_vertex_pos_buffer_.allocate(aux_cube_vertexs_pos_.data(), static_cast<int>(aux_cube_vertexs_pos_.size() * sizeof(float)));
	auto vv_pos_loc = volume_rendering_program_.attributeLocation("a_pos");
	volume_rendering_program_.enableAttributeArray(vv_pos_loc);
	volume_rendering_program_.setAttributeBuffer(vv_pos_loc, GL_FLOAT, 0, 3);
	volume_aux_cube_vertex_pos_buffer_.release();
	
	volume_vao_.release();
	volume_rendering_program_.release();
}
void CSceneVolumeObject::SetBackTexcoordBufferFromElements()
{
	back_texcoord_rendering_program_.bind();
	if (back_texcoord_vao_.isCreated())
		back_texcoord_vao_.destroy();
	back_texcoord_vao_.create();
	back_texcoord_vao_.bind();
	if (back_texcoord_aux_cube_vertex_pos_buffer_.isCreated())
		back_texcoord_aux_cube_vertex_pos_buffer_.destroy();
	back_texcoord_aux_cube_vertex_pos_buffer_.create();
	back_texcoord_aux_cube_vertex_pos_buffer_.bind();
	back_texcoord_aux_cube_vertex_pos_buffer_.allocate(aux_cube_vertexs_pos_.data(), static_cast<int>(aux_cube_vertexs_pos_.size() * sizeof(float)));
	auto bv_pos_loc = back_texcoord_rendering_program_.attributeLocation("a_pos");
	back_texcoord_rendering_program_.enableAttributeArray(bv_pos_loc);
	back_texcoord_rendering_program_.setAttributeBuffer(bv_pos_loc, GL_FLOAT, 0, 3);
	back_texcoord_aux_cube_vertex_pos_buffer_.release();
	if (back_texcoord_aux_cube_tex_coords_buffer_.isCreated())
		back_texcoord_aux_cube_tex_coords_buffer_.destroy();
	back_texcoord_aux_cube_tex_coords_buffer_.create();
	back_texcoord_aux_cube_tex_coords_buffer_.bind();
	back_texcoord_aux_cube_tex_coords_buffer_.allocate(aux_cube_tex_coords_.data(), static_cast<int>(aux_cube_tex_coords_.size() * sizeof(float)));
	auto btex_coords_loc = back_texcoord_rendering_program_.attributeLocation("a_texcoords");
	back_texcoord_rendering_program_.enableAttributeArray(btex_coords_loc);
	back_texcoord_rendering_program_.setAttributeBuffer(btex_coords_loc, GL_FLOAT, 0, 3);
	back_texcoord_aux_cube_tex_coords_buffer_.release();
	back_texcoord_vao_.release();
	back_texcoord_rendering_program_.release();
}
void CSceneVolumeObject::SetFrontTexcoordBufferFromElements()
{
	front_texcoord_rendering_program_.bind();
	if (front_cube_texcoord_vao_.isCreated())
		front_cube_texcoord_vao_.destroy();
	front_cube_texcoord_vao_.create();
	front_cube_texcoord_vao_.bind();
	if (front_texcoord_aux_cube_vertex_pos_buffer_.isCreated())
		front_texcoord_aux_cube_vertex_pos_buffer_.destroy();
	front_texcoord_aux_cube_vertex_pos_buffer_.create();
	front_texcoord_aux_cube_vertex_pos_buffer_.bind();
	front_texcoord_aux_cube_vertex_pos_buffer_.allocate(aux_cube_vertexs_pos_.data(), static_cast<int>(aux_cube_vertexs_pos_.size() * sizeof(float)));
	auto fv_pos_loc = front_texcoord_rendering_program_.attributeLocation("a_pos");
	front_texcoord_rendering_program_.enableAttributeArray(fv_pos_loc);
	front_texcoord_rendering_program_.setAttributeBuffer(fv_pos_loc, GL_FLOAT, 0, 3);
	front_texcoord_aux_cube_vertex_pos_buffer_.release();
	if (front_texcoord_aux_cube_tex_coords_buffer_.isCreated())
		front_texcoord_aux_cube_tex_coords_buffer_.destroy();
	front_texcoord_aux_cube_tex_coords_buffer_.create();
	front_texcoord_aux_cube_tex_coords_buffer_.bind();
	front_texcoord_aux_cube_tex_coords_buffer_.allocate(aux_cube_tex_coords_.data(), static_cast<int>(aux_cube_tex_coords_.size() * sizeof(float)));
	auto ftex_coords_loc = front_texcoord_rendering_program_.attributeLocation("a_texcoords");
	front_texcoord_rendering_program_.enableAttributeArray(ftex_coords_loc);
	front_texcoord_rendering_program_.setAttributeBuffer(ftex_coords_loc, GL_FLOAT, 0, 3);
	front_texcoord_aux_cube_tex_coords_buffer_.release();
	front_cube_texcoord_vao_.release();



	if (front_plane_texcoord_vao_.isCreated())
		front_plane_texcoord_vao_.destroy();
	front_plane_texcoord_vao_.create();
	front_plane_texcoord_vao_.bind();
	if (front_texcoord_aux_plane_vertex_pos_buffer_.isCreated())
		front_texcoord_aux_plane_vertex_pos_buffer_.destroy();
	front_texcoord_aux_plane_vertex_pos_buffer_.create();
	front_texcoord_aux_plane_vertex_pos_buffer_.bind();
	front_texcoord_aux_plane_vertex_pos_buffer_.allocate(aux_plane_vertexs_pos_.data(), static_cast<int>(aux_plane_vertexs_pos_.size() * sizeof(float)));
	fv_pos_loc = front_texcoord_rendering_program_.attributeLocation("a_pos");
	front_texcoord_rendering_program_.enableAttributeArray(fv_pos_loc);
	front_texcoord_rendering_program_.setAttributeBuffer(fv_pos_loc, GL_FLOAT, 0, 3);
	front_texcoord_aux_plane_vertex_pos_buffer_.release();
	if (front_texcoord_aux_plane_tex_coords_buffer_.isCreated())
		front_texcoord_aux_plane_tex_coords_buffer_.destroy();
	front_texcoord_aux_plane_tex_coords_buffer_.create();
	front_texcoord_aux_plane_tex_coords_buffer_.bind();
	front_texcoord_aux_plane_tex_coords_buffer_.allocate(aux_plane_tex_coords_.data(), static_cast<int>(aux_plane_tex_coords_.size() * sizeof(float)));
	ftex_coords_loc = front_texcoord_rendering_program_.attributeLocation("a_texcoords");
	front_texcoord_rendering_program_.enableAttributeArray(ftex_coords_loc);
	front_texcoord_rendering_program_.setAttributeBuffer(ftex_coords_loc, GL_FLOAT, 0, 3);
	front_texcoord_aux_plane_tex_coords_buffer_.release();
	front_plane_texcoord_vao_.release();
	front_texcoord_rendering_program_.release();
}
void CSceneVolumeObject::SetGLBufferDataFromElements()
{
	SetBackTexcoordBufferFromElements();
	SetFrontTexcoordBufferFromElements();
	SetVolumeRenderBufferFromElements();
	
}
void CSceneVolumeObject::InitShader()
{

	std::ifstream vback_texcoord_ifstream("back_texcoord_render.vert");
	std::string vback_texcoord_shader_source((std::istreambuf_iterator<char>(vback_texcoord_ifstream)), std::istreambuf_iterator<char>());
	vback_texcoord_ifstream.close();
	QOpenGLShader *vback_coord_shader = new QOpenGLShader(QOpenGLShader::Vertex);
	if (!vback_coord_shader->compileSourceCode(vback_texcoord_shader_source.c_str()))
	{
		std::cerr << "Compiling back_texcoord_render.vert FAILED" << std::endl;
	}

	std::ifstream fback_texcoord_ifstream("back_texcoord_render.frag");
	std::string fback_texcoord_shader_source((std::istreambuf_iterator<char>(fback_texcoord_ifstream)), std::istreambuf_iterator<char>());
	fback_texcoord_ifstream.close();
	QOpenGLShader *fback_coord_shader = new QOpenGLShader(QOpenGLShader::Fragment);
	if (!fback_coord_shader->compileSourceCode(fback_texcoord_shader_source.c_str()))
	{
		std::cerr << "Compiling back_texcoord_render.frag FAILED" << std::endl;
	}
	if (!back_texcoord_rendering_program_.addShader(vback_coord_shader))
	{
		std::cerr << "adding vback_coord_shader FAILED" << std::endl;
	}
	if (!back_texcoord_rendering_program_.addShader(fback_coord_shader))
	{
		std::cerr << "adding fback_coord_shader FAILED" << std::endl;
	}
	if (!back_texcoord_rendering_program_.link())
	{
		std::cerr << "linking Program FAILED" << std::endl;
	}

////////////////////////////////////////////////////////////
//	std::ifstream vdepth_ifstream("depth_render.vert");
//	std::string vdepth_shader_source((std::istreambuf_iterator<char>(vdepth_ifstream)), std::istreambuf_iterator<char>());
//	vdepth_ifstream.close();
//	QOpenGLShader *vdepth_shader = new QOpenGLShader(QOpenGLShader::Vertex);
//	if (!vdepth_shader->compileSourceCode(vdepth_shader_source.c_str()))
//	{
//		std::cerr << "Compiling depth_render.vert FAILED" << std::endl;
//	}
//
//	std::ifstream fdepth_ifstream("depth_render.frag");
//	std::string fdepth_shader_source((std::istreambuf_iterator<char>(fdepth_ifstream)), std::istreambuf_iterator<char>());
//	fdepth_ifstream.close();
//	QOpenGLShader *fdepth_shader = new QOpenGLShader(QOpenGLShader::Fragment);
//	if (!fdepth_shader->compileSourceCode(fdepth_shader_source.c_str()))
//	{
//		std::cerr << "Compiling depth_render.frag FAILED" << std::endl;
//	}
//
//	if (!depth_rendering_program_.addShader(vdepth_shader))
//	{
//		std::cerr << "adding vdepth_shader FAILED" << std::endl;
//	}
//	if (!depth_rendering_program_.addShader(fdepth_shader))
//	{
//		std::cerr << "adding fdepth_shader FAILED" << std::endl;
//	}
//
//
//////////////////////////////////////////////////////////
	std::ifstream vfront_texcoord_ifstream("front_texcoord_render.vert");
	std::string vfront_texcoord_shader_source((std::istreambuf_iterator<char>(vfront_texcoord_ifstream)), std::istreambuf_iterator<char>());
	vfront_texcoord_ifstream.close();
	QOpenGLShader *vfront_texcoord_shader = new QOpenGLShader(QOpenGLShader::Vertex);
	if (!vfront_texcoord_shader->compileSourceCode(vfront_texcoord_shader_source.c_str()))
	{
		std::cerr << "Compiling front_texcoord_shader.vert FAILED" << std::endl;
	}

	std::ifstream ffront_texcoord_ifstream("front_texcoord_render.frag");
	std::string ffront_texcoord_shader_source((std::istreambuf_iterator<char>(ffront_texcoord_ifstream)), std::istreambuf_iterator<char>());
	ffront_texcoord_ifstream.close();
	QOpenGLShader *ffront_texcoord_shader = new QOpenGLShader(QOpenGLShader::Fragment);
	if (!ffront_texcoord_shader->compileSourceCode(ffront_texcoord_shader_source.c_str()))
	{
		std::cerr << "Compiling back_texcoord_render.frag FAILED" << std::endl;
	}
	if (!front_texcoord_rendering_program_.addShader(vfront_texcoord_shader))
	{
		std::cerr << "adding vfront_coord_shader FAILED" << std::endl;
	}
	if (!front_texcoord_rendering_program_.addShader(ffront_texcoord_shader))
	{
		std::cerr << "adding ffront_coord_shader FAILED" << std::endl;
	}
	if (!front_texcoord_rendering_program_.link())
	{
		std::cerr << "linking Program FAILED" << std::endl;
	}

/////////////////////////////////////////////////////////////////
	std::ifstream vvolume_ifstream("volumeRender.vert");
	std::string vvolume_shader_source((std::istreambuf_iterator<char>(vvolume_ifstream)), std::istreambuf_iterator<char>());
	vvolume_ifstream.close();
	QOpenGLShader *vvolume_shader = new QOpenGLShader(QOpenGLShader::Vertex);
	if (!vvolume_shader->compileSourceCode(vvolume_shader_source.c_str()))
	{
		std::cerr << "Compiling volumeRender.vert FAILED" << std::endl;
	}

	std::ifstream fvolume_ifstream("volumeRender.frag");
	std::string fvolume_shader_source((std::istreambuf_iterator<char>(fvolume_ifstream)), std::istreambuf_iterator<char>());
	fvolume_ifstream.close();
	QOpenGLShader *fvolume_shader = new QOpenGLShader(QOpenGLShader::Fragment);
	if (!fvolume_shader->compileSourceCode(fvolume_shader_source.c_str()))
	{
		std::cerr << "Compiling volumeRender.frag FAILED" << std::endl;
	}
	if (!volume_rendering_program_.addShader(vvolume_shader))
	{
		std::cerr << "adding volumeRender FAILED" << std::endl;
	}
	if (!volume_rendering_program_.addShader(fvolume_shader))
	{
		std::cerr << "adding volumeRender FAILED" << std::endl;
	}
	if (!volume_rendering_program_.link())
	{
		std::cerr << "linking Program FAILED" << std::endl;
	}


	
}

CSceneVolumeObject::CSceneVolumeObject(std::weak_ptr<CVolumeDataObject>  volumedata)
{
	offscreen_frame_size_ = QSize(960, 960);
	back_texcoord_fbo_ = NULL;

	front_texcoord_fbo_ = NULL;
	InitShader();
	GenBackTexcoordFrameBuffer();
	GenFrontTexcoordFrameBuffer();
	
	this->initializeOpenGLFunctions();
	volume_data_ = volumedata;
	volume_data_.lock().get()->SetChanged(true);
	UpdateRenderInfo();
	volume_data_.lock().get()->SetChanged(false);

}
CSceneVolumeObject::~CSceneVolumeObject()
{
	if (back_texcoord_fbo_ != NULL)
	{
		back_texcoord_fbo_->release();
		delete back_texcoord_fbo_;
	}
	if (front_texcoord_fbo_ != NULL)
	{
		front_texcoord_fbo_->release();
		delete front_texcoord_fbo_;
	}
	
}
void CSceneVolumeObject::UpdateRenderInfo()// if the is_changed_ of mesh is seted, the opengl buffer will be updated
{

	
	if (volume_data_.lock().get()->IsChanged())
	{
	
		SetVolumeDataTextureFromVolumeData();
		
		
		volume_data_.lock().get()->SetChanged(false);

	}
}

void CSceneVolumeObject::RenderBackTexcoord(CCamera& camera)
{
	glEnableClientState(GL_VERTEX_ARRAY);
	glEnableClientState(GL_TEXTURE_COORD_ARRAY);
	back_texcoord_fbo_->bind();
	back_texcoord_vao_.bind();
	back_texcoord_rendering_program_.bind();
	//glClearDepth(0);
	glClear( GL_DEPTH_BUFFER_BIT);
	SetCommonRenderAttributes(back_texcoord_rendering_program_, camera);
	//(GL_GREATER);
	glEnable(GL_CULL_FACE);
	glCullFace(GL_FRONT);
	
	glDrawArrays(GL_TRIANGLES, 0, static_cast<GLsizei>(aux_cube_vertexs_pos_.size() / 3));
	glFlush();
	glDisable(GL_CULL_FACE);
	//glDepthFunc(GL_LESS);

	//glClearDepth(0xffffff);

	back_texcoord_rendering_program_.release();
	back_texcoord_vao_.release();
	back_texcoord_fbo_->bindDefault();
	//auto img=back_texcoord_fbo_->toImage();
	//img.save("testimg.bmp");
	
	
}
void CSceneVolumeObject::GenBackTexcoordFrameBuffer()
{
	QOpenGLFramebufferObjectFormat format;
	format.setAttachment(QOpenGLFramebufferObject::CombinedDepthStencil);
	back_texcoord_fbo_=new  QOpenGLFramebufferObject(offscreen_frame_size_, format);
	back_texcoord_fbo_->bindDefault();
}
void CSceneVolumeObject::GenFrontTexcoordFrameBuffer()
{
	QOpenGLFramebufferObjectFormat format;
	format.setAttachment(QOpenGLFramebufferObject::CombinedDepthStencil);
	front_texcoord_fbo_ = new  QOpenGLFramebufferObject(offscreen_frame_size_, format);
	front_texcoord_fbo_->bindDefault();
}
void CSceneVolumeObject::RenderVolumeData(CCamera&camera)
{
	volume_vao_.bind();
	volume_rendering_program_.bind();
	SetCommonRenderAttributes(volume_rendering_program_, camera);

	glActiveTexture(GL_TEXTURE1);
	glBindTexture(GL_TEXTURE_2D, front_texcoord_fbo_->texture());


	glActiveTexture(GL_TEXTURE2);
	glBindTexture(GL_TEXTURE_2D, back_texcoord_fbo_->texture());

	glActiveTexture(GL_TEXTURE3);
	glBindTexture(GL_TEXTURE_2D, volume_data_texture_id_);

	volume_rendering_program_.setUniformValue("u_front_texcoord_map",1);
	volume_rendering_program_.setUniformValue("u_back_texcoord_map",2);
	volume_rendering_program_.setUniformValue("u_volume_data", 3);

	volume_rendering_program_.setUniformValue("u_number_of_slices", number_of_slices_*1.0f);
	volume_rendering_program_.setUniformValue("u_slices_over_x", slices_over_x_*1.0f);
	volume_rendering_program_.setUniformValue("u_slices_over_y", slices_over_y_*1.0f);

	volume_rendering_program_.setUniformValue("u_render_as_image",0);
	volume_rendering_program_.setUniformValue("u_with_env_depth", 0);
	glDisable(GL_CULL_FACE);
	glDrawArrays(GL_TRIANGLES, 0, static_cast<GLsizei>(aux_cube_vertexs_pos_.size() / 3));
	glFlush();
	

	
	
	volume_rendering_program_.release();
	volume_vao_.release();

}
void CSceneVolumeObject::RenderFrontTexcoord(CCamera& camera)
{
	front_texcoord_fbo_->bind();
	front_texcoord_rendering_program_.bind();

	front_cube_texcoord_vao_.bind();
	glClearDepth(0);
	glClear(GL_DEPTH_BUFFER_BIT);
	glDepthFunc(GL_GREATER);

	SetCommonRenderAttributes(front_texcoord_rendering_program_, camera);
	glEnable(GL_CULL_FACE);
	glCullFace(GL_BACK);
	glDrawArrays(GL_TRIANGLES, 0, static_cast<GLsizei>(aux_cube_vertexs_pos_.size() / 3));
	glFlush();
	glDisable(GL_CULL_FACE);
	front_cube_texcoord_vao_.release();
	
	front_plane_texcoord_vao_.bind();
	SetCommonRenderAttributes(front_texcoord_rendering_program_, camera);
	glDepthFunc(GL_GREATER);
	//glClear(GL_DEPTH_BUFFER_BIT);
	glDrawArrays(GL_TRIANGLES, 0, static_cast<GLsizei>(aux_plane_vertexs_pos_.size() / 3));
	glFlush();
	glDepthFunc(GL_LESS);
	
	front_plane_texcoord_vao_.release();
	glClearDepth(0xfffff);
	//

	//glDepthFunc(GL_GREATER);
	//glEnable(GL_CULL_FACE);
	//glCullFace(GL_BACK);
	//glDisable(GL_CULL_FACE);
	//glDepthFunc(GL_LESS);


	//glClearDepth(0xffffff);
	front_texcoord_rendering_program_.release();
	front_texcoord_fbo_->bindDefault();
}
void CSceneVolumeObject::SetCommonRenderAttributes(QOpenGLShaderProgram &program, CCamera& camera)
{
	//compute attributes
	QMatrix4x4 mvpMatrix;
	QMatrix4x4 mvMatrix, frame_matrix, scale_matrix;
	double mat[16];
	camera.getModelViewProjectionMatrix(mat);

	auto local_mat = volume_data_.lock().get()->GetMatrix();
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




	int mvp_mat_loc = program.uniformLocation("mvp_matrix");
	int mv_mat_loc = program.uniformLocation("mv_matrix");


	program.setUniformValue(mvp_mat_loc, mvpMatrix);
	program.setUniformValue(mv_mat_loc, mvMatrix);

}
void CSceneVolumeObject::Render(CCamera camera)
{
	ComputeRenderingElements(camera);
	SetGLBufferDataFromElements();
	UpdateRenderInfo();
	
	RenderBackTexcoord(camera);
	RenderFrontTexcoord(camera);

	RenderVolumeData(camera);
	
	
	
	
}