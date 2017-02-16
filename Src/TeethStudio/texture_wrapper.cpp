#include"texture_wrapper.h"
int TextureWraper::max_texture_id_ = 0;
void TextureWraper::UpdateTexture()
{
	//if (texture_ != NULL)
	//	texture_->release();
	//else
	if (texture_ == NULL)
		texture_ = new QOpenGLTexture(QOpenGLTexture::Target::Target2D);
	/*texture_->setda
	texture_->destroy();*/
	texture_->setData(img_.mirrored()); 
	//texture_->setWrapMode(QOpenGLTexture::WrapMode::ClampToEdge);

	is_updated_ = true;
}
void TextureWraper::SetImage(QImage* img)
{
	img_ = *img;
	//	if (texture_ != NULL)
	//texture_->release();
	//delete texture_;
	//texture_ = NULL; 
	is_updated_ = false;
};


int TextureWraper::GetTextureId()
{
	return texture_id_;
}
void TextureWraper::ReleaseTexture()
{
	if (texture_ != NULL)
		texture_->release();
	delete texture_;
	texture_ = NULL;
}