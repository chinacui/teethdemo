#ifndef CTEXTURE_WRAPPER_H
#define CTEXTURE_WRAPPER_H
#include <QOpenGLFunctions_2_1>
#include<qopengltexture.h>
class TextureWraper :public QOpenGLFunctions_2_1
{
	QImage img_;
	QOpenGLTexture* texture_ = NULL;
	bool is_updated_ = false;
	int texture_id_ = 0;
	static int max_texture_id_;
public:

	QImage & GetImage() { return img_; };
	QOpenGLTexture *GetTexture() { return texture_; };
	void ReleaseTexture();
	TextureWraper(QImage &img) :img_(img) { is_updated_ = false; texture_id_ = max_texture_id_++; };
	void SetImage(QImage* img);
	void UpdateTexture();
	int GetTextureId();
	bool IsUpdated() { return is_updated_; }
};
#endif
