// Copyright 2016_9 by ChenNenglun
#include "ui_context.h"
#include "scene.h"
#include "imface_window.h"
#include "texture_wrapper.h"
CScene* CUIContext::scene_ = NULL;
CMainWindow* CUIContext::main_window_ = NULL;
int CUIContext::selected_mesh_object_id_ = -1;
CMorphSkelDentalMeshSeg* CUIContext::msdm_seg_ = NULL;
std::string CUIContext::cm_current_adjust_param_= "Curvature";
int CUIContext::color_stripe_texture_id_ = -1;
int CUIContext::color_bar_texture_id_ = -1;
std::map<int, TextureWraper*>CUIContext::textures_ = std::map<int, TextureWraper*>();
int CUIContext::max_textures_id_ = 0;

CScene* CUIContext::GetScene()
{
	return scene_;
}

void CUIContext::SetSelectedMeshObjectId(int id)
{
	selected_mesh_object_id_ = id;
}

int CUIContext::GetSelectedMeshObjectId()
{
	if (DataPool::GetMeshObject(selected_mesh_object_id_) == NULL)
	{
		selected_mesh_object_id_ = -1;
	}
	return selected_mesh_object_id_;
}

void CUIContext::Init()
{
	scene_ = new CScene();
	main_window_ = new CMainWindow();
	main_window_->show();

	/*QImage color_stripe("stripe2.bmp");
	color_stripe_texture_id_ = CUIContext::AddTexture(new QImage(color_stripe));

	QImage color_bar("color.bmp");
	color_bar_texture_id_ = CUIContext::AddTexture(new QImage(color_bar));*/
}

int CUIContext::AddTexture(QImage *img)
{
	TextureWraper *tw = new TextureWraper(*img);

	int id = tw->GetTextureId();
	textures_.insert(std::make_pair(id, tw));

	return id;
}

bool CUIContext::DeleteTextures(int id)
{
	auto iter = textures_.find(id);
	if (iter != textures_.end())
	{
		iter->second->ReleaseTexture();
		textures_.erase(iter);
		return true;
	}
	return false;
}

TextureWraper* CUIContext::GetTexture(int id)
{
	auto iter = textures_.find(id);
	if (iter != textures_.end())
	{
		return iter->second;
	}
	return NULL;
}
