// Copyright 2016_9 by ChenNenglun
#include"ui_context.h"
#include"scene.h"
#include"imface_window.h"
CScene* CUIContext::scene_=NULL;
CImfaceWindow* CUIContext::main_window_ = NULL;
int CUIContext::selected_mesh_object_id_ = -1;
CMorphSkelDentalMeshSeg* CUIContext::msdm_seg_ = NULL;
std::string CUIContext::cm_current_adjust_param_= "Curvature";
CScene* CUIContext::GetScene()
{
	return scene_;
}
void CUIContext::SetSelectedMeshObjectId(int id)
{
	selected_mesh_object_id_=id;
}
int CUIContext::GetSelectedMeshObjectId()
{
	return selected_mesh_object_id_;
}
void CUIContext::Init()
{
	scene_ = new CScene();
	main_window_ = new CImfaceWindow();
	main_window_->show();
}