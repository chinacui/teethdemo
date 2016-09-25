#include"ui_context.h"
#include"scene.h"
#include"imface_window.h"
CScene* CUIContext::scene_=NULL;
CImfaceWindow* CUIContext::main_window_ = NULL;
CScene* CUIContext::GetScene()
{
	return scene_;
}
void CUIContext::Init()
{
	scene_ = new CScene();
	main_window_ = new CImfaceWindow();
	main_window_->show();
}