// Copyright 2016_9 by ChenNenglun
#ifndef UI_CONTEXT_H
#define UI_CONTEXT_H
#include"../AlgColle/morph_skel_dental_mesh_seg.h"
class CScene;
class CImfaceWindow;
class TextureWraper;
class QImage;
class CUIContext
{
protected:
	static CScene* scene_;
	static CImfaceWindow* main_window_;
	static int selected_mesh_object_id_;
	static int color_bar_texture_id_;
public:
	static CScene* GetScene();
	static void Init();
	static int GetSelectedMeshObjectId();
	static void SetSelectedMeshObjectId(int id);
	static int &ColorBarTextureId() { return color_bar_texture_id_; }



	static std::map<int, TextureWraper*>textures_;
	static int max_textures_id_;
	static int AddTexture(QImage *texture);
	static bool DeleteTextures(int id);
	static TextureWraper* GetTexture(int id);

	/////////////for debug
	static CMorphSkelDentalMeshSeg* msdm_seg_;
	static std::string cm_current_adjust_param_;

};
#endif