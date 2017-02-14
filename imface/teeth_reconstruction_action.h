#ifndef CTEETH_RECONSTRUCTION_ACTION_H
#define CTEETH_RECONSTRUCTION_ACTION_H
#include"action_base.h"
#include<map>
#include"../DataColle/mesh_object.h"
#include"../AlgColle/arap_deform.h"
#include"../DataColle/teeth_template_object.h"
#include"../AlgColle/arap_deform.h"
#include"../AlgColle/non_rigid_icp.h"
class CTeethReconstructionAction :public CActionBase
{
protected:
	std::vector<int>crown_ids_;
	std::map<std::string, CTeethTemplateObject*>template_meshes_;
	std::map<int, std::string>crown2temptype_map_;
	std::map<int, CTeethTemplateObject*>crown_temp_map_;
	std::map<int,CARAPDeformation*> crown_temp_arap_ ;
	std::map<int,CNonRigidICP*> non_rigid_icp_ ;
	std::map<int, bool>is_temp_matching_finished_;

	bool is_selecting_teeth_ = false;
	std::set<int>sel_teeth_ids_;
	void Init();
	void MousePressEvent(QMouseEvent *e);
	void MouseMoveEvent(QMouseEvent *e);
	void MouseReleaseEvent(QMouseEvent *e);
	void KeyPressEvent(QKeyEvent *e);
	void KeyReleaseEvent(QKeyEvent *e);
	void QuitAction();
public:
	CTeethReconstructionAction();

};

#endif