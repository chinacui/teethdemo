#ifndef MANIPULATION_ACTION_H
#define MANIPULATION_ACTION_H
#include"action_base.h"
#include "../DataColle/custom_openmesh_type.h"
class CManipulationAction : public CActionBase
{
protected:
	std::vector<int>sel_mesh_ids_;
	bool is_moving_mesh_ = false;
	OpenMesh::Vec3d pre_move_pos_;
	OpenMesh::Vec3d sel_mesh_center_;
protected:

	void MousePressEvent(QMouseEvent *e);
	void MouseMoveEvent(QMouseEvent *e);
	void MouseReleaseEvent(QMouseEvent *e);
	void KeyPressEvent(QKeyEvent *e);
	void KeyReleaseEvent(QKeyEvent *e);
	void WheelEvent(QWheelEvent *);
public:
	CManipulationAction() { type_ = Common; }
};
#endif