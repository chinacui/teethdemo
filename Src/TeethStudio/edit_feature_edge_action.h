#ifndef EDIT_FEATURE_EDGE_TOOL_H
#define EDIT_FEATURE_EDGE_TOOL_H
#include"action_base.h"
#include "../DataColle/custom_openmesh_type.h"
class CEditFeatureEdgeAction:public CActionBase
{
protected:
	std::vector<COpenMeshT::VertexHandle>picked_vhs_;
	bool is_pick_ = false;
protected:

	void MousePressEvent(QMouseEvent *e);
	void MouseMoveEvent(QMouseEvent *e);
	void MouseReleaseEvent(QMouseEvent *e);
	void KeyPressEvent(QKeyEvent *e);
	void KeyReleaseEvent(QKeyEvent *e);
public:
	CEditFeatureEdgeAction() { type_ = EditFeatureEdge; }
};
#endif