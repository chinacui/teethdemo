#ifndef CACTION_MANAGER_H
#define CACTION_MANAGER_H
#include<qevent.h>
#include"action_base.h"
class CActionManager
{
protected:

	std::vector<CActionBase*>actions_;

	CActionType current_action_type_ = Common;
	CModelViewer* viewer_=NULL;

	
	
public:

	virtual void MousePressEvent(QMouseEvent *e);
	virtual void MouseMoveEvent(QMouseEvent *e);
	virtual void MouseReleaseEvent(QMouseEvent *e);
	virtual void KeyPressEvent(QKeyEvent *e);
	virtual void KeyReleaseEvent(QKeyEvent *e);
	virtual void WheelEvent(QWheelEvent *);

	virtual void MouseDoubleClickEvent(QMouseEvent * e);

	CActionType GetCurrentActionType() { return current_action_type_; }
	void SetCurrentActionType(CActionType t);
	void Init(CModelViewer* viewer);
	CActionManager() {};


	
};
#endif