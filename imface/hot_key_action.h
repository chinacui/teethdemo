#ifndef CHOT_KEY_ACTION_H
#define CHOT_KEY_ACTION_H
#include"action_base.h"
class CModelViewer;
class CHotKeyAction:public CActionBase
{

protected:
	
	void KeyPressEvent(QKeyEvent *e);
	void KeyReleaseEvent(QKeyEvent *e);
public:
	CHotKeyAction() { type_ = Common; }
};
#endif