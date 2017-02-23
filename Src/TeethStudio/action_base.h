#ifndef CACTION_BASE_H
#define CACTION_BASE_H
#include<iostream>

#include<qevent.h>

#include"qobject.h"
#include"qmetaobject.h"
class CModelViewer;
class CActionManager;
enum CActionType { Common,EditFeatureEdge,HarmonicFieldSegmentation,VolumeDataSegmentation, CTeethReconstructionTest, CTeethReconstruction,CPanoramicSimulationTest};
class CActionBase :public QObject
{
	Q_OBJECT
protected:
	bool is_done_ = true;
	CModelViewer *viewer_;
	CActionManager*manager_;
	CActionType type_;
	bool shielding_key_modifiers_ = true;
	bool coexist_with_common_tool_ = false;
public:

	virtual void SetViewer(CModelViewer * v)
	{
		viewer_ = v;
	}
	virtual void SetManager(CActionManager *m)
	{
		manager_ = m;
	}
	virtual CActionType GetType() { return type_; }
	

	
	virtual void MousePressEvent(QMouseEvent *e) { return; }
	virtual void MouseMoveEvent(QMouseEvent *e) { return; }
	virtual void MouseReleaseEvent(QMouseEvent *e) { return; }
	virtual void KeyPressEvent(QKeyEvent *e) { return; }
	virtual void KeyReleaseEvent(QKeyEvent *e) { return; }
	virtual void WheelEvent(QWheelEvent *) { return; }
	
	virtual void MouseDoubleClickEvent(QMouseEvent * e) { return ; }
	virtual void Init() {  }



	

};
#endif

