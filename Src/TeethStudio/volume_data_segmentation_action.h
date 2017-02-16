#ifndef CVOLUME_DATA_SEGMENTATION_ACTION_H
#define CVOLUME_DATA_SEGMENTATION_ACTION_H
#include"action_base.h"
#include "../DataColle/custom_openmesh_type.h"
#include"../DataColle/mesh_object.h"
#include"../DataColle/curve_object.h"
#include"../DataColle/volume_data_object.h"
class CVolumeDataSegmentationAction :public CActionBase
{
protected:

protected:
	void Init();
	void MousePressEvent(QMouseEvent *e);
	void MouseMoveEvent(QMouseEvent *e);
	void MouseReleaseEvent(QMouseEvent *e);
	void KeyPressEvent(QKeyEvent *e);
	void KeyReleaseEvent(QKeyEvent *e);


public:
	CVolumeDataSegmentationAction();
};

#endif