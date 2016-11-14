#ifndef CHARMONIC_FIELD_SEGMENTATION_ACTION_H
#define CHARMONIC_FIELD_SEGMENTATION_ACTION_H
#include"action_base.h"
#include "../DataColle/custom_openmesh_type.h"
#include"../DataColle/mesh_object.h"
#include"../DataColle/curve_object.h"
class CHarmonicFieldSegmentation :public CActionBase
{
protected:
	std::vector<COpenMeshT::VertexHandle>picked_vhs_fore_,picked_vhs_back_;
	std::set<COpenMeshT::VertexHandle>picked_vhs_fore_mark_, picked_vhs_back_mark_;

	std::set<OpenMesh::Vec3d>rand_color_set_;
	std::vector<int>teeth_seg_mark_;//mark of each vertex
	std::vector<int>gingiva_seg_mark_;
	int teeth_seg_mark_id_;
	std::map<int,int>teeth_seg_count_;
	std::map<int,OpenMesh::Vec3d>teeth_seg_color_;//color of each segmented teeth
	bool is_pick_fore_ = true;
	CMeshObject *p_mesh_obj_;
	CCurveObject *p_curve_obj_=NULL;
	bool is_eliminating_feature_ = false;
	bool is_seg_teeth_from_gingiva_ = false;
	bool is_drawing_ = false;
	OpenMesh::Vec3d pca_mean_;
	std::vector<OpenMesh::Vec3d>pca_frame_;
	OpenMesh::Vec3d GetRandColor();
	void RenderFeature();
protected:
	void Init();
	void MousePressEvent(QMouseEvent *e);
	void MouseMoveEvent(QMouseEvent *e);
	void MouseReleaseEvent(QMouseEvent *e);
	void KeyPressEvent(QKeyEvent *e);
	void KeyReleaseEvent(QKeyEvent *e);
public:
	CHarmonicFieldSegmentation();
};
#endif