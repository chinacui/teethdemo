#ifndef CTEETH_RECONSTRUCTIONTEST_ACTION_H
#define CTEETH_RECONSTRUCTIONTEST_ACTION_H
#include"action_base.h"
#include<map>
#include"../DataColle/mesh_object.h"
#include"../AlgColle/arap_deform.h"
class CNonRigidICP;
class CTeethReconstructionTestAction:public CActionBase
{
protected:
	bool is_picking_;
	bool is_picking_fa_point_;
	std::map<int,CMeshObject*>crowns_;
	std::map<int,CMeshObject*>template_tooth_;
	std::map<int,CMeshObject*>crown_of_template_tooth_;
	std::map<int,std::map<COpenMeshT::FaceHandle, COpenMeshT::FaceHandle>>crown_to_template_fmap_;
	std::map<int, std::map<COpenMeshT::VertexHandle, COpenMeshT::VertexHandle>>crown_to_template_vmap_;
	CMeshObject *segmented_crown_=NULL;
	CMeshObject* orig_tooth_template_=NULL;
	int sel_crown_id_, sel_temp_id_,sel_curve_id_,sel_temp_crown_id_;
	CARAPDeformation* template_arap_=NULL;
	CNonRigidICP* non_rigid_icp_=NULL;
	std::vector<std::pair<COpenMeshT::FaceHandle, OpenMesh::Vec3d>>test_sel_ps;
	std::map<int, COpenMeshT::VertexHandle>fa_point_map_;

	std::map<int, OpenMesh::Vec3d>long_axis_;
	std::map<int, OpenMesh::Vec3d>crown_centers_;
	std::map<int, Eigen::Matrix4d>crown_mats_;
	bool IsTemplateTeeth(int id);
	void Init();
	void MousePressEvent(QMouseEvent *e);
	void MouseMoveEvent(QMouseEvent *e);
	void MouseReleaseEvent(QMouseEvent *e);
	void KeyPressEvent(QKeyEvent *e);
	void KeyReleaseEvent(QKeyEvent *e);
public:
	CTeethReconstructionTestAction();

};
#endif 
