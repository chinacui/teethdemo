#ifndef CPANORAMIC_SIMULATION_TEST_H
#define CPANORAMIC_SIMULATION_TEST_H
#include"action_base.h"
#include<map>
#include"../DataColle/mesh_object.h"
#include"../AlgColle/arap_deform.h"
#include"../TeethRootRecoAlg/panoramic_fitting.h"



class CPanoramicSimulationTestAction :public CActionBase
{
protected:
	std::map<int, CMeshObject*>crowns_;
	std::vector<CMeshObject *>segmented_jaws_ ;
	std::map<int, std::vector<int>>jaw_crowns_;
	OpenMesh::Vec3d  center_;
	std::vector<OpenMesh::Vec3d>frame_;
	std::vector<OpenMesh::Vec2d>pick_pts_on_panorama_;
	std::vector<std::pair<int,COpenMeshT::VertexHandle>>pick_vhs_on_crown_;

	std::map<int, std::vector<std::pair<COpenMeshT::VertexHandle,OpenMesh::Vec2d>>>proj_coords_;
	OpenMesh::Vec3d param_rendering_offset_;
	
	double rot_degree_step_;
	double center_move_step_;
	double radius_step_;
	double radius_;
	int panoramic_pic_id_;
	bool is_selecting_pair_;
	void Init();
	void MousePressEvent(QMouseEvent *e);
	void MouseMoveEvent(QMouseEvent *e);
	void MouseReleaseEvent(QMouseEvent *e);
	void KeyPressEvent(QKeyEvent *e);
	void KeyReleaseEvent(QKeyEvent *e);

	//size of frame equal to 3
	void RotateFrame(std::vector<OpenMesh::Vec3d>&frame,OpenMesh::Vec3d axis,double degree);
	void RenderAuxilliaryShape();
	void ComputeFrameAndCenter();
	void GenPanorama();

	void SavePoints();
public:
	CPanoramicSimulationTestAction();

};
#endif 
