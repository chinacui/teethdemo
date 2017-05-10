#ifndef SINGLE_TEETH_PROJECTION_ACTION
#define SINGLE_TEETH_PROJECTION_ACTION
#include"action_base.h"
#include "../DataColle/mesh_object.h"
#include <Eigen\Dense>
#include <iostream>
#include "ui_imface.h"
#include <QWidget>
#include<qpushbutton.h>
#include "cmodelviewer.h"
#include "../TeethRootRecoAlg/teethProjection.h"
class QLineEdit;


class CSingleTeethProjectionAction:public CActionBase
{
	Q_OBJECT
public:
	CSingleTeethProjectionAction();
protected:
	void Init();
	void KeyPressEvent(QKeyEvent *e);
	void calProjectDir(int teeth_select_id, std::vector<OpenMesh::Vec3d> teethCenter, std::vector<OpenMesh::Vec3d> poly_teeth_cruve);
	void calProjectPlane(int teeth_select_id, std::vector<OpenMesh::Vec3d> teethCenter, std::vector<int> tagTeeth_num, COpenMeshT mesh, std::map<int, int> tags);
	void calToothEdge(std::map<int, std::vector<bool>>&judge_volume_edge_temp, std::vector<OpenMesh::Vec2d> projection_plane, int teeth_select_id, std::vector<std::vector<int>>vids);
	void getCrownRootPoint(int i);
	void bubblesort(int count);
	bool LoadRawData(std::string fname, int i);
	void segTeethToIndividual(COpenMeshT &mesh, std::vector<int>&v_tag);
	void segTemplateCrownRoot(std::map<int, COpenMeshT*>&temp_crown_root_mesh, std::map<int, std::map<COpenMeshT::VertexHandle, COpenMeshT::VertexHandle>>&temp_crown_root_map); //segment the template mesh to crown mesh and root mesh
	void crownRegistration(int id);
	void tempPositionInit(COpenMeshT &template_mesh, COpenMeshT& tooth_mesh, int id);
	void topPointHarmonic(COpenMeshT &template_mesh_crown, std::map<int, std::map<COpenMeshT::VertexHandle, COpenMeshT::VertexHandle>>&temp_crown_root_map, int id);
	bool TempToMeshDisplacement(COpenMeshT &template_mesh, COpenMeshT& tooth_mesh, int id);
	bool meshSmooth(COpenMeshT &template_mesh);
	bool handlePointBoundary(COpenMeshT &template_mesh, int teeth_id);
public:
	void templateIdNum(int id_1_, int id_2_, int id_3_, int id_4_, int id_5_, int id_6_, int id_7_);


public:
	static CSingleTeethProjectionAction* GetInstance();
	static void DeleteInstance();
	std::vector<OpenMesh::Vec2d>projection_Plane;
	std::map<int, std::vector<bool>> judge_volume_edge;
	//std::vector<bool> judge_plane_edge;
	std::vector<OpenMesh::Vec3d> edge_Point;
	OpenMesh::Vec3d single_Teeth_center;
	OpenMesh::Vec3d projection_dir;
	std::vector<std::vector<int>>vids;

	std::vector<OpenMesh::Vec3d> project_dir; //the projection direction of one tooth
	std::vector<OpenMesh::Vec3d> project_tangent;
	std::vector<OpenMesh::Vec3d> project_z;
	OpenMesh::Vec3d project_dir_temp;
	OpenMesh::Vec3d project_tangent_temp;
	OpenMesh::Vec3d project_z_temp;
	
	std::vector<OpenMesh::Vec3d> teethCenter;
	std::vector<OpenMesh::Vec3d> teethPoints;
	std::vector<int> tagTeeth_num;
	std::vector<int> tagTeeth;

	std::vector<std::map<int, COpenMeshT>> teeth_data_;
	std::vector<int> tag_mesh;
	std::map<int, COpenMeshT*>sep_meshes;
	std::vector<OpenMesh::Vec2d>projection_Plane_front;
	std::map<int, std::vector<bool>> judge_volume_edge_front;
	std::vector<OpenMesh::Vec2d>projection_Plane_back;
	std::map<int, std::vector<bool>> judge_volume_edge_back;
	std::map<int, OpenMesh::Vec3d> edge_crown_root_top_point;
	CSingleTeethProjectionAction *SingleTeethProjectionAction;
	CMeshObject *mesh_template = new CMeshObject();
	std::map<int, std::vector<OpenMesh::Vec3d>> id_edge_crown_root_Point;
	int template_id;
	std::map<int, CMeshObject*> template_map;
	std::map<int, std::vector<OpenMesh::Vec2d>> projection_teeth_image_;
	std::map<int, OpenMesh::Vec3d> template_crown_displacement_;
private:
	static CSingleTeethProjectionAction* instance_;
    teethProjection * teeth_projection_;
	static int id_1_num_;
	static int id_2_num_;
	static int id_3_num_;
	static int id_4_num_;
	static int id_5_num_;
	static int id_6_num_;
	static int id_7_num_;
	//CMainWindow * aasa;
};
#endif
/*#ifndef SINGLE_TEETH_PROJECTION_ACTION
#define SINGLE_TEETH_PROJECTION_ACTION
#include"action_base.h"
#include "../DataColle/mesh_object.h"
#include <Eigen\Dense>
#include <iostream>

class CSingleTeethProjectionAction:public CActionBase
{
protected:
	void Init();
	void KeyPressEvent(QKeyEvent *e);
	void calProjectDir(int teeth_select_id, std::vector<OpenMesh::Vec3d> teethCenter, std::vector<OpenMesh::Vec3d> poly_teeth_cruve);
	void calProjectPlane(int teeth_select_id, std::vector<OpenMesh::Vec3d> teethCenter, std::vector<int> tagTeeth_num, COpenMeshT mesh, std::map<int, int> tags);
	void calToothEdge(std::map<int, std::vector<bool>>&judge_volume_edge_temp, std::vector<OpenMesh::Vec2d> projection_plane, int teeth_select_id, std::vector<std::vector<int>>vids);
	void getCrownRootPoint(int i);
	void bubblesort(int count);
	bool LoadRawData(std::string fname);
	void crownRegistration(int id);
	void tempPositionInit(COpenMeshT &template_mesh, COpenMeshT& tooth_mesh, int id);
	void topPointHarmonic(COpenMeshT &template_mesh_crown, std::map<int, std::map<COpenMeshT::VertexHandle, COpenMeshT::VertexHandle>>&temp_crown_root_map, int id);
	void handlePointBoundary(COpenMeshT &template_mesh_crown,int teeth_id);
	void boundPointHarmonic(COpenMeshT &template_mesh);
	void segTemplateCrownRoot(std::map<int, COpenMeshT*>&temp_crown_root_mesh, std::map<int, std::map<COpenMeshT::VertexHandle, COpenMeshT::VertexHandle>>&temp_crown_root_map); //segment the template mesh to crown mesh and root mesh
	void segTeethToIndividual(COpenMeshT &mesh, std::vector<int>&v_tag);
public:
	CSingleTeethProjectionAction();
	std::vector<OpenMesh::Vec2d>projection_Plane_front;
	std::vector<OpenMesh::Vec2d>projection_Plane_back;
	std::map<int, std::vector<bool>> judge_volume_edge;
	std::map<int, std::vector<bool>> judge_volume_edge_front;
	std::map<int, std::vector<bool>> judge_volume_edge_back;
	//std::vector<bool> judge_plane_edge;
	std::vector<OpenMesh::Vec3d> edge_Point;
	OpenMesh::Vec3d single_Teeth_center;
	OpenMesh::Vec3d projection_dir;
	std::vector<std::vector<int>>vids;

	std::vector<std::vector<int>>temp_vids;

	std::vector<OpenMesh::Vec3d> project_dir; //the projection direction of one tooth
	std::vector<OpenMesh::Vec3d> project_tangent;
	std::vector<OpenMesh::Vec3d> project_z;
	OpenMesh::Vec3d project_dir_temp;
	OpenMesh::Vec3d project_tangent_temp;
	OpenMesh::Vec3d project_z_temp;
	
	std::vector<OpenMesh::Vec3d> teethCenter;
	std::vector<OpenMesh::Vec3d> teethPoints;
	std::vector<int> tagTeeth_num;
	std::vector<int> tagTeeth;
	std::map<int, std::vector<OpenMesh::Vec3d>> id_edge_crown_root_Point;
	std::map<int, OpenMesh::Vec3d> handle_point_displacement;

	std::vector<std::map<int, COpenMeshT>> teeth_data_;
	/*harmonic deformation
	//Eigen::MatrixXd V, U, V_bc, U_bc;
CSingleTeethProjectionAction *SingleTeethProjectionAction;

CMeshObject *mesh_template = new CMeshObject();
std::vector<int> tag_mesh;

std::map<int, OpenMesh::Vec3d> edge_crown_root_top_point;

std::map<int, COpenMeshT*>sep_meshes;
};
#endif*/