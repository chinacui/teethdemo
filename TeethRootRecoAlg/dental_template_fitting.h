#ifndef CDENTAL_TEMPLATE_FITTING_H
#define CDENTAL_TEMPLATE_FITTING_H


#include"prereq.h"
#include"../DataColle/mesh_object.h"
#include"../DataColle/teeth_template_object.h"
#include<map>
class TEETHROOTRECOALG_CLASS CDentalTemplateFitting
{
public:
	static void RefineFittingTemplate(COpenMeshT &mesh_template, COpenMeshT &mesh_target);
	static void RefineCrownBoundaryAfterFitting(CMeshObject &mesh_crown, CMeshObject &teeth_template, std::map<COpenMeshT::VertexHandle, COpenMeshT::FaceHandle>&temp_crown_map);
	static void ReplaceTemplateCrown(CMeshObject &teeth_template, CMeshObject &mesh_crown);
	static void AlignTemplate2Crown(CMeshObject &crown, OpenMesh::Vec3d front_dir,OpenMesh::Vec3d up_dir,CTeethTemplateObject& teeth_template);
	static void MergeBoundaryRegisteredMeshes(COpenMeshT &mesh_a, Eigen::Matrix4d mat_a, COpenMeshT&mesh_b, Eigen::Matrix4d mat_b,std::map<COpenMeshT::VertexHandle,COpenMeshT::VertexHandle>&bound_map_a2b,COpenMeshT&res_mesh,std::vector<COpenMeshT::VertexHandle>&res_seam_vhs=std::vector<COpenMeshT::VertexHandle>());
	static void LoadTemplates(std::string dirname, std::map<std::string, CTeethTemplateObject*>&res_templates);
	static void ComputeStretchCrowns2LineMatrix(std::vector<CMeshObject*>&crowns, std::vector<Eigen::Matrix4d>&res_mats);
	static void ComputeCrownFrontDirs(std::vector<CMeshObject*>&crowns, std::vector<OpenMesh::Vec3d>&res_front_dirs);
	static void OrderCrowns(std::vector<CMeshObject*>&crowns,OpenMesh::Vec3d crown_updir, bool from_left2right);
	static void ComputeCrownsType(std::vector<CMeshObject*>crowns,std::map<CMeshObject*,int>&fixed_id, std::map<CMeshObject *, std::string>&res_crown_type);
};
#endif