#ifndef CCNON_RIGID_ICP_H
#define CCNON_RIGID_ICP_h
#include"prereq.h"
#include <Eigen/Sparse>
#include<ANN/ANN.h>
#include"../DataColle/mesh_object.h"
 /* Options : structured object with fields:
%       gamm : real valued, weights differences in the rotational and skew
%           part of the deformation against the translational part.
%       epsilon : real values, tolerence for change in transformation.
%       lambda : If using the bi-directional distance metric this weights
%           the contribution of the target -> source term.
%       alphaSet : decreasing vector of real-valued stiffness parameters. 
%           High stiffness parameters force global transformations whereas 
%           low values allow for local deformations.
%       biDirectional : logical, specifies that a bi-directional distance 
%           is used.
%       useNormals : logical, specifies that surface normals are to be used
%           to project the source onto the target surface. If this term is
%           used then the Source input should contain a normals field.
%       plot : logical, specifies that the transformations should be
%           plotted.
%       rigidInit : logical, specifies that rigid ICP should be performed
%           first before allowing non-rigid and non-global deformations.*/
class ALGCOLLE_CLASS CNonRigidICPParams
{
public:
	bool use_normals_=true;
	bool rigid_init_ = true;
	bool normal_weighting_ = false;
	bool bi_directional_ = false;//
	double gama_ = 0.1;//weight differences in the rotational and skew part of the deformation against the translational part of the deformation
	std::vector<double>stiffness_set_;//different stiffness
	std::vector<double>iters_each_stiffness_;
	double epsilon_ = 1e-4;
	bool project_to_tgt_;
	CNonRigidICPParams();

};
class ALGCOLLE_CLASS CNonRigidICP
{
protected:
	COpenMeshT* src_mesh_,* tgt_mesh_;
	CNonRigidICPParams params_;
	Eigen::MatrixXd src_vertexs_;
	
	Eigen::MatrixXd src_vertex_normals_;
	Eigen::MatrixXi src_faces_;
	Eigen::MatrixXd tgt_vertexs_;
	std::vector<bool>tgt_is_border_vertexs_;
	Eigen::MatrixXi tgt_faces_;
	Eigen::MatrixXd tgt_vertex_normals_;
	Eigen::MatrixXd tgt_sample_vert_;
	std::vector<double>stiffness_;//
	int current_iter_ = 0;
	Eigen::Matrix4d G_;//matrix in equ3, weighting matrix of stiffness term
	Eigen::SparseMatrix<double>M_;//arc_node incidence matrix
	Eigen::SparseMatrix<double>kron_M_G_;//kroneckerProduct of M_ and G_ in equ 10 
	Eigen::SparseMatrix<double>D_;//the sparse matrixDmapping the 4n¡Á3 matrix of unknowns X onto displaced source vertices, equ 8
	Eigen::SparseMatrix<double>W_;//weight matrix
	Eigen::SparseMatrix<double,Eigen::RowMajor>A_;
	Eigen::MatrixXd B_;
	Eigen::MatrixXd U_;
	Eigen::VectorXd W_diag_;//diag of W
	Eigen::MatrixXd X_;//transform matrix
	Eigen::MatrixXd transformed_;
	ANNpointArray tgt_ann_pts_,src_ann_pts_;
	ANNkd_tree* tgt_kd_tree_;
	ANNkd_tree*src_kd_tree_;
	std::map<COpenMeshT::VertexHandle, COpenMeshT::FaceHandle>vh_fh_map_;//tgt to src map;
	void SampleVerts(Eigen::MatrixXd& vertexs, Eigen::MatrixXi & faces, double radius,Eigen::MatrixXd &res_sample_verts);
	void ComputeArcNodeMatrix(Eigen::MatrixXd& vertexs, Eigen::MatrixXi & faces, Eigen::SparseMatrix<double> &res_incident_mat);
	void InitKdTree(Eigen::MatrixXd &pts, ANNpointArray* res_pts, ANNkd_tree** res_kd_tree);
	int NNSearch(ANNkd_tree *kd_tree, Eigen::Vector3d p, double &res_d);
	void UpdateSrcMesh();
	

	bool CNonRigidICP::WriteDMAT(
		 std::string file_name,
		 Eigen::MatrixXd & W);
	bool CNonRigidICP::WriteDMAT(
		std::string file_name,
		Eigen::Matrix4d & W);
	bool CNonRigidICP::WriteDMAT(
		std::string file_name,
		Eigen::MatrixXi & W);
	void Init();
public:
	//CNonRigidICP(Eigen::MatrixXd& src_vertexs, Eigen::MatrixXi& src_faces, Eigen::MatrixXd& tgt_vertexs, Eigen::MatrixXi& tgt_faces, CNonRigidICPParams params);
	CNonRigidICP(COpenMeshT *src_mesh, COpenMeshT* tgt_mesh, CNonRigidICPParams params=CNonRigidICPParams());
	bool CTemplateToDentalmesh(COpenMeshT *src_mesh, COpenMeshT* tgt_mesh);
	void GetTgt2SrcMap(std::map<COpenMeshT::VertexHandle, COpenMeshT::FaceHandle>&vh_fh_map);
	bool Run(int step=-1);//return true if finished
	void GetTransformedSourceVertexs(Eigen::MatrixXd &res_transformed);
};
#endif