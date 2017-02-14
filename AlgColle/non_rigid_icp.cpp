#include"non_rigid_icp.h"
#include<igl/per_vertex_normals.h>
#include<time.h>
#include<unsupported/Eigen/src/KroneckerProduct/KroneckerTensorProduct.h>
#include <ANN/ANN.h>
#include<math.h>
#include"../AdditionalLibs/libicp/libicp_wrapper.h"
#include"../DataColle/cgal_igl_converter.h"
#include<igl/writeDMAT.h>
#include"../AlgColle/linear_algebra_alg.h"
#include<Eigen/IterativeLinearSolvers>
#include"../AdditionalLibs/obb/cpqp_obb_wrapper.h"
CNonRigidICPParams::CNonRigidICPParams()
{
	stiffness_set_.resize(0);
	iters_each_stiffness_.resize(0);
	double step = 10;
	for (int i = 41; i >= 3; i -= step)
	{
		stiffness_set_.push_back(i);
		iters_each_stiffness_.push_back(6);
		step--;
	}
	iters_each_stiffness_.back() = 5;
	project_to_tgt_ = true;
}
//CNonRigidICP::CNonRigidICP(Eigen::MatrixXd& src_vertexs, Eigen::MatrixXi& src_faces, Eigen::MatrixXd& tgt_vertexs, Eigen::MatrixXi& tgt_faces, CNonRigidICPParams params)
//{
//	src_vertexs_ = src_vertexs;
//	src_faces_ = src_faces;
//	tgt_vertexs_ = tgt_vertexs;
//	tgt_faces_ = tgt_faces;
//	params_ = params;
//	if (params_.normal_weighting_)
//	{
//		igl::per_vertex_normals(src_vertexs_, src_faces_, src_vertex_normals_);
//		igl::per_vertex_normals(tgt_vertexs_, tgt_faces_, tgt_vertex_normals_);
//	}
//	InitKdTree(tgt_vertexs_, &tgt_ann_pts_, &tgt_kd_tree_);
//	
//}
CNonRigidICP::CNonRigidICP(COpenMeshT *src_mesh, COpenMeshT* tgt_mesh, CNonRigidICPParams params )
{
	src_mesh_ = src_mesh;
	tgt_mesh_ = tgt_mesh;
	CConverter::ConvertFromOpenMeshToIGL(*src_mesh, src_vertexs_, src_faces_);
	CConverter::ConvertFromOpenMeshToIGL(*tgt_mesh, tgt_vertexs_, tgt_faces_);
	tgt_is_border_vertexs_.resize(tgt_mesh->n_vertices(), 0);
	for (auto viter = tgt_mesh->vertices_begin(); viter != tgt_mesh->vertices_end(); viter++)
	{
		tgt_is_border_vertexs_[viter->idx()] = tgt_mesh->is_boundary(viter);
		tgt_is_border_vertexs_[viter->idx()] = false;
	}
	params_ = params;
	if (params_.normal_weighting_)
	{
		igl::per_vertex_normals(src_vertexs_, src_faces_, src_vertex_normals_);
		igl::per_vertex_normals(tgt_vertexs_, tgt_faces_, tgt_vertex_normals_);
	} 
	tgt_ann_pts_ = NULL;
	src_ann_pts_ = NULL;
	tgt_kd_tree_ = NULL;
	src_kd_tree_ = NULL;
	InitKdTree(tgt_vertexs_, &tgt_ann_pts_, &tgt_kd_tree_);
	

	Init();
}
void CNonRigidICP::InitKdTree(Eigen::MatrixXd &pts, ANNpointArray* res_pts, ANNkd_tree** res_kd_tree)
{
	if(*res_pts==NULL)
	*res_pts = annAllocPts(pts.rows(), 3);
	for (int i = 0; i < pts.rows(); i++)
	{
		(*res_pts)[i][0] = pts(i, 0);
		(*res_pts)[i][1] = pts(i, 1);
		(*res_pts)[i][2] = pts(i, 2);
	}

	if (*res_kd_tree != NULL)
	{
		delete *res_kd_tree;
	}
	*res_kd_tree =new ANNkd_tree(*res_pts, pts.rows(), 3);
	
//	double d;
	//NNSearch(*res_kd_tree, Eigen::Vector3d(0, 0, 0), d);

}

void CNonRigidICP::UpdateSrcMesh()
{
	for (auto viter = src_mesh_->vertices_begin(); viter != src_mesh_->vertices_end(); viter++)
	{
		int id = viter->idx();
		OpenMesh::Vec3d v(transformed_(id, 0), transformed_(id, 1), transformed_(id, 2));
		src_mesh_->set_point(viter, v);
	}
}
int CNonRigidICP::NNSearch(ANNkd_tree *kd_tree, Eigen::Vector3d p, double &res_d)
{
	ANNpoint ap = annAllocPt(3);
	ap[0] = p[0];
	ap[1] = p[1];
	ap[2] = p[2];
	ANNidxArray res_ids= new ANNidx[1];
	ANNdistArray res_diss = new ANNdist[1];
	kd_tree->annkSearch(ap, 1, res_ids, res_diss);
	res_d = res_diss[0];
	return res_ids[0];
}
void CNonRigidICP::GetTransformedSourceVertexs(Eigen::MatrixXd &res_transformed)
{
	res_transformed=transformed_ ;
}

bool CNonRigidICP::WriteDMAT(
	 std::string file_name,
	 Eigen::MatrixXd & W)
{
	FILE * fp = fopen(file_name.c_str(), "w");
	if (fp == NULL)
	{
		fprintf(stderr, "IOError: writeDMAT() could not open %s...", file_name.c_str());
		return false;
	}

		// first line contains number of rows and number of columns
		//fprintf(fp, "%d %d\n", (int)W.cols(), (int)W.rows());
		// Loop over columns slowly
		for (int i = 0; i < W.rows(); i++)
		{
			// loop over rows (down columns) quickly
			for (int j = 0; j < W.cols(); j++)
			{
				fprintf(fp, "%0.17lg", (double)W(i, j));
				if (j != W.cols() - 1)
					fprintf(fp, " ");
			}
			fprintf(fp, "\n");
		}


	fclose(fp);
	return true;
}
bool CNonRigidICP::WriteDMAT(
	std::string file_name,
	Eigen::MatrixXi & W)
{
	FILE * fp = fopen(file_name.c_str(), "w");
	if (fp == NULL)
	{
		fprintf(stderr, "IOError: writeDMAT() could not open %s...", file_name.c_str());
		return false;
	}

	// first line contains number of rows and number of columns
	//fprintf(fp, "%d %d\n", (int)W.cols(), (int)W.rows());
	// Loop over columns slowly
	for (int i = 0; i < W.rows(); i++)
	{
		// loop over rows (down columns) quickly
		for (int j = 0; j < W.cols(); j++)
		{
			fprintf(fp, "%0.17lg", (double)W(i, j));
			if (j != W.cols() - 1)
				fprintf(fp, " ");
		}
		fprintf(fp, "\n");
	}


	fclose(fp);
	return true;
}
bool CNonRigidICP::WriteDMAT(
	std::string file_name,
	Eigen::Matrix4d & W)
{
	FILE * fp = fopen(file_name.c_str(), "w");
	if (fp == NULL)
	{
		fprintf(stderr, "IOError: writeDMAT() could not open %s...", file_name.c_str());
		return false;
	}

	// first line contains number of rows and number of columns
	//fprintf(fp, "%d %d\n", (int)W.cols(), (int)W.rows());
	// Loop over columns slowly
	for (int i = 0; i < W.rows(); i++)
	{
		// loop over rows (down columns) quickly
		for (int j = 0; j < W.cols(); j++)
		{
			fprintf(fp, "%0.17lg", (double)W(i, j));
			if (j != W.cols() - 1)
				fprintf(fp, " ");
		}
		fprintf(fp, "\n");
	}


	fclose(fp);
	return true;
}
void CNonRigidICP::Init()
{
	if (params_.bi_directional_)
	{
		SampleVerts(tgt_vertexs_, tgt_faces_, 15, tgt_sample_vert_);
	}
	G_.setIdentity();
	G_(3, 3) = params_.gama_;
	ComputeArcNodeMatrix(src_vertexs_, src_faces_, M_);
	kron_M_G_ = Eigen::kroneckerProduct(M_, G_);
	std::cerr << "m size " << M_.rows() << " " << M_.cols() << std::endl;
	std::cerr << "g size " << G_.rows() << " " << G_.cols() << std::endl;
	std::cerr << "kron mg size " << kron_M_G_.rows() << " " << kron_M_G_.cols() << std::endl;
	D_.resize(src_vertexs_.rows(), src_vertexs_.rows() * 4);
	std::vector < Eigen::Triplet<double>>trip_d(0);

	for (int i = 0; i < src_vertexs_.rows(); i++)
	{
		for (int j = 0; j < 3; j++)
		{
			trip_d.push_back(Eigen::Triplet<double>(i, i * 4 + j, src_vertexs_(i, j)));
		}
		trip_d.push_back(Eigen::Triplet<double>(i, i * 4 + 3, 1));
	}


	D_.setFromTriplets(trip_d.begin(), trip_d.end());

	W_diag_.resize(src_vertexs_.rows());
	W_diag_.setOnes();

	if (params_.rigid_init_)
	{
		//perform icp
		CIcpPoint2Point icp(tgt_vertexs_);
		Eigen::Matrix3d rot;
		Eigen::Vector3d trans;
		icp.Fit(src_vertexs_, rot, trans);
		X_.resize(src_vertexs_.rows() * 4, 3);
		X_.setZero();
		rot.transposeInPlace();
		for (int i = 0; i < src_vertexs_.rows(); i++)
		{
			for (int j = 0; j < 3; j++)
			{
				X_(i * 4 + j, 0) = rot(j, 0);
				X_(i * 4 + j, 1) = rot(j, 1);
				X_(i * 4 + j, 2) = rot(j, 2);
			}
			X_(i * 4 + 3, 0) = trans(0);
			X_(i * 4 + 3, 1) = trans(1);
			X_(i * 4 + 3, 2) = trans(2);
		}
		transformed_ = (D_*X_).eval();
		//UpdateSrcMesh();
		//return;

		/*std::string fname("X_rigid");

		WriteDMAT(fname, X_);*/
	}
	else
	{
		X_.resize(src_vertexs_.rows() * 4, 3);
		X_.setZero();

		for (int i = 0; i < src_vertexs_.rows(); i++)
		{

			X_(i * 4, 0) = 1;
			X_(i * 4 + 1, 1) = 1;
			X_(i * 4 + 2, 2) = 1;
		}
		//std::string fname("X_init");
		//WriteDMAT(fname, X_);
		transformed_ = (D_*X_).eval();
	}
	stiffness_.clear();
	for (int i = 0; i < params_.stiffness_set_.size(); i++)
	{
		for (int j = 0; j < params_.iters_each_stiffness_[i]; j++)
		{
			stiffness_.push_back(params_.stiffness_set_[i]);
		}
		
	}
	UpdateSrcMesh();
}
bool CNonRigidICP::Run(int step)
{


	//perform nonrigid icp
	//for (int i = 0; i < params_.stiffness_set_.size(); i++)
	for(int i= current_iter_;i<current_iter_+step;i++)
	{
	//	double cstiffness = params_.stiffness_set_[i];
		if (i == stiffness_.size())
		{
			current_iter_ = i;
			break;
		}
		double cstiffness = stiffness_[i];
		std::cerr <<"iter "<<i<< " stiffness " << cstiffness <<" "<<i<< std::endl;
		Eigen::MatrixXd old_X = (10.0*X_).eval();
		/*int count = 0;
		double error = ((X_ - old_X).eval()).norm();
		while (error> params_.epsilon_)
		{*/

	
			
			std::map<int, double>src_tgt_dis;
			//transformed_ = (D_*X_).eval();
			Eigen::VectorXi tids(transformed_.rows());
			U_.resize(transformed_.rows(), 3);

			InitKdTree(transformed_, &src_ann_pts_, &src_kd_tree_);
			for (int j = 0; j < transformed_.rows(); j++)
			{
				Eigen::Vector3d p(transformed_(j, 0), transformed_(j, 1), transformed_(j, 2));
				double dis;
				tids(j)=NNSearch(tgt_kd_tree_, p, dis);
				bool flag = false;
				if ((src_tgt_dis.find(tids(j)) == src_tgt_dis.end()|| src_tgt_dis[tids(j)] > dis)&&(tgt_is_border_vertexs_[tids(j)] == false))
				{


					
					Eigen::Vector3d tgt_p(tgt_vertexs_(tids(j), 0), tgt_vertexs_(tids(j), 1), tgt_vertexs_(tids(j), 2));
					double inv_dis;
					NNSearch(src_kd_tree_, tgt_p, inv_dis);

					if (inv_dis > dis / 2)
					//if(true)
					{
						src_tgt_dis[tids(j)] = dis;
						U_(j, 0) = tgt_vertexs_(tids(j), 0);
						U_(j, 1) = tgt_vertexs_(tids(j), 1);
						U_(j, 2) = tgt_vertexs_(tids(j), 2);
					}
					else
					{
						flag = true;
					}
					
				}
				else
				{
					flag = true;
				}
				
				if(flag)
				{
					U_(j, 0) = transformed_(j, 0);
					U_(j, 1) = transformed_(j, 1);
					U_(j, 2) = transformed_(j, 2);
				}
				
			
			}

			if (params_.normal_weighting_)
			{
				Eigen::SparseMatrix<double> N;
				N.resize(src_vertexs_.rows(), src_vertexs_.rows() * 4);
				std::vector<Eigen::Triplet<double>>trip_N(0);
				for (int j = 0; j < src_vertexs_.rows(); j++)
				{
					for (int k = 0; k < 3; k++)
					{
						trip_N.push_back(Eigen::Triplet<double>(j, 4 * j + k, src_vertex_normals_(j, k)));
					}
					trip_N.push_back(Eigen::Triplet<double>(j, 4 * j + 3, 1));
				}
				N.setFromTriplets(trip_N.begin(), trip_N.end());
				Eigen::MatrixXd normals_transformed = N*X_;
				Eigen::Vector3d cor_normals_target;
				
				
				for (int j = 0; j < tids.size();j++)
				{
					if (tgt_is_border_vertexs_[tids(j)])
					{
						W_diag_(j) =0.3;
						continue;
					}
					cor_normals_target(0) = tgt_vertex_normals_(tids[j], 0);
					cor_normals_target(1) = tgt_vertex_normals_(tids[j], 1);
					cor_normals_target(2) = tgt_vertex_normals_(tids[j], 2);
					Eigen::Vector3d nt(normals_transformed(j, 0), normals_transformed(j, 1), normals_transformed(j, 2));
					Eigen::Vector3d tmpcross = cor_normals_target.cross(nt);
					double cross_norm = tmpcross.norm();
					double dot_normal = cor_normals_target.dot(nt);
					double angle = std::atan2(cross_norm, dot_normal);
					W_diag_(j) = W_diag_(j)*(angle < M_PI / 4 ? 1 : 0.5);
				}
			
				
			}
		
			std::vector<Eigen::Triplet<double>>trip_W(W_diag_.size());

			for (int j = 0; j < trip_W.size(); j++)
			{
				trip_W[j]=(Eigen::Triplet<double>(j, j, W_diag_(j)));
			}
			W_.resize(W_diag_.size(), W_diag_.size());
			W_.setFromTriplets(trip_W.begin(),trip_W.end());

			Eigen::SparseMatrix<double>W_D = (W_*D_).eval();
			Eigen::SparseMatrix<double>tmp_kron_M_G = (cstiffness*kron_M_G_).eval();
			
			A_.resize(tmp_kron_M_G.rows() + W_D.rows(), tmp_kron_M_G.cols());
			A_.reserve(W_D.nonZeros() + tmp_kron_M_G.nonZeros());
			A_.topRows(tmp_kron_M_G.rows()) = tmp_kron_M_G;
			A_.bottomRows(W_D.rows()) = W_D.eval();
		
		
		//	std::string fname("X_");
			//fname = fname + std::to_string(i) + " " + std::to_string(count);
			//WriteDMAT(fname, X_);

			Eigen::MatrixXd W_U = (W_*U_).eval();
			B_.resize(M_.rows()*G_.rows()+ W_U.rows(), 3);
			B_.setZero();
			B_.bottomRows(W_U.rows()) = W_U.eval();
		
			old_X = X_;
			Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double, Eigen::ColMajor>> solver;
			//Eigen::SparseLU<Eigen::SparseMatrix<double, Eigen::ColMajor>> solver;
			Eigen::SparseMatrix<double>AT = A_.transpose();
			Eigen::SparseMatrix<double>ATA = (AT*A_);
			ATA.makeCompressed();
		
			solver.compute(A_);
			

			Eigen::MatrixXd BT = (AT*B_);
			
			X_ = solver.solve(B_).eval();
			
			

			double trans_diff_error = ((X_ - old_X).eval()).norm();
			std::cerr << "trans diff " << trans_diff_error << std::endl;
		/*	count++;
			if (count >10)
			{
				break;
			}*/
		

		//}
			transformed_ = (D_*X_).eval();
			double error = ((transformed_ - U_).eval()).norm();
			std::cerr << "error " << error << std::endl;
	}
	
	if (current_iter_>= stiffness_.size()&&params_.project_to_tgt_)
	{
		CPqpObb obbsrc(transformed_, src_faces_);
		CPqpObb obbtgt(tgt_vertexs_, tgt_faces_);
		//InitKdTree(transformed_, &src_ann_pts_, &src_kd_tree_);
		vh_fh_map_.clear();
		for (int i = 0; i < transformed_.rows(); i++)
		{
		
				Eigen::Vector3d p(transformed_(i, 0), transformed_(i, 1), transformed_(i, 2));
				double dis;
				CPqpObb::QueryResult tgt_clos_res=obbtgt.QueryClosestPoint(p);
				dis=tgt_clos_res.distance;
				Eigen::Vector3d tgt_p=tgt_clos_res.closest_pnt_;
				CPqpObb::QueryResult src_clos_res = obbsrc.QueryClosestPoint(tgt_p);
				double inv_dis= src_clos_res.distance;
				
				if (inv_dis > dis / 3.0*2.0)
				{
					transformed_.row(i) = tgt_p;
					//std::cerr << "project " << transformed_.row(i) << " " << tgt_p << std::endl;
					vh_fh_map_[src_mesh_->vertex_handle(i)] = tgt_mesh_->face_handle(tgt_clos_res.fid_);	
				}

		}
		UpdateSrcMesh();
		return true;
	}
	else
	{
		current_iter_ = current_iter_ + step;
		UpdateSrcMesh();
		return false;
	}
	
	


	
}
void CNonRigidICP::GetTgt2SrcMap(std::map<COpenMeshT::VertexHandle, COpenMeshT::FaceHandle>&vh_fh_map)
{
	vh_fh_map = vh_fh_map_;
}
void CNonRigidICP::ComputeArcNodeMatrix(Eigen::MatrixXd& vertexs, Eigen::MatrixXi & faces, Eigen::SparseMatrix<double> &res_incident_mat)
{
	
	std::vector<Eigen::Triplet<double>>trips(0);
	int eid = 0;
	std::set<long long>Eset;
	for (int i = 0; i <faces.rows(); i++)
	{
		for (int j = 0; j < 3; j++)
		{
			int va = faces(i, j);
			int vb = faces(i, (j + 1) % 3);
			if (va > vb)
			{
				int tmp = va;
				va = vb;
				vb = tmp;
			}
			long long tmpid = va*vertexs.rows() + vb;
			if (Eset.find(tmpid) != Eset.end())
				continue;
			Eset.insert(tmpid);
			
			trips.push_back(Eigen::Triplet<double>(eid, va, -1));
			trips.push_back(Eigen::Triplet<double>(eid, vb, 1));
			eid++;
			
		}
	}
	res_incident_mat.resize(eid, vertexs.rows());
	res_incident_mat.setFromTriplets(trips.begin(),trips.end());

}
void CNonRigidICP::SampleVerts(Eigen::MatrixXd& vertexs, Eigen::MatrixXi & faces, double radius,Eigen::MatrixXd &res_sample_verts)
{
	srand((unsigned)time(NULL));

	std::vector<bool>is_sampled(vertexs.rows(),false);
	int num_unsampled_vert = vertexs.rows();
	std::vector<int>sample_vids(0);
	std::vector<int>unsampled_vids(vertexs.rows());
	//this sampling method maybe slow, to see if there exists other methods
	for (int i = 0; i < unsampled_vids.size(); i++)
		unsampled_vids[i] = i;
	while (num_unsampled_vert > 0)
	{
		int sample_vid= unsampled_vids[rand() % unsampled_vids.size()];
		sample_vids.push_back(sample_vid);
		Eigen::Vector3d pv(vertexs(sample_vid, 0), vertexs(sample_vid, 1), vertexs(sample_vid, 2));
		std::vector<int>to_be_removed(0);
		for (int i = 0; i < unsampled_vids.size(); i++)
		{
			Eigen::Vector3d tv(vertexs(unsampled_vids[i], 0), vertexs(unsampled_vids[i], 1), vertexs(unsampled_vids[i], 2));
			double dis=(pv - tv).norm();
			if (dis < radius)
			{
				to_be_removed.push_back(i);
			}
		}
		for (int i = 0; i < to_be_removed.size(); i++)
		{
			int l = to_be_removed[i];
			int r = (i == to_be_removed.size() - 1 ? unsampled_vids.size() - 1 : to_be_removed[i + 1]-1);
			for (int j = l + 1; j <= r; j++)
			{
				unsampled_vids[j - 1] = unsampled_vids[j];
			}
		}
		int count = to_be_removed.size();
		while (count--)
		{
			unsampled_vids.pop_back();
		}
		
	}
	res_sample_verts.resize(sample_vids.size(),3);
	for (int i = 0; i < sample_vids.size(); i++)
	{
		res_sample_verts(i, 0) = vertexs(sample_vids[i], 0);
		res_sample_verts(i, 1) = vertexs(sample_vids[i], 1);
		res_sample_verts(i, 2) = vertexs(sample_vids[i], 2);
	}


}