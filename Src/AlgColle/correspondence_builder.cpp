#include "correspondence_builder.h"
#include"geo_base_alg.h"
#include"geo_alg.h"
#include "../DataColle/curve_object.h"

#include "curve_base_alg.h"
#include "hmm.h"

CCorrespondenceBuilder::CCorrespondenceBuilder()
{
	sigma_p_ = 0.25;
	sigma_n_ = 0.8;
	sigma_c_ = 0.20;
}

std::map<COpenMeshT::VertexHandle, OpenMesh::Vec3d> CCorrespondenceBuilder::FindCorrespondenceMap(CMeshObject* p_polyobj, std::vector<OpenMesh::Vec3d> p_curveobj, std::vector<COpenMeshT::VertexHandle> roi_region)
{
	p_polyobj->ApplyTransform();
	p_polyobj->SetChanged();

	p_polyobj_ = p_polyobj;
	p_curveobj_ = p_curveobj;
	if (roi_region.empty())
	{
		CGeoBaseAlg::ComputePerVertexNormal(p_polyobj->GetMesh(), roi_vertices_normal_map_);
		
	}
	else
	{
		std::map<COpenMeshT::VertexHandle,OpenMesh::Vec3d> poly_normal_map;
		CGeoBaseAlg::ComputePerVertexNormal(p_polyobj->GetMesh(), poly_normal_map);
		roi_vertices_normal_map_.clear();
		for (auto p = roi_region.begin(); p != roi_region.end(); p++)
		{
			roi_vertices_normal_map_[*p] =poly_normal_map[*p];
		}
	}
	
	std::vector<OpenMesh::Vec3d> backup_curve = p_curveobj_;
	//std::cout << "size of curve" << p_curveobj_->curve_.size() << std::endl;

	int sample_method = 1;
	if (sample_method == 1)
	{
		//////////////////////////////////////////////////////////////////////////

		double cl=CCurveBaseAlg::ComputeLenOfCurve(p_curveobj_, true);
	
		double l = CGeoBaseAlg::ComputeAverageEdgeLength(p_polyobj_->GetMesh());
		std::cerr << "ave mesh edge len " << l << std::endl;
		int sample_num = cl / l;
		std::cout << "Arc length sampling; num of sample: " << sample_num << std::endl;

		CCurveBaseAlg::ResampleCurve(p_curveobj_, 2*sample_num, true);
		std::cerr << "size of p_curveobject_ " << p_curveobj_.size() << std::endl;
		//////////////////////////////////////////////////////////////////////////
	}
	/*else
	{
		CCurveBaseAlg::CurveResampling(p_curveobj_);
		CCurveBaseAlg::SubdivideCurve(p_curveobj_->curve_, CPolyhedronBaseAlg::AverageEdgeLengthOfPolyhedron(p_polyobj_->poly_));
		
	}*/
	//std::cout << "size of curve" << p_curveobj_->curve_.size() << std::endl;




	CCurveBaseAlg::ComputeCurveNormals(p_curveobj_, curve_normal_vec_);
	
	ComputeEmissionProbabilities();
	ComputeTransitionProbabilities();
	std::cerr << "start hmm matching " << std::endl;
	if (!HMMMatching())
	{
		std::cout << "HMM Matching error!" << std::endl;
	}
	std::cerr << "hmm finished" << std::endl;
	COpenMeshT&mesh = p_polyobj_->GetMesh();
	for (auto p = deform_handle_map_.begin(); p != deform_handle_map_.end(); p++)
	{
		p->second = OpenMesh::Vec3d(p->second[0], p->second[1], mesh.point(p->first)[2]);
	}

	p_curveobj_.clear();
	p_curveobj_ = backup_curve;

	return deform_handle_map_;

}

std::map<COpenMeshT::VertexHandle, OpenMesh::Vec3d> CCorrespondenceBuilder::FindCorrespondenceMapNormal(CMeshObject* p_polyobj, std::vector<OpenMesh::Vec3d>& p_curveobj, std::vector<COpenMeshT::VertexHandle> roi_region)
{
	p_polyobj->ApplyTransform();
	p_polyobj->SetChanged();
	p_polyobj_ = p_polyobj;
	p_curveobj_ = p_curveobj;
	if (roi_region.empty())
	{
		CGeoBaseAlg::ComputePerVertexNormal(p_polyobj->GetMesh(), roi_vertices_normal_map_);
	}
	else
	{

		std::map<COpenMeshT::VertexHandle, OpenMesh::Vec3d> poly_normal_map;
		CGeoBaseAlg::ComputePerVertexNormal(p_polyobj->GetMesh(), poly_normal_map);
		roi_vertices_normal_map_.clear();
		for (auto p = roi_region.begin(); p != roi_region.end(); p++)
		{
			roi_vertices_normal_map_[*p] = poly_normal_map[*p];
		}
	}
	CCurveBaseAlg::ComputeCurveNormals(p_curveobj_, curve_normal_vec_);

	std::cout << "size of normal in curve" << curve_normal_vec_.size() << std::endl;

	//for (int i = 0; i < curve_normal_vec_.size(); i++)
	//{
	//	std::cout << curve_normal_vec_[i] << std::endl;
	//}
	ComputeEmissionProbabilities();
	ComputeTransitionProbabilities();
	if (!HMMMatching())
	{
		std::cout << "HMM Matching error!" << std::endl;
	}
	COpenMeshT&mesh = p_polyobj_->GetMesh();
	for (auto p = deform_handle_map_.begin(); p != deform_handle_map_.end(); p++)
	{
		p->second = OpenMesh::Vec3d(p->second[0], p->second[1], mesh.point(p->first)[2]);
	}
	return deform_handle_map_;
}

void CCorrespondenceBuilder::ComputeEmissionProbabilities()
{
	size_t nv = roi_vertices_normal_map_.size();
	size_t np = p_curveobj_.size();
	/*emission*/
	ep_matrix_.clear();
	ep_matrix_.resize(nv);
	//Polyhedron::Vertex_iterator viter = p_polyobject_->poly_.vertices_begin();
	auto mapiter = roi_vertices_normal_map_.begin();
	COpenMeshT &mesh = p_polyobj_->GetMesh();
	for (int vi = 0; vi < nv; ++vi, ++mapiter) //for every vertex in candidate
	{
		ep_matrix_[vi].clear();
		ep_matrix_[vi].resize(np);
		for (int pi = 0; pi < np; ++pi) // for every curve point
		{
			//Default z-axis is zero, for generalization, a plane projection should be added here!!!
			OpenMesh::Vec3d point_v(mesh.point(mapiter->first)[0], mesh.point(mapiter->first)[1], 0);
			OpenMesh::Vec3d point_p(p_curveobj_[pi][0], p_curveobj_[pi][1], 0);
			double dp = (point_v - point_p).sqrnorm();
			double dn = OpenMesh::dot((mapiter->second) ,curve_normal_vec_[pi]);

			// YL modification
			OpenMesh::Vec3d v_vp = (point_p - point_v) / sqrt(dp);
			OpenMesh::Vec3d n_xy(mapiter->second[0], mapiter->second[1], 0);
			n_xy = n_xy / n_xy.length();
			double dn2 = std::fabs(OpenMesh::dot(n_xy , v_vp));
			
			//double exp_pjvi = dp / sigma_p_ / sigma_p_
			//+(dn - 1.0)*(dn - 1.0) / sigma_n_ / sigma_n_;

			// the first term is reduced to dp (rather than dp * dp)
		/*	double exp_pjvi = dp / sigma_p_ / sigma_p_
				+ (dn2 - 1.0)*(dn2 - 1.0) / sigma_n_ / sigma_n_;*/
			
			// Original 
			double exp_pjvi = dp * dp / sigma_p_ / sigma_p_ + (dn - 1.0)*(dn - 1.0) / sigma_n_ / sigma_n_; 
			ep_matrix_[vi][pi] = -exp_pjvi;
		}
	}
}

void CCorrespondenceBuilder::ComputeTransitionProbabilities()
{
	size_t nv = roi_vertices_normal_map_.size();
	size_t np = p_curveobj_.size();
	/*transition*/
	vdist_matrix_.clear();
	vdist_matrix_.resize(nv);
	auto mapiter = roi_vertices_normal_map_.begin();
	COpenMeshT &mesh = p_polyobj_->GetMesh();
	for (int vi1 = 0; vi1 < nv; vi1++, mapiter++) //for every vertex in candidate
	{
		vdist_matrix_[vi1].clear(); vdist_matrix_[vi1].resize(nv);
		OpenMesh::Vec3d p1, p2;
		//for (int j = 0; j<3; j++)
		//	v1[j] = _mesh->vertices[candIdx[vi1]][j];
		p1 = mesh.point(mapiter->first);
		auto mapiter2 = roi_vertices_normal_map_.begin();
		for (int vi2 = 0; vi2 < nv; ++vi2, ++mapiter2) //for every vertex in candidate
		{
			p2 = mesh.point(mapiter2->first);
			vdist_matrix_[vi1][vi2] = (p1 - p2).sqrnorm();
		}
	}
	pdist_vec_.clear();
	for (int i = 0; i < np - 1; i++)
	{
		pdist_vec_.push_back((p_curveobj_[i + 1] - p_curveobj_[i]).sqrnorm());
	}
	//pdist_vec_.push_back((curve_point_vec_[np - 1] - curve_point_vec_[np - 2]).squared_length());
}

bool CCorrespondenceBuilder::HMMMatching()
{
	COpenMeshT&mesh = p_polyobj_->GetMesh();
	size_t nv = roi_vertices_normal_map_.size();
	size_t np = p_curveobj_.size();
	//init HMM
	HMM<int> hmm;
	std::vector<double> init_probs(nv);  //init probs has same size as hidden state, or say vertex
	for (int i = 0; i < nv; ++i)
	{
		init_probs[i] = ep_matrix_[i][0];
	}
	hmm.Init((int)nv, (int)np, init_probs, &vdist_matrix_, &ep_matrix_);

	///init Viterbi triple
	viterbi_triple_t<int> vtb;

	std::vector<int> ob_set(np);
	for (int i = 0; i < np; ++i)
	{
		ob_set[i] = i;
	}
	std::vector<int> &ob_sequence = ob_set;
	std::vector<double> ob_evolve(np);
	ob_evolve[0] = 1.0; //euqal prob for the init evolving
	std::copy(pdist_vec_.begin(), pdist_vec_.end(), ++ob_evolve.begin());

	hmm.ViterbiEvolve_Modified(vtb, ob_sequence, ob_evolve, sigma_c_);

	if (vtb.vpath.empty())
	{
		return false;
	}
	deform_handle_map_.clear();
	OpenMesh::Vec3d aa;

	for (int i = 0; i < vtb.vpath.size(); i++)
	{
		COpenMeshT::VertexHandle vd = std::next(roi_vertices_normal_map_.begin(), vtb.vpath[i])->first;
		//Here we filtered an one to one map
		if (deform_handle_map_.find(vd) != deform_handle_map_.end())
		{
			if ((deform_handle_map_[vd] - mesh.point(vd)).sqrnorm() >
				(p_curveobj_[i] - mesh.point(vd)).sqrnorm())
				deform_handle_map_[vd] = p_curveobj_[i];
		}
		else
		{
			deform_handle_map_[vd] = p_curveobj_[i];
		}
	}
	return true;
}