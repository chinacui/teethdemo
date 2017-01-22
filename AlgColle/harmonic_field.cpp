#include"harmonic_field.h"
#include <igl/cotmatrix.h>
#include"../DataColle/cgal_igl_converter.h"
#include <igl/gaussian_curvature.h>
#include <igl/massmatrix.h>
#include <igl/invert_diag.h>
#include <Eigen/Sparse>
#include<Eigen/SparseLU>
#include<igl/writeDMAT.h>
#include"numerical_base_alg.h"
#include <igl/grad.h>
#include"morphlogic_operation.h"
#include"geo_base_alg.h"
#include"geo_alg.h"
bool CHarmonicFieldSeg::IsConcave(COpenMeshT &mesh, COpenMeshT::VertexHandle vh)
{
	double theta = 0.01;
	OpenMesh::Vec3d v = mesh.point(vh);
	OpenMesh::Vec3d vnorm = mesh.normal(vh);
	//std::cerr << vnorm << std::endl;
	for (auto vviter = mesh.vv_begin(vh); vviter != mesh.vv_end(vh); vviter++)
	{
		OpenMesh::Vec3d vv = mesh.point(vviter);
		OpenMesh::Vec3d dif_norm=(v - vv).normalized();
		OpenMesh::Vec3d vvnorm = mesh.normal(vviter);
		OpenMesh::Vec3d norm_dif = vvnorm - vnorm;
		if (OpenMesh::dot(dif_norm,norm_dif) > theta)
			return true;
	}
	return false;
}
void CHarmonicFieldSeg::RefineSegTwoTooth(COpenMeshT &mesh, std::vector<COpenMeshT::VertexHandle>&vhs0, std::vector<COpenMeshT::VertexHandle>&vhs1)
{


	std::vector<bool>teeth_mark0(mesh.n_vertices(), false);
	std::vector<bool>teeth_mark1(mesh.n_vertices(), false);
	for (int i = 0; i < vhs0.size(); i++)
	{
		teeth_mark0[vhs0[i].idx()] = true;
	}
	for (int i = 0; i < vhs1.size(); i++)
	{
		teeth_mark1[vhs1[i].idx()] = true;
	}

	int ecount = 6;
	while (ecount--)
	{
		CMorphlogicOperation::Erode(mesh, teeth_mark0);
	}
	ecount = 6;
	while (ecount--)
	{
		CMorphlogicOperation::Erode(mesh, teeth_mark1);
	}
	std::vector<COpenMeshT::VertexHandle>tmp_orig_teeth_vhs0(0);
	std::vector<COpenMeshT::VertexHandle>tmp_orig_teeth_vhs1(0);
	for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); viter++)
	{
		if (teeth_mark0[viter->idx()] == true)
		{
			tmp_orig_teeth_vhs0.push_back(viter);
		}
		if (teeth_mark1[viter->idx()] == true)
		{
			tmp_orig_teeth_vhs1.push_back(viter);
		}
		mesh.set_color(viter, OpenMesh::Vec3d(1, 0, 0));
	}
	SegTwoTooth(mesh, tmp_orig_teeth_vhs0, tmp_orig_teeth_vhs1,vhs0,vhs1);
	/*for (int i = 0; i < vhs0.size(); i++)
	{
		mesh.set_color(vhs0[i], OpenMesh::Vec3d(0,1, 0));
	}
	for (int i = 0; i < vhs1.size(); i++)
	{
		mesh.set_color(vhs1[i], OpenMesh::Vec3d(0, 0, 1));
	}*/

}
void CHarmonicFieldSeg::RefineTeethGingivaSeg(COpenMeshT &mesh, std::vector<int>&teeth_mark)
{
	std::map<int, std::vector<COpenMeshT::VertexHandle>>teeth_vhs_map;
	teeth_vhs_map.clear();
	for(auto viter=mesh.vertices_begin();viter!=mesh.vertices_end();viter++)
	{
		int vid = viter->idx();
		int tid = teeth_mark[vid];
		if (tid == -1)
			continue;
		if (teeth_vhs_map.find(tid) == teeth_vhs_map.end())
		{
			teeth_vhs_map[tid] = std::vector<COpenMeshT::VertexHandle>();
		}
		teeth_vhs_map[tid].push_back(viter);
	}
	int minmum_teeth_num = 10;
	std::vector<int>tobe_del;
	for (auto iter = teeth_vhs_map.begin(); iter != teeth_vhs_map.end(); iter++)
	{
		if (iter->second.size() < minmum_teeth_num)
			tobe_del.push_back(iter->first);
	}
	for (int i = 0; i < tobe_del.size(); i++)
	{
		teeth_vhs_map.erase(tobe_del[i]);
	}
	teeth_mark.resize(mesh.n_vertices(), -1);
	auto iter_pre = teeth_vhs_map.begin();
	auto iter = iter_pre;
	iter++;
	for (; iter != teeth_vhs_map.end()&&iter_pre!=teeth_vhs_map.end(); iter++)
	{
		RefineSegTwoTooth(mesh, teeth_vhs_map[iter_pre->first], teeth_vhs_map[iter->first]);

		for (int i = 0; i < teeth_vhs_map[iter_pre->first].size(); i++)
		{
			mesh.set_color(teeth_vhs_map[iter_pre->first][i], OpenMesh::Vec3d(0, 1, 0));
		}
		for (int i = 0; i < teeth_vhs_map[iter->first].size(); i++)
		{
			mesh.set_color(teeth_vhs_map[iter->first][i], OpenMesh::Vec3d(0, 0, 1));
		}
		iter++;
		iter_pre = iter;
	}
	if (iter == teeth_vhs_map.end()&&iter_pre!=teeth_vhs_map.end()&&iter_pre!=teeth_vhs_map.begin())
	{
		RefineSegTwoTooth(mesh, teeth_vhs_map[iter_pre->first], teeth_vhs_map[teeth_vhs_map.begin()->first]);
	
	}
	for (int i = 0; i < teeth_mark.size(); i++)
	{
		teeth_mark[i] = -1;
	}
	for (auto iter = teeth_vhs_map.begin(); iter != teeth_vhs_map.end(); iter++)
	{
		std::vector<COpenMeshT::VertexHandle>&vhs = iter->second;
		for (int i = 0; i < vhs.size(); i++)
		{
			teeth_mark[vhs[i].idx()] = iter->first;
		}
	}
}
void CHarmonicFieldSeg::RefineTeethGingivaSeg(COpenMeshT&mesh, std::vector<COpenMeshT::VertexHandle>&orig_teeth_vhs, std::vector<COpenMeshT::VertexHandle>&res_teeth)
{
	std::vector<COpenMeshT::VertexHandle>tmp_orig_teeth_vhs(0); 
	std::vector<COpenMeshT::VertexHandle>tmp_orig_gingiva_vhs(0);

	std::vector<bool>teeth_mark(mesh.n_vertices(), false);
	std::vector<bool>gingiva_mark(mesh.n_vertices(), true);
	for (int i = 0; i < orig_teeth_vhs.size(); i++)
	{
		teeth_mark[orig_teeth_vhs[i].idx()] = true;
		gingiva_mark[orig_teeth_vhs[i].idx()] = false;
	}

	int ecount = 6;
	while (ecount--)
	{
		CMorphlogicOperation::Erode(mesh, teeth_mark);
	}
	ecount = 4;
	while (ecount--)
	{
		CMorphlogicOperation::Erode(mesh, gingiva_mark);
	}
	

	for (int i = 0; i < orig_teeth_vhs.size(); i++)
	{
		if (teeth_mark[orig_teeth_vhs[i].idx()])
		{
			tmp_orig_teeth_vhs.push_back(orig_teeth_vhs[i]);
		}
	}
	for(auto viter=mesh.vertices_begin();viter!=mesh.vertices_end();viter++)
	{
		if (gingiva_mark[viter->idx()] == true)
		{
			tmp_orig_gingiva_vhs.push_back(viter);
		}
	}

	std::vector<COpenMeshT::VertexHandle>res_gingiva_vhs;
	SegTwoSet(mesh, tmp_orig_teeth_vhs, tmp_orig_gingiva_vhs, res_teeth, res_gingiva_vhs);
	
}
void CHarmonicFieldSeg::GetConcavityAwareLaplacianMatrix(COpenMeshT &mesh, std::vector<Eigen::Triplet<double>>&triplets)
{
	triplets.clear();
	Eigen::MatrixXd vertexs;
	Eigen::MatrixXi faces;
	
	CConverter::ConvertFromOpenMeshToIGL(mesh, vertexs, faces);
	Eigen::VectorXd gaussian_curv;
	igl::gaussian_curvature(vertexs, faces, gaussian_curv);
//	igl::writeDMAT("gaussian.dmat", gaussian_curv);
	Eigen::SparseMatrix<double> M, Minv;
	igl::massmatrix(vertexs, faces, igl::MASSMATRIX_TYPE_DEFAULT, M);
	igl::invert_diag(M, Minv);
	gaussian_curv = (Minv*gaussian_curv).eval();

	double gama = 0.0001;
	double beta = 0.01;

	std::vector<bool>is_concave(mesh.n_vertices(), false);
	for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); viter++)
	{
		is_concave[viter->idx()] = IsConcave(mesh, viter);
	}

	for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); viter++)
	{
		int vid = viter->idx();
		OpenMesh::Vec3d pv = mesh.point(viter);
		for (auto hiter = mesh.voh_begin(viter); hiter != mesh.voh_end(viter); hiter++)
		{
			auto vviter=mesh.to_vertex_handle(hiter);
			OpenMesh::Vec3d pvv = mesh.point(vviter);
			int vvid = vviter.idx();
			double elen = (pvv - pv).length();
			double w = elen / (std::abs(gaussian_curv[vid] + gaussian_curv[vvid]) + gama);

			if (is_concave[vid] || is_concave[vvid])
			{
				w *= beta;
				//mesh.data(*hiter).SetUV(OpenMesh::Vec2f(0.8, 0.8));
			}
			else
			{
				//w *= 10;
			}
		//	else
				//mesh.data(*hiter).SetUV(OpenMesh::Vec2f(0.1, 0.1));
			triplets.push_back(Eigen::Triplet<double>(vid, vvid, -w));
		//	std::cerr << w << std::endl;
			mesh.data(*hiter).SetUV(OpenMesh::Vec2f(w,w));
		}
		/*for (auto vviter = mesh.vv_begin(viter); vviter != mesh.vv_end(viter); vviter++)
		{
			OpenMesh::Vec3d pvv = mesh.point(vviter);
			int vvid = vviter->idx();
			double elen = (pvv - pv).length();
			double w=elen/(std::abs(gaussian_curv[vid]+gaussian_curv[vvid])+gama);
			
			if (is_concave[vid] || is_concave[vvid])
			{
				w *= beta;
			}
			triplets.push_back(Eigen::Triplet<double>(vid,vvid,-w));
		}*/
	}
	std::vector<double>diag(mesh.n_vertices(), 0);
	for (int i = 0; i < triplets.size(); i++)
	{
		diag[triplets[i].row()] -= triplets[i].value();
	}
	for (int i = 0; i < diag.size(); i++)
	{
		triplets.push_back(Eigen::Triplet<double>(i, i, diag[i]));
	}
}
void CHarmonicFieldSeg::SegOneTeeth(COpenMeshT&mesh, std::vector<COpenMeshT::VertexHandle>&vhs, std::vector<COpenMeshT::VertexHandle>&res_teeth)
{
	res_teeth.clear();
	std::vector<std::pair<COpenMeshT::VertexHandle, double>>cons;
	for (int i = 0; i < vhs.size(); i++)
	{
		for (auto vviter = mesh.vv_begin(vhs[i]); vviter != mesh.vv_end(vhs[i]); vviter++)
			cons.push_back(std::make_pair(vviter, 1));
	}
	std::set<COpenMeshT::VertexHandle>vhset;
	for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); viter++)
	{
		if (mesh.is_boundary(viter))
		{
			cons.push_back(std::make_pair(viter, 0));
			for (auto vviter = mesh.vv_begin(viter); vviter != mesh.vv_end(viter); vviter++)
			{
				if (mesh.is_boundary(vviter)==false&&vhset.find(vviter) == vhset.end())
				{
					vhset.insert(vviter);
					cons.push_back(std::make_pair(vviter,0));
				}
			}
		}
	}
	Eigen::VectorXd harmonic_field;
	ComputeConcavityAwareHarmonicField(mesh, cons, harmonic_field);
	CNumericalBaseAlg::NormalizeScalarField(harmonic_field);
	std::vector<std::pair<double, double>>bins;
	std::vector<std::vector<int>>bin_vids;
	double bin_size = 0.0006;
	CNumericalBaseAlg::ComputeHistgram(harmonic_field, bin_size, bins, bin_vids);
	Eigen::MatrixXd vertexs;
	Eigen::MatrixXi faces;

	CConverter::ConvertFromOpenMeshToIGL(mesh, vertexs, faces);
	Eigen::SparseMatrix<double> G;
	igl::grad(vertexs, faces, G);
	// Compute gradient of U
	Eigen::MatrixXd fgrad_harmonic = Eigen::Map<const Eigen::MatrixXd>((G*harmonic_field).eval().data(), faces.rows(), 3);
	Eigen::VectorXd fgrad_harmonic_mag = fgrad_harmonic.rowwise().norm();

	std::vector<double>ave_bin_grad(bin_vids.size(), 0);
	std::vector<int>bin_count(bin_vids.size(), 0);
	for (auto fiter = mesh.faces_begin(); fiter != mesh.faces_end(); fiter++)
	{
		double fharmonic_min = std::numeric_limits<double>::max();
		double fharmonic_max = std::numeric_limits<double>::min();
		int fid = fiter->idx();
		for (auto fviter = mesh.fv_begin(fiter); fviter != mesh.fv_end(fiter); fviter++)
		{
			double value = harmonic_field(fviter->idx());
			if (value > fharmonic_max)
				fharmonic_max = value;
			if (value < fharmonic_min)
				fharmonic_min = value;
		}
		int sbin_id = fharmonic_min / bin_size;
		int ebin_id = fharmonic_max / bin_size;
		for (int i = sbin_id; i <= ebin_id; i++)
		{
			ave_bin_grad[i] += fgrad_harmonic_mag(fid);
			if (mesh.is_boundary(fiter))
			{
				ave_bin_grad[i] = 0;
			}
			bin_count[i]++;
		}

	}
	for (int i = 0; i < bin_count.size(); i++)
	{
		ave_bin_grad[i] /= bin_count[i];
	}
	int max_bin_id = 0;
	for (int i = 0; i < ave_bin_grad.size(); i++)
	{
		if (i <= max_bin_id + 2)
		{
			if (ave_bin_grad[i] > ave_bin_grad[max_bin_id])
			{
				max_bin_id = i;
			}
		}
	}
	double threshold = bin_size*(max_bin_id + 1);
	for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); viter++)
	{
		int vid = viter->idx();
		if (harmonic_field(vid) >= threshold)
		{
			res_teeth.push_back(viter);
		}
	}

}
void CHarmonicFieldSeg::RefineToothGingivaSeg(COpenMeshT&mesh, std::vector<COpenMeshT::VertexHandle>&tooth_vhs, std::vector<COpenMeshT::VertexHandle>&res_tooth)
{
	std::vector<bool>tags_teeth(mesh.n_vertices(),false), tags_gingiva(mesh.n_vertices(),false);
	for (int i = 0; i < tooth_vhs.size(); i++)
	{
		tags_teeth[tooth_vhs[i].idx()] = true;
	}
	for (int i = 0; i < tags_teeth.size(); i++)
	{
		if (tags_teeth[i] == false)
			tags_gingiva[i] = true;
	}
	
	int mcount = 3;
	while (mcount--)
	{
		CMorphlogicOperation::Erode(mesh, tags_teeth);

	}
	mcount = 8;
	while (mcount--)
	{
		CMorphlogicOperation::Erode(mesh, tags_gingiva);
	}

	std::vector<OpenMesh::VertexHandle>teeth_edges, gingiva_edges;
	CGeoBaseAlg::GetEdgeVertexs(mesh, tags_teeth, teeth_edges);
	CGeoBaseAlg::GetEdgeVertexs(mesh, tags_gingiva, gingiva_edges);
	std::vector<std::pair<COpenMeshT::VertexHandle, double>>cons;
	for (int i = 0; i < teeth_edges.size(); i++)
	{
		cons.push_back(std::make_pair(teeth_edges[i], 0.01));
	}
	for (int i = 0; i < gingiva_edges.size(); i++)
	{
		cons.push_back(std::make_pair(gingiva_edges[i], 0.99));
	}
	Eigen::VectorXd harmonic_field;
	CHarmonicFieldSeg harmonicSeg;
	harmonicSeg.ComputeConcavityAwareHarmonicField(mesh, cons, harmonic_field);

	CNumericalBaseAlg::NormalizeScalarField(harmonic_field);

	std::vector<std::pair<double, double>>bins;
	std::vector<std::vector<int>>bin_vids;
	double bin_size = 0.0006;
	CNumericalBaseAlg::ComputeHistgram(harmonic_field, bin_size, bins, bin_vids);
	double bin_start = bins[0].first;


	double max_threshold = 0.8;

	Eigen::MatrixXd vertexs;
	Eigen::MatrixXi faces;

	CConverter::ConvertFromOpenMeshToIGL(mesh, vertexs, faces);
	Eigen::SparseMatrix<double> G;
	igl::grad(vertexs, faces, G);
	// Compute gradient of U
	Eigen::MatrixXd fgrad_harmonic = Eigen::Map<const Eigen::MatrixXd>((G*harmonic_field).eval().data(), faces.rows(), 3);
	Eigen::VectorXd fgrad_harmonic_mag = fgrad_harmonic.rowwise().norm();

	std::vector<double>ave_bin_grad(bin_vids.size(), 0);
	std::vector<int>bin_count(bin_vids.size(), 0);

	for (auto fiter = mesh.faces_begin(); fiter != mesh.faces_end(); fiter++)
	{
		double fharmonic_min = std::numeric_limits<double>::max();
		double fharmonic_max = std::numeric_limits<double>::min();
		int fid = fiter->idx();
		for (auto fviter = mesh.fv_begin(fiter); fviter != mesh.fv_end(fiter); fviter++)
		{
			double value = harmonic_field(fviter->idx());
			if (value > fharmonic_max)
				fharmonic_max = value;
			if (value < fharmonic_min)
				fharmonic_min = value;
		}
		int sbin_id = (fharmonic_min - bin_start) / bin_size;
		int ebin_id = (fharmonic_max - bin_start) / bin_size;
		for (int i = sbin_id; i <= ebin_id; i++)
		{
			ave_bin_grad[i] += fgrad_harmonic_mag(fid);
			if (mesh.is_boundary(fiter))
			{
				ave_bin_grad[i] = -99999;
			}
			bin_count[i]++;
		}

	}
	for (int i = 0; i < bin_count.size(); i++)
	{
		ave_bin_grad[i] /= bin_count[i];
	}
	int max_bin_id = 0;
	for (int i = 0; i < ave_bin_grad.size(); i++)
	{
		if ((i + 1)*bin_size + bin_start > max_threshold)
		{
			if (ave_bin_grad[i] > ave_bin_grad[max_bin_id])
			{
				max_bin_id = i;
			}
		}

	}
	double threshold = bin_size*(max_bin_id + 1) + bin_start;
	res_tooth.clear();
	for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); viter++)
	{
		int vid = viter->idx();

		if (harmonic_field(vid) < threshold)
		{
			res_tooth.push_back(viter);
		}

	}
	std::cerr << res_tooth.size() << std::endl;
}
void CHarmonicFieldSeg::SegToothGingiva(COpenMeshT&mesh, std::vector<COpenMeshT::VertexHandle>&vhs_tooth, std::vector<COpenMeshT::VertexHandle>&vhs_gingiva, std::vector<COpenMeshT::VertexHandle>&res_tooth, std::vector<COpenMeshT::VertexHandle>&res_gingiva)
{

	res_tooth.clear();
	std::vector<std::pair<COpenMeshT::VertexHandle, double>>cons;
	std::set<COpenMeshT::VertexHandle>tooth_vhs_set,gingiva_vhs_set;
	for (int i = 0; i < vhs_tooth.size(); i++)
	{
		std::vector<COpenMeshT::VertexHandle>neis;
		CGeoAlg::ExtractNRing(mesh, vhs_tooth[i], 2, neis);
		if (tooth_vhs_set.find(vhs_tooth[i]) == tooth_vhs_set.end())
		{
			tooth_vhs_set.insert(vhs_tooth[i]);
			cons.push_back(std::make_pair(vhs_tooth[i], 0.001));
		}
		for (int j = 0; j < neis.size(); j++)
		{
			if (tooth_vhs_set.find(neis[j]) == tooth_vhs_set.end())
			{
				tooth_vhs_set.insert(neis[j]);
				cons.push_back(std::make_pair(neis[j], 0.001));
			}
		}
		
		
	}
	for (int i = 0; i < vhs_gingiva.size(); i++)
	{

		if (gingiva_vhs_set.find(vhs_gingiva[i]) == gingiva_vhs_set.end())
		{
			gingiva_vhs_set.insert(vhs_gingiva[i]);
			cons.push_back(std::make_pair(vhs_gingiva[i], 0.999));
		}
	}
	Eigen::VectorXd harmonic_field;
	ComputeConcavityAwareHarmonicField(mesh, cons, harmonic_field);
	//for (auto hiter = mesh.halfedges_begin(); hiter != mesh.halfedges_end(); hiter++)
	//{
	//auto vh = mesh.to_vertex_handle(hiter);
	//mesh.data(hiter).SetUV(OpenMesh::Vec2f(harmonic_field(vh.idx()), harmonic_field(vh.idx())));
	//}

	//return;
	CNumericalBaseAlg::NormalizeScalarField(harmonic_field);

	std::vector<std::pair<double, double>>bins;
	std::vector<std::vector<int>>bin_vids;
	double bin_size = 0.0006;
	CNumericalBaseAlg::ComputeHistgram(harmonic_field, bin_size, bins, bin_vids);
	double bin_start = bins[0].first;
	

	double mid_threshold =0.3,max_threshold=0.7;

	Eigen::MatrixXd vertexs;
	Eigen::MatrixXi faces;

	CConverter::ConvertFromOpenMeshToIGL(mesh, vertexs, faces);
	Eigen::SparseMatrix<double> G;
	igl::grad(vertexs, faces, G);
	// Compute gradient of U
	Eigen::MatrixXd fgrad_harmonic = Eigen::Map<const Eigen::MatrixXd>((G*harmonic_field).eval().data(), faces.rows(), 3);
	Eigen::VectorXd fgrad_harmonic_mag = fgrad_harmonic.rowwise().norm();

	std::vector<double>ave_bin_grad(bin_vids.size(), 0);
	std::vector<int>bin_count(bin_vids.size(), 0);

	for (auto fiter = mesh.faces_begin(); fiter != mesh.faces_end(); fiter++)
	{
		double fharmonic_min = std::numeric_limits<double>::max();
		double fharmonic_max = std::numeric_limits<double>::min();
		int fid = fiter->idx();
		for (auto fviter = mesh.fv_begin(fiter); fviter != mesh.fv_end(fiter); fviter++)
		{
			double value = harmonic_field(fviter->idx());
			if (value > fharmonic_max)
				fharmonic_max = value;
			if (value < fharmonic_min)
				fharmonic_min = value;
		}
		int sbin_id = (fharmonic_min - bin_start) / bin_size;
		int ebin_id = (fharmonic_max - bin_start) / bin_size;
		for (int i = sbin_id; i <= ebin_id; i++)
		{
			ave_bin_grad[i] += fgrad_harmonic_mag(fid);
			if (mesh.is_boundary(fiter))
			{
				ave_bin_grad[i] = -99999;
			}
			bin_count[i]++;
		}

	}
	for (int i = 0; i < bin_count.size(); i++)
	{
		ave_bin_grad[i] /= bin_count[i];
	}
	int max_bin_id0=0,max_bin_id1=0;
	for (int i = 0; i < ave_bin_grad.size(); i++)
	{
		if ((i+1)*bin_size+bin_start<mid_threshold)
		{
			if (ave_bin_grad[i] > ave_bin_grad[max_bin_id0])
			{
				max_bin_id0 = i;
			}
		}
		if (i*bin_size + bin_start > max_threshold)
		{
			if (ave_bin_grad[i] > ave_bin_grad[max_bin_id1])
			{
				max_bin_id1 = i;
			}
		}
		
	}
	double threshold0 = bin_size*(max_bin_id0 + 1) + bin_start;
	double threshold1 = bin_size*(max_bin_id1 + 1) + bin_start;
	
	for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); viter++)
	{
		int vid = viter->idx();

		if (harmonic_field(vid) <= threshold0)
		{
			res_tooth.push_back(viter);
		}
		else if (harmonic_field(vid) >= threshold1)
		{
			res_gingiva.push_back(viter);
		}
		else
		{
			int c = 0;
			for (auto vviter = mesh.vv_begin(viter); vviter != mesh.vv_end(viter); vviter++)
			{
				if (harmonic_field(vviter->idx()) < threshold1)
				{
					c++;
				}
			}
			if (c <= 2)
			{
				res_gingiva.push_back(viter);
			}
		}
		
	}
	std::cerr << res_tooth.size() << std::endl;
}

void CHarmonicFieldSeg::SegTwoSet(COpenMeshT&mesh, std::vector<COpenMeshT::VertexHandle>&vhs0, std::vector<COpenMeshT::VertexHandle>&vhs1, std::vector<COpenMeshT::VertexHandle>&res_vhs0, std::vector<COpenMeshT::VertexHandle>&res_vhs1)
{
	SegTwoTooth(mesh, vhs0, vhs1, res_vhs0, res_vhs1);
}
void CHarmonicFieldSeg::SegTwoTooth(COpenMeshT&mesh, std::vector<COpenMeshT::VertexHandle>&vhs0, std::vector<COpenMeshT::VertexHandle>&vhs1, std::vector<COpenMeshT::VertexHandle>&res_teeth0, std::vector<COpenMeshT::VertexHandle>&res_teeth1)
{
	res_teeth1.clear();
	res_teeth0.clear();
	std::cerr << "seg two tooth" << std::endl;
	std::cerr << vhs0.size() << std::endl;
	std::cerr << vhs1.size() << std::endl;
	std::vector<std::pair<COpenMeshT::VertexHandle, double>>cons;
	for (int i = 0; i < vhs0.size(); i++)
	{
		//for(auto vviter=mesh.vv_begin(vhs0[i]);vviter!=mesh.vv_end(vhs0[i]);vviter++)
		cons.push_back(std::make_pair(vhs0[i], 0.1));
	}
	for (int i = 0; i < vhs1.size(); i++)
	{
		//for (auto vviter = mesh.vv_begin(vhs1[i]); vviter != mesh.vv_end(vhs1[i]); vviter++)
		cons.push_back(std::make_pair(vhs1[i], 0.9));
	}
	Eigen::VectorXd harmonic_field;
	ComputeConcavityAwareHarmonicField(mesh, cons, harmonic_field);
	for (auto hiter = mesh.halfedges_begin(); hiter != mesh.halfedges_end(); hiter++)
	{
		auto vh = mesh.to_vertex_handle(hiter);
		mesh.data(hiter).SetUV(OpenMesh::Vec2f(harmonic_field(vh.idx()), harmonic_field(vh.idx())));
	}
	
	CNumericalBaseAlg::NormalizeScalarField(harmonic_field);
	
	
	std::vector<std::pair<double, double>>bins;
	std::vector<std::vector<int>>bin_vids;
	double bin_size = 0.0006;
	CNumericalBaseAlg::ComputeHistgram(harmonic_field, bin_size, bins, bin_vids);
	double bin_start = bins[0].first;
	int max_bin = -1;
	int max_bin_id=0;

	for (int i = 0; i < bin_vids.size(); i++)
	{
	
		int bin_s = bin_vids[i].size();
		if (max_bin <bin_s)
		{
			max_bin = bin_s;
			max_bin_id = i;
			
		}
	}

	double mid_threshold = bin_size*(max_bin_id + 0.5)+bin_start;

	Eigen::MatrixXd vertexs;
	Eigen::MatrixXi faces;

	CConverter::ConvertFromOpenMeshToIGL(mesh, vertexs, faces);
	Eigen::SparseMatrix<double> G;
	igl::grad(vertexs, faces, G);
	// Compute gradient of U
	Eigen::MatrixXd fgrad_harmonic = Eigen::Map<const Eigen::MatrixXd>((G*harmonic_field).eval().data(), faces.rows(), 3);
	Eigen::VectorXd fgrad_harmonic_mag = fgrad_harmonic.rowwise().norm();

	std::vector<double>ave_bin_grad(bin_vids.size(),0);
	std::vector<int>bin_count(bin_vids.size(), 0);
	
	for (auto fiter = mesh.faces_begin(); fiter != mesh.faces_end(); fiter++)
	{
		double fharmonic_min = std::numeric_limits<double>::max();
		double fharmonic_max = std::numeric_limits<double>::min();
		int fid = fiter->idx();
		for (auto fviter = mesh.fv_begin(fiter); fviter != mesh.fv_end(fiter); fviter++)
		{
			double value=harmonic_field(fviter->idx());
			if (value > fharmonic_max)
				fharmonic_max = value;
			if (value < fharmonic_min)
				fharmonic_min = value;
		}
		int sbin_id = (fharmonic_min-bin_start) / bin_size;
		int ebin_id = (fharmonic_max-bin_start) / bin_size;
		for (int i = sbin_id; i <=ebin_id; i++)
		{
			ave_bin_grad[i] += fgrad_harmonic_mag(fid);
			if (mesh.is_boundary(fiter))
			{
				ave_bin_grad[i] = -99999;
			}
			bin_count[i]++;
		}
		
	}
	for (int i = 0; i < bin_count.size(); i++)
	{
		ave_bin_grad[i] /= bin_count[i];
	}
	int max_bin_id0=0, max_bin_id1 = 0;
	for (int i = 0; i < ave_bin_grad.size(); i++)
	{
		if (i <=max_bin_id+2)
		{
			if (ave_bin_grad[i] > ave_bin_grad[max_bin_id0])
			{
				max_bin_id0 = i;
			}
		}
		if(i>=max_bin_id-1)
		{
			if (ave_bin_grad[i] > ave_bin_grad[max_bin_id1])
			{
				max_bin_id1 = i;
			}
		}
	}
	double threshold0 = bin_size*(max_bin_id0+1)+bin_start;
	double threshold1 = bin_size*max_bin_id1+bin_start;
	for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); viter++)
	{
		int vid = viter->idx();
		//mesh.set_color(viter, OpenMesh::Vec3d(1, 1, 1));
		if (harmonic_field(vid) <= threshold0)
		{
			res_teeth0.push_back(viter);
			//mesh.set_color(viter, OpenMesh::Vec3d(1, 0, 0));
		}
		if (harmonic_field(vid) >= threshold1)
		{
			res_teeth1.push_back(viter);
			//mesh.set_color(viter, OpenMesh::Vec3d(0, 0, 1));
		}

	}
	/*std::cerr << "threshold " << threshold0 << " " << threshold1 << std::endl;
	std::cerr << "res_t" << res_teeth0.size() << " " << res_teeth1.size() << std::endl;*/
	/*for (auto hiter = mesh.halfedges_begin(); hiter != mesh.halfedges_end(); hiter++)
	{
		int vid = mesh.to_vertex_handle(hiter).idx();
		if (harmonic_field(vid) > mid_threshold)
		{
			mesh.data(hiter).SetUV(OpenMesh::Vec2f(0.1, 0.1));
		}
		else
		{
			mesh.data(hiter).SetUV(OpenMesh::Vec2f(0.8, 0.8));
		}
	}*/
	/*std::vector<std::vector<COpenMeshT::VertexHandle>>bin_vhs;
	bin_vhs.resize(bin_vids.size());
	for (int i = 0; i < bin_vids.size(); i++)
	{
		bin_vids[i].clear();
		for (int j = 0; j < bin_vids[i].size(); j++)
		{
			bin_vhs[i].push_back(mesh.vertex_handle(bin_vids[i][j]));
		}
	}*/



	
	
}
void CHarmonicFieldSeg::ComputeConcavityAwareHarmonicField(COpenMeshT&mesh, std::vector<std::pair<COpenMeshT::VertexHandle, double>>&cons, Eigen::VectorXd &res_u)
{
	std::vector<Eigen::Triplet<double>>ca_lp_trip;
	GetConcavityAwareLaplacianMatrix(mesh, ca_lp_trip);
	Eigen::VectorXd b;
	b.resize(mesh.n_vertices());

	b.setZero();
	double alpha = 1e8;
	for (int i = 0; i < cons.size(); i++)
	{
	
		int vid = cons[i].first.idx();

		ca_lp_trip.push_back(Eigen::Triplet<double>(vid, vid, alpha));
	
	
		b(vid) = cons[i].second*alpha;
	}
	Eigen::SparseMatrix<double>L(mesh.n_vertices(),mesh.n_vertices());
	L.setFromTriplets(ca_lp_trip.begin(),ca_lp_trip.end());
	Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>>ldlt;
	//auto Lt = L*L.transpose();
	ldlt.compute(L);
	res_u=ldlt.solve(b).eval();
	igl::writeDMAT("res_u.dmat", res_u);

}