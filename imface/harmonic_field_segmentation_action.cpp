#include"harmonic_field_segmentation_action.h"
#include"../AlgColle/geo_alg.h"
#include"ui_context.h"
#include"../DataColle/data_pool.h"
#include"cmodelviewer.h"
#include"camera.h"
#include"../AlgColle/geo_base_alg.h"
#include<qimage.h>
#include"../AlgColle/harmonic_field.h"
#include <OpenMesh/Core/IO/MeshIO.hh>
#include<igl/writeDMAT.h>
#include"../DataColle/cgal_igl_converter.h"
#include<igl/readDMAT.h>
#include<igl/writeSTL.h>
#include<igl/writeOBJ.h>
#include<set>
#include<time.h>
#include"../AlgColle/dental_base_alg.h"
#include"../DataColle/aux_geo_utils.h"
#include"../AlgColle/curve_base_alg.h"
#include"../AlgColle/image_base_alg.h"
#include<opencv2/opencv.hpp>
#include"../AlgColle/numerical_base_alg.h"
#include"../AlgColle/morphlogic_operation.h"
#include "qfiledialog.h"
#include"../DataColle/data_io.h"
#include<sstream>
void CHarmonicFieldSegmentation::MousePressEvent(QMouseEvent *e)
{
	if (is_drawing_)
	{
		
		CMeshObject *meshobj = DataPool::GetMeshObject(dental_mesh_id_);
		if (meshobj != NULL)
		{
			OpenMesh::Vec3d  curve_color;
		
			if (is_pick_fore_)
			{
				curve_color = OpenMesh::Vec3d(0, 0, 1);
				if (!is_seg_teeth_from_gingiva_)
				{
					picked_vhs_fore_.clear();
					picked_vhs_fore_mark_.clear();
				}
				
			}
			else
			{
				curve_color = OpenMesh::Vec3d(1, 0, 0);
				if (!is_seg_teeth_from_gingiva_)
				{
					picked_vhs_back_.clear();
					picked_vhs_back_mark_.clear();
				}
			}


			auto camera = viewer_->GetCamera();
			OpenMesh::Vec3d orig, dir;
			camera.ConvertClickToLine(e->pos(), orig, dir);

			COpenMeshT::VertexHandle vh;
			if (CGeoAlg::RayMeshIntersection(orig, dir, *meshobj, vh))
			{
				COpenMeshT &mesh = meshobj->GetMesh();
				if (is_pick_fore_)
				{
					if (picked_vhs_fore_mark_.find(vh) == picked_vhs_fore_mark_.end())
					{
						picked_vhs_fore_.push_back(vh);
						picked_vhs_fore_mark_.insert(vh);
					}
						
				}	
				else
				{
					if (picked_vhs_back_mark_.find(vh) == picked_vhs_back_mark_.end())
					{
						picked_vhs_back_.push_back(vh);
						picked_vhs_back_mark_.insert(vh);
					}
				
				}
				//RenderFeature();
			}

			orig = orig + dir*0.01;
			if (p_curve_obj_ != NULL)
				DataPool::DeleteCurveObject(p_curve_obj_->GetId());
			p_curve_obj_ = new CCurveObject();
		
			p_curve_obj_->GetCurve().push_back(orig);
			p_curve_obj_->SetChanged();
			p_curve_obj_->SetColor(curve_color);
			DataPool::AddCurveObject(p_curve_obj_);
			
		}
	}
	else if (is_eliminating_feature_)
	{
		if (e->button() == Qt::LeftButton)
		{
		
			CMeshObject *meshobj = DataPool::GetMeshObject(dental_mesh_id_);
			if (meshobj != NULL)
			{
				auto camera = viewer_->GetCamera();
				OpenMesh::Vec3d orig, dir;
				camera.ConvertClickToLine(e->pos(), orig, dir);
				COpenMeshT::VertexHandle vh;
				std::cerr << "eliminating_feature0" << std::endl;
				if (CGeoAlg::RayMeshIntersection(orig, dir, *meshobj, vh))
				{
					std::cerr << "eliminating_feature" << std::endl;
					COpenMeshT&mesh = meshobj->GetMesh();
					OpenMesh::Vec3d pickedp = mesh.point(vh);
					double min_dis = std::numeric_limits<double>::max();
					int mi;
					for (int i = 0; i < picked_vhs_fore_.size(); i++)
					{
						OpenMesh::Vec3d p = mesh.point(picked_vhs_fore_[i]);
						double dis = (p - pickedp).length();
						if (dis < min_dis)
						{
							min_dis = dis;
							mi = i;
						}
					}

					picked_vhs_fore_mark_.erase(picked_vhs_fore_[mi]);
					picked_vhs_fore_.clear();
					for (auto iter = picked_vhs_fore_mark_.begin(); iter != picked_vhs_fore_mark_.end(); iter++)
					{
						picked_vhs_fore_.push_back(*iter);
					}
					RenderFeature();
				}
			}
		}
	}
	else
	{
		if (e->button() == Qt::LeftButton)
		{
		
			CMeshObject *meshobj = DataPool::GetMeshObject(dental_mesh_id_);
			if (meshobj != NULL)
			{
				auto camera = viewer_->GetCamera();
				OpenMesh::Vec3d orig, dir;
				camera.ConvertClickToLine(e->pos(), orig, dir);
				COpenMeshT::VertexHandle vh;
				if (CGeoAlg::RayMeshIntersection(orig, dir, *meshobj, vh))
				{
					COpenMeshT &mesh = meshobj->GetMesh();
					if (is_pick_fore_)
					{
						if (picked_vhs_fore_mark_.find(vh) == picked_vhs_fore_mark_.end())
						{
							picked_vhs_fore_.push_back(vh);
							picked_vhs_fore_mark_.insert(vh);
						}
						
					}
					else
					{
						if (picked_vhs_back_mark_.find(vh) == picked_vhs_back_mark_.end())
						{
							picked_vhs_back_.push_back(vh);
							picked_vhs_back_mark_.insert(vh);
						}
							
					}
						

					std::cerr << "pick " << vh.idx() << std::endl;

				}

			}
		}
		else
		{
		
			CMeshObject *meshobj = DataPool::GetMeshObject(dental_mesh_id_);
			if (meshobj != NULL)
			{
				auto camera = viewer_->GetCamera();
				OpenMesh::Vec3d orig, dir;
				camera.ConvertClickToLine(e->pos(), orig, dir);
				COpenMeshT::VertexHandle vh;
				if (CGeoAlg::RayMeshIntersection(orig, dir, *meshobj, vh))
				{
					//CUIContext::msdm_seg_->SwitchGingivaAndTeeth(vh);
					//CUIContext::msdm_seg_->RemoveSmallIsolateTeethRegion();
				}
			}
		}
	}


}
void CHarmonicFieldSegmentation::RenderFeature()
{

	CMeshObject *meshobj = DataPool::GetMeshObject(dental_mesh_id_);
	if (meshobj != NULL)
	{
		COpenMeshT &mesh = meshobj->GetMesh();
		for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); viter++)
		{
			mesh.set_color(viter, OpenMesh::Vec3d(0.8, 0.8, 0.8));
		}
		for (int i = 0; i < picked_vhs_fore_.size(); i++)
		{
			mesh.set_color(picked_vhs_fore_[i], OpenMesh::Vec3d(1, 0, 0));
			for (auto vviter = mesh.vv_begin(picked_vhs_fore_[i]); vviter != mesh.vv_end(picked_vhs_fore_[i]); vviter++)
			{
				mesh.set_color(vviter, OpenMesh::Vec3d(1, 0, 0));
			}
	
		}

		for (int i = 0; i < picked_vhs_back_.size(); i++)
		{
			mesh.set_color(picked_vhs_back_[i], OpenMesh::Vec3d(1, 1, 0));
		}
		meshobj->SetAttrChanged();
		meshobj->UseTexture() = false;
	}
}
void CHarmonicFieldSegmentation::MouseMoveEvent(QMouseEvent *e)
{

	if (is_drawing_)
	{
		
		CMeshObject *meshobj = DataPool::GetMeshObject(dental_mesh_id_);
		if (meshobj != NULL)
		{
			auto camera = viewer_->GetCamera();
			OpenMesh::Vec3d orig, dir;
			camera.ConvertClickToLine(e->pos(), orig, dir);

			COpenMeshT::VertexHandle vh;
			if (CGeoAlg::RayMeshIntersection(orig, dir, *meshobj, vh))
			{
				COpenMeshT &mesh = meshobj->GetMesh();
				if (is_pick_fore_)
				{
					if (picked_vhs_fore_mark_.find(vh) == picked_vhs_fore_mark_.end())
					{
						picked_vhs_fore_.push_back(vh);
						picked_vhs_fore_mark_.insert(vh);
					}
				
					
					//std::cerr << "picked fore " << vh.idx()<< std::endl;
				}	
				else
				{
					if (picked_vhs_back_mark_.find(vh) == picked_vhs_back_mark_.end())
					{
						picked_vhs_back_.push_back(vh);
						picked_vhs_back_mark_.insert(vh);
					}
					
				//	std::cerr << "picked back "<<vh.idx() << std::endl;
				}
				
				
			}

			orig = orig + dir*0.01;
		
			p_curve_obj_->GetCurve().push_back(orig);
			p_curve_obj_->SetChanged();
			//std::cerr << porig << std::endl;
		}
	}
}
void CHarmonicFieldSegmentation::MouseReleaseEvent(QMouseEvent *e)
{
	if (is_drawing_&&(!is_seg_teeth_from_gingiva_))
	{
		is_pick_fore_ = !is_pick_fore_;
	}
}
CHarmonicFieldSegmentation::CHarmonicFieldSegmentation()
{
	type_ = HarmonicFieldSegmentation;

}
OpenMesh::Vec3d CHarmonicFieldSegmentation::GetRandColor()
{

	double r = (rand() % 10000)/11000.0;
	double g = (rand() % 10000) / 11000.0;
	double b = (rand() % 10000) / 11000.0;
	OpenMesh::Vec3d color(r, g, b);
	return color;
}
void CHarmonicFieldSegmentation::Init()
{
	srand((unsigned)time(NULL));
	picked_vhs_fore_.clear();
	picked_vhs_fore_mark_.clear();
	picked_vhs_back_.clear();
	picked_vhs_back_mark_.clear();

	teeth_seg_mark_.clear();
	teeth_seg_color_.clear();
	rand_color_set_.clear();
	teeth_seg_count_.clear();
	teeth_seg_mark_id_ = 0;

	dental_mesh_id_ = CUIContext::GetSelectedMeshObjectId();
	CMeshObject *meshobj = DataPool::GetMeshObject(dental_mesh_id_);
	if (meshobj != NULL)
	{
		meshobj->SetMeshColor(OpenMesh::Vec3d(0.8, 0.8, 0.8));
		meshobj->SetAttrChanged();
	}
	seg_tooth_mesh_id_ = -1;
}
void CHarmonicFieldSegmentation::KeyPressEvent(QKeyEvent *e)
{
	switch (e->key())
	{

	case Qt::Key_W:
	{
		is_drawing_ = true;

		
		break;
	}
	case Qt::Key_Q:
	{
		std::cerr << "switch to common action" << std::endl;
		manager_->SetCurrentActionType(CActionType::Common);
		break;
	}
	case Qt::Key_O:
	{
		QString path = QFileDialog::getSaveFileName((QWidget*)CUIContext::GetMainWindow(), "save seg result");

		if (path.length() == 0)
		{
			std::cerr << "unable to open path\n" << std::endl;
			return;
		}

		p_mesh_obj_ = DataPool::GetMeshObject(dental_mesh_id_);
		p_mesh_obj_->RecoverCurrentVPos();
		if (p_mesh_obj_ != NULL)
		{

			
			std::string fname,ftagname;
			fname = path.toStdString();
			if (fname.length() - 4>=0&&fname[fname.length() - 4] != '.')
			{
				ftagname = fname + ".dmat";
				fname += ".obj";
			}
			else
			{
				ftagname = fname.substr(0,fname.length()-4) + ".dmat";
			}
			CDataIO::WriteMesh(fname, *p_mesh_obj_);
			COpenMeshT &mesh = p_mesh_obj_->GetMesh();
			Eigen::VectorXi tags(mesh.n_vertices());
			std::vector<int>vtags_vec;
			std::map<int, int>&vtags = p_mesh_obj_->GetVertexTags();
			vtags_vec.resize(vtags.size());
			for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); viter++)
			{

				tags(viter->idx()) =vtags[viter->idx()];
				vtags_vec[viter->idx()] = vtags[viter->idx()];
			//	std::cerr << vtags[viter] << std::endl;
			}

			//for (int i = 0; i < tags.size(); i++)
			//{
			//	std::cerr << tags(i) << std::endl;
			//}
			igl::writeDMAT(ftagname, tags);
			std::cerr << "save " << fname << "successfully" << std::endl;

			bool write_sep_crowns = false;

			if (write_sep_crowns)
			{
				std::map<int, COpenMeshT*>crowns;
				CGeoAlg::SeparateMeshByVertexTag(p_mesh_obj_->GetMesh(), vtags_vec, crowns, std::map<int, std::map<OpenMesh::VertexHandle, OpenMesh::VertexHandle>>());
				int cid = 0;
				for (auto iter = crowns.begin(); iter != crowns.end(); iter++)
				{
					std::stringstream sstr;
					if (fname.length() - 4 >= 0 && fname[fname.length() - 4] != '.')
					{
						sstr << fname;
					}
					else
					{
						sstr << fname.substr(0, fname.length() - 4);
					}
					sstr << "_" << cid++ << ".stl";
					Eigen::MatrixXd v;
					Eigen::MatrixXi f;
					CConverter::ConvertFromOpenMeshToIGL(*(iter->second), v, f);
					//igl::writeSTL(sstr.str(), v, f);
					//igl::writeOBJ(fname, v, f);
					if (!OpenMesh::IO::write_mesh(*(iter->second), sstr.str()))
					{
						std::cerr << "write error\n";

					}
				}
			}
			
			
		}
	
		break;
	}
	case Qt::Key_Z:
	{
		is_seg_teeth_from_gingiva_ = true;
		is_drawing_ = true;
		is_pick_fore_ = false;
		break;
	}
	case Qt::Key_A:
	{
		is_seg_teeth_from_gingiva_ = true;
		is_drawing_ = true;
		is_pick_fore_ = true;
		break;
	}
	case  Qt::Key_C:
	{
		Init();
		DataPool::DeleteAllCurveObjects();
		std::cerr << "clear" << std::endl;
		break;
	}
	case Qt::Key_T:
	{
		is_pick_fore_ = !is_pick_fore_;
		std::cerr << "is pick fore " << is_pick_fore_ << std::endl;
		break;
	}
	case Qt::Key_Space:
	{
		e->setAccepted(false);
		break;
	}

	case Qt::Key_S:
	{
	
		CMeshObject *meshobj = DataPool::GetMeshObject(dental_mesh_id_);
		if (meshobj != NULL)
		{
			
			COpenMeshT&mesh = meshobj->GetMesh();
			std::cerr << "seg gingiva from teeth" << std::endl;
			

			std::map<int, COpenMeshT*>sep_meshes;
			std::map<int, std::map<COpenMeshT::VertexHandle, COpenMeshT::VertexHandle>> vid_org_map_;
			//CGeoAlg::SeparateMeshByVertexTag(mesh, teeth_seg_mark_, sep_meshes, vid_org_map_);
			CGeoAlg::SeparateMeshByVertexTag(mesh, gingiva_seg_mark_, sep_meshes, vid_org_map_);
			all_seg_tooth_orig_vhs_map_.clear();
	
			//all_seg_tooth_orig_vhs_map_ = vid_org_map_[-1];
			for(auto iter=sep_meshes.begin();iter!=sep_meshes.end();iter++)
			{
				int tid = iter->first;
				if (tid == -1)
				{
					CMeshObject* p_seg_tooth_mesh = new CMeshObject();
				
				//	CGeoBaseAlg::RemoveNonManifold(*iter->second);
					p_seg_tooth_mesh->GetMesh() = *iter->second;
				
					for (auto viter = vid_org_map_[-1].begin(); viter != vid_org_map_[-1].end(); viter++)
					{
						all_seg_tooth_orig_vhs_map_[p_seg_tooth_mesh->GetMesh().vertex_handle(viter->first.idx())] = viter->second;
					}
					//CGeoAlg::RefineMeshBoundary(p_seg_tooth_mesh->GetMesh());
				//	CGeoAlg::SimplifyMesh(p_mesh_obj->GetMesh(), 10000);
					//CGeoAlg::FillHoles(p_mesh_obj->GetMesh());
					delete sep_meshes[tid];
				
					p_seg_tooth_mesh->SetMeshColor(teeth_seg_color_[tid]);
					p_seg_tooth_mesh->SetChanged();
					p_seg_tooth_mesh->SetMeshColor(OpenMesh::Vec3d(0.8, 0.8, 0.8));

					seg_tooth_mesh_id_=DataPool::AddMeshObject(p_seg_tooth_mesh);
				}
				

			}
			meshobj->IsVisiable() = false;
			//DataPool::DeleteMeshObject(meshobj->GetId());
		}
		break;
	}
	case Qt::Key_U:
	{
		
		CMeshObject *meshobj = DataPool::GetMeshObject(dental_mesh_id_);
		if (meshobj != NULL)
		{
			meshobj->UseTexture() = !meshobj->UseTexture();
		}
		break;
	}
	case Qt::Key_P:
	{
		
		CMeshObject *meshobj = DataPool::GetMeshObject(dental_mesh_id_);
		if (meshobj != NULL)
		{
			COpenMeshT&mesh = meshobj->GetMesh();
			std::vector<bool>tags_teeth(teeth_seg_mark_.size(),false),tags_gingiva(teeth_seg_mark_.size(),false);
			for (int i = 0; i < teeth_seg_mark_.size(); i++)
			{
				if (teeth_seg_mark_[i] != -1)
				{
					tags_teeth[i] = true;
				}
					
			}
			for (int i = 0; i < gingiva_seg_mark_.size(); i++)
			{
				if (gingiva_seg_mark_[i] != -1)
				{
					tags_gingiva[i] = true;
					
				}
			}
			int mcount = 8;
			while (mcount--)
			{
				CMorphlogicOperation::Erode(mesh, tags_teeth);
			
			}
			mcount = 1;
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
			/*tags_gingiva.resize(mesh.n_vertices(), false);
			for (int i = 0; i < picked_vhs_back_.size(); i++)
			{
				cons.push_back(std::make_pair(picked_vhs_back_[i], 0.99));
				tags_gingiva[picked_vhs_back_[i].idx()] = true;
			}*/
			Eigen::VectorXd harmonic_u;
			CHarmonicFieldSeg harmonicSeg;
			harmonicSeg.ComputeConcavityAwareHarmonicField(mesh, cons, harmonic_u);
			static bool render_grad = false;
		//	render_grad = !render_grad;
			if (!render_grad)
			{
				for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); viter++)
				{
					OpenMesh::Vec3d color(1, 1, 1);
					if (tags_teeth[viter->idx()] == true  )
						color = OpenMesh::Vec3d(0, 0, 1);
					else if(tags_gingiva[viter->idx()] == true)
						color = OpenMesh::Vec3d(1, 1, 0);
					mesh.set_color(viter, color);

					for (auto hiter = mesh.vih_begin(viter); hiter != mesh.vih_end(viter); hiter++)
					{
						mesh.data(*hiter).SetUV(OpenMesh::Vec2f(harmonic_u(viter->idx()), harmonic_u(viter->idx())));
					}

				}
				meshobj->UseTexture() = true;
			}
			else
			{
				Eigen::MatrixXd u_grad;
				CGeoAlg::ComputeGradientOfScalarField(mesh, harmonic_u, u_grad);
				Eigen::VectorXd grad_mag=u_grad.rowwise().norm();
				CNumericalBaseAlg::NormalizeScalarField(grad_mag);
				for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); viter++)
				{
					OpenMesh::Vec3d color(1, 1, 1);
					if (tags_teeth[viter->idx()] == true || tags_gingiva[viter->idx()] == true)
						color = OpenMesh::Vec3d(0, 0, 1);
					mesh.set_color(viter, color);

					for (auto hiter = mesh.vih_begin(viter); hiter != mesh.vih_end(viter); hiter++)
					{
						mesh.data(*hiter).SetUV(OpenMesh::Vec2f(grad_mag(viter->idx()), grad_mag(viter->idx())));
					}

				}
				meshobj->UseTexture() = true;
			}
			
			/*for (int i = 0; i < teeth_edges.size(); i++)
			{
				mesh.set_color(teeth_edges[i], OpenMesh::Vec3d(1, 0, 0));
			}
			for (int i = 0; i < gingiva_edges.size(); i++)
			{
				mesh.set_color(gingiva_edges[i], OpenMesh::Vec3d(0, 1, 1));
			}*/
		
			meshobj->SetChanged();
		}
		break;
	}
	case Qt::Key_K:
	{
		
		CMeshObject *meshobj = DataPool::GetMeshObject(dental_mesh_id_);
		if (meshobj != NULL)
		{
			COpenMeshT&mesh = meshobj->GetMesh();
			CHarmonicFieldSeg harmonicSeg;
			std::vector<OpenMesh::VertexHandle>tooth_vhs,res_tooth,gingiva_vhs,res_gingiva;


			std::vector<bool>tags_teeth(teeth_seg_mark_.size(), false), tags_gingiva(teeth_seg_mark_.size(), false);
			for (int i = 0; i < teeth_seg_mark_.size(); i++)
			{
				if (teeth_seg_mark_[i] != -1)
				{
					tags_teeth[i] = true;
				}

			}
			for (int i = 0; i < gingiva_seg_mark_.size(); i++)
			{
				if (gingiva_seg_mark_[i] != -1)
				{
					tags_gingiva[i] = true;

				}
			}
			int mcount = 5;
			while (mcount--)
			{
				CMorphlogicOperation::Erode(mesh, tags_teeth);

			}
			mcount = 2;
			while (mcount--)
			{
				CMorphlogicOperation::Erode(mesh, tags_gingiva);
			}


			CGeoBaseAlg::GetEdgeVertexs(mesh, tags_teeth, tooth_vhs);
			CGeoBaseAlg::GetEdgeVertexs(mesh, tags_gingiva, gingiva_vhs);

			/*for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); viter++)
			{
				if (teeth_seg_mark_[viter->idx()] != -1)
				{
					tooth_vhs.push_back(viter);
				}
			}
			for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); viter++)
			{
				if (gingiva_seg_mark_[viter->idx()] != -1)
				{
					gingiva_vhs.push_back(viter);
				}
			}*/
			harmonicSeg.SegToothGingiva(mesh, tooth_vhs, gingiva_vhs, res_tooth, res_gingiva);
			teeth_seg_mark_.resize(mesh.n_vertices(), -1);
			for (int i = 0; i < res_tooth.size(); i++)
			{
				teeth_seg_mark_[res_tooth[i].idx()] = 1;
			}
			gingiva_seg_mark_.resize(mesh.n_vertices(), -1);
			for (int i = 0; i < res_gingiva.size(); i++)
			{
				gingiva_seg_mark_[res_gingiva[i].idx()] = 1;
			}

			for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); viter++)
			{
				int seg_id = teeth_seg_mark_[viter->idx()];

				OpenMesh::Vec3d color(1, 1, 1);
				if (seg_id != -1)
					color = OpenMesh::Vec3d(0, 0, 1);
				else if(gingiva_seg_mark_[viter->idx()]!=-1)
					color = OpenMesh::Vec3d(1, 1, 0);
				mesh.set_color(viter, color);
			}
			/*harmonicSeg.RefineToothGingivaSeg(mesh, tooth_vhs, res_tooth);
			teeth_seg_mark_.resize(mesh.n_vertices(), -1);
			for (int i = 0; i < res_tooth.size(); i++)
			{
				teeth_seg_mark_[res_tooth[i].idx()] = 1;
			}
			for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); viter++)
			{
				int seg_id = teeth_seg_mark_[viter->idx()];

				OpenMesh::Vec3d color(1, 1, 1);
				if (seg_id != -1)
					color = OpenMesh::Vec3d(0, 0, 1);
				mesh.set_color(viter, color);
			}*/
			meshobj->UseTexture() = false;
			meshobj->SetAttrChanged();
		}
		break;
	}
	case Qt::Key_Y://smooth
	{
		
		CMeshObject *meshobj = DataPool::GetMeshObject(dental_mesh_id_);
		if (meshobj != NULL)
		{
			COpenMeshT&mesh = meshobj->GetMesh();
			CGeoAlg::LaplacianSmooth(mesh, 5, 0.5);
			meshobj->SetChanged();
		}
		break;
	}
	case Qt::Key_F:
	{
		
		CMeshObject *meshobj = DataPool::GetMeshObject(dental_mesh_id_);
		if (meshobj != NULL)
		{
			std::cerr << "compute feature points" << std::endl;
			COpenMeshT&mesh = meshobj->GetMesh();
			std::vector<OpenMesh::VertexHandle>fvhs;
			//CDentalBaseAlg::ComputeTeethFeaturePoints(mesh, fvhs);
			CDentalBaseAlg::ComputeTeethFeaturePointsUsingSmoothedMesh(mesh, fvhs);
	
			picked_vhs_fore_.clear();
			picked_vhs_fore_mark_.clear();
		
			for (int i = 0; i < fvhs.size(); i++)
			{
				auto vh = fvhs[i];
			
				
				if (picked_vhs_fore_mark_.find(vh) == picked_vhs_fore_mark_.end())
				{
					picked_vhs_fore_.push_back(vh);
					picked_vhs_fore_mark_.insert(vh);
					if (picked_vhs_back_mark_.find(vh) != picked_vhs_back_mark_.end())
					{
						picked_vhs_back_mark_.erase(vh);
					}
				}
			}
			
			picked_vhs_back_.clear();
			for (auto iter = picked_vhs_back_mark_.begin(); iter != picked_vhs_back_mark_.end(); iter++)
			{
				picked_vhs_back_.push_back(*iter);
			}
			RenderFeature();
			meshobj->TextureId() = CUIContext::ColorBarTextureId();
			meshobj->UseTexture() = false;
			meshobj->SetChanged();
		
		
		}
		break;
	}

	case Qt::Key_V:
	{
		
		CMeshObject *meshobj = DataPool::GetMeshObject(dental_mesh_id_);
		if (meshobj != NULL)
		{
			COpenMeshT&mesh = meshobj->GetMesh();
			std::vector<OpenMesh::VertexHandle> gvhs;
		

			CDentalBaseAlg::ComputeGingivaVhs(mesh, gvhs);

		
			picked_vhs_back_.clear();
			picked_vhs_back_mark_.clear();

			
			for (int i = 0; i < gvhs.size(); i++)
			{
				auto vh = gvhs[i];
				
				if (picked_vhs_back_mark_.find(vh) == picked_vhs_back_mark_.end())
				{
					picked_vhs_back_.push_back(vh);
					picked_vhs_back_mark_.insert(vh);
					if (picked_vhs_fore_mark_.find(vh) != picked_vhs_fore_mark_.end())
					{
						picked_vhs_fore_mark_.erase(vh);
					}
				}
			}
			picked_vhs_fore_.clear();
			for (auto iter = picked_vhs_fore_mark_.begin(); iter != picked_vhs_fore_mark_.end(); iter++)
			{
				picked_vhs_fore_.push_back(*iter);
			}

			meshobj->TextureId() = CUIContext::ColorStripeTextureId();
			meshobj->UseTexture() = false;
			meshobj->SetChanged();
			RenderFeature();

		}
		break;
	}
	/*case Qt::Key_L:
	{
		QString path = QFileDialog::getOpenFileName(NULL, "load dental mesh", ".", "Files(*.stl *.off *.obj)");

		if (path.length() == 0)
		{
			std::cerr << "unable to load mesh\n" << std::endl;
			return;
		}
		CMeshObject *meshobj = new CMeshObject();
		COpenMeshT &mesh = meshobj->GetMesh();
		if (!CDataIO::ReadMesh(path.toStdString(), *meshobj))
		{
			std::cerr << "unable to load mesh\n" << std::endl;
		}
		path[path.length() - 3] = 'd';
		path[path.length() - 2] = 'm';
		path[path.length() - 1] = 'a';
		path = path + 't';
		Eigen::VectorXd scalars;
		igl::readDMAT(path.toStdString(), scalars);
		CNumericalBaseAlg::NormalizeScalarField(scalars);
		for (auto hiter = mesh.halfedges_begin(); hiter != mesh.halfedges_end(); hiter++)
		{
			auto vh = mesh.to_vertex_handle(hiter);
			mesh.data(*hiter).SetUV(OpenMesh::Vec2f(scalars[vh.idx()], scalars[vh.idx()]));
		}
		std::vector<COpenMeshT::VertexHandle>vhs;
		CGeoAlg::ComputeLocalExtremum(mesh, scalars, 3, vhs);
		std::vector<COpenMeshT::VertexHandle> features;
		for (int i = 0; i < vhs.size(); i++)
		{
			if ((!CGeoBaseAlg::IsConcavity(mesh, vhs[i])) && (!mesh.is_boundary(vhs[i])))
			{
				bool flag = true;
				std::vector<COpenMeshT::VertexHandle>neivhs;
				CGeoAlg::ExtractNRing(mesh, vhs[i], 1, neivhs);
				for (int j=0;j<neivhs.size();j++)
				{
					if (CGeoBaseAlg::IsConcavity(mesh, neivhs[j]))
					{
						flag = false;
						break;
					}
					
				}
				if (flag)
				{
					features.push_back(vhs[i]);
				}
			
			}
			
		}
		meshobj->SetMeshColor(OpenMesh::Vec3d(0.8, 0.8, 0.8));
		picked_vhs_fore_.clear();
		picked_vhs_fore_mark_.clear();
		OpenMesh::Vec3d mean;
		std::vector<OpenMesh::Vec3d>frame;
		CDentalBaseAlg::ComputePCAFrameFromHighCurvaturePoints(mesh, 10, mean, frame);
		std::vector<double>feature_z(features.size());
		OpenMesh::Vec3d bound_mean, dir;
		CDentalBaseAlg::ComputeHoleVertMean(mesh, bound_mean);
		bound_mean = bound_mean - mean;
		dir = frame[2];
		if (dir[0] * bound_mean[0] + dir[1] * bound_mean[1] + dir[2] * bound_mean[2] > 0)
			dir = -dir;
		dir = dir.normalize();
		for (int i = 0; i < features.size(); i++)
		{
			auto p = mesh.point(features[i]);
			auto dirp = p - mean;
			feature_z[i] = dirp[0] * dir[0] + dirp[1] * dir[1] + dirp[2] * dir[2];
		}
		int max_id;
		CNumericalBaseAlg::GetMaxValue(feature_z, max_id);
		for (int i = 0; i < features.size(); i++)
		{
			if (feature_z[max_id] - feature_z[i] < 0.1)
			{
				mesh.set_color(features[i], OpenMesh::Vec3d(1, 0, 0));
				for (auto vviter = mesh.vv_begin(features[i]); vviter != mesh.vv_end(features[i]); vviter++)
				{
					mesh.set_color(vviter, OpenMesh::Vec3d(1, 0, 0));
				}
				picked_vhs_fore_.push_back(features[i]);
				picked_vhs_fore_mark_.insert(features[i]);
			}
			
		}
		meshobj->UseTexture() = false;
		meshobj->TextureId() = CUIContext::ColorStripeTextureId();
		int id=DataPool::AddMeshObject(meshobj);
		CUIContext::SetSelectedMeshObjectId(id);
		break;
	}*/

		case Qt::Key_L:
		{
			//	std::cerr << "l" << std::endl;
			QString path = QFileDialog::getOpenFileName(NULL, "load dental mesh", ".", "Stl Files(*.stl)");

			if (path.length() == 0)
			{
				std::cerr << "unable to load mesh\n" << std::endl;
				return;
			}
			CMeshObject *meshobj = new CMeshObject();
			if (!CDataIO::ReadMesh(path.toStdString(), *meshobj))
			{
				std::cerr << "unable to load mesh\n" << std::endl;
			}
			/*int mid = DataPool::AddMeshObject(meshobj);
			CUIContext::SetSelectedMeshObjectId(mid);
			break;*/
			COpenMeshT &mesh = meshobj->GetMesh();
			std::vector<COpenMeshT::FaceHandle>res_fhs;
			std::vector<OpenMesh::Vec3d>res_bary_coords;



			std::vector<COpenMeshT*>resmeshes;
			std::cerr << "separate disconnected parts" << std::endl;
			CGeoAlg::SeparateDisconnectedParts(mesh, resmeshes, std::vector<std::map<COpenMeshT::VertexHandle, COpenMeshT::VertexHandle>>());
			std::cerr << "separate disconnected parts finish" << std::endl;
			int max_vnum = -1;
			int mi;
			for (int i = 0; i < resmeshes.size(); i++)
			{
				int vnum = resmeshes[i]->n_vertices();
				//std::cerr << "num " << resmeshes[i]->n_vertices() <<" max"<<max_vnum<<" flag"<<(max_vnum<resmeshes[i]->n_vertices())<< std::endl;
				if (max_vnum <vnum)
				{

					max_vnum = vnum;
					mi = i;
				}
			}

			mesh = *resmeshes[mi];
			std::cerr << mesh.n_vertices() << std::endl;
			for (int i = 0; i < resmeshes.size(); i++)
			{
				delete resmeshes[i];
			}


			//CGeoAlg::FillHoles(mesh);
			CGeoBaseAlg::RemoveNonManifold(mesh);
			CGeoAlg::SelfIntersectionRemoval(mesh);
			CGeoAlg::FillHoles(mesh, true);
			CGeoAlg::SimplifyMesh(mesh, 160000);

			//CGeoAlg::LaplacianSmooth(mesh, 20, 0.5);
			std::cerr << "vnum " << mesh.n_vertices() << std::endl;
			std::cerr << "fnum " << mesh.n_faces() << std::endl;
			std::cerr << "enum " << mesh.n_edges() << std::endl;
			for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); viter++)
			{
				mesh.set_color(viter, OpenMesh::Vec3d(0.8, 0.8, 0.8));
			}
			meshobj->RestoreCurrentVPos();
			CDentalBaseAlg::PCABasedOrientationCorrection(meshobj->GetMesh());
			CGeoBaseAlg::NormalizeMeshSize(mesh);
			meshobj->SetChanged();
			int id = DataPool::AddMeshObject(meshobj);
			CUIContext::SetSelectedMeshObjectId(id);

			dental_mesh_id_ = CUIContext::GetSelectedMeshObjectId();
			
			seg_tooth_mesh_id_ = -1;

			break;
		}

	case Qt::Key_I:
	{
	
		std::map<int, int>&vtags = p_mesh_obj_->GetVertexTags();

		CMeshObject *dental_mesh_obj = DataPool::GetMeshObject(dental_mesh_id_);
		if (dental_mesh_obj != NULL)
		{
			std::cerr << "refine teeth seg" << std::endl;
			COpenMeshT&mesh = dental_mesh_obj->GetMesh();
			CHarmonicFieldSeg harmonic_seg;
			harmonic_seg.RefineTeethGingivaSeg(mesh, teeth_seg_mark_);
			for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); viter++)
			{
				int vid = viter->idx();
				
				int tid = teeth_seg_mark_[vid];
				if(tid!=-1)
				mesh.set_color(viter, teeth_seg_color_[tid]);
				else
				{
					mesh.set_color(viter,OpenMesh::Vec3d(0.8,0.8,0.8));
				}
				vtags[viter->idx()] = tid;
			
			}
			
			dental_mesh_obj->SetAttrChanged();

			
		}
		break;
	}
	case Qt::Key_G:
	{
		std::map<int, int>&vtags = p_mesh_obj_->GetVertexTags();

		CMeshObject *tooth_mesh_obj = DataPool::GetMeshObject(seg_tooth_mesh_id_);
		if (tooth_mesh_obj != NULL)
		{
			COpenMeshT&tooth_mesh = tooth_mesh_obj->GetMesh();
		
			std::vector<std::vector<COpenMeshT::VertexHandle>>inside_vhs, outside_vhs;
			CDentalBaseAlg::ComputeBoundCuttingPointsOfToothMesh(*tooth_mesh_obj, inside_vhs, outside_vhs);
			tooth_mesh_obj->UseTexture() = false;
			tooth_mesh_obj->SetAttrChanged();
		
			CDentalBaseAlg::MergeCuttingPointsByDis(tooth_mesh, inside_vhs, outside_vhs, 0.02);
			std::vector<CDentalBaseAlg::CCuttingPath>cutting_paths;
			CDentalBaseAlg::ComputeCuttingPath(*tooth_mesh_obj, inside_vhs, outside_vhs, cutting_paths);
			for (int i = 0; i < cutting_paths.size(); i++)
			{
				CCurveObject *ccobj = new CCurveObject();
				cutting_paths[i].GetPathPoints(ccobj->GetCurve());
				ccobj->SetChanged();
				ccobj->SetColor(OpenMesh::Vec3d(1, 0, 0));
				DataPool::AddCurveObject(ccobj);
			}
			std::vector<int>tooth_tags;
			
			int tagnum=CDentalBaseAlg::TagToothByCuttingPath(tooth_mesh, cutting_paths, tooth_tags);
		

			CMeshObject *dental_mesh_obj = DataPool::GetMeshObject(dental_mesh_id_);
			COpenMeshT &dental_mesh = dental_mesh_obj->GetMesh();
			dental_mesh_obj->SetMeshColor(OpenMesh::Vec3d(0.8, 0.8, 0.8));

			if (teeth_seg_mark_.size() != dental_mesh.n_vertices())
			{
				teeth_seg_mark_.resize(dental_mesh.n_vertices(), -1);
				teeth_seg_color_.clear();
				rand_color_set_.clear();
				teeth_seg_count_.clear();
				teeth_seg_mark_id_ = 0;
			}
		
			for (auto viter = tooth_mesh.vertices_begin(); viter != tooth_mesh.vertices_end(); viter++)
			{
				int tid = tooth_tags[viter->idx()];
				vtags[all_seg_tooth_orig_vhs_map_[viter].idx()] = tid;
				teeth_seg_mark_[all_seg_tooth_orig_vhs_map_[viter].idx()] = tid;
				//std::cerr << all_seg_tooth_orig_vhs_map_[viter].idx()<<" "<<teeth_seg_mark_[all_seg_tooth_orig_vhs_map_[viter].idx()] << std::endl;
				if (teeth_seg_mark_id_ < tid)
					teeth_seg_mark_id_ = tid;

				if (teeth_seg_count_.find(tid) == teeth_seg_count_.end())
					teeth_seg_count_[tid] = 0;
				else
					teeth_seg_count_[tid]++;


				if (teeth_seg_color_.find(tid) == teeth_seg_color_.end())
				{
					auto c = GetRandColor();
					if (rand_color_set_.find(c) != rand_color_set_.end())
						c = GetRandColor();
					teeth_seg_color_[tid] = c;
					rand_color_set_.insert(c);
				}

				if (tid == -1)
				{
					tooth_mesh.set_color(viter, OpenMesh::Vec3d(0.8, 0.8, 0.8));
					dental_mesh.set_color(all_seg_tooth_orig_vhs_map_[viter], OpenMesh::Vec3d(0.8, 0.8, 0.8));
				}
				else
				{
					
						
					
					tooth_mesh.set_color(viter, teeth_seg_color_[tid]);
					dental_mesh.set_color(all_seg_tooth_orig_vhs_map_[viter], teeth_seg_color_[tid]);
				}
			}
		
			
			dental_mesh_obj->IsVisiable() = true;
			dental_mesh_obj->UseTexture() = false;
			dental_mesh_obj->SetAttrChanged();

			tooth_mesh_obj->IsVisiable() = false;
			tooth_mesh_obj->UseTexture() =false;
			tooth_mesh_obj->SetAttrChanged();
		}
		
		
		break;
	}
	case Qt::Key_X:
	{
		p_mesh_obj_ = DataPool::GetMeshObject(dental_mesh_id_);
		//p_mesh_obj_->TextureId() = CUIContext::ColorBarTextureId();
		//if (p_mesh_obj_->TextureId() == -1)
		//{
		//	std::cerr << "color bar texture is NULL" << std::endl;
		//}
		//else
		//{
		//	p_mesh_obj_->UseTexture() = true;
		//}
		COpenMeshT &mesh = p_mesh_obj_->GetMesh();
	
		teeth_seg_mark_.resize(mesh.n_vertices(), -1);
		gingiva_seg_mark_.resize(mesh.n_vertices(), -1);

		teeth_seg_mark_id_ = 0;
	
		CHarmonicFieldSeg harmonicSeg;
		std::vector<COpenMeshT::VertexHandle>res_vhs_tooth,res_vhs_gingiva;
		harmonicSeg.SegToothGingiva(mesh, picked_vhs_fore_, picked_vhs_back_, res_vhs_tooth, res_vhs_gingiva);

		for (int i = 0; i < teeth_seg_mark_.size(); i++)
		{
			teeth_seg_mark_[i] = -1;
			gingiva_seg_mark_[i] = -1;
		}
			
		for (int i = 0; i < res_vhs_tooth.size(); i++)
		{
			int vid = res_vhs_tooth[i].idx();
			teeth_seg_mark_[vid] = 1;
		}
		for (int i = 0; i < res_vhs_gingiva.size(); i++)
		{
			int vid = res_vhs_gingiva[i].idx();
			gingiva_seg_mark_[vid] = 1;
		}
		for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); viter++)
		{
			int seg_id = teeth_seg_mark_[viter->idx()];

			OpenMesh::Vec3d color(0.8,0.8,0.8);
			if (seg_id != -1)
				color = OpenMesh::Vec3d(0,0,1);
			else if(gingiva_seg_mark_[viter->idx()]!=-1)
				color= OpenMesh::Vec3d(1, 1, 0);
			mesh.set_color(viter, color);
		}
		

		p_mesh_obj_->UseTexture() = false;
		p_mesh_obj_->SetAttrChanged();
		break;
	}
	case Qt::Key_R:
	{
		p_mesh_obj_ = DataPool::GetMeshObject(dental_mesh_id_);


		/*p_mesh_obj_->TextureId() = CUIContext::ColorBarTextureId();
		if (p_mesh_obj_->TextureId() == -1)
		{
			std::cerr << "color bar texture is NULL" << std::endl;
		}
		else
		{
			p_mesh_obj_->UseTexture() = true;
		}*/
		std::map<int, int>&vtags = p_mesh_obj_->GetVertexTags();

	
		COpenMeshT &mesh = p_mesh_obj_->GetMesh();
		if (teeth_seg_mark_.size() != mesh.n_vertices())
		{
			teeth_seg_mark_.resize(mesh.n_vertices(), -1);
			teeth_seg_color_.clear();
			rand_color_set_.clear();
			teeth_seg_count_.clear();
			teeth_seg_mark_id_ = 0;
		}
		CHarmonicFieldSeg harmonicSeg;
		std::vector<COpenMeshT::VertexHandle>res_vhs0, res_vhs1;
		harmonicSeg.SegTwoTooth(mesh, picked_vhs_fore_, picked_vhs_back_, res_vhs0, res_vhs1);
		for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); viter++)
		{
			for (auto hiter = mesh.vih_begin(viter); hiter != mesh.vih_end(viter); hiter++)
			{
				mesh.data(hiter).SetUV(OpenMesh::Vec2f(0, 0));
			}
		}
		OpenMesh::Vec3d color0, color1;
		std::map<int, int>overlap_mark_count;
		do {
			color0 = GetRandColor();
		} while (rand_color_set_.find(color0) != rand_color_set_.end());
		rand_color_set_.insert(color0);
		teeth_seg_color_[teeth_seg_mark_id_] = color0;
		teeth_seg_count_[teeth_seg_mark_id_] = 0;
		for (int i = 0; i < res_vhs0.size(); i++)
		{
			auto vh = res_vhs0[i];
			int pre_mark = teeth_seg_mark_[vh.idx()];
			if (pre_mark != -1)
			{
				if (overlap_mark_count.find(pre_mark) == overlap_mark_count.end())
				{
					overlap_mark_count[pre_mark] = 1;
				}
				else
				{
					overlap_mark_count[pre_mark] = overlap_mark_count[pre_mark] + 1;
				}
			}
			teeth_seg_mark_[vh.idx()] = teeth_seg_mark_id_;
			vtags[vh.idx()] = teeth_seg_mark_id_;
			teeth_seg_count_[teeth_seg_mark_id_]= teeth_seg_count_[teeth_seg_mark_id_]+1;
		}

		teeth_seg_mark_id_++;
		teeth_seg_count_[teeth_seg_mark_id_] = 0;
		do {
			color1 = GetRandColor();
		} while (rand_color_set_.find(color1) != rand_color_set_.end());
		rand_color_set_.insert(color1);
		teeth_seg_color_[teeth_seg_mark_id_] = color1;
		for (int i = 0; i < res_vhs1.size(); i++)
		{
			auto vh = res_vhs1[i];
			int pre_mark = teeth_seg_mark_[vh.idx()];
			if (pre_mark != -1)
			{
				if (overlap_mark_count.find(pre_mark) == overlap_mark_count.end())
				{
					overlap_mark_count[pre_mark] = 1;
				}
				else
				{
					overlap_mark_count[pre_mark]= overlap_mark_count[pre_mark]+1;
				}
			}
			teeth_seg_mark_[vh.idx()] = teeth_seg_mark_id_;
			vtags[vh.idx()] = teeth_seg_mark_id_;
			teeth_seg_count_[teeth_seg_mark_id_] = teeth_seg_count_[teeth_seg_mark_id_] + 1;
		}
		teeth_seg_mark_id_++;
	
		for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); viter++)
		{
			int seg_id = teeth_seg_mark_[viter->idx()];
			if (overlap_mark_count.find(seg_id)!=overlap_mark_count.end())
			{
			//	std::cerr << overlap_mark_count[seg_id]<<" "<< teeth_seg_count_[seg_id] << std::endl;
				if (overlap_mark_count[seg_id]*1.0 / teeth_seg_count_[seg_id] > 0.5)
				{
					teeth_seg_count_[seg_id] = 0;
					seg_id = -1;
					teeth_seg_mark_[viter->idx()] = -1;
					vtags[viter->idx()] = -1;
				}
			
			}
			OpenMesh::Vec3d color(1,1,1);
			if(seg_id !=-1)
				color = teeth_seg_color_[seg_id];
			mesh.set_color(viter, color);
		}
		/*for (int i = 0; i < res_vhs0.size(); i++)
		{
			for (auto hiter = mesh.vih_begin(res_vhs0[i]); hiter != mesh.vih_end(res_vhs0[i]); hiter++)
			{
				mesh.data(hiter).SetUV(OpenMesh::Vec2f(0.5, 0.5));
			}
		}
		for (int i = 0; i < res_vhs1.size(); i++)
		{
			for (auto hiter = mesh.vih_begin(res_vhs1[i]); hiter != mesh.vih_end(res_vhs1[i]); hiter++)
			{
				mesh.data(hiter).SetUV(OpenMesh::Vec2f(0.8, 0.8));
			}
		}*/
		
		p_mesh_obj_->UseTexture() = false;
		p_mesh_obj_->SetAttrChanged();
		break;
	}
	/*case Qt::Key_F:
	{
		p_mesh_obj_ = DataPool::GetMeshObject(dental_mesh_id_);


		p_mesh_obj_->TextureId() = CUIContext::ColorBarTextureId();
		if (p_mesh_obj_->TextureId() == -1)
		{
			std::cerr << "color bar texture is NULL" << std::endl;
		}
		else
		{
			p_mesh_obj_->UseTexture() = true;
		}
		COpenMeshT &mesh = p_mesh_obj_->GetMesh();
		CHarmonicFieldSeg harmonicSeg;
		std::vector<COpenMeshT::VertexHandle>res_vhs;
		harmonicSeg.SegOneTeeth(mesh, picked_vhs_fore_, res_vhs);
		for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); viter++)
		{
			for (auto hiter = mesh.vih_begin(viter); hiter != mesh.vih_end(viter); hiter++)
			{
				mesh.data(hiter).SetUV(OpenMesh::Vec2f(0, 0));
			}
		}
		for (int i = 0; i < res_vhs.size(); i++)
		{
			for (auto hiter = mesh.vih_begin(res_vhs[i]); hiter != mesh.vih_end(res_vhs[i]); hiter++)
			{
				mesh.data(hiter).SetUV(OpenMesh::Vec2f(0.5, 0.5));
			}
		}
		
		p_mesh_obj_->SetAttrChanged();
		break;
	}*/
	case Qt::Key_N:
	{
		is_eliminating_feature_ = true;
		break;
	}
	//case Qt::Key_D:
	//{

	//	p_mesh_obj_ = DataPool::GetMeshObject(dental_mesh_id_);


	//	p_mesh_obj_->TextureId() = CUIContext::ColorStripeTextureId();
	//	if (p_mesh_obj_->TextureId() == -1)
	//	{
	//		std::cerr << "color bar texture is NULL" << std::endl;
	//	}
	//	else
	//	{
	//		p_mesh_obj_->UseTexture() = true;
	//	}
	//	COpenMeshT &mesh = p_mesh_obj_->GetMesh();
	//	CHarmonicFieldSeg harmonicSeg;
	//	std::vector<std::pair<COpenMeshT::VertexHandle, double>>cons;
	//	for (int i = 0; i < picked_vhs_fore_.size(); i++)
	//	{
	//		cons.push_back(std::pair<COpenMeshT::VertexHandle, double>(picked_vhs_fore_[i], 1));
	//	}
	//
	//	std::set<COpenMeshT::VertexHandle>vhset;
	//	for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); viter++)
	//	{
	//	if (mesh.is_boundary(viter))
	//	{
	//	cons.push_back(std::make_pair(viter, 0));
	//	for (auto vviter = mesh.vv_begin(viter); vviter != mesh.vv_end(viter); vviter++)
	//	{
	//	if (mesh.is_boundary(vviter)==false&&vhset.find(vviter) == vhset.end())
	//	{
	//	vhset.insert(vviter);
	//	cons.push_back(std::make_pair(vviter,0));
	//	}
	//	}
	//	}
	//	}
	//	Eigen::VectorXd res_u;
	//	harmonicSeg.ComputeConcavityAwareHarmonicField(mesh, cons, res_u);
	//	//p_mesh_obj_->NormalizeUV();
	//	//igl::writeDMAT("res_u.dmat", res_u);
	//	for (auto hiter = mesh.halfedges_begin(); hiter != mesh.halfedges_end(); hiter++)
	//	{
	//		auto vh = mesh.to_vertex_handle(hiter);
	//		mesh.data(*hiter).SetUV(OpenMesh::Vec2f(res_u(vh.idx()), res_u(vh.idx())));
	//	}
	//	p_mesh_obj_->SetAttrChanged();
	//	break;
	//}
	case Qt::Key_E:
	{
		
		p_mesh_obj_ = DataPool::GetMeshObject(dental_mesh_id_);
		

		p_mesh_obj_->TextureId() = CUIContext::ColorStripeTextureId();
		if (p_mesh_obj_->TextureId() == -1)
		{
			std::cerr << "color bar texture is NULL" << std::endl;
		}
		else
		{
			p_mesh_obj_->UseTexture() = true;
		}
		COpenMeshT &mesh = p_mesh_obj_->GetMesh();
		CHarmonicFieldSeg harmonicSeg;
		std::vector<bool>mmark(mesh.n_vertices(),false);
		std::vector<std::pair<COpenMeshT::VertexHandle, double>>cons;
		for (int i = 0; i < picked_vhs_fore_.size(); i++)
		{
			mmark[picked_vhs_fore_[i].idx()] = true;
			cons.push_back(std::pair<COpenMeshT::VertexHandle, double>(picked_vhs_fore_[i], 0.1));
			std::vector<COpenMeshT::VertexHandle>neis;
			CGeoAlg::ExtractNRing(mesh, picked_vhs_fore_[i], 2, neis);
			for (int j = 0; j < neis.size(); j++)
			{
				if (mmark[neis[j].idx()])
					continue;
				else
				{
					mmark[neis[j].idx()] = true;
					cons.push_back(std::pair<COpenMeshT::VertexHandle, double>(neis[j], 0.1));
				}
			}
		}
		std::cerr << "fore size " << picked_vhs_fore_.size()<<" "<< picked_vhs_back_.size() << std::endl;
		for (int i = 0; i < picked_vhs_back_.size(); i++)
		{
			if (mmark[picked_vhs_back_[i].idx()])
				continue;
			cons.push_back(std::pair<COpenMeshT::VertexHandle, double>(picked_vhs_back_[i], 0.9));
		}
		/*std::set<COpenMeshT::VertexHandle>vhset;
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
		}*/
		Eigen::VectorXd res_u;
		harmonicSeg.ComputeConcavityAwareHarmonicField(mesh, cons, res_u);
		//p_mesh_obj_->NormalizeUV();
		//igl::writeDMAT("res_u.dmat", res_u);
		for (auto hiter = mesh.halfedges_begin(); hiter != mesh.halfedges_end(); hiter++)
		{
			auto vh = mesh.to_vertex_handle(hiter);
			mesh.data(*hiter).SetUV(OpenMesh::Vec2f(res_u(vh.idx()), res_u(vh.idx())));
		}
		p_mesh_obj_->SetAttrChanged();
		break;
	}
	}
}
void CHarmonicFieldSegmentation::KeyReleaseEvent(QKeyEvent *e)
{
	switch (e->key())
	{

		case Qt::Key_W:
		{
			is_drawing_ = false;
			break;
		}
		case Qt::Key_Z:
		{
			is_seg_teeth_from_gingiva_ = false;
			is_drawing_ = false;
			break;
		}
		case Qt::Key_N:
		{
			is_eliminating_feature_ = false;
			break;
		}
		case Qt::Key_A:
		{
			is_seg_teeth_from_gingiva_ = false;
			is_drawing_ = false;
			break;
		}
	}
}