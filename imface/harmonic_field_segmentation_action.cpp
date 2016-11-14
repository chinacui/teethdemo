#include"harmonic_field_segmentation_action.h"
#include"../AlgColle/geo_alg.h"
#include"ui_context.h"
#include"../DataColle/data_pool.h"
#include"cmodelviewer.h"
#include"camera.h"
#include"../AlgColle/geo_base_alg.h"
#include<qimage.h>
#include"../AlgColle/harmonic_field.h"
#include<igl/writeDMAT.h>
#include<igl/readDMAT.h>
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
void CHarmonicFieldSegmentation::MousePressEvent(QMouseEvent *e)
{
	if (is_drawing_)
	{
		int mid = CUIContext::GetSelectedMeshObjectId();
		CMeshObject *meshobj = DataPool::GetMeshObject(mid);
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
			int mid = CUIContext::GetSelectedMeshObjectId();
			CMeshObject *meshobj = DataPool::GetMeshObject(mid);
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
			int mid = CUIContext::GetSelectedMeshObjectId();
			CMeshObject *meshobj = DataPool::GetMeshObject(mid);
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
			int mid = CUIContext::GetSelectedMeshObjectId();
			CMeshObject *meshobj = DataPool::GetMeshObject(mid);
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
	int mid = CUIContext::GetSelectedMeshObjectId();
	CMeshObject *meshobj = DataPool::GetMeshObject(mid);
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
		int mid = CUIContext::GetSelectedMeshObjectId();
		CMeshObject *meshobj = DataPool::GetMeshObject(mid);
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

	int mid = CUIContext::GetSelectedMeshObjectId();
	CMeshObject *meshobj = DataPool::GetMeshObject(mid);
	if (meshobj != NULL)
	{
		meshobj->SetMeshColor(OpenMesh::Vec3d(0.8, 0.8, 0.8));
		meshobj->SetAttrChanged();
	}
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
		manager_->SetCurrentActionType(CActionType::Common);
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
		std::cerr << "clear" << std::endl;
		break;
	}
	case Qt::Key_T:
	{
		is_pick_fore_ = !is_pick_fore_;
		std::cerr << "is pick fore " << is_pick_fore_ << std::endl;
		break;
	}
	case Qt::Key_S:
	{
		int mid = CUIContext::GetSelectedMeshObjectId();
		CMeshObject *meshobj = DataPool::GetMeshObject(mid);
		if (meshobj != NULL)
		{
			COpenMeshT&mesh = meshobj->GetMesh();
		
			CDentalBaseAlg::ComputePCAFrameFromHighCurvaturePoints(mesh,10, pca_mean_,pca_frame_);
			//OpenMesh::Vec3d dir = pca_frame_[2];
			//CPlane cplane(pca_mean_, dir);
			//CMeshObject *plane_obj = new CMeshObject();
			//CAuxGeoUtils::GetPlaneMeshFromPointAndAxis(cplane.p(), pca_frame_[0], pca_frame_[1], pca_frame_[2], 1.2, plane_obj->GetMesh());
			//plane_obj->SetChanged();
			//DataPool::AddMeshObject(plane_obj);

			std::map<int, COpenMeshT*>sep_meshes;
			std::map<int, std::map<COpenMeshT::VertexHandle, COpenMeshT::VertexHandle>>vid_orig;
			//CGeoAlg::SeparateMeshByVertexTag(mesh, teeth_seg_mark_, sep_meshes, vid_orig);
			CGeoAlg::SeparateMeshByVertexTag(mesh, gingiva_seg_mark_, sep_meshes, vid_orig);
			for(auto iter=sep_meshes.begin();iter!=sep_meshes.end();iter++)
			{
				int tid = iter->first;
				if (tid == -1)
				{
					CMeshObject* p_mesh_obj = new CMeshObject();
				/*	std::vector<COpenMeshT*>res_meshes;
					CGeoAlg::SeparateDisconnectedParts(*sep_meshes[tid], res_meshes, std::vector<std::map<COpenMeshT::VertexHandle, COpenMeshT::VertexHandle>>());
					int max_i=0;
					for (int i = 0; i < res_meshes.size(); i++)
					{
						if (res_meshes[i]->n_vertices() > res_meshes[max_i]->n_vertices())
						{
							max_i = i;
						}
					}
					p_mesh_obj->GetMesh() = *res_meshes[max_i];*/
					CGeoBaseAlg::RemoveNonManifold(*iter->second);
					p_mesh_obj->GetMesh() = *iter->second;
				//	CGeoAlg::SimplifyMesh(p_mesh_obj->GetMesh(), 10000);
					//CGeoAlg::FillHoles(p_mesh_obj->GetMesh());
					delete sep_meshes[tid];
					p_mesh_obj->SetMeshColor(teeth_seg_color_[tid]);
					p_mesh_obj->SetChanged();
					p_mesh_obj->SetMeshColor(OpenMesh::Vec3d(0.8, 0.8, 0.8));
					CUIContext::SetSelectedMeshObjectId(DataPool::AddMeshObject(p_mesh_obj));
				}
				

			}
			DataPool::DeleteMeshObject(meshobj->GetId());
		}
		break;
	}
	case Qt::Key_U:
	{
		int mid = CUIContext::GetSelectedMeshObjectId();
		CMeshObject *meshobj = DataPool::GetMeshObject(mid);
		if (meshobj != NULL)
		{
			meshobj->UseTexture() = !meshobj->UseTexture();
		}
		break;
	}
	case Qt::Key_P:
	{
		int mid = CUIContext::GetSelectedMeshObjectId();
		CMeshObject *meshobj = DataPool::GetMeshObject(mid);
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
		int mid = CUIContext::GetSelectedMeshObjectId();
		CMeshObject *meshobj = DataPool::GetMeshObject(mid);
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
		int mid = CUIContext::GetSelectedMeshObjectId();
		CMeshObject *meshobj = DataPool::GetMeshObject(mid);
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
		int mid = CUIContext::GetSelectedMeshObjectId();
		CMeshObject *meshobj = DataPool::GetMeshObject(mid);
		if (meshobj != NULL)
		{
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
		int mid = CUIContext::GetSelectedMeshObjectId();
		CMeshObject *meshobj = DataPool::GetMeshObject(mid);
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
	case Qt::Key_L:
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
	}
	case Qt::Key_G:
	{
		int mid = CUIContext::GetSelectedMeshObjectId();
		CMeshObject *meshobj = DataPool::GetMeshObject(mid);
		if (meshobj != NULL)
		{
			COpenMeshT&mesh = meshobj->GetMesh();
			std::vector<std::vector<COpenMeshT::VertexHandle>>inside_vhs, outside_vhs;
			CDentalBaseAlg::ComputeBoundCuttingPointsOfToothMesh(*meshobj, inside_vhs, outside_vhs);
			meshobj->UseTexture() = false;
			meshobj->SetAttrChanged();
		
			CDentalBaseAlg::MergeCuttingPointsByDis(mesh, inside_vhs, outside_vhs, 0.02);
			std::vector<CDentalBaseAlg::CCuttingPath>cutting_paths;
			CDentalBaseAlg::ComputeCuttingPath(*meshobj, inside_vhs, outside_vhs, cutting_paths);
			for (int i = 0; i < cutting_paths.size(); i++)
			{
				CCurveObject *ccobj = new CCurveObject();
				cutting_paths[i].GetPathPoints(ccobj->GetCurve());
				ccobj->SetChanged();
				ccobj->SetColor(OpenMesh::Vec3d(1, 0, 0));
				DataPool::AddCurveObject(ccobj);
			}
			std::vector<int>tooth_tags;
			
			int tagnum=CDentalBaseAlg::TagToothByCuttingPath(mesh, cutting_paths, tooth_tags);
			std::map<int,OpenMesh::Vec3d>tag_colors;
			for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); viter++)
			{
				int tid = tooth_tags[viter->idx()];
				if (tid == -1)
				{
					mesh.set_color(viter, OpenMesh::Vec3d(0.8, 0.8, 0.8));
				}
				else
				{
					if (tag_colors.find(tid) == tag_colors.end())
						tag_colors[tid] = GetRandColor();
					mesh.set_color(viter, tag_colors[tid]);
				}
			}
		


			meshobj->UseTexture() =false;
			meshobj->SetAttrChanged();
		}
		
		
		break;
	}
	case Qt::Key_X:
	{
		p_mesh_obj_ = DataPool::GetMeshObject(CUIContext::GetSelectedMeshObjectId());
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
		p_mesh_obj_ = DataPool::GetMeshObject(CUIContext::GetSelectedMeshObjectId());


		/*p_mesh_obj_->TextureId() = CUIContext::ColorBarTextureId();
		if (p_mesh_obj_->TextureId() == -1)
		{
			std::cerr << "color bar texture is NULL" << std::endl;
		}
		else
		{
			p_mesh_obj_->UseTexture() = true;
		}*/
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
		p_mesh_obj_ = DataPool::GetMeshObject(CUIContext::GetSelectedMeshObjectId());


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
	case Qt::Key_D:
	{

		p_mesh_obj_ = DataPool::GetMeshObject(CUIContext::GetSelectedMeshObjectId());


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
		std::vector<std::pair<COpenMeshT::VertexHandle, double>>cons;
		for (int i = 0; i < picked_vhs_fore_.size(); i++)
		{
			cons.push_back(std::pair<COpenMeshT::VertexHandle, double>(picked_vhs_fore_[i], 1));
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
	case Qt::Key_E:
	{
		
		p_mesh_obj_ = DataPool::GetMeshObject(CUIContext::GetSelectedMeshObjectId());
		

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