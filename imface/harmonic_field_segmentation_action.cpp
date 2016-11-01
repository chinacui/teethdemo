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
#include<set>
#include<time.h>
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
				picked_vhs_fore_.clear();
				picked_vhs_fore_mark_.clear();
			}
			else
			{
				curve_color = OpenMesh::Vec3d(1, 0, 0);
				picked_vhs_back_.clear();
				picked_vhs_back_mark_.clear();
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
	if (is_drawing_)
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
	teeth_seg_mark_id_ = 0;
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
		for (int i = 0; i < res_vhs0.size(); i++)
		{
			auto vh = res_vhs0[i];
			int pre_mark=
			teeth_seg_mark_[vh.idx()] = teeth_seg_mark_id_;
		}

		teeth_seg_mark_id_++;

		do {
			color1 = GetRandColor();
		} while (rand_color_set_.find(color1) != rand_color_set_.end());
		rand_color_set_.insert(color1);
		teeth_seg_color_[teeth_seg_mark_id_] = color1;
		for (int i = 0; i < res_vhs1.size(); i++)
		{
			auto vh = res_vhs1[i];
			teeth_seg_mark_[vh.idx()] = teeth_seg_mark_id_;
		}
		teeth_seg_mark_id_++;
		
		for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); viter++)
		{
			int mark_id = teeth_seg_mark_[viter->idx()];
			OpenMesh::Vec3d color(1,1,1);
			if(mark_id!=-1)
				color = teeth_seg_color_[mark_id];
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
	case Qt::Key_F:
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
	}
	case Qt::Key_D:
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
		std::vector<std::pair<COpenMeshT::VertexHandle, double>>cons;
		for (int i = 0; i < picked_vhs_fore_.size(); i++)
		{
			cons.push_back(std::pair<COpenMeshT::VertexHandle, double>(picked_vhs_fore_[i], 0.1));
		}
		std::cerr << "fore size " << picked_vhs_fore_.size()<<" "<< picked_vhs_back_.size() << std::endl;
		for (int i = 0; i < picked_vhs_back_.size(); i++)
		{
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
	}
}