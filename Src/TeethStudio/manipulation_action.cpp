#include"manipulation_action.h"
#include"../AlgColle/geo_alg.h"
#include"../AlgColle/geo_base_alg.h"
#include"cmodelviewer.h"
#include"action_manager.h"
#include"ui_context.h"
#include"../DataColle/mesh_object.h"
#include"ui_context.h"
#include"imface_window.h"
void CManipulationAction::MousePressEvent(QMouseEvent *e)
{
	
	if (e->modifiers() == Qt::Modifier::CTRL)
	{
		std::cerr << "picking" << std::endl;
		auto camera = viewer_->GetCamera();
		OpenMesh::Vec3d orig, dir;

		camera.ConvertClickToLine(e->pos(), orig, dir);
		auto data_pool = DataPool::GetMeshObjectPool();

		int meshid = CGeoAlg::PickMesh(orig, dir, DataPool::GetMeshObjectPool(), true);
		CMeshObject * mesh_obj = DataPool::GetMeshObject(meshid);
		auto datapool = DataPool::GetMeshObjectPool();
		for (auto iter = datapool.begin(); iter != datapool.end(); iter++)
		{
			iter->second.get()->IsShinning() = false;
		}
		sel_mesh_ids_.clear();
		if (mesh_obj != NULL)
		{
			mesh_obj->IsShinning() = true;
			sel_mesh_center_=CGeoBaseAlg::ComputeMeshCenter(mesh_obj->GetMesh());
			sel_mesh_center_=mesh_obj->TransformPointByLocalMatrix(sel_mesh_center_);
			sel_mesh_ids_.push_back(meshid);
		}
	}
	else if (e->modifiers() == Qt::Modifier::ALT)
	{
		if (e->button() == Qt::MouseButton::LeftButton)
		{
			is_moving_mesh_ = true;
			OpenMesh::Vec3d orig, dir;
			viewer_->GetCamera().ConvertClickToLine(e->pos(), orig, dir);
			pre_move_pos_ = OpenMesh::Vec3d(orig[0], orig[1], orig[2]);
		}
		
	}
}
void CManipulationAction::MouseMoveEvent(QMouseEvent *e)
{
	if (e->modifiers() == Qt::Modifier::ALT)
	{
		if (is_moving_mesh_)
		{
			OpenMesh::Vec3d orig, dir;
			viewer_->GetCamera().ConvertClickToLine(e->pos(), orig, dir);
			OpenMesh::Vec3d cpos = OpenMesh::Vec3d(orig[0], orig[1], orig[2]);
			OpenMesh::Vec3d trans = cpos - pre_move_pos_;
			for (int i = 0; i < sel_mesh_ids_.size(); i++)
			{
				CMeshObject *mesh_obj = DataPool::GetMeshObject(sel_mesh_ids_[i]);
				mesh_obj->Transform(trans);
				mesh_obj->SetChanged();

			}
			pre_move_pos_ = cpos;
		}
	}
	else
	{
		is_moving_mesh_ = false;
	}
}
void CManipulationAction::MouseReleaseEvent(QMouseEvent *e)
{
	if (is_moving_mesh_)
	{
		for (int i = 0; i < sel_mesh_ids_.size(); i++)
		{
			CMeshObject *mesh_obj = DataPool::GetMeshObject(sel_mesh_ids_[i]);
			mesh_obj->ApplyTransform();
			mesh_obj->SetChanged();
			int count_vn = 0;
			for (int i = 0; i < sel_mesh_ids_.size(); i++)
			{
				CMeshObject *mesh_obj = DataPool::GetMeshObject(sel_mesh_ids_[i]);
				OpenMesh::Vec3d c = CGeoBaseAlg::ComputeMeshCenter(mesh_obj->GetMesh());
				c = mesh_obj->TransformPointByLocalMatrix(c);
				c = c*mesh_obj->GetMesh().n_vertices();
				count_vn += mesh_obj->GetMesh().n_vertices();
				sel_mesh_center_ += c;
			


			}
			sel_mesh_center_ /= count_vn;

		}
		is_moving_mesh_ = false;
	}
	
}
void CManipulationAction::WheelEvent(QWheelEvent *e)
{
	if (e->modifiers() == Qt::Modifier::ALT)
	{
		double scale = 1;
		if (e->delta() > 0)
		{
			scale = 1.1;
		}
		else
		{
			scale = 0.9;

		}
		for (int i = 0; i < sel_mesh_ids_.size(); i++)
		{
			CMeshObject* mesh_obj = DataPool::GetMeshObject(sel_mesh_ids_[i]);
			if (mesh_obj != NULL)
			{
				COpenMeshT &mesh = mesh_obj->GetMesh();
				for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); viter++)
				{
					OpenMesh::Vec3d p = mesh.point(viter);
					p=sel_mesh_center_+(p - sel_mesh_center_)*scale;
					mesh.set_point(viter, p);
				}
				mesh_obj->ApplyTransform();
				mesh_obj->SetChanged();
			}
		}
	}
}
void CManipulationAction::KeyPressEvent(QKeyEvent *e)
{
	if (e->modifiers() == Qt::Modifier::CTRL)
	{
		if (e->key() == Qt::Key_A)
		{
			sel_mesh_ids_.clear();
			auto datapool = DataPool::GetMeshObjectPool();
			sel_mesh_center_ = OpenMesh::Vec3d(0, 0, 0);
			int count_vn = 0;
			for (auto iter = datapool.begin(); iter != datapool.end(); iter++)
			{
				if (iter->second.get()->IsVisiable()&& iter->second.get()->IsPickAble())
				{
					iter->second.get()->IsShinning() = true;
					sel_mesh_ids_.push_back(iter->second.get()->GetId());
					COpenMeshT&mesh = iter->second.get()->GetMesh();
					OpenMesh::Vec3d c=CGeoBaseAlg::ComputeMeshCenter(mesh);
					c=iter->second.get()->TransformPointByLocalMatrix(c);
					c = c*mesh.n_vertices();
					count_vn += mesh.n_vertices();
					sel_mesh_center_ += c;
					iter->second.get()->IsShinning() = true;
				
				}
			}
			sel_mesh_center_ /= count_vn;
		}
		
	}
	else
	{
		if (e->key() == Qt::Key_Space)
		{
			CUIContext::GetMainWindow()->GetModelViewer()->InitCamera();

		}
	}

}
void CManipulationAction::KeyReleaseEvent(QKeyEvent *e)
{

}