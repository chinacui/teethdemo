#include"edit_feature_edge_action.h"
#include"../AlgColle/geo_alg.h"
#include"ui_context.h"
#include"../DataColle/data_pool.h"
#include"cmodelviewer.h"
#include"camera.h"
#include"../AlgColle/geo_base_alg.h"
void CEditFeatureEdgeAction::MousePressEvent(QMouseEvent *e)
{
	if(e->button()==Qt::LeftButton)
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
				
				picked_vhs_.push_back(vh);

				if (picked_vhs_.size() == 2)
				{
					std::vector<COpenMeshT::VertexHandle>path;
					CGeoAlg::ComputeShortestPathAlongEdge(mesh, picked_vhs_[0], picked_vhs_[1], path);
					std::vector<CMorphSkelDentalMeshSeg::CTag>tags(path.size());
					for (int i = 0; i < path.size(); i++)
					{
						tags[i] = CMorphSkelDentalMeshSeg::CTag::Feature;
					}
					CUIContext::msdm_seg_->SetVertexTags(path, tags);
					CUIContext::msdm_seg_->TestRender();
					picked_vhs_.clear();
				}

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
				CUIContext::msdm_seg_->SwitchGingivaAndTeeth(vh);
				CUIContext::msdm_seg_->RemoveSmallIsolateTeethRegion();
			}
		}
	}
	

}
void CEditFeatureEdgeAction::MouseMoveEvent(QMouseEvent *e)
{

}
void CEditFeatureEdgeAction::MouseReleaseEvent(QMouseEvent *e)
{

}
void CEditFeatureEdgeAction::KeyPressEvent(QKeyEvent *e)
{
	switch (e->key())
	{

	case Qt::Key_Q:
	{
		manager_->SetCurrentActionType(CActionType::Common);
		break;
	}
	}
}
void CEditFeatureEdgeAction::KeyReleaseEvent(QKeyEvent *e)
{

}