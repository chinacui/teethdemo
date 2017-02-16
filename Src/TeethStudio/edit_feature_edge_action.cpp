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
	case Qt::Key_A:
	{
		is_pick_ = true;
		break;
	}
	case Qt::Key_T:
	{
		int mid = CUIContext::GetSelectedMeshObjectId();
		auto p_mesh_object = DataPool::GetMeshObject(mid);
		std::map<int, COpenMeshT*>meshes;
		std::vector<int>tags;
		for (int i = 0; i < CUIContext::msdm_seg_->tags_.size(); i++)
		{
			if (CUIContext::msdm_seg_->tags_[i] == CMorphSkelDentalMeshSeg::CTag::Teeth || CUIContext::msdm_seg_->tags_[i] == CMorphSkelDentalMeshSeg::CTag::Feature)
				tags.push_back(1);
			else
				tags.push_back(0);

		}
		CGeoAlg::SeparateMeshByVertexTag(p_mesh_object->GetMesh(), tags, meshes, std::map<int, std::map<COpenMeshT::VertexHandle, COpenMeshT::VertexHandle>>());
		CMeshObject m;
		m.GetMesh() = *meshes[1];
		//CDataIO::WriteMesh("res0.obj", m);
		break;
	}
	case 16777237:
	{
		if (CUIContext::cm_current_adjust_param_ == "Curvature")
		{
			CUIContext::msdm_seg_->TestGetCurvatureThreshold() -= 1;
			CUIContext::msdm_seg_->TestCurvature();
			std::cerr << "current curvature threshold " << CUIContext::msdm_seg_->TestGetCurvatureThreshold() << std::endl;


		}
		else if (CUIContext::cm_current_adjust_param_ == "SmallRegion")
		{
			CUIContext::msdm_seg_->TestGetSmallRegionThreshold() -= 1;

			std::cerr << "current small region threshold " << CUIContext::msdm_seg_->TestGetSmallRegionThreshold() << std::endl;

		}
		break;
	}
	case Qt::Key_E:
	{
		CUIContext::msdm_seg_->TestErode();

		break;
	}
	case Qt::Key_R:
	{
		CUIContext::msdm_seg_->TestRemoveSmallFeatureRegions();

		break;
	}

	case Qt::Key_P:
	{
		std::cerr << "please input the prarm you want to adjust" << std::endl;
		std::cerr << "0: curvature" << std::endl;
		std::cerr << "1: small region num" << std::endl;
		int id;
		std::cin >> id;
		if (id == 0)
		{
			CUIContext::cm_current_adjust_param_ = "Curvature";
		}
		else if (id == 1)
		{
			CUIContext::cm_current_adjust_param_ = "SmallRegion";
		}
		break;
	}
	case Qt::Key_D:
	{
		CUIContext::msdm_seg_->TestDilate();

		break;
	}

	case Qt::Key_G:
	{
		CUIContext::msdm_seg_->TestTagGingiva();

		break;
	}
	case 16777235:
	{
		if (CUIContext::cm_current_adjust_param_ == "Curvature")
		{
			CUIContext::msdm_seg_->TestGetCurvatureThreshold() += 1;
			CUIContext::msdm_seg_->TestCurvature();
			std::cerr << "current curvature threshold " << CUIContext::msdm_seg_->TestGetCurvatureThreshold() << std::endl;


		}
		else if (CUIContext::cm_current_adjust_param_ == "SmallRegion")
		{
			CUIContext::msdm_seg_->TestGetSmallRegionThreshold() += 1;

			std::cerr << "current small region threshold " << CUIContext::msdm_seg_->TestGetSmallRegionThreshold() << std::endl;

		}
		break;
	}
	case Qt::Key_F:
	{
		int mid = CUIContext::GetSelectedMeshObjectId();
		auto p_mesh_object = DataPool::GetMeshObject(mid);
		if (CUIContext::msdm_seg_ == NULL)
			CUIContext::msdm_seg_ = new CMorphSkelDentalMeshSeg(*p_mesh_object);
		CUIContext::msdm_seg_->TestCurvature();

		break;
	}
	case Qt::Key_S:
	{
		CUIContext::msdm_seg_->TestSkeletonize();

		break;
	}
	case Qt::Key_C://compute mean curvature
	{
		std::cout << "compute mean curvature" << std::endl;
		int mid = CUIContext::GetSelectedMeshObjectId();
		auto p_mesh_object = DataPool::GetMeshObject(mid);
		if (CUIContext::msdm_seg_ == NULL)
			CUIContext::msdm_seg_ = new CMorphSkelDentalMeshSeg(*p_mesh_object);

		CUIContext::msdm_seg_->ComputeSegmentation(true);

		break;
	}

	}
}
void CEditFeatureEdgeAction::KeyReleaseEvent(QKeyEvent *e)
{
	switch (e->key())
	{


	case Qt::Key_A:
	{
		is_pick_ = false;
		break;
	}
	}
}