#include "hot_key_action.h"
#include"ui_context.h"
#include <cmath>
#include <float.h>
#include"../DataColle/data_io.h"
#include"../AlgColle/geo_base_alg.h"
#include"../DataColle/data_pool.h"
#include"../AlgColle/geo_alg.h"
#include"../DataColle/aux_geo_utils.h"
#include"../AlgColle/morph_skel_dental_mesh_seg.h"
#include "qfiledialog.h"
#include"cmodelviewer.h"
#include"action_manager.h"

void CHotKeyAction::KeyPressEvent(QKeyEvent *e)
{

	switch (e->key())
	{

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
		CGeoBaseAlg::NormalizeMeshSize(mesh);
		//CGeoAlg::FillHoles(mesh);
		CGeoBaseAlg::RemoveNonManifold(mesh);
		CGeoAlg::SelfIntersectionRemoval(mesh);
		CGeoAlg::FillHoles(mesh);
		for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); viter++)
		{
			mesh.set_color(viter, OpenMesh::Vec3d(0.8, 0.8, 0.8));
		}
		meshobj->SetChanged();
		int id = DataPool::AddMeshObject(meshobj);
		CUIContext::SetSelectedMeshObjectId(id);
		if (CUIContext::msdm_seg_ == NULL)
			CUIContext::msdm_seg_ = new CMorphSkelDentalMeshSeg(*meshobj);
		break;
	}
	case Qt::Key_D:
	{
		CUIContext::msdm_seg_->TestDilate();
	
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
	case Qt::Key_M:
	{
		manager_->SetCurrentActionType(EditFeatureEdge);
		std::cerr << "switch to edit feature edge action" << std::endl;
		break;
	}
	case Qt::Key_H:
	{
		manager_->SetCurrentActionType(HarmonicFieldSegmentation);
		std::cerr << "switch to harmonic action" << std::endl;
		break;
		break;
	}
	case Qt::Key_T:
	{
		int mid = CUIContext::GetSelectedMeshObjectId();
		auto p_mesh_object = DataPool::GetMeshObject(mid);
		std::vector<COpenMeshT*>meshes;
		std::vector<int>tags;
		for (int i = 0; i < CUIContext::msdm_seg_->tags_.size(); i++)
		{
			if (CUIContext::msdm_seg_->tags_[i] == CMorphSkelDentalMeshSeg::CTag::Teeth || CUIContext::msdm_seg_->tags_[i] == CMorphSkelDentalMeshSeg::CTag::Feature)
				tags.push_back(1);
			else
				tags.push_back(0);
			
		}
		CGeoAlg::SeparateMeshByVertexTag(p_mesh_object->GetMesh(), tags, meshes, std::vector<std::map<COpenMeshT::VertexHandle, COpenMeshT::VertexHandle>>());
		CMeshObject m;
		m.GetMesh() = *meshes[1];
		CDataIO::WriteMesh("res0.obj", m);
		break;
	}
	case Qt::Key_O://
	{
		std::cerr << "save" << std::endl;
		int mid = CUIContext::GetSelectedMeshObjectId();
		auto p_mesh_object = DataPool::GetMeshObject(mid);
		CDataIO::WriteMesh("Resources\\testmodel.obj", *p_mesh_object);
	}
	}
}
void CHotKeyAction::KeyReleaseEvent(QKeyEvent *e)
{

}