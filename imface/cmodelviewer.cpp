// Copyright 2016_9_25 by ChenNenglun
#include "cmodelviewer.h"
#include"ui_context.h"
#include <cmath>
#include <float.h>
#include"../DataColle/data_io.h"
#include <QKeyEvent>
#include"../AlgColle/geo_base_alg.h"
#include"../AlgColle/geo_alg.h"
#include"../DataColle/aux_geo_utils.h"
#include"../AlgColle/morph_skel_dental_mesh_seg.h"
using namespace std;


void CModelViewer::draw()
{
	glClearColor(background_color_(0), background_color_(1), background_color_(2),1);
	if(scene_!=NULL)
	scene_->Render( camera_);

}
void CModelViewer::SetBackgroundColor(Eigen::Vector3d background_color)
{
	background_color_ = background_color;
}
void CModelViewer::initializeGL()
{


	QGLViewer::initializeGL();

	initializeOpenGLFunctions();



}
CModelViewer::CModelViewer(QWidget *parent)
{
	background_color_ = Eigen::Vector3d(1, 1, 1);
	scene_ = NULL;
	
	
}
void CModelViewer::InitCamera()
{
	
	camera_.setAspectRatio(1);
	camera_.setPosition(qglviewer::Vec(0, 0, 1));
	camera_.lookAt(qglviewer::Vec(0, 0, 0));

}
void CModelViewer::keyPressEvent(QKeyEvent *e)
{
	
	switch (e->key())
	{
	case Qt::Key_M:
	{
		int mid = CUIContext::GetSelectedMeshObjectId();
		auto p_mesh_object = DataPool::GetMeshObject(mid);
		CGeoBaseAlg::RemoveNonManifold(p_mesh_object->GetMesh());
		p_mesh_object->SetChanged();
		break;
	}
	case Qt::Key_L:
	{
		CMeshObject *meshobj = new CMeshObject();
		if (!CDataIO::ReadMesh("Resources\\dentalMesh\\UpperJaw3.stl", *meshobj))
		{
			std::cerr << "unable to load mesh\n" << std::endl;
		}
		
		COpenMeshT &mesh = meshobj->GetMesh();
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
				std::cerr << "a" << std::endl;
				max_vnum = vnum;
				mi = i;
			}
		}
		std::cerr << "mi " << mi << std::endl;
		mesh = *resmeshes[mi];
		std::cerr << mesh.n_vertices() << std::endl;
		/*for (int i = 0; i < resmeshes.size(); i++)
		{
			delete resmeshes[i];
		}*/

		CGeoAlg::FillHoles(mesh);
		//CGeoBaseAlg::RemoveNonManifold(mesh);
		//CGeoAlg::FillHoles(mesh);
		for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); viter++)
		{
			mesh.set_color(viter, OpenMesh::Vec3d(0.8, 0.8, 0.8));
		}
		meshobj->SetChanged();
		int mid = DataPool::AddMeshObject(meshobj);
		CUIContext::SetSelectedMeshObjectId(mid);
		this->update();
		break;
	}
	case Qt::Key_D:
	{
		CUIContext::msdm_seg_->TestDilate();
		this->update();
		break;
	}
	case Qt::Key_E:
	{
		CUIContext::msdm_seg_->TestErode();
		this->update();
		break;
	}
	case Qt::Key_R:
	{
		CUIContext::msdm_seg_->TestRemoveSmallFeatureRegions();
		this->update();
		break;
	}
	case Qt::Key_P:
	{
		std::cerr << "please input the prarm you want to adjust" << std::endl;
		std::cerr << "0: curvature" << std::endl;
		std::cerr << "1: small region num" << std::endl;
		int id;
		std::cin>>id;
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
			this->update();

		}
		else if (CUIContext::cm_current_adjust_param_ == "SmallRegion")
		{
			CUIContext::msdm_seg_->TestGetSmallRegionThreshold() -= 1;
			
			std::cerr << "current small region threshold " << CUIContext::msdm_seg_->TestGetSmallRegionThreshold() << std::endl;
			this->update();
		}
		break;
	}
	case Qt::Key_G:
	{
		CUIContext::msdm_seg_->TestTagGingiva();
		this->update();
		break;
	}
	case 16777235:
	{
		if (CUIContext::cm_current_adjust_param_ == "Curvature")
		{
			CUIContext::msdm_seg_->TestGetCurvatureThreshold() += 1;
			CUIContext::msdm_seg_->TestCurvature();
			std::cerr << "current curvature threshold " << CUIContext::msdm_seg_->TestGetCurvatureThreshold() << std::endl;
			this->update();

		}
		else if (CUIContext::cm_current_adjust_param_ == "SmallRegion")
		{
			CUIContext::msdm_seg_->TestGetSmallRegionThreshold() += 1;
		
			std::cerr << "current small region threshold " << CUIContext::msdm_seg_->TestGetSmallRegionThreshold() << std::endl;
			this->update();
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
		this->update();
		break;
	}
	case Qt::Key_S:
	{
		CUIContext::msdm_seg_->TestSkeletonize();
		this->update();
		break;
	}
	case Qt::Key_C://compute mean curvature
	{
		std::cout << "compute mean curvature" << std::endl;
		int mid = CUIContext::GetSelectedMeshObjectId();
		auto p_mesh_object = DataPool::GetMeshObject(mid);
		if(CUIContext::msdm_seg_==NULL)
		CUIContext::msdm_seg_ = new CMorphSkelDentalMeshSeg(*p_mesh_object);
		
		CUIContext::msdm_seg_->ComputeSegmentation(true);
	
		this->update();
		break;
	}
	case Qt::Key_O://
	{
		std::cerr << "save" << std::endl;
		int mid = CUIContext::GetSelectedMeshObjectId();
		auto p_mesh_object = DataPool::GetMeshObject(mid);
		CDataIO::WriteMesh("Resources\\LowerJawCut.obj", *p_mesh_object);
	}
	}
}
void CModelViewer::SetScene(CScene* scene)
{
	scene_ = scene;
}
CModelViewer::~CModelViewer()
{

}
void CModelViewer::init()
{
	InitCamera();
	setCamera(&camera_);
	setMouseBinding(Qt::ControlModifier, Qt::LeftButton, QGLViewer::CAMERA, QGLViewer::ROTATE);
	setMouseBinding(Qt::ControlModifier, Qt::RightButton, QGLViewer::CAMERA, QGLViewer::TRANSLATE);
	setMouseBinding(Qt::ControlModifier, Qt::MidButton, QGLViewer::CAMERA, QGLViewer::ZOOM);
	setWheelBinding(Qt::ControlModifier, QGLViewer::CAMERA, QGLViewer::ZOOM);
}
