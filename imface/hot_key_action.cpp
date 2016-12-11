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
#include"../AlgColle/dental_base_alg.h"
#include"imface_window.h"
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
		CDentalBaseAlg::PCABasedOrientationCorrection(meshobj->GetMesh());
		CGeoBaseAlg::NormalizeMeshSize(mesh);
		meshobj->SetChanged();
		int id = DataPool::AddMeshObject(meshobj);
		CUIContext::SetSelectedMeshObjectId(id);
		if (CUIContext::msdm_seg_ == NULL)
			CUIContext::msdm_seg_ = new CMorphSkelDentalMeshSeg(*meshobj);

		
		break;
	}


	case Qt::Key_Space:
	{
		CUIContext::GetMainWindow()->GetModelViewer()->InitCamera();
		break;
	}
	case Qt::Key_M:
	{
		manager_->SetCurrentActionType(EditFeatureEdge);
		std::cerr << "switch to edit feature edge action" << std::endl;
		break;
	}
	case Qt::Key_J:
	{
		manager_->SetCurrentActionType(VolumeDataSegmentation);
		std::cerr << "switch to volume segmentation action" << std::endl;
		break;
	
	}
	case Qt::Key_H:
	{
		manager_->SetCurrentActionType(HarmonicFieldSegmentation);
		std::cerr << "switch to harmonic action" << std::endl;
		break;
	
	}
	case Qt::Key_C:
	{
		DataPool::DeleteAllCurveObjects();
	
		break;
	}
	case Qt::Key_P:
	{
		manager_->SetCurrentActionType(CPanoramaMeshRegistration);
	std::cerr << "switch to panoramia registration action" << std::endl;
		break;
	}
	case Qt::Key_O://
	{
		std::cerr << "save" << std::endl;
		int mid = CUIContext::GetSelectedMeshObjectId();
		auto p_mesh_object = DataPool::GetMeshObject(mid);
		CDataIO::WriteMesh("Resources\\out.obj", *p_mesh_object);
	}
	}
}
void CHotKeyAction::KeyReleaseEvent(QKeyEvent *e)
{

}