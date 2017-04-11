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
#include<igl\writeOBJ.h>
#include "panoramic_image_registration.h"

void CHotKeyAction::KeyPressEvent(QKeyEvent *e)
{

	switch (e->key())
	{



	case Qt::Key_V:
	{
		manager_->SetCurrentActionType(VolumeDataSegmentation);
		std::cerr << "switch to volume segmentation action" << std::endl;
		break;

	}
	case Qt::Key_T://test panoramic simulation
	{
		manager_->SetCurrentActionType(CPanoramicSimulationTest);
		std::cerr << "switch to volume CPanoramicSimulationTest action" << std::endl;
		break;
	}
	case Qt::Key_H:
	{
		manager_->SetCurrentActionType(HarmonicFieldSegmentation);
		std::cerr << "switch to harmonic action" << std::endl;
		break;

	}
	case Qt::Key_P:
	{
		manager_->SetCurrentActionType(CTeethReconstructionTest);
		std::cerr << "switch to panoramia registration test action" << std::endl;
		break;
	}
	case Qt::Key_K:
	{
		manager_->SetCurrentActionType(CTeethReconstruction);
		std::cerr << "switch to panoramia registration action" << std::endl;
		break;
	}
	case Qt::Key_S:
	{
		PanoramicImageRegistration* panoramic_image_registration =
			new PanoramicImageRegistration();
		panoramic_image_registration->show();
	}

	case Qt::Key_1:
	{
		manager_->SetCurrentActionType(CSingleTeethProjection);
		std::cerr << "switch to single teeth projection!" << std::endl;

	}
	//case Qt::Key_M:
	//{
	//	manager_->SetCurrentActionType(EditFeatureEdge);
	//	std::cerr << "switch to edit feature edge action" << std::endl;
	//	break;
	//}

	//case Qt::Key_C:
	//{
	//	DataPool::DeleteAllCurveObjects();
	//
	//	break;
	//}

	/*case Qt::Key_R:
	{
		std::cerr << "test tri" << std::endl;
		std::vector<Eigen::Vector3d>bounds,inner;
		bounds.push_back(Eigen::Vector3d(0, 0,0));
		bounds.push_back(Eigen::Vector3d( 0.00898122,-1.0842e-19,0));
		bounds.push_back(Eigen::Vector3d(0.00776125, 0.00372481, 0));

		bounds.push_back(Eigen::Vector3d(0.00833388,0.00197644, 0));
		bounds.push_back(Eigen::Vector3d(0.00490583 ,0.000803329, 0));
		Eigen::MatrixXd res_v;
		Eigen::MatrixXi res_f;
		CGeoAlg::DelaunayTriangulation3d(bounds, res_v, res_f);
	
		igl::writeOBJ("test_tri.obj", res_v, res_f);
		break;
	}*/
	//case Qt::Key_O://
	//{
	//	std::cerr << "save" << std::endl;
	//	int mid = CUIContext::GetSelectedMeshObjectId();
	//	auto p_mesh_object = DataPool::GetMeshObject(mid);
	//	CDataIO::WriteMesh("Resources\\out.obj", *p_mesh_object);
	//}
	}
}
void CHotKeyAction::KeyReleaseEvent(QKeyEvent *e)
{

}