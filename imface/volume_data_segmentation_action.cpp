#include"volume_data_segmentation_action.h"
#include"ui_context.h"
#include "qfiledialog.h"

#include"cmodelviewer.h"
#include"camera.h"
#include<qimage.h>
#include"../DataColle/volume_data_object.h"
#include"../DataColle/data_pool.h"
#include"../DataColle/aux_geo_utils.h"
#include"../DataColle/data_io.h"
void CVolumeDataSegmentationAction::Init()
{

}
void CVolumeDataSegmentationAction::MousePressEvent(QMouseEvent *e)
{
	
}
void CVolumeDataSegmentationAction::MouseMoveEvent(QMouseEvent *e)
{

}
void CVolumeDataSegmentationAction::MouseReleaseEvent(QMouseEvent *e)
{

}
void CVolumeDataSegmentationAction::KeyPressEvent(QKeyEvent *e)
{
	switch (e->key())
	{

	//	case Qt::Key_L:
	//	{
	///*		CMeshObject *mesh_obj = new CMeshObject();
	//	
	//		CAuxGeoUtils::GetVolumeRenderAuxCube(1, 1, 1, mesh_obj->GetMesh());
	//		DataPool::AddMeshObject(mesh_obj);
	//		mesh_obj->SetMeshColor(OpenMesh::Vec3d(1, 0, 0));
	//		mesh_obj->SetChanged();*/
	//		std::cerr << "load volume data" << std::endl;
	//		CVolumeDataObject *volume_data_obj = new CVolumeDataObject();
	//		DataPool::AddVolumeDataObject(volume_data_obj);
	//		break;
	//	}
		case Qt::Key_L:
		{
			QString path = QFileDialog::getExistingDirectory(NULL, "load volume data");

			if (path.length() == 0)
			{
				std::cerr << "unable to load volume data\n" << std::endl;
				return;
			}
			CVolumeDataObject *volume_data_obj = new CVolumeDataObject();
			CDataIO::ReadVolumeDataObjFromDICOMSeries(path.toStdString(), *volume_data_obj);
			DataPool::AddVolumeDataObject(volume_data_obj);
			volume_data_obj->SetChanged();
			break;
		}
		case Qt::Key_Q:
		{
			std::cerr << "switch to common action" << std::endl;
			manager_->SetCurrentActionType(CActionType::Common);
			break;
		}
	}
}
void CVolumeDataSegmentationAction::KeyReleaseEvent(QKeyEvent *e)
{

}

CVolumeDataSegmentationAction::CVolumeDataSegmentationAction()
{
	type_ = VolumeDataSegmentation;
}