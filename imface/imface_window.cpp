// Copyright 2016_9 by ChenNenglun
#include "imface_window.h"
#include"ui_context.h"
#include"../DataColle/data_io.h"
#include <QKeyEvent>
#include"../AlgColle/geo_base_alg.h"
CImfaceWindow::CImfaceWindow(QWidget *parent)
	: QMainWindow(parent)
{
	ui.setupUi(this);
	/*connect(ui.button_load_model, SIGNAL(clicked()), this, SLOT(OnClickButtonLoadData()));*/

	ui.model_viewer->SetScene(CUIContext::GetScene());
	//ui.base_cut_slilder->hide();
	//ui.region_threshold_slilder->hide();
	this->connect(ui.base_cut_slilder, SIGNAL(valueChanged(int)), this, SLOT(AdjustBaseCuttingPlane(int)));
	this->connect(ui.region_threshold_slilder, SIGNAL(valueChanged(int)), this, SLOT(AdjustSmallRegionThreshold(int)));
	
}
void CImfaceWindow::AdjustSmallRegionThreshold(int v)
{
	double percent = v / 100.0*0.4;
	CUIContext::msdm_seg_->AdjustSmallRegionThreshold(percent);
	CUIContext::msdm_seg_->RemoveSmallIsolateTeethRegion();
}
void CImfaceWindow::AdjustBaseCuttingPlane(int v)
{
	
	static int prev = 0;
	double l = (v - prev) / 100.0;
	prev = v;
	CUIContext::msdm_seg_->AdjustBaseCuttingPlane(l/5.0);
	
	this->UpdateRequest();

}
//void CImfaceWindow::OnClickButtonLoadData()
//{
//	std::cerr << "load data clicked" << std::endl;
//	CMeshObject *meshobj = new CMeshObject();
//	CDataIO::ReadPly("groundtruth_model.ply", *meshobj);
//	DataPool::AddMeshObject(meshobj);
//
//	CMeshObject *meshobj2 = new CMeshObject();
//	CDataIO::ReadPly("model_result_rgb.ply", *meshobj2);
//	DataPool::AddMeshObject(meshobj2);
//}

void CImfaceWindow::UpdateRequest()
{
	this->ui.model_viewer->update();
	this->update();	
}
CImfaceWindow::~CImfaceWindow()
{
	
}
