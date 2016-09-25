#include "imface_window.h"
#include"ui_context.h"
#include"../DataColle/data_io.h"
CImfaceWindow::CImfaceWindow(QWidget *parent)
	: QMainWindow(parent)
{
	ui.setupUi(this);
	connect(ui.button_load_model, SIGNAL(clicked()), this, SLOT(OnClick_Button_Load_Data()));

	ui.model_viewer->SetScene(CUIContext::GetScene());
	
}
void CImfaceWindow::OnClick_Button_Load_Data()
{
	std::cerr << "load data clicked" << std::endl;
	CMeshObject *meshobj = new CMeshObject();
	CDataIO::ReadPly("groundtruth_model.ply", *meshobj);
	DataPool::AddMeshObject(meshobj);
}

CImfaceWindow::~CImfaceWindow()
{
	
}
