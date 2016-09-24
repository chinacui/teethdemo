#include "imface_window.h"

CImfaceWindow::CImfaceWindow(QWidget *parent)
	: QMainWindow(parent)
{
	ui.setupUi(this);
	connect(ui.button_load_model, SIGNAL(clicked()), this, SLOT(OnClick_Button_Load_Data()));
	
}
void CImfaceWindow::OnClick_Button_Load_Data()
{
	std::cerr << "load data clicked" << std::endl;
}
CImfaceWindow::~CImfaceWindow()
{
	
}
