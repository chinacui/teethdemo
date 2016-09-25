#include "imface_window.h"
#include <QtWidgets/QApplication>
#include"ui_context.h"
int main(int argc, char *argv[])
{
	QApplication a(argc, argv);

	
	CUIContext::Init();
	
	return a.exec();
}
