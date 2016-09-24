#include "imface_window.h"
#include <QtWidgets/QApplication>

int main(int argc, char *argv[])
{
	QApplication a(argc, argv);

	CImfaceWindow w;

	//CImface m;
	//w.setCentralWidget(&m);
	w.show();

	return a.exec();
}
