#ifndef IMFACE_WINDOW_H
#define IMFACE_WINDOW_H

#include <QtWidgets/QMainWindow>
#include "ui_imface.h"
#include"cmodelviewer.h"
class CImfaceWindow : public QMainWindow
{
	Q_OBJECT

public:
	CImfaceWindow(QWidget *parent = 0);
	~CImfaceWindow();


protected:
	Ui::ImfaceWindow ui;
protected:
	 
private slots:
	void OnClick_Button_Load_Data();
};

#endif // MY3DMM_H
