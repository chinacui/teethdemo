// Copyright 2016_9 by ChenNenglun
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
	CModelViewer* GetModelViewer() { return ui.model_viewer; }

protected:
	Ui::ImfaceWindow ui;
protected:
	
	 void UpdateRequest();
private slots:
	//void OnClickButtonLoadData();
void AdjustBaseCuttingPlane(int v);
void AdjustSmallRegionThreshold(int);
};

#endif // MY3DMM_H
