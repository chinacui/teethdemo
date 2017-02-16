// Copyright 2016_9 by ChenNenglun
#ifndef IMFACE_WINDOW_H
#define IMFACE_WINDOW_H

#include <QtWidgets/QMainWindow>
#include "ui_imface.h"
#include"cmodelviewer.h"
class CMainWindow : public QMainWindow
{
	Q_OBJECT

public:
	CMainWindow(QWidget *parent = 0);
	~CMainWindow();
	CModelViewer* GetModelViewer() { return ui.model_viewer; }

protected:
	Ui::ImfaceWindow ui;
protected:
	
	 void UpdateRequest();
private slots:
	//void OnClickButtonLoadData();
//void AdjustBaseCuttingPlane(int v);
//void AdjustSmallRegionThreshold(int);
};

#endif // MY3DMM_H
