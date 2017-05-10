// Copyright 2016_9 by ChenNenglun
#ifndef IMFACE_WINDOW_H
#define IMFACE_WINDOW_H

#include <QWidget>
#include <QtWidgets/QMainWindow>
#include "ui_imface.h"
#include"cmodelviewer.h"
#include<qpushbutton.h>
#include <single_Teeth_Projection_Action.h>
#include <panoramic_image_registration.h>
#include "../TeethRootRecoAlg/teethProjection.h"

class QImage;
class NDTRegistration;

class CMainWindow : public QMainWindow
{
	Q_OBJECT

public:
	CMainWindow(QWidget *parent = 0);
	~CMainWindow();

	CModelViewer* GetModelViewer();
	void SetConnections();
	void SetComponents();
	void loadCrownData();
	int test();

	QLineEdit *id_7_;
	QLineEdit *id_6_;
	QLineEdit *id_5_;
	QLineEdit *id_4_;
	QLineEdit *id_3_;
	QLineEdit *id_2_;
	QLineEdit *id_1_;

	double image_crown_minW_, image_crown_maxW_;
	double projectiom_crown_minW_, projection_crown_maxW_;
	Eigen::Vector2d image_crown_center_;
	Eigen::Vector2d projection_crown_center_;
	std::map<int, std::vector<OpenMesh::Vec2d>> projection_tooth_crown_;
	std::vector<Eigen::Vector2d> projection_teeth_crown_;
	std::map<int, Eigen::Vector2d> teeth_root_tip_;
	std::map<int, std::map<int, std::map<int, double>>> map_toothId_map_meshId_tipLen_;
	std::vector<std::map<int, int>> map_3dvertice_2dpoints_;
protected:
	Ui::ImfaceWindow ui;
	QPushButton *push_template_kind;
	QPushButton *push_load_image_;
protected:
	 void UpdateRequest();

private slots:
	void pushTemplateKind();
	void loadXImage();
private:
	CSingleTeethProjectionAction * single_teeth_projection_action_;
	PanoramicImageRegistration *panoramic_image_registration_;
	teethProjection * teeth_projection_;
	NDTRegistration *ndt_registration_;
	QImage *panoramic_image_;
};

#endif // MY3DMM_H
