#pragma once
#ifndef PANORAMIC_IMAGE_REGISTRATION_H_
#define PANORAMIC_IMAGE_REGISTRATION_H_
#include <QWidget>
#include "ui_panoramic_image_registration.h"
#include <vector>
class QPushButton;
class QLabel;
class QImage;
class QProgressBar;
class QLineEdit;
class NDTRegistration;
using std::vector;

class PanoramicImageRegistration : public QWidget
{
	Q_OBJECT

public:
	PanoramicImageRegistration(QWidget *parent = Q_NULLPTR);
	~PanoramicImageRegistration();
	void SetComponents();
	void SetConnections();
	void ShowBoundaryPointsInImage();	

public slots:
	void OnLoadPrimiticeImage(void);
	void OnLoadLabeledImage(void);
	void OnLoadOBJPoints();
	void OnTransformAndRotateBoundaries(void);
	void OnRegisterWithNDTAlgorithm(void);

private:
	Ui::PanoramicImageRegistration ui;

private:
	QPushButton* button_load_primitive_image_;
	QPushButton* button_labeled_image_;
	QPushButton* button_load_obj_points_;
	QPushButton* button_boundaries_up_;
	QPushButton* button_boundaries_down_;
	QPushButton* button_boundaries_left_;
	QPushButton* button_boundaries_right_;
	QLineEdit* lineedit_cell_num_w_;
	QLineEdit* lineedit_cell_num_h_;
	QLineEdit* lineedit_cell_num_addition_w_;
	QLineEdit* lineedit_cell_num_addition_h_;
	QLineEdit* lineedit_cell_num_count_;
	QLineEdit* lineedit_max_iterator_times_;
	QPushButton* button_registering_with_NDT_algorithm_;
	QLabel* label_image_panel1_;
	QLabel* label_image_panel2_;

private:
	NDTRegistration* ndt_registration_;
	QImage* primitive_image_;
	QImage* labeled_image_;
};

#endif