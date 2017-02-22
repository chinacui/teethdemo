#pragma once
#ifndef PANORAMIC_IMAGE_REGISTRATION_H_
#define PANORAMIC_IMAGE_REGISTRATION_H_
#include <QWidget>
#include "ui_panoramic_image_registration.h"
#include <vector>
class QPushButton;
class QLabel;
class QPointF;
class QImage;
class QProgressBar;
class QLineEdit;
class NDTRegisterThread;
using std::vector;

class PanoramicImageRegistration : public QWidget
{
	Q_OBJECT

public:
	~PanoramicImageRegistration();
	static PanoramicImageRegistration* GetInstance();
	static void DeleteInstance();
	void SetComponents();
	void SetConnections();
	void GetBoundaryPointsFromLabels();
	void ComputeGradient();
	void ComputeCellMeanPoints();
	void ComputeCellConvarianceMatrix();
	void ComputeOneNDTRegister();

private:
	PanoramicImageRegistration(QWidget *parent = Q_NULLPTR);

public slots:
	void OnLoadPrimiticeImage(void);
	void OnLoadLabeledImage(void);
	void OnComputeCrownBoundaries(void);
	void OnTransformAndRotateBoundaries(void);
	void OnRegisterWithNDTAlgorithm(void);
	void OnUpdateImageOnce(int iterator_time, int max_times, 
		double cell_size_w, double cell_size_h);
	void OnLoadOBJPoints();

private:
	Ui::PanoramicImageRegistration ui;
	static PanoramicImageRegistration* instance_;

private:
	QPushButton* button_load_primitive_image_;
	QPushButton* button_labeled_image_;
	QPushButton* button_compute_crown_boundaries_;
	QPushButton* button_boundaries_up_;
	QPushButton* button_boundaries_down_;
	QPushButton* button_boundaries_left_;
	QPushButton* button_boundaries_right_;
	QPushButton* button_registering_with_NDT_algorithm_;
	QLabel* label_image_panel1_;
	QLabel* label_image_panel2_;
	QProgressBar* progress_bar_;
	QLineEdit* lineedit_cell_num_w_;
	QLineEdit* lineedit_cell_num_h_;
	QLineEdit* lineedit_cell_num_addition_w_;
	QLineEdit* lineedit_cell_num_addition_h_;
	QLineEdit* lineedit_cell_num_count_;
	QLineEdit* lineedit_max_iterator_times_;
	QLabel* label_info_;
	QPushButton* button_load_obj_points_;

private:
	QImage* primitive_image_;
	QImage* labeled_image_;
	QImage* crown_boundary_image_;
	vector<QPointF> crown_boundary_points_;
	vector<QPointF> crown_boundary_points_pre_;
	vector<double> image_gradient_;
	vector<vector<bool>> boundary_label_;
	vector<double> cell_weights_;
	vector<QPointF> cell_mean_points_;
	vector<double> cell_convariance_matrix_;
	vector<double> cell_inverse_convariance_matrix_;
	int image_width_;
	int image_height_;
	int cell_num_w_;
	int cell_num_h_;
	NDTRegisterThread* ndt_register_thread_;
};

#endif