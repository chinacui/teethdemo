#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include "QObject"
#include "QFileDialog"
#include "QProgressBar"
#include "QLineEdit"
#include "QLabel"
#include "QRgb"
#include <Eigen\Dense>
#include "panoramic_image_registration.h"
#include "ndt_registration.h"

PanoramicImageRegistration::PanoramicImageRegistration(QWidget *parent)
	: QWidget(parent)
{
	ui.setupUi(this);

	this->SetComponents();
	this->SetConnections();
}

PanoramicImageRegistration::~PanoramicImageRegistration()
{
	this->ndt_registration_->DeleteInstance();
}

void PanoramicImageRegistration::SetComponents() {
	this->button_load_primitive_image_ = this->ui.button_load_primitive_image;
	this->button_labeled_image_ = this->ui.button_labeled_image;
	this->button_boundaries_up_ = this->ui.button_boundaries_up;
	this->button_boundaries_down_ = this->ui.button_boundaries_down;
	this->button_boundaries_left_ = this->ui.button_boundaries_left;
	this->button_boundaries_right_ = this->ui.button_boundaries_right;
	this->button_registering_with_NDT_algorithm_ = this->ui.button_registering_with_NDT_algorithm;
	this->label_image_panel1_ = this->ui.label_image_panel1;
	this->label_image_panel2_ = this->ui.label_image_panel2;
	this->lineedit_cell_num_w_ = this->ui.lineedit_cell_num_w;
	this->lineedit_cell_num_h_ = this->ui.lineedit_cell_num_h;
	this->lineedit_cell_num_addition_w_ = this->ui.lineedit_cell_num_addition_w;
	this->lineedit_cell_num_addition_h_ = this->ui.lineedit_cell_num_addition_h;
	this->lineedit_cell_num_count_ = this->ui.lineedit_cell_num_count;
	this->lineedit_max_iterator_times_ = this->ui.lineedit_max_iterator_times;
	this->button_load_obj_points_ = this->ui.button_load_obj_points;
	this->ndt_registration_ = NDTRegistration::GetInstance();
}

void PanoramicImageRegistration::SetConnections() {
	connect(this->button_load_primitive_image_, SIGNAL(clicked(void)),
		this, SLOT(OnLoadPrimiticeImage(void)));
	connect(this->button_labeled_image_, SIGNAL(clicked(void)),
		this, SLOT(OnLoadLabeledImage(void)));
	connect(this->button_boundaries_up_, SIGNAL(clicked(void)),
		this, SLOT(OnTransformAndRotateBoundaries(void)));
	connect(this->button_boundaries_down_, SIGNAL(clicked(void)),
		this, SLOT(OnTransformAndRotateBoundaries(void)));
	connect(this->button_boundaries_left_, SIGNAL(clicked(void)),
		this, SLOT(OnTransformAndRotateBoundaries(void)));
	connect(this->button_boundaries_right_, SIGNAL(clicked(void)),
		this, SLOT(OnTransformAndRotateBoundaries(void)));
	connect(this->button_registering_with_NDT_algorithm_, SIGNAL(clicked(void)),
		this, SLOT(OnRegisterWithNDTAlgorithm(void)));
	connect(this->button_load_obj_points_, SIGNAL(clicked(void)),
		this, SLOT(OnLoadOBJPoints(void)));
}

void PanoramicImageRegistration::OnLoadPrimiticeImage(void) {
	QString path = QFileDialog::getOpenFileName(NULL, "load primitive image", ".", 
		"Files(*.png *.jpg *.bmp )");
	if (0 == path.size()) return;
	this->primitive_image_ = new QImage(path);
	int image_width = this->primitive_image_->width();
	int image_height = this->primitive_image_->height();
	vector<double> image_data(image_width * image_height, 0.0);
	for (int w = 0; w < image_width; ++w) {
		for (int h = 0; h < image_height; ++h) {
			image_data[w * image_height + h] = QColor(this->primitive_image_->pixel(w, h)).red();
		}
	}
	this->ndt_registration_->SetPrimitiveImage(image_data, image_width, image_height);
	this->label_image_panel1_->setPixmap(QPixmap::fromImage(*this->primitive_image_));
}

void PanoramicImageRegistration::OnLoadLabeledImage(void) {
	QString path = QFileDialog::getOpenFileName(NULL, "load labeled image", ".",
		"Files(*.png *.jpg *.bmp )");
	if (0 == path.size()) return;
	this->labeled_image_ = new QImage(path);
	int image_width = this->labeled_image_->width();
	int image_height = this->labeled_image_->height();
	vector<Eigen::Vector2d> boundary_points;
	for (int w = 0; w < image_width; ++w) {
		for (int h = 0; h < image_height; ++h) {
			if (QColor(this->labeled_image_->pixel(w, h)).red() > 200 &&
				QColor(this->labeled_image_->pixel(w, h)).green() < 100 &&
				QColor(this->labeled_image_->pixel(w, h)).blue() < 100) {
				boundary_points.push_back(Eigen::Vector2d(w, h));
			}
		}
	}
	this->ndt_registration_->SetBoundaryPoints(boundary_points);
	this->label_image_panel2_->setPixmap(QPixmap::fromImage(*this->labeled_image_));
}

void PanoramicImageRegistration::OnLoadOBJPoints(void) {
	QString path = QFileDialog::getOpenFileName(NULL, "load obj points", ".",
		"Files(*.obj)");
	std::ifstream in(path.toLocal8Bit());
	double minw = 1e30, maxw = -1e30;
	double minh = 1e30, maxh = -1e30;
	double valuew, valueh;
	vector<Eigen::Vector2d> boundary_points;
	while (in >> valuew >> valueh) {
		valuew *= 3;
		valueh *= 3;
		valueh = -valueh;
		minw = std::min(minw, valuew);
		maxw = std::max(maxw, valuew);
		minh = std::min(minh, valueh);
		maxh = std::max(maxh, valueh);
		boundary_points.push_back(Eigen::Vector2d(valuew, valueh));
	}
	int points_count = (int)boundary_points.size();
	int image_width = this->ndt_registration_->GetImageWidth();
	int image_height = this->ndt_registration_->GetImageHeight();
	QImage image(*primitive_image_);
	for (int i = 0; i < points_count; ++i) {
		boundary_points[i](0) -= minw;
		boundary_points[i](1) -= minh;
		int idx = (int)boundary_points[i](0);
		int idy = (int)boundary_points[i](1);
		if (idx >= 0 && idx < image_width && idy >= 0 && idy < image_height) {
			image.setPixel(idx, idy, QColor(255, 0, 0).rgb());
		}
	}
	this->ndt_registration_->SetBoundaryPoints(boundary_points);
	this->label_image_panel2_->setPixmap(QPixmap::fromImage(image));
}

void PanoramicImageRegistration::OnTransformAndRotateBoundaries(void) {
	QPushButton* sender = qobject_cast<QPushButton*>(this->sender());
	if (sender == this->button_boundaries_up_) {
		this->ndt_registration_->TransformAndRotateBoundaryPoints(1, 3);
	}
	else if (sender == this->button_boundaries_down_) {
		this->ndt_registration_->TransformAndRotateBoundaryPoints(2, 3);
	}
	else if (sender == this->button_boundaries_left_) {
		this->ndt_registration_->TransformAndRotateBoundaryPoints(3, 3);
	}
	else if (sender == this->button_boundaries_right_) {
		this->ndt_registration_->TransformAndRotateBoundaryPoints(4, 3);
	}
	this->ShowBoundaryPointsInImage();
}

void PanoramicImageRegistration::OnRegisterWithNDTAlgorithm(void) {
	int cell_init_num_w = this->lineedit_cell_num_w_->text().toInt();
	int cell_init_num_h = this->lineedit_cell_num_h_->text().toInt();
	int cell_num_addition_w = this->lineedit_cell_num_addition_w_->text().toInt();
	int cell_num_addition_h = this->lineedit_cell_num_addition_h_->text().toInt();
	int cell_num_count = this->lineedit_cell_num_count_->text().toInt();
	int max_iterator_times = this->lineedit_max_iterator_times_->text().toInt();
	this->ndt_registration_->SetCellPartitionParameters(cell_init_num_w,
		cell_init_num_h, cell_num_addition_w, cell_num_addition_h, cell_num_count);
	this->ndt_registration_->SetMaxIteratorTimes(max_iterator_times);
	this->ndt_registration_->ExecuteNDTRegister();
	this->ShowBoundaryPointsInImage();
}

void PanoramicImageRegistration::ShowBoundaryPointsInImage() {
	QImage image(*primitive_image_);
	vector<Eigen::Vector2d>& boundary_points = this->ndt_registration_->GetBoundaryPoints();
	int image_width = this->ndt_registration_->GetImageWidth();
	int image_height = this->ndt_registration_->GetImageHeight();
	int boundary_points_count = boundary_points.size();
	for (int p = 0; p < boundary_points_count; ++p) {
		int idx = (int)boundary_points[p](0);
		int idy = (int)boundary_points[p](1);
		if (idx >= 0 && idx < image_width && idy >= 0 && idy < image_height) {
			image.setPixel(idx, idy, QColor(255, 0, 0).rgb());
		}
	}
	this->label_image_panel1_->setPixmap(QPixmap::fromImage(image));
}