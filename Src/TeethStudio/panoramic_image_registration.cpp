#include <iostream>
#include "panoramic_image_registration.h"
#include "ndt_register_thread.h"
#include "QObject"
#include "QFileDialog"
#include "QProgressBar"
#include "QLineEdit"
#include "QLabel"
#include "QRgb"
#include <Eigen\Dense>
#include <fstream>
#include <string>
#include <algorithm>

PanoramicImageRegistration* PanoramicImageRegistration::instance_ = NULL;

PanoramicImageRegistration::PanoramicImageRegistration(QWidget *parent)
	: QWidget(parent)
{
	ui.setupUi(this);

	this->SetComponents();
	this->SetConnections();
}

PanoramicImageRegistration::~PanoramicImageRegistration()
{
}

PanoramicImageRegistration* PanoramicImageRegistration::GetInstance() {
	if (instance_ == NULL) {
		instance_ = new PanoramicImageRegistration();
	}
	return instance_;
}

void PanoramicImageRegistration::DeleteInstance() {
	if (NULL != instance_) {
		delete instance_;
		instance_ = NULL;
	}
}

void PanoramicImageRegistration::SetComponents() {
	this->button_load_primitive_image_ = this->ui.button_load_primitive_image;
	this->button_labeled_image_ = this->ui.button_labeled_image;
	this->button_compute_crown_boundaries_ = this->ui.button_compute_crown_boundaries;
	this->button_boundaries_up_ = this->ui.button_boundaries_up;
	this->button_boundaries_down_ = this->ui.button_boundaries_down;
	this->button_boundaries_left_ = this->ui.button_boundaries_left;
	this->button_boundaries_right_ = this->ui.button_boundaries_right;
	this->button_registering_with_NDT_algorithm_ = this->ui.button_registering_with_NDT_algorithm;
	this->label_image_panel1_ = this->ui.label_image_panel1;
	this->label_image_panel2_ = this->ui.label_image_panel2;
	this->progress_bar_ = this->ui.progress_bar;
	this->lineedit_cell_num_w_ = this->ui.lineedit_cell_num_w;
	this->lineedit_cell_num_h_ = this->ui.lineedit_cell_num_h;
	this->lineedit_cell_num_addition_w_ = this->ui.lineedit_cell_num_addition_w;
	this->lineedit_cell_num_addition_h_ = this->ui.lineedit_cell_num_addition_h;
	this->lineedit_cell_num_count_ = this->ui.lineedit_cell_num_count;
	this->label_info_ = this->ui.label_info;
	this->lineedit_max_iterator_times_ = this->ui.lineedit_max_iterator_times;
	this->button_load_obj_points_ = this->ui.button_load_obj_points;
	this->ndt_register_thread_ = new NDTRegisterThread;
}

void PanoramicImageRegistration::SetConnections() {
	connect(this->button_load_primitive_image_, SIGNAL(clicked(void)),
		this, SLOT(OnLoadPrimiticeImage(void)));
	connect(this->button_labeled_image_, SIGNAL(clicked(void)),
		this, SLOT(OnLoadLabeledImage(void)));
	connect(this->button_compute_crown_boundaries_, SIGNAL(clicked(void)),
		this, SLOT(OnComputeCrownBoundaries(void)));
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
	connect(this->ndt_register_thread_, SIGNAL(UpdateOnce(int, int, double, double)),
		this, SLOT(OnUpdateImageOnce(int, int, double, double)));
	connect(this->button_load_obj_points_, SIGNAL(clicked(void)),
		this, SLOT(OnLoadOBJPoints(void)));
}

void PanoramicImageRegistration::OnLoadPrimiticeImage(void) {
	QString path = QFileDialog::getOpenFileName(NULL, "load primitive image", ".", 
		"Files(*.png *.jpg *.bmp )");
	if (0 == path.size()) return;
	this->primitive_image_ = new QImage(path);
	this->image_width_ = this->primitive_image_->width();
	this->image_height_ = this->primitive_image_->height();
	this->label_image_panel1_->setPixmap(QPixmap::fromImage(*this->primitive_image_));
}

void PanoramicImageRegistration::OnLoadLabeledImage(void) {
	QString path = QFileDialog::getOpenFileName(NULL, "load labeled image", ".",
		"Files(*.png *.jpg *.bmp )");
	if (0 == path.size()) return;
	this->labeled_image_ = new QImage(path);
	this->label_image_panel2_->setPixmap(QPixmap::fromImage(*this->labeled_image_));
}

void PanoramicImageRegistration::OnLoadOBJPoints(void) {
	QString path = QFileDialog::getOpenFileName(NULL, "load obj points", ".",
		"Files(*.obj)");
	std::ifstream in(path.toLocal8Bit());
	double minw = 1e30, maxw = -1e30;
	double minh = 1e30, maxh = -1e30;
	double valuew, valueh;
	this->crown_boundary_points_.clear();
	while (in >> valuew >> valueh) {
		valuew *= 3;
		valueh *= 3;
		valueh = -valueh;
		minw = std::min(minw, valuew);
		maxw = std::max(maxw, valuew);
		minh = std::min(minh, valueh);
		maxh = std::max(maxh, valueh);
		this->crown_boundary_points_.push_back(QPointF(valuew, valueh));
	}
	this->labeled_image_ = new QImage(*this->primitive_image_);
	int points_count = (int)this->crown_boundary_points_.size();
	for (int i = 0; i < points_count; ++i) {
		this->crown_boundary_points_[i].rx() -= minw;
		this->crown_boundary_points_[i].ry() -= minh;
		int idx = (int)this->crown_boundary_points_[i].rx();
		int idy = (int)this->crown_boundary_points_[i].ry();
		if (idx >= 0 && idx < image_width_ && idy >= 0 && idy < image_height_) {
			this->labeled_image_->setPixel(idx, idy, QColor(255, 0, 0).rgb());
		}
	}
	this->label_image_panel2_->setPixmap(QPixmap::fromImage(*this->labeled_image_));
}

void PanoramicImageRegistration::OnComputeCrownBoundaries(void) {
	int image_width = this->labeled_image_->width();
	int image_height = this->labeled_image_->height();
	this->boundary_label_.clear();
	this->boundary_label_.resize(image_width);
	for (auto i = 0; i < image_width; ++i) {
		this->boundary_label_[i].resize(image_height, false);
	}

	this->crown_boundary_image_ = new QImage(*this->primitive_image_);
	for (int w = 0; w < image_width; ++w) {
		for (int h = 0; h < image_height; ++h) {
			if (QColor(this->labeled_image_->pixel(w, h)).red() > 200 && 
				QColor(this->labeled_image_->pixel(w, h)).green() < 100 && 
				QColor(this->labeled_image_->pixel(w, h)).blue() < 100) {
				this->crown_boundary_image_->setPixel(w, h, QColor(255, 255, 0).rgb());
				this->boundary_label_[w][h] = true;
			}
		}
	}
	this->label_image_panel1_->setPixmap(QPixmap::fromImage(*this->crown_boundary_image_));
	// compute the crown boundary points according to 2D label array
	this->GetBoundaryPointsFromLabels();
}

void PanoramicImageRegistration::OnTransformAndRotateBoundaries(void) {
	QPushButton* sender = qobject_cast<QPushButton*>(this->sender());
	int image_width = this->labeled_image_->width();
	int image_height = this->labeled_image_->height();
	vector<vector<bool>> boundary_label_2(this->boundary_label_);
	for (int w = 0; w < image_width; ++w) {
		for (int h = 0; h < image_height; ++h) {
			this->boundary_label_[w][h] = false;
			this->crown_boundary_image_->setPixel(w, h, this->primitive_image_->pixel(w, h));
		}
	}
	for (int w = 0; w < image_width; ++w) {
		for (int h = 0; h < image_height; ++h) {
			if (boundary_label_2[w][h]) {
				if (sender == this->button_boundaries_up_ && h > 3) {
					this->crown_boundary_image_->setPixel(w, h - 3, QColor(255, 255, 0).rgb());
					this->boundary_label_[w][h - 3] = true;
				}
				else if (sender == this->button_boundaries_down_ && h < image_height - 3) {
					this->crown_boundary_image_->setPixel(w, h + 3, QColor(255, 255, 0).rgb());
					this->boundary_label_[w][h + 3] = true;
				}
				else if (sender == this->button_boundaries_left_ && w > 3) {
					this->crown_boundary_image_->setPixel(w - 3, h, QColor(255, 255, 0).rgb());
					this->boundary_label_[w - 3][h] = true;
				}
				else if (sender == this->button_boundaries_right_ && w < image_width - 3) {
					this->crown_boundary_image_->setPixel(w + 3, h, QColor(255, 255, 0).rgb());
					this->boundary_label_[w + 3][h] = true;
				}
			}
		}
	}
	this->label_image_panel1_->setPixmap(QPixmap::fromImage(*this->crown_boundary_image_));
	// compute the crown boundary points according to 2D label array
	this->GetBoundaryPointsFromLabels();
}

void PanoramicImageRegistration::GetBoundaryPointsFromLabels() {
	int image_width = this->labeled_image_->width();
	int image_height = this->labeled_image_->height();
	this->crown_boundary_points_.clear();
	for (int w = 0; w < this->image_width_; ++w) {
		for (int h = 0; h < this->image_height_; ++h) {
			if (this->boundary_label_[w][h]) {
				this->crown_boundary_points_.push_back(QPointF(w, h));
			}
		}
	}
}

void PanoramicImageRegistration::ComputeGradient() {
	this->image_gradient_.resize(this->image_width_ * this->image_height_);
	fill(this->image_gradient_.begin(), this->image_gradient_.end(), 0.0);
	for (int w = 0; w < this->image_width_; ++w) {
		for (int h = 0; h < this->image_height_; ++h) {
			int upw = w;
			int uph = (h - 1 >= 0) ? h - 1 : 0;
			int downw = w;
			int downh = (h + 1 <= this->image_height_ - 1) ? h + 1 : this->image_height_ - 1;
			int leftw = (w - 1 >= 0) ? w - 1 : 0;
			int lefth = h;
			int rightw = (w + 1 <= this->image_width_ - 1) ? w + 1 : this->image_width_ - 1;
			int righth = h;
			double gradientw = (QColor(this->primitive_image_->pixel(leftw, lefth)).red() -
				QColor(this->primitive_image_->pixel(rightw, righth)).red()) / 2.0;
			double gradienth = (QColor(this->primitive_image_->pixel(upw, uph)).red() -
				QColor(this->primitive_image_->pixel(downw, downh)).red()) / 2.0;
			double gradient = std::sqrt(gradientw * gradientw + gradienth * gradienth);
			this->image_gradient_[w * this->image_height_ + h] = gradient;
		}
	}
}

void PanoramicImageRegistration::ComputeCellMeanPoints() {
}

void PanoramicImageRegistration::ComputeCellConvarianceMatrix() {
}

void PanoramicImageRegistration::ComputeOneNDTRegister() {
}

void PanoramicImageRegistration::OnRegisterWithNDTAlgorithm(void) {
	// compute image gradient
	this->ComputeGradient();
	// compute mean value and variance for all cells
	this->cell_num_w_ = this->lineedit_cell_num_w_->text().toInt();
	this->cell_num_h_ = this->lineedit_cell_num_h_->text().toInt();
	int cell_num_addition_w = this->lineedit_cell_num_addition_w_->text().toInt();
	int cell_num_addition_h = this->lineedit_cell_num_addition_h_->text().toInt();
	int cell_num_count = this->lineedit_cell_num_count_->text().toInt();
	int max_iterator_times = this->lineedit_max_iterator_times_->text().toInt();
	this->ndt_register_thread_->SetData(this->image_width_, this->image_height_,
		this->cell_num_w_, this->cell_num_h_, cell_num_addition_w, cell_num_addition_h,
		cell_num_count, &this->crown_boundary_points_,
		&this->crown_boundary_points_pre_, &this->cell_mean_points_,
		&this->cell_convariance_matrix_, &this->cell_inverse_convariance_matrix_,
		&this->image_gradient_, &this->cell_weights_, max_iterator_times);
	this->ndt_register_thread_->start();
}

void PanoramicImageRegistration::OnUpdateImageOnce(int iterator_time, int max_times,
	double cell_size_w, double cell_size_h) {
	this->label_info_->setText(QString("Current Cell Size: ") +
		QString::number(cell_size_w) + QString("*") + QString::number(cell_size_h));
	this->progress_bar_->setRange(0, max_times);
	this->progress_bar_->setValue(iterator_time);
	if (iterator_time % 5 != 0) return;
	int number_crown_boundary_points = crown_boundary_points_.size();
	for (int w = 0; w < image_width_; ++w) {
		for (int h = 0; h < image_height_; ++h) {
			this->crown_boundary_image_->setPixel(w, h, this->primitive_image_->pixel(w, h));
			this->boundary_label_[w][h] = false;
		}
	}
	for (int id = 0; id < number_crown_boundary_points; ++id) {
		int idx = (int)crown_boundary_points_[id].rx();
		int idy = (int)crown_boundary_points_[id].ry();
		idx = (idx < 0) ? 0 : idx;
		idx = (idx >= image_width_) ? image_width_ - 1 : idx;
		idy = (idy < 0) ? 0 : idy;
		idy = (idy >= image_height_) ? image_height_ - 1 : idy;
		this->crown_boundary_image_->setPixel(idx, idy, QColor(255, 0, 0).rgb());
		this->boundary_label_[idx][idy] = true;
	}
	this->label_image_panel1_->setPixmap(QPixmap::fromImage(*this->crown_boundary_image_));
}