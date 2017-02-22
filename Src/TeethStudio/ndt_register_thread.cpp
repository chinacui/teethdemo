#include "ndt_register_thread.h"
#include <Eigen\Dense>
#include <iostream>
#include "QPointF"

NDTRegisterThread::NDTRegisterThread()
{
}


NDTRegisterThread::~NDTRegisterThread()
{
}

void NDTRegisterThread::SetData(int image_width, int image_height, int cell_num_w, int cell_num_h,
	int cell_num_addition_w, int cell_num_addition_h, int cell_num_count,
	vector<QPointF>* crown_boundary_points, vector<QPointF>* crown_boundary_points_pre,
	vector<QPointF>* cell_mean_points, vector<double>* cell_convariance_matrix,
	vector<double>* cell_inverse_convariance_matrix, vector<double>* image_gradient,
	vector<double>* cell_weights, int max_iterator_times) {
	this->image_width_ = image_width;
	this->image_height_ = image_height;
	this->cell_num_w_ = cell_num_w;
	this->cell_num_h_ = cell_num_h;
	this->cell_num_addition_w_ = cell_num_addition_w;
	this->cell_num_addition_h_ = cell_num_addition_h;
	this->cell_num_count_ = cell_num_count;
	this->crown_boundary_points_ = crown_boundary_points;
	this->crown_boundary_points_pre_ = crown_boundary_points_pre;
	this->cell_mean_points_ = cell_mean_points;
	this->cell_convariance_matrix_ = cell_convariance_matrix;
	this->cell_inverse_convariance_matrix_ = cell_inverse_convariance_matrix;
	this->image_gradient_ = image_gradient;
	this->cell_weights_ = cell_weights;
	this->max_iterator_times_ = max_iterator_times;
}

void NDTRegisterThread::ComputeCellMeanPoints() {
	(*cell_mean_points_).resize(cell_num_w_ * cell_num_h_);
	(*cell_weights_).resize(cell_num_w_ * cell_num_h_, 0.0);
	fill((*cell_weights_).begin(), (*cell_weights_).end(), 0.0);
	double cell_size_w = (image_width_ + 1.0) / cell_num_w_;
	double cell_size_h = (image_height_ + 1.0) / cell_num_h_;
	for (int w = 0; w < cell_num_w_; ++w) {
		for (int h = 0; h < cell_num_h_; ++h) {
			int cell_id = w * cell_num_h_ + h;
			(*cell_mean_points_)[cell_id].rx() = 0.0;
			(*cell_mean_points_)[cell_id].ry() = 0.0;
			int w_start = ceilf(w * cell_size_w);
			int w_end = floorf((w + 1) * cell_size_w);
			int h_start = ceilf(h * cell_size_h);
			int h_end = floorf((h + 1) * cell_size_h);
			if (w_start < 0) w_start = 0;
			if (w_end >= image_width_) w_end = image_width_ - 1;
			if (h_start < 0) h_start = 0;
			if (h_end >= image_height_) h_end = image_height_ - 1;
			for (int x = w_start; x <= w_end; ++x) {
				for (int y = h_start; y <= h_end; ++y) {
					int pixel_id = x * image_height_ + y;
					(*cell_mean_points_)[cell_id].rx() += (*image_gradient_)[pixel_id] * x;
					(*cell_mean_points_)[cell_id].ry() += (*image_gradient_)[pixel_id] * y;
					(*cell_weights_)[cell_id] += (*image_gradient_)[pixel_id];
				}
			}
		}
	}
	for (int w = 0; w < cell_num_w_; ++w) {
		for (int h = 0; h < cell_num_h_; ++h) {
			int cell_id = w * cell_num_h_ + h;
			if (0 == (*cell_weights_)[cell_id]) continue;
			(*cell_mean_points_)[cell_id].rx() /= (*cell_weights_)[cell_id];
			(*cell_mean_points_)[cell_id].ry() /= (*cell_weights_)[cell_id];
		}
	}
}

void NDTRegisterThread::ComputeCellConvarianceMatrix() {
	(*cell_convariance_matrix_).resize(cell_num_w_ * cell_num_h_ * 4, 0.0);
	fill((*cell_convariance_matrix_).begin(), (*cell_convariance_matrix_).end(), 0.0);
	(*cell_inverse_convariance_matrix_).resize(cell_num_w_ * cell_num_h_ * 4, 0.0);
	fill((*cell_inverse_convariance_matrix_).begin(), (*cell_inverse_convariance_matrix_).end(), 0.0);
	double cell_size_w = (image_width_ + 1.0) / cell_num_w_;
	double cell_size_h = (image_height_ + 1.0) / cell_num_h_;
	for (int w = 0; w < cell_num_w_; ++w) {
		for (int h = 0; h < cell_num_h_; ++h) {
			int cell_id = w * cell_num_h_ + h;
			int w_start = ceilf(w * cell_size_w);
			int w_end = floorf((w + 1) * cell_size_w);
			int h_start = ceilf(h * cell_size_h);
			int h_end = floorf((h + 1) * cell_size_h);
			if (w_start < 0) w_start = 0;
			if (w_end >= image_width_) w_end = image_width_ - 1;
			if (h_start < 0) h_start = 0;
			if (h_end >= image_height_) h_end = image_height_ - 1;
			for (int x = w_start; x <= w_end; ++x) {
				for (int y = h_start; y <= h_end; ++y) {
					int pixel_id = x * image_height_ + y;
					double tempx = (x - (*cell_mean_points_)[cell_id].rx());
					double tempy = (y - (*cell_mean_points_)[cell_id].ry());
					(*cell_convariance_matrix_)[cell_id * 4] += tempx * tempx * (*image_gradient_)[pixel_id];
					(*cell_convariance_matrix_)[cell_id * 4 + 1] += tempx * tempy * (*image_gradient_)[pixel_id];
					(*cell_convariance_matrix_)[cell_id * 4 + 2] += tempy * tempx * (*image_gradient_)[pixel_id];
					(*cell_convariance_matrix_)[cell_id * 4 + 3] += tempy * tempy * (*image_gradient_)[pixel_id];
				}
			}
		}
	}
	// calculate inverse covariance matrixs
	Eigen::Matrix2f cov_mat;
	Eigen::Matrix2f inv_cov_mat;
	Eigen::SelfAdjointEigenSolver<Eigen::Matrix2f> eigensolver;
	Eigen::Matrix2f eigen_val;
	Eigen::Matrix2f eigen_vecs;
	for (int w = 0; w < cell_num_w_; ++w) {
		for (int h = 0; h < cell_num_h_; ++h) {
			int cell_id = w * cell_num_h_ + h;
			if (1 == (*cell_weights_)[cell_id]) continue;
			(*cell_convariance_matrix_)[cell_id * 4] /= ((*cell_weights_)[cell_id] - 1);
			(*cell_convariance_matrix_)[cell_id * 4 + 1] /= ((*cell_weights_)[cell_id] - 1);
			(*cell_convariance_matrix_)[cell_id * 4 + 2] /= ((*cell_weights_)[cell_id] - 1);
			(*cell_convariance_matrix_)[cell_id * 4 + 3] /= ((*cell_weights_)[cell_id] - 1);
			cov_mat(0, 0) = (*cell_convariance_matrix_)[cell_id * 4];
			cov_mat(0, 1) = (*cell_convariance_matrix_)[cell_id * 4 + 1];
			cov_mat(1, 0) = (*cell_convariance_matrix_)[cell_id * 4 + 2];
			cov_mat(1, 1) = (*cell_convariance_matrix_)[cell_id * 4 + 3];
			//Normalize Eigen Val such that max no more than 100x min.
			eigensolver.compute(cov_mat);
			eigen_val = eigensolver.eigenvalues().asDiagonal();
			eigen_vecs = eigensolver.eigenvectors();

			float min_covar_eigvalue = 0.01 * eigen_val(1, 1);
			if (eigen_val(0, 0) < min_covar_eigvalue)
			{
				eigen_val(0, 0) = min_covar_eigvalue;
				cov_mat = eigen_vecs * eigen_val * eigen_vecs.inverse();
			}
			inv_cov_mat = cov_mat.inverse();
			(*cell_inverse_convariance_matrix_)[cell_id * 4] = inv_cov_mat(0, 0);
			(*cell_inverse_convariance_matrix_)[cell_id * 4 + 1] = inv_cov_mat(0, 1);
			(*cell_inverse_convariance_matrix_)[cell_id * 4 + 2] = inv_cov_mat(1, 0);
			(*cell_inverse_convariance_matrix_)[cell_id * 4 + 3] = inv_cov_mat(1, 1);
		}
	}
}

void NDTRegisterThread::ComputeOneNDTRegister() {
	double outlier_ratio, resolution, gauss_c1, gauss_c2, gauss_d1, gauss_d2, gauss_d3;
	// Initializes the guassian fitting parameters (eq. 6.8) [Magnusson 2009]
	outlier_ratio = 0.55;
	resolution = (image_width_ + 1.0) / (float)cell_num_w_;
	gauss_c1 = 10 * (1 - outlier_ratio);
	gauss_c2 = outlier_ratio / pow(resolution, 3);
	gauss_d3 = -log(gauss_c2);
	gauss_d1 = -log(gauss_c1 + gauss_c2) - gauss_d3;
	gauss_d2 = -2 * log((-log(gauss_c1 * exp(-0.5) + gauss_c2) - gauss_d3) / gauss_d1);

	Eigen::Matrix<double, 3, 1> p, delta_p, score_gradient;
	p.setZero();
	delta_p.setZero();
	Eigen::Matrix<double, 3, 3> hessian;

	// Initialize Point Gradient and Hessian
	Eigen::Matrix<double, 2, 3> point_gradient;
	point_gradient.setZero();
	point_gradient.block<2, 2>(0, 0).setIdentity();
	Eigen::Matrix<double, 6, 3> point_hessian;
	point_hessian.setZero();
	int number_crown_boundary_points = (*crown_boundary_points_).size();
	double cell_size_w = (image_width_ + 1.0) / cell_num_w_;
	double cell_size_h = (image_height_ + 1.0) / cell_num_h_;
	Eigen::Matrix2d c_inv;
	Eigen::Vector2d x, x_trans;
	crown_boundary_points_pre_ = crown_boundary_points_;
	for (int it = 0; it < this->max_iterator_times_; ++it) {
		double score = 0.0;
		hessian.setZero();
		score_gradient.setZero();
		for (int id = 0; id < number_crown_boundary_points; ++id) {
			x = Eigen::Vector2d((*crown_boundary_points_pre_)[id].rx(),
				(*crown_boundary_points_pre_)[id].ry());
			x_trans = Eigen::Vector2d((*crown_boundary_points_)[id].rx(),
				(*crown_boundary_points_)[id].ry());

			int w = (int)(x_trans(0) / cell_size_w);
			int h = (int)(x_trans(1) / cell_size_h);
			w = (w >= cell_num_w_) ? cell_num_w_ - 1 : w;
			h = (h >= cell_num_h_) ? cell_num_h_ - 1 : h;
			w = (w < 0) ? 0 : w;
			h = (h < 0) ? 0 : h;

			int cell_id = w * cell_num_h_ + h;
			x_trans(0) -= (*cell_mean_points_)[cell_id].rx();
			x_trans(1) -= (*cell_mean_points_)[cell_id].ry();
			c_inv << (*this->cell_inverse_convariance_matrix_)[cell_id * 4],
				(*this->cell_inverse_convariance_matrix_)[cell_id * 4 + 1],
				(*this->cell_inverse_convariance_matrix_)[cell_id * 4 + 2],
				(*this->cell_inverse_convariance_matrix_)[cell_id * 4 + 3];

			point_gradient(0, 2) = -x(0) * sin(p(2)) - x(1) * cos(p(2));
			point_gradient(1, 2) = x(0) * cos(p(2)) - x(1) * sin(p(2));
			point_hessian(4, 2) = -x(0) * cos(p(2)) + x(1) * sin(p(2));
			point_hessian(5, 2) = -x(0) * sin(p(2)) - x(1) * cos(p(2));

			Eigen::Vector2d cov_dxd_pi;
			// e^(-d_2/2 * (x_k - mu_k)^T Sigma_k^-1 (x_k - mu_k)) Equation 6.9 [Magnusson 2009]
			double e_x_cov_x = exp(-gauss_d2 * x_trans.dot(c_inv * x_trans) / 2);
			// Calculate probability of transtormed points existance, Equation 6.9 [Magnusson 2009]
			double score_inc = -gauss_d1 * e_x_cov_x;
			score += score_inc;

			e_x_cov_x = gauss_d2 * e_x_cov_x;

			// Error checking for invalid values.
			if (e_x_cov_x > 1 || e_x_cov_x < 0 || e_x_cov_x != e_x_cov_x) {
				continue;
			}

			// Reusable portion of Equation 6.12 and 6.13 [Magnusson 2009]
			e_x_cov_x *= gauss_d1;

			for (int i = 0; i < 3; i++) {
				// Sigma_k^-1 d(T(x,p))/dpi, Reusable portion of Equation 6.12 and 6.13 [Magnusson 2009]
				cov_dxd_pi = c_inv * point_gradient.col(i);

				// Update gradient, Equation 6.12 [Magnusson 2009]
				score_gradient(i) += x_trans.dot(cov_dxd_pi) * e_x_cov_x;

				for (int j = 0; j < hessian.cols(); j++) {
					// Update hessian, Equation 6.13 [Magnusson 2009]
					hessian(i, j) += e_x_cov_x * (-gauss_d2 * x_trans.dot(cov_dxd_pi) *
						x_trans.dot(c_inv * point_gradient.col(j)) +
						x_trans.dot(c_inv * point_hessian.block<2, 1>(2 * i, j)) +
						point_gradient.col(j).dot(cov_dxd_pi));
				}
			}
		}
		// Solve for decent direction using newton method, line 23 in Algorithm 2 [Magnusson 2009]
		Eigen::JacobiSVD<Eigen::Matrix<double, 3, 3> > sv(hessian, Eigen::ComputeFullU | Eigen::ComputeFullV);
		// Negative for maximization as opposed to minimization
		delta_p = sv.solve(-score_gradient);
		delta_p.normalize();
		p = p + 0.01 * delta_p;

		for (int id = 0; id < number_crown_boundary_points; ++id) {
			double x = (*crown_boundary_points_)[id].rx();
			double y = (*crown_boundary_points_)[id].ry();
			(*crown_boundary_points_)[id].rx() = x * cos(p(2)) - y * sin(p(2)) + p(0);
			(*crown_boundary_points_)[id].ry() = x * sin(p(2)) + y * cos(p(2)) + p(1);
		}
		std::cout << it << "th iterator, score sum is: " << score << std::endl;
		emit(UpdateOnce(it + 1, this->max_iterator_times_, cell_size_w, cell_size_h));
	}
	emit(UpdateOnce(this->max_iterator_times_, this->max_iterator_times_, cell_size_w, cell_size_h));
}

void NDTRegisterThread::run() {
	for (int i = 0; i < this->cell_num_count_; ++i) {
		this->ComputeCellMeanPoints();
		this->ComputeCellConvarianceMatrix();
		this->ComputeOneNDTRegister();

		this->cell_num_w_ += this->cell_num_addition_w_;
		this->cell_num_h_ += this->cell_num_addition_h_;
	}
}