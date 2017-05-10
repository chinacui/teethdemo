#include "ndt_registration.h"


NDTRegistration* NDTRegistration::instance_ = nullptr;

NDTRegistration::NDTRegistration()
{
}

NDTRegistration::~NDTRegistration()
{
}

NDTRegistration* NDTRegistration::GetInstance() {
	if (nullptr == instance_) {
		instance_ = new NDTRegistration();
	}
	return instance_;
}

void NDTRegistration::DeleteInstance() {
	if (nullptr != instance_) {
		delete instance_;
		instance_ = nullptr;
	}
}

bool NDTRegistration::SetPrimitiveImage(vector<double>& image_data,
	const int image_width, const int image_height) {
	this->primitive_image_ = image_data;
	this->image_width_ = image_width;
	this->image_height_ = image_height;
	return true;
}

int NDTRegistration::GetImageWidth() {
	return this->image_width_;
}

int NDTRegistration::GetImageHeight() {
	return this->image_height_;
}

vector<double>& NDTRegistration::GetPrimitiveImage() {
	return this->primitive_image_;
}

bool NDTRegistration::SetBoundaryPoints(vector<Eigen::Vector2d>& points) {
	this->boundary_points_ = points;
	return true;
}
bool NDTRegistration::setBoundaryRootPoint(vector<Eigen::Vector2d>& boundary_points_root)
{
	this->boundary_points_root = boundary_points_root;
	return true;
}

bool NDTRegistration::SetMaxIteratorTimes(const int max_iterator_times) {
	this->max_iterator_times_ = max_iterator_times;
	return true;
}

vector<Eigen::Vector2d>& NDTRegistration::GetBoundaryPoints() {
	return this->boundary_points_;
}

void NDTRegistration::SetBoundaryCenterPoints(Eigen::Vector2d & boundaryCenterPoints) {
	boundary_center_points_ = boundaryCenterPoints;
}

bool NDTRegistration::SetCellPartitionParameters(const int cell_init_num_w,
	const int cell_init_num_h, const int cell_add_num_w,
	const int cell_add_num_h, const int cell_add_count) {
	this->cell_init_num_w_ = cell_init_num_w;
	this->cell_init_num_h_ = cell_init_num_h;
	this->cell_add_num_w_ = cell_add_num_w;
	this->cell_add_num_h_ = cell_add_num_h;
	this->cell_add_count_ = cell_add_count;
	return true;
}

bool NDTRegistration::ComputeImageGradient() {
	this->image_gradient_.resize(this->image_width_ * this->image_height_);
	fill(this->image_gradient_.begin(), this->image_gradient_.end(), 0.0);
	for (int w = 0; w < this->image_width_; ++w) {
		for (int h = 0; h < this->image_height_; ++h) {
			double gradientw, gradienth, gradient;
			int upw = w;
			int uph = (h - 1 >= 0) ? h - 1 : 0;
			int downw = w;
			int downh = (h + 1 <= this->image_height_ - 1) ? h + 1 : this->image_height_ - 1;
			int leftw = (w - 1 >= 0) ? w - 1 : 0;
			int lefth = h;
			int rightw = (w + 1 <= this->image_width_ - 1) ? w + 1 : this->image_width_ - 1;
			int righth = h;
		    gradientw = abs(this->primitive_image_[rightw * this->image_height_ + righth] -
					this->primitive_image_[leftw * this->image_height_ + lefth]) / 2.0;
		    gradienth = abs(this->primitive_image_[upw * this->image_height_ + uph] -
					this->primitive_image_[downw * this->image_height_ + downh]) / 2.0;
			gradient = std::sqrt(gradientw * gradientw + gradienth * gradienth);
			this->image_gradient_[w * this->image_height_ + h] = gradient;
		}
	}
	return true;
}

bool NDTRegistration::ComputeCellMeanPoints() {
	cell_mean_points_.resize(cell_num_w_ * cell_num_h_);
	cell_weights_.resize(cell_num_w_ * cell_num_h_, 0.0);
	fill(cell_weights_.begin(), cell_weights_.end(), 0.0);
	double cell_size_w = (image_width_ + 1.0) / cell_num_w_;
	double cell_size_h = (image_height_ + 1.0) / cell_num_h_;
	for (int w = 0; w < cell_num_w_; ++w) {
		for (int h = 0; h < cell_num_h_; ++h) {
			int cell_id = w * cell_num_h_ + h;
			cell_mean_points_[cell_id](0) = 0.0;
			cell_mean_points_[cell_id](1) = 0.0;
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
					cell_mean_points_[cell_id](0) += image_gradient_[pixel_id] * x;
					cell_mean_points_[cell_id](1) += image_gradient_[pixel_id] * y;
					cell_weights_[cell_id] += image_gradient_[pixel_id];
				}
			}
		}
	}
	for (int w = 0; w < cell_num_w_; ++w) {
		for (int h = 0; h < cell_num_h_; ++h) {
			int cell_id = w * cell_num_h_ + h;
			if (0 == cell_weights_[cell_id]) continue;
			cell_mean_points_[cell_id](0) /= cell_weights_[cell_id];
			cell_mean_points_[cell_id](1) /= cell_weights_[cell_id];
		}
	}
	return true;
}

bool NDTRegistration::ComputeCellCovarianceMatrix() {
	cell_covariance_matrix_.resize(cell_num_w_ * cell_num_h_ * 4, 0.0);
	fill(cell_covariance_matrix_.begin(), cell_covariance_matrix_.end(), 0.0);
	cell_inverse_covariance_matrix_.resize(cell_num_w_ * cell_num_h_ * 4, 0.0);
	fill(cell_inverse_covariance_matrix_.begin(), cell_inverse_covariance_matrix_.end(), 0.0);
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
					double tempx = (x - cell_mean_points_[cell_id](0));
					double tempy = (y - cell_mean_points_[cell_id](1));
					cell_covariance_matrix_[cell_id * 4] += tempx * tempx * image_gradient_[pixel_id];
					cell_covariance_matrix_[cell_id * 4 + 1] += tempx * tempy * image_gradient_[pixel_id];
					cell_covariance_matrix_[cell_id * 4 + 2] += tempy * tempx * image_gradient_[pixel_id];
					cell_covariance_matrix_[cell_id * 4 + 3] += tempy * tempy * image_gradient_[pixel_id];
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
			if (1 == cell_weights_[cell_id]) continue;
			cell_covariance_matrix_[cell_id * 4] /= (cell_weights_[cell_id] - 1);
			cell_covariance_matrix_[cell_id * 4 + 1] /= (cell_weights_[cell_id] - 1);
			cell_covariance_matrix_[cell_id * 4 + 2] /= (cell_weights_[cell_id] - 1);
			cell_covariance_matrix_[cell_id * 4 + 3] /= (cell_weights_[cell_id] - 1);
			cov_mat(0, 0) = cell_covariance_matrix_[cell_id * 4];
			cov_mat(0, 1) = cell_covariance_matrix_[cell_id * 4 + 1];
			cov_mat(1, 0) = cell_covariance_matrix_[cell_id * 4 + 2];
			cov_mat(1, 1) = cell_covariance_matrix_[cell_id * 4 + 3];
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
			cell_inverse_covariance_matrix_[cell_id * 4] = inv_cov_mat(0, 0);
			cell_inverse_covariance_matrix_[cell_id * 4 + 1] = inv_cov_mat(0, 1);
			cell_inverse_covariance_matrix_[cell_id * 4 + 2] = inv_cov_mat(1, 0);
			cell_inverse_covariance_matrix_[cell_id * 4 + 3] = inv_cov_mat(1, 1);
		}
	}
	return true;
}

bool NDTRegistration::ComputeOneNDTRegister() {
	double outlier_ratio, resolution, gauss_c1, gauss_c2, gauss_d1, gauss_d2, gauss_d3;
	// Initializes the guass fitting parameters (eq. 6.8) [Magnusson 2009]
	outlier_ratio = 0.55;
	resolution = (image_width_ + 1.0) / (float)cell_num_w_;
	gauss_c1 = 10 * (1 - outlier_ratio);
	gauss_c2 = outlier_ratio / pow(resolution, 3);
	gauss_d3 = -log(gauss_c2);
	gauss_d1 = -log(gauss_c1 + gauss_c2) - gauss_d3;
	gauss_d2 = -2 * log((-log(gauss_c1 * exp(-0.5) + gauss_c2) - gauss_d3) / gauss_d1);


	Eigen::Matrix<double, 4, 1> p, delta_p, score_gradient;
	p.setZero();
	delta_p.setZero();
	Eigen::Matrix<double, 4, 4> hessian;

	// Initialize Point Gradient and Hessian
	Eigen::Matrix<double, 2, 4> point_gradient;
	point_gradient.setZero();
	point_gradient.block<2, 2>(0, 0).setIdentity();
	Eigen::Matrix<double, 8, 4> point_hessian;
	point_hessian.setZero();
	int number_crown_boundary_points = boundary_points_.size();
	double cell_size_w = (image_width_ + 1.0) / cell_num_w_;
	double cell_size_h = (image_height_ + 1.0) / cell_num_h_;
	Eigen::Matrix2d c_inv;
	Eigen::Vector2d x, x_trans;
	vector<Eigen::Vector2d> boundary_points_pre = boundary_points_;
	for (int it = 0; it < this->max_iterator_times_; ++it) {
		double score = 0.0;
		hessian.setZero();
		score_gradient.setZero();
		for (int id = 0; id < number_crown_boundary_points; ++id) {
			x = Eigen::Vector2d(boundary_points_pre[id](0),
				boundary_points_pre[id](1));
			x_trans = Eigen::Vector2d(boundary_points_[id](0),
				boundary_points_[id](1));

			int w = (int)(x_trans(0) / cell_size_w);
			int h = (int)(x_trans(1) / cell_size_h);
			w = (w >= cell_num_w_) ? cell_num_w_ - 1 : w;
			h = (h >= cell_num_h_) ? cell_num_h_ - 1 : h;
			w = (w < 0) ? 0 : w;
			h = (h < 0) ? 0 : h;
			int cell_id = w * cell_num_h_ + h;
			x_trans(0) -= cell_mean_points_[cell_id](0);
			x_trans(1) -= cell_mean_points_[cell_id](1);
			c_inv << this->cell_inverse_covariance_matrix_[cell_id * 4],
				this->cell_inverse_covariance_matrix_[cell_id * 4 + 1],
				this->cell_inverse_covariance_matrix_[cell_id * 4 + 2],
				this->cell_inverse_covariance_matrix_[cell_id * 4 + 3];

			point_gradient(0, 2) = -x(0) * sin(p(2)) - x(1) * cos(p(2));
			point_gradient(1, 2) = x(0) * cos(p(2)) - x(1) * sin(p(2));
			point_gradient(0, 3) = x(0) - boundary_center_points_(0);
			point_gradient(1, 3) = x(1) - boundary_center_points_(1);
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
			for (int i = 0; i < 4; i++) {
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
		
		//score = score / number_crown_boundary_points;
		// Solve for decent direction using newton method, line 23 in Algorithm 2 [Magnusson 2009]
		Eigen::JacobiSVD<Eigen::Matrix<double, 4, 4> > sv(hessian, Eigen::ComputeFullU | Eigen::ComputeFullV);
		// Negative for maximization as opposed to minimization
		delta_p = sv.solve(-score_gradient);
		delta_p.normalize();
		p = p + 0.01 * delta_p;
		p(3) = 0;
		if (it == 1 && cell_add_count_end_ == 0)
		{
			image_to_model_ro(0, 0) = cos(p(2)) + p(3);
			image_to_model_ro(0, 1) = -1 * sin(p(2));
			image_to_model_ro(1, 0) = sin(p(2));
			image_to_model_ro(1, 1) = cos(p(2)) + p(3);

			image_to_model_tra(0, 0) = p(0) - p(3)*boundary_center_points_(0);
			image_to_model_tra(1, 0) = p(1) - p(3)*boundary_center_points_(1);
		}
		else
		{
			Eigen::Matrix<double, 2, 2> temp_ro;
			Eigen::Matrix<double, 2, 1> temp_tra;
			temp_ro(0, 0) = cos(p(2)) + p(3);
			temp_ro(0, 1) = -1 * sin(p(2));
			temp_ro(1, 0) = sin(p(2));
			temp_ro(1, 1) = cos(p(2)) + p(3);

			temp_tra(0, 0) = p(0) - p(3)*boundary_center_points_(0);
			temp_tra(1, 0) = p(1) - p(3)*boundary_center_points_(1);

			image_to_model_ro = temp_ro*image_to_model_ro;
			image_to_model_tra = temp_ro*image_to_model_tra + temp_tra;
		}
		for (int i = 0; i < number_crown_boundary_points; ++i) {
			double x = boundary_points_[i](0);
			double y = boundary_points_[i](1);
			boundary_points_[i](0) = x * cos(p(2)) - y * sin(p(2)) + p(0) + p(3) * (x - boundary_center_points_(0));
			boundary_points_[i](1) = x * sin(p(2)) + y * cos(p(2)) + p(1) + p(3) * (y - boundary_center_points_(1));
		}
	}

	for (int i = 0; i < number_crown_boundary_points; ++i)
	{
		Eigen::Vector2d temp;
		temp(0) = image_to_model_ro.inverse()(0, 0)*(boundary_points_[i](0) - image_to_model_tra(0, 0)) + image_to_model_ro.inverse()(0, 1)*(boundary_points_[i](1) - image_to_model_tra(1, 0));
		temp(1) = image_to_model_ro.inverse()(1, 0)*(boundary_points_[i](0) - image_to_model_tra(0, 0)) + image_to_model_ro.inverse()(1, 1)*(boundary_points_[i](1) - image_to_model_tra(1, 0));
		if (cell_add_count_end_ == cell_add_count_ - 1)
		{
			boundary_points_image_to_model.push_back(temp);

		}
	}
		for (int i = 0; i < boundary_points_root.size(); i++)
		{
			Eigen::Vector2d temp;
			temp(0) = image_to_model_ro.inverse()(0, 0)*(boundary_points_root[i](0) - image_to_model_tra(0, 0)) + image_to_model_ro.inverse()(0, 1)*(boundary_points_root[i](1) - image_to_model_tra(1, 0));
			temp(1) = image_to_model_ro.inverse()(1, 0)*(boundary_points_root[i](0) - image_to_model_tra(0, 0)) + image_to_model_ro.inverse()(1, 1)*(boundary_points_root[i](1) - image_to_model_tra(1, 0));
			if (cell_add_count_end_ == cell_add_count_ - 1)
			{
				boundary_points_image_to_model.push_back(temp);
			}
		}
	
	return true;
}

bool NDTRegistration::ExecuteNDTRegister() {
	// compute image gradient
	this->ComputeImageGradient();
	this->cell_num_w_ = this->cell_init_num_w_;
	this->cell_num_h_ = this->cell_init_num_h_;
	for (int i = 0; i < this->cell_add_count_; ++i) {
		this->cell_add_count_end_ = i;
		this->ComputeCellMeanPoints();
		this->ComputeCellCovarianceMatrix();
		this->ComputeOneNDTRegister();

		this->cell_num_w_ += this->cell_add_num_w_;
		this->cell_num_h_ += this->cell_add_num_h_;
	}
	return true;
}

bool NDTRegistration::TransformAndRotateBoundaryPoints(int type, double distance) {
	// 1: up, 2: down, 3: left, 4: right
	int boundary_points_count = this->boundary_points_.size();
	for (int p = 0; p < boundary_points_count; ++p) {
		if (1 == type) {
			this->boundary_points_[p](1) -= distance;
		}
		else if (2 == type) {
			this->boundary_points_[p](1) += distance;
		}
		else if (3 == type) {
			this->boundary_points_[p](0) -= distance;
		}
		else if (4 == type) {
			this->boundary_points_[p](0) += distance;
		}
	}
	return true;
}