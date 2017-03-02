#include "panoramic_fitting.h"
#include"../AlgColle/geo_alg.h"
CNDTEnergyComputer::CNDTEnergyComputer(cv::Mat& panoramic_img)
{
	panoramic_img_ = panoramic_img;
	cell_num_w_ = 40;
	cell_num_h_ = 40;
}
void CNDTEnergyComputer::ComputeGradient()
{
	cv::Mat gradientx, gradienty;
	cv::Sobel(panoramic_img_, gradientx, CV_64FC1, 1, 0, 3);
	cv::Sobel(panoramic_img_, gradienty, CV_64FC1, 0, 1, 3);
	pan_grad_ = cv::Mat(gradientx.size(), CV_64FC1);
	for (int i = 0; i < pan_grad_.rows; i++)
	{
		for (int j = 0; j < pan_grad_.cols; j++)
		{
			pan_grad_.at<double>(i, j) = std::sqrt(std::pow(gradientx.at<double>(i, j), 2) + std::pow(gradienty.at<double>(i, j), 2));
		}
	}
	/*cv::imshow("aa", gradient);
	std::cerr << "grad " << gradient.rows << " " << gradient.cols << std::endl;
	std::cerr << "img " << img.rows << " " << img.cols << std::endl;*/

	ComputeCellMeanPoints();
	ComputeCellCovarianceMatrix();
}
bool CNDTEnergyComputer::ComputeCellMeanPoints()
{
	cell_mean_points_.resize(cell_num_w_ * cell_num_h_);
	cell_weights_.resize(cell_num_w_ * cell_num_h_, 0.0);
	fill(cell_weights_.begin(), cell_weights_.end(), 0.0);
	double cell_size_w = (panoramic_img_.cols + 1.0) / cell_num_w_;
	double cell_size_h = (panoramic_img_.rows + 1.0) / cell_num_h_;
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
			if (w_end >= panoramic_img_.cols) w_end = panoramic_img_.cols - 1;
			if (h_start < 0) h_start = 0;
			if (h_end >= panoramic_img_.rows) h_end = panoramic_img_.rows - 1;
			for (int x = w_start; x <= w_end; ++x) {
				for (int y = h_start; y <= h_end; ++y) {
					int pixel_id = x * panoramic_img_.rows + y;
					cell_mean_points_[cell_id](0) += pan_grad_.at<double>(y, x) * x;
					cell_mean_points_[cell_id](1) += pan_grad_.at<double>(y, x) * y;
					cell_weights_[cell_id] += pan_grad_.at<double>(y, x);
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
}
void CNDTEnergyComputer::ChangeCellNum(int cell_num_w, int cell_num_h)
{
	cell_num_w_ = cell_num_w;
	cell_num_h_ = cell_num_h;
	ComputeCellMeanPoints();
	ComputeCellCovarianceMatrix();
}
bool CNDTEnergyComputer::ComputeCellCovarianceMatrix() {
	cell_covariance_matrix_.resize(cell_num_w_ * cell_num_h_ * 4, 0.0);
	fill(cell_covariance_matrix_.begin(), cell_covariance_matrix_.end(), 0.0);
	cell_inverse_covariance_matrix_.resize(cell_num_w_ * cell_num_h_ * 4, 0.0);
	fill(cell_inverse_covariance_matrix_.begin(), cell_inverse_covariance_matrix_.end(), 0.0);
	double cell_size_w = (panoramic_img_.cols + 1.0) / cell_num_w_;
	double cell_size_h = (panoramic_img_.rows + 1.0) / cell_num_h_;
	for (int w = 0; w < cell_num_w_; ++w) {
		for (int h = 0; h < cell_num_h_; ++h) {
			int cell_id = w * cell_num_h_ + h;
			int w_start = ceilf(w * cell_size_w);
			int w_end = floorf((w + 1) * cell_size_w);
			int h_start = ceilf(h * cell_size_h);
			int h_end = floorf((h + 1) * cell_size_h);
			if (w_start < 0) w_start = 0;
			if (w_end >= panoramic_img_.cols) w_end = panoramic_img_.cols - 1;
			if (h_start < 0) h_start = 0;
			if (h_end >= panoramic_img_.rows) h_end = panoramic_img_.rows - 1;
			for (int x = w_start; x <= w_end; ++x) {
				for (int y = h_start; y <= h_end; ++y) {
					int pixel_id = x * panoramic_img_.rows + y;
					double tempx = (x - cell_mean_points_[cell_id](0));
					double tempy = (y - cell_mean_points_[cell_id](1));
					cell_covariance_matrix_[cell_id * 4] += tempx * tempx * pan_grad_.at<double>(y, x);
					cell_covariance_matrix_[cell_id * 4 + 1] += tempx * tempy * pan_grad_.at<double>(y, x);
					cell_covariance_matrix_[cell_id * 4 + 2] += tempy * tempx * pan_grad_.at<double>(y, x);
					cell_covariance_matrix_[cell_id * 4 + 3] += tempy * tempy * pan_grad_.at<double>(y, x);
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

double CNDTEnergyComputer::ComputeNDTEnergy( std::vector<Eigen::Vector2d>&pts)
{
	
	double outlier_ratio, resolution, gauss_c1, gauss_c2, gauss_d1, gauss_d2, gauss_d3;
	// Initializes the guass fitting parameters (eq. 6.8) [Magnusson 2009]
	outlier_ratio = 0.55;
	resolution = (panoramic_img_.cols + 1.0) / (float)cell_num_w_;
	gauss_c1 = 10 * (1 - outlier_ratio);
	gauss_c2 = outlier_ratio / pow(resolution, 3);
	gauss_d3 = -log(gauss_c2);
	gauss_d1 = -log(gauss_c1 + gauss_c2) - gauss_d3;
	gauss_d2 = -2 * log((-log(gauss_c1 * exp(-0.5) + gauss_c2) - gauss_d3) / gauss_d1);




	int number_crown_boundary_points = pts.size();
	double cell_size_w = (panoramic_img_.cols + 1.0) / cell_num_w_;
	double cell_size_h = (panoramic_img_.rows + 1.0) / cell_num_h_;
	Eigen::Matrix2d c_inv;
	Eigen::Vector2d x, x_trans;
	std::vector<Eigen::Vector2d> boundary_points_pre = pts;
	
	double score = 0.0;

	for (int id = 0; id < number_crown_boundary_points; ++id) {
		x = Eigen::Vector2d(boundary_points_pre[id](0),
			boundary_points_pre[id](1));
		x_trans = Eigen::Vector2d(pts[id](0),
			pts[id](1));

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


		Eigen::Vector2d cov_dxd_pi;
		// e^(-d_2/2 * (x_k - mu_k)^T Sigma_k^-1 (x_k - mu_k)) Equation 6.9 [Magnusson 2009]
		double e_x_cov_x = exp(-gauss_d2 * x_trans.dot(c_inv * x_trans) / 2);
		// Calculate probability of transtormed points existance, Equation 6.9 [Magnusson 2009]
		double score_inc = -gauss_d1 * e_x_cov_x;
		score += score_inc;

	

	
	}
	
	
	return score;


}
double CPanoramicFittingOptimization::ComputeError(std::vector<double>&params)
{
	if (fitting_stage_ == InitialFitting)
	{
		projector_.SetParams(params);
		std::vector<OpenMesh::Vec2d>projected_ps;
		projected_ps.resize(points_.size());
		for (int i = 0; i < points_.size(); i++)
		{
			projected_ps[i] = projector_.ComputeProjParam(points_[i]);
		}
		double error = 0;
		for (int i = 1; i < points_.size(); i++)
		{
			double ratio = (projected_ps[i][1] - projected_ps[0][1]) / (projected_ps[i][0] - projected_ps[0][0] + 1e-7);
			double e = std::pow(ratio - ratio_tgt_proj_params_[i], 2);
			error += e;
		}
		return error;
	}
	else if(fitting_stage_== NTDRefine)
	{
		projector_.SetParams(params);
		std::vector<OpenMesh::Vec2d>projected_ps;
		projected_ps.resize(points_.size());
		for (int i = 0; i < points_.size(); i++)
		{
			projected_ps[i] = projector_.ComputeProjParam(points_[i]);
		}
		
		double errora = 0,errorb=0;
		for (int i = 1; i < points_.size(); i++)
		{
			double ratio = (projected_ps[i][1] - projected_ps[0][1]) / (projected_ps[i][0] - projected_ps[0][0] + 1e-7);
			double e = std::pow(ratio - ratio_tgt_proj_params_[i], 2);
			errora += e;
		}
		return errora+errorb;
	}
	
	
}
void CPanoramicFittingOptimization::GetPanoScaleAndTrans(double &res_scale, OpenMesh::Vec2d &res_trans)
{
	res_scale = pano_scale_;
	res_trans = pano_trans_;
}
void CPanoramicFittingOptimization::RunFitting()
{
	run();
	std::vector<OpenMesh::Vec2d>pano_coords;
	projector_.ComputeProjParams(points_, pano_coords);
	
	CGeoAlg::ComputeScaleAndTransformFrom2SetOfRegisteredPoints(tgt_proj_params_, pano_coords, pano_scale_, pano_trans_);


}
lbfgsfloatval_t CPanoramicFittingOptimization::evaluate(
		const lbfgsfloatval_t *x,
		lbfgsfloatval_t *g,
		const int n,
		const lbfgsfloatval_t step
		) 
{
	std::cerr << "evaluate " << std::endl;
	auto tx = const_cast<lbfgsfloatval_t*>(x);
	std::vector<double>params;
	for (int i = 0; i < m_variables.size(); i++)
	{
		params.push_back(tx[i]);
	}
	
	double error=ComputeError(params);
	for (int i = 0; i < params.size(); i++)
	{
		tx[i]=params[i];
	}

	double eps = 1e-8;
	std::cerr << "params ";
	for (int i = 0; i < params.size(); i++)
	{

		params[i] += eps;
		double e = ComputeError(params);
		params[i] -= eps;
		g[i]=(e - error) / eps;
		std::cerr << params[i]<<" ";
	}
	std::cerr<<std::endl;
	projector_.SetParams(params);



	
	return error;
}
