#ifndef CPANORAMIC_FITTING_H
#define CPANORAMIC_FITTING_H
#include"panoramic_simulation.h"
#include "../AdditionalLibs/lbfgs/non_linear_optimization.h"
#include <vector>
#include <iostream>
#include"prereq.h"
#include"opencv2\opencv.hpp"
class CNDTEnergyComputer
{
protected:
	cv::Mat panoramic_img_;
	cv::Mat pan_grad_;
	int cell_num_w_;
	int cell_num_h_;
	std::vector<double>cell_weights_;

	std::vector<Eigen::Vector2d> cell_mean_points_;
	// covariance matrix for all cells
	std::vector<double> cell_covariance_matrix_;
	// inverse matrix of covariance matrix for all cells
	std::vector<double> cell_inverse_covariance_matrix_;
	void ComputeGradient();
	bool ComputeCellCovarianceMatrix();
	bool ComputeCellMeanPoints();
public:
	CNDTEnergyComputer(cv::Mat& panoramic_img);
	double ComputeNDTEnergy(std::vector<Eigen::Vector2d>&pts);

	void ChangeCellNum(int cell_num_w, int cell_num_h);
};


class TEETHROOTRECOALG_CLASS CPanoramicFittingOptimization :public Non_Linear_Optimization
{
protected:
	enum FittingStage { InitialFitting, NTDRefine };
	FittingStage fitting_stage_= InitialFitting;
	std::vector<OpenMesh::Vec3d> points_;
	std::vector<OpenMesh::Vec2d>tgt_proj_params_;
	std::vector<double>ratio_tgt_proj_params_;
	CPanoramicProjectorBase &projector_;

	int proj_param_num_;
	double ComputeError(std::vector<double>&params);
	CNDTEnergyComputer ndt_energy_computer_;

	double pano_size_w_;
	double pano_size_h_;


	double pano_scale_;
	OpenMesh::Vec2d pano_trans_;
public:

	
	CPanoramicFittingOptimization(CPanoramicProjectorBase& projector, std::vector<OpenMesh::Vec3d>&points, std::vector<OpenMesh::Vec2d>&tgt_proj_params,cv::Mat &panoramic_img,double pano_size_w,double pano_size_h) :Non_Linear_Optimization(),projector_(projector), ndt_energy_computer_(panoramic_img)
	{
		pano_size_w_ = pano_size_w;
		pano_size_h_ = pano_size_h;
		m_variables.resize(6);

		projector.GetParams(m_variables);
		points_ = points;
		tgt_proj_params_ = tgt_proj_params;
		m_parameters.epsilon = 1e-6;
		ratio_tgt_proj_params_.resize(tgt_proj_params.size(), 0);
		for (int i = 1; i < tgt_proj_params.size(); i++)
		{
			ratio_tgt_proj_params_[i] = (tgt_proj_params[i][1] - tgt_proj_params[0][1]) / (tgt_proj_params[i][0] - tgt_proj_params[0][0] + 1e-7);
		}
		

	}
	void GetResProjector(CPanoramicProjectorBase&res_projector)
	{
		res_projector = projector_;
	}
	lbfgsfloatval_t evaluate(
		const lbfgsfloatval_t *x,
		lbfgsfloatval_t *g,
		const int n,
		const lbfgsfloatval_t step
		);
	std::vector<double> GetVariables() const
	{
		return m_variables;
	}
	void GetPanoScaleAndTrans(double &res_scale, OpenMesh::Vec2d &res_trans);
	void RunFitting();
};

#endif

