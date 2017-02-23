#pragma once
#ifndef NDT_REGISTRATION_H_
#define NDT_REGISTRATION_H_
#include <iostream>
#include <vector>
#include <Eigen\Dense>
#include "prereq.h"
using std::vector;

class TEETHROOTRECOALG_CLASS NDTRegistration
{
public:
	~NDTRegistration();
	static NDTRegistration* GetInstance();
	static void DeleteInstance();
	bool SetPrimitiveImage(vector<double>& image_data, 
		const int image_width, const int image_height);
	int GetImageWidth();
	int GetImageHeight();
	vector<double>& GetPrimitiveImage();
	bool SetBoundaryPoints(vector<Eigen::Vector2d>& points);
	vector<Eigen::Vector2d>& GetBoundaryPoints();
	bool SetCellPartitionParameters(const int cell_init_num_w,
		const int cell_init_num_h, const int cell_add_num_w,
		const int cell_add_num_h, const int cell_add_count);
	bool ComputeImageGradient();
	bool ComputeCellMeanPoints();
	bool ComputeCellCovarianceMatrix();
	bool ComputeOneNDTRegister();
	bool ExecuteNDTRegister();
	bool SetMaxIteratorTimes(const int max_iterator_times);
	bool TransformAndRotateBoundaryPoints(int type, double distance);

private:
	NDTRegistration();

private:
	static NDTRegistration* instance_;
	// the image which points will be registered with
	vector<double> primitive_image_;
	// the points which will be registered with primitive image
	vector<Eigen::Vector2d> boundary_points_;
	// image gradient for primitive image
	vector<double> image_gradient_;
	// the sum of gradient in all cells
	vector<double> cell_weights_;
	// the mean points in all cells
	vector<Eigen::Vector2d> cell_mean_points_;
	// covariance matrix for all cells
	vector<double> cell_covariance_matrix_;
	// inverse matrix of covariance matrix for all cells
	vector<double> cell_inverse_covariance_matrix_;
	// the image width and height for primitive image
	int image_width_;
	int image_height_;
	// the cell num in width and height direction
	// for example, the width and height of primitive image is
	// 1000 and 800, and the cell num of width and height is 
	// 100 and 100. So the size of every cell is 10 * 8
	int cell_num_w_;
	int cell_num_h_;
	// the initial cell num, for example, 40 * 40
	int cell_init_num_w_;
	int cell_init_num_h_;
	// the increment of cell num in one step, generally it's 10 * 10
	int cell_add_num_w_;
	int cell_add_num_h_;
	// the number of increment steps, for example, 6
	int cell_add_count_;
	// the maximum iterator times in one step
	int max_iterator_times_;
};

#endif