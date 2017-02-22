#pragma once
#ifndef NDT_REGISTER_THREAD_H_
#define NDT_REGISTER_THREAD_H_
#include <QThread>
#include <iostream>
using std::vector;
class QPointF;

class NDTRegisterThread : public QThread
{
	Q_OBJECT

public:
	NDTRegisterThread();
	~NDTRegisterThread();
	void SetData(int image_width, int image_height, int cell_num_w, int cell_num_h,
		int cell_num_addition_w, int cell_num_addition_h, int cell_num_count,
		vector<QPointF>* crown_boundary_points, vector<QPointF>* crown_boundary_points_pre,
		vector<QPointF>* cell_mean_points, vector<double>* cell_convariance_matrix,
		vector<double>* cell_inverse_convariance_matrix, vector<double>* image_gradient,
		vector<double>* cell_weights, int max_iterator_times);
	void ComputeCellMeanPoints();
	void ComputeCellConvarianceMatrix();
	void ComputeOneNDTRegister();

signals:
	void UpdateOnce(int iterator_time, int max_times, double cell_size_w, double cell_size_h);

protected:
	virtual void run();

private:
	int image_width_;
	int image_height_;
	int cell_num_w_;
	int cell_num_h_;
	int cell_num_addition_w_;
	int cell_num_addition_h_;
	int cell_num_count_;
	vector<QPointF>* crown_boundary_points_;
	vector<QPointF>* crown_boundary_points_pre_;
	vector<QPointF>* cell_mean_points_;
	vector<double>* cell_convariance_matrix_;
	vector<double>* cell_inverse_convariance_matrix_;
	vector<double>* image_gradient_;
	vector<double>* cell_weights_;
	int max_iterator_times_;
};

#endif
