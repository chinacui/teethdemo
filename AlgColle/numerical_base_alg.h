#ifndef NUMERICAL_BASE_ALG_H
#define NUMERICAL_BASE_ALG_H
#include"prereq.h"
#include<iostream>
#include<vector>
#include<Eigen/Core>
#include"../DataColle/custom_openmesh_type.h"
class ALGCOLLE_CLASS CNumericalBaseAlg
{
protected:
	static double PI;
public:
	static double ComputeGaussian(double mean, double std_devi, double x);
	static double ComputeMean(std::vector<double>&values);
	static double ComputeStdDeviation(std::vector<double>&values);
	static void GetPlaneFromMeanAndDir(OpenMesh::Vec3d mean, OpenMesh::Vec3d dir,double & res_a,double &res_b,double &res_c,double &res_d);
	static void ComputeHistgram(Eigen::VectorXd &data, double bin_size, std::vector<std::pair<double, double>>&bins, std::vector<std::vector<int>>&bin_eles);
	static void NormalizeScalarField(Eigen::VectorXd &data);
	static void NormalizeScalarField(std::vector<double>&data);
	template<typename T>//T should have <
	static void ComputeOrder(std::vector<T>&data, std::vector<int>&res_order)
	{
		class CT
		{
		public:
			T data_;
			int id_;
			bool operator<(const CT&b)
			{
				return  data_ < b.data_;
			}
		};
		std::vector<CT> tmp_data;
		tmp_data.resize(data.size());
		for (int i = 0; i < data.size(); i++)
		{
			tmp_data[i].data_ = data[i];
			tmp_data[i].id_ = i;
		}
		std::sort(tmp_data.begin(), tmp_data.end());
		res_order.resize(data.size());
		for (int i = 0; i < res_order.size(); i++)
		{
			res_order[tmp_data[i].id_] = i;
		}
	}
	static void GetMaxValue(std::vector<double>&data, int &res_id);
	static void GetMinValue(std::vector<double>&data, int &res_id);
	static double Sigmoid(double x);
};
#endif