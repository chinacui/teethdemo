#include"numerical_base_alg.h"
#include<math.h>
#include<iostream>
#include<Eigen\Core>
double CNumericalBaseAlg::PI = 3.14159265359;
double CNumericalBaseAlg::ComputeGaussian(double mean, double std_devi, double x)
{	
	double res=(1.0 / (std_devi*std::sqrt(2*PI)))*std::exp(-std::pow(x-mean,2)/(2*std_devi*std_devi));
	return res;
}
void CNumericalBaseAlg::GetMaxValue(std::vector<double>&data, int &res_id)
{
	res_id = 0;
	for (int i = 0; i < data.size(); i++)
	{
		if (data[res_id] < data[i])
		{
			res_id = i;
		}
	}
}
void CNumericalBaseAlg::GetMinValue(std::vector<double>&data, int &res_id)
{
	res_id = 0;
	for (int i = 0; i < data.size(); i++)
	{
		if (data[res_id] > data[i])
		{
			res_id = i;
		}
	}
}
double CNumericalBaseAlg::Sigmoid(double x)
{
	return 1.0 / (1.0 + std::exp(-x));
}

void CNumericalBaseAlg::NormalizeScalarField(Eigen::VectorXd &data)
{
	double mmin = std::numeric_limits<double>::max();
	double mmax = std::numeric_limits<double>::min();
	int maxi, mini;
	for (int i = 0; i < data.size(); i++)
	{
		if (mmin > data(i))
		{
			mmin = data(i);
			mini = i;
		}
			
		if (mmax < data(i))
		{
			mmax = data(i);
			maxi = i;
		}
		
	}
	double len = mmax - mmin;
	
	if (len != 0)
	{
		int data_size = data.size();
		for (int i = 0; i<data_size;i++)
		{
			data(i) = (data(i) - mmin) / len;
		}
	}

}
void CNumericalBaseAlg::NormalizeScalarField(std::vector<double> &data)
{
	double mmin = std::numeric_limits<double>::max();
	double mmax = std::numeric_limits<double>::min();
	for (int i = 0; i < data.size(); i++)
	{
		if (mmin > data[i])
			mmin = data[i];
		if (mmax < data[i])
			mmax = data[i];
	}
	double len = mmax - mmin;

	if (len != 0)
	{
		int data_size = data.size();
		for (int i = 0; i < data_size; i++)
		{
			data[i] = (data[i] - mmin) / len;
		}
	}

}
void CNumericalBaseAlg::ComputeHistgram(Eigen::VectorXd &data, double bin_size, std::vector<std::pair<double, double>>&bins, std::vector<std::vector<int>>&bin_eles)
{
	double mmin = std::numeric_limits<double>::max();
	double mmax = std::numeric_limits<double>::min();
	for (int i = 0; i < data.size(); i++)
	{
		if (data(i) > mmax)
			mmax = data(i);
		if (data(i) < mmin)
			mmin = data(i);
	}
	bins.clear();
	bin_eles.clear();
	int bin_num = std::ceil((mmax - mmin) / bin_size);

	if ((mmax - mmin) / bin_size == bin_num*1.0)
		bin_num++;
	
	for (int i = 0; i < bin_num; i++)
	{
		double s, t;
		s = i*bin_size+mmin;
		t = (i + 1)*bin_size+mmin;
		bins.push_back(std::make_pair(s, t));
		bin_eles.push_back(std::vector<int>());
		bin_eles.back().clear();
	}
	
	for (int i = 0; i < data.size(); i++)
	{
		int bin_id=(data(i)- mmin) / bin_size;
		bin_eles[bin_id].push_back(i);

	}
}
double CNumericalBaseAlg::ComputeStdDeviation(std::vector<double>&values)
{
	double mean = ComputeMean(values);
	double res = 0;
	for (int i = 0; i < values.size(); i++)
	{
		res += std::pow(values[i] - mean, 2);
	}
	res /= values.size();//Uncorrected sample standard deviation
	res = std::sqrt(res);
	return res;
}
double CNumericalBaseAlg::ComputeMean(std::vector<double>&values)
{
	double sum = 0;
	for (int i = 0; i < values.size(); i++)
	{
		sum += values[i];
	}
	double ave=sum / values.size();
	return ave;
}
void CNumericalBaseAlg::GetPlaneFromMeanAndDir(OpenMesh::Vec3d mean, OpenMesh::Vec3d dir, double & res_a, double &res_b, double &res_c, double &res_d)
{
	res_a = dir[0];
	res_b = dir[1];
	res_c = dir[2];
	res_d = -dir[0]*mean[0] - dir[1]*mean[1] - dir[2]*mean[2];
}