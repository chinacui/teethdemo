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