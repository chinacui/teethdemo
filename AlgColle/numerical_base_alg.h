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
};
#endif