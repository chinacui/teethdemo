#ifndef EXTREMITY_POINT_UTILS_H
#define EXTREMITY_POINT_UTILS_H
#include "Xin_Wang.h"
class ExtremityPointUtils
{
protected:
	CRichModel &m_model;
public:
	ExtremityPointUtils(CRichModel &model);
	double ComputeMeanDistanceAtSourcePoint(int source, double radius);
	void ComputeMeanDistances(double radius, std::vector<double>& meanDis, int num_thread = 4);
};
#endif