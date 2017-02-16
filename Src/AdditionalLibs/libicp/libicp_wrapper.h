#ifndef CLIB_ICP_WRAPPER_H
#define CLIB_ICP_WRAPPER_H
#include"../prereq.h"
#include <Eigen/Sparse>

class IcpPointToPoint;
class ADDITIONALLIBS_CLASS CIcpPoint2Point
{
protected:
	IcpPointToPoint *icp_p2p_=NULL;
	std::vector<double>target_pts_, src_pts_;
public:
	CIcpPoint2Point(Eigen::MatrixXd target_pts);
	void Fit(Eigen::MatrixXd src_pts,Eigen::Matrix3d& res_rot,Eigen::Vector3d& trans);

	~CIcpPoint2Point();

};

#endif
