#include"libicp_wrapper.h"
#include "icpPointToPoint.h"


CIcpPoint2Point::CIcpPoint2Point(Eigen::MatrixXd target_pts)
{
	int dim = target_pts.cols();
	target_pts_.resize(target_pts.rows() * dim);
	for (int i = 0; i < target_pts.rows(); i++)
	{
		for(int j=0;j<dim;j++)
		target_pts_[i * dim +j] = target_pts(i, j);
	}
	icp_p2p_ = new IcpPointToPoint(&target_pts_[0], target_pts.rows(), dim);
}
void CIcpPoint2Point::Fit(Eigen::MatrixXd src_pts, Eigen::Matrix3d& res_rot, Eigen::Vector3d& trans)
{
	int dim = src_pts.cols();
	src_pts_.resize(src_pts.rows() * dim);
	for (int i = 0; i < src_pts.rows(); i++)
	{
		for (int j = 0; j<dim; j++)
			src_pts_[i * dim + j] = src_pts(i, j);
	}
	ICPMatrix R = ICPMatrix::eye(3);
	ICPMatrix t(3, 1);
	icp_p2p_->fit(&src_pts_[0],src_pts.rows(),R,t,-1);
	
	FLOAT rdata[9];
	R.getData(rdata);
	FLOAT tdata[3];
	t.getData(tdata);

	
	trans(0)= tdata[0];
	trans(1) = tdata[1];
	trans(2) = tdata[2];

	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			res_rot(i, j) = rdata[i * 3 + j];
		}
	}


}
CIcpPoint2Point::~CIcpPoint2Point()
{
	if(icp_p2p_!=NULL)
	delete icp_p2p_;
}