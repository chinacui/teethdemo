#ifndef MESH_OBJECT_H
#define MESH_OBJECT_H
#include <Eigen/Dense>
class CMeshObject 
{
protected:

public:
	Eigen::MatrixXd vertexs_;
	Eigen::MatrixXi faces_;

};
#endif