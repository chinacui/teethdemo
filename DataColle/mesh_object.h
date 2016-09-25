#ifndef MESH_OBJECT_H
#define MESH_OBJECT_H
#include"prereq.h"
#include<Eigen/Dense>


class DATACOLLE_CLASS CMeshObject
{
protected:
	int mesh_id_;
	bool is_changed_;
public:
	Eigen::MatrixXd vertexs_;
	Eigen::MatrixXi faces_;
	Eigen::MatrixXi vertex_colors_;//rgb 0-255
	int GetId();
	void SetId(int id);
	bool IsChanged();
	void SetChanged(bool is_changed=true);
	CMeshObject();
	CMeshObject(CMeshObject &b);
	~CMeshObject();
};
#endif