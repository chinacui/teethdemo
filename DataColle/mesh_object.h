// Copyright 2016_9 by ChenNenglun
#ifndef MESH_OBJECT_H
#define MESH_OBJECT_H
#include"prereq.h"
#include<Eigen/Dense>
#include"base_object.h"

class DATACOLLE_CLASS CMeshObject:public CBaseObject
{
protected:
	int mesh_id_;// id of mesh object
	bool is_changed_;//if the geometry of mesh is changed , it should be true, thus the opengl render will re-computeds
public:
	Eigen::MatrixXd vertexs_;//n x 3 matrix,n is vertex number
	Eigen::MatrixXi faces_;//m x 3 matrix, m is the face number
	Eigen::MatrixXd vertex_colors_;//n x 3 matrix ,(rgb)
	int GetId();//get id of mesh object
	void SetId(int id);//set id of mesh object
	bool IsChanged();
	void SetChanged(bool is_changed=true); //if the geometry of mesh is changed, it should be true, and the opengl render will re - computeds
	CMeshObject();
	CMeshObject(CMeshObject &b);
	~CMeshObject();
};
#endif