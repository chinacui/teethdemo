// Copyright 2016_9 by ChenNenglun
#ifndef MESH_OBJECT_H
#define MESH_OBJECT_H
#include"prereq.h"
#include<Eigen/Dense>
#include"base_object.h"
#include<vector>
#include<set>
#include"custom_openmesh_type.h"
class DATACOLLE_CLASS CMeshObject:public CBaseObject
{
protected:
	int mesh_id_;// id of mesh object
	bool is_changed_;//if the geometry of mesh is changed , it should be true, thus the opengl render will re-computeds



protected:
	COpenMeshT mesh_;
public:

	int GetId();//get id of mesh object
	void SetId(int id);//set id of mesh object
	bool IsChanged();
	void SetChanged(bool is_changed=true); //if the geometry of mesh is changed, it should be true, and the opengl render will re - computes
	
	void CopyFrom(CMeshObject& b);
	COpenMeshT& GetMesh() { return mesh_; }

	CMeshObject();
	CMeshObject(CMeshObject &b);
	~CMeshObject();
};


#endif