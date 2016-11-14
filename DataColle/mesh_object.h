// Copyright 2016_9 by ChenNenglun
#ifndef MESH_OBJECT_H
#define MESH_OBJECT_H
#include"prereq.h"
#include<Eigen/Dense>
#include"base_object.h"
#include<vector>
#include<set>
#include"custom_openmesh_type.h"
class CXWGeodesic;
class CAABBTree;
class DATACOLLE_CLASS CMeshObject:public CBaseObject
{
protected:
	int mesh_id_;// id of mesh object
	bool is_changed_;//if the geometry of mesh is changed , it should be true, thus the opengl render will re-computeds
	double tot_surface_area_;
	int texture_id_;
	bool use_texture_;

protected:
	COpenMeshT mesh_;


	CAABBTree* aabb_tree_=NULL;
	CXWGeodesic *geodesic_model_ = NULL;


public:

	int GetId();//get id of mesh object
	void SetId(int id);//set id of mesh object
	bool IsChanged();
	double GetTotSurfaceArea();
	void SetChanged(bool is_changed=true); //if the geometry of mesh is changed, it should be true, and the opengl render will re - computes
	void SetAttrChanged(bool is_attrchanged = true);
	void CopyFrom(CMeshObject& b);
	bool& UseTexture() { return use_texture_; }
	COpenMeshT& GetMesh() { return mesh_; }
	void SetMeshColor(OpenMesh::Vec3d color);

	CMeshObject();
	CMeshObject(CMeshObject &b);
	~CMeshObject();



	void NormalizeUV();
	CAABBTree* GetAABBTree() { return aabb_tree_; }
	void SetAABBTree(CAABBTree* t) { aabb_tree_ = t; }

	CXWGeodesic*GetGeodesicModel() { return geodesic_model_; }
	void SetGeodesicModel(CXWGeodesic *model) { geodesic_model_ = model; }
	

	int& TextureId() { return texture_id_; }
};


#endif