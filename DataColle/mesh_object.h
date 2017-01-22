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
	bool is_visiable_;
	std::map<int,int>vtag_;
	bool is_shinning_=false;
	bool is_pickable_ = true;
	std::map<COpenMeshT::VertexHandle, OpenMesh::Vec3d>restore_pos_;
protected:
	COpenMeshT mesh_;


	CAABBTree* aabb_tree_=NULL;
	CXWGeodesic *geodesic_model_ = NULL;


public:
	void RestoreCurrentVPos();
	void RecoverCurrentVPos();
	std::map<int,int>&GetVertexTags() { return vtag_; }
	bool& IsVisiable() { return is_visiable_; }
	bool &IsShinning() { return is_shinning_; }
	int GetId();//get id of mesh object
	void SetId(int id);//set id of mesh object
	bool IsChanged();
	bool& IsPickAble() { return is_pickable_; }
	double GetTotSurfaceArea();
	void SetChanged(bool is_changed=true); //if the geometry of mesh is changed, it should be true, and the opengl render will re - computes
	void SetAttrChanged(bool is_attrchanged = true);
	OpenMesh::Vec3d TransformPointByLocalMatrix(OpenMesh::Vec3d p);
	void CopyFrom(CMeshObject& b);
	bool& UseTexture() { return use_texture_; }
	COpenMeshT& GetMesh() { return mesh_; }
	void SetMeshColor(OpenMesh::Vec3d color);
	void ApplyTransform();
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