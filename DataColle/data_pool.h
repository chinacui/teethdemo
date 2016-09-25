// Copyright 2016_9 by ChenNenglun
#ifndef DATA_POOL_H
#define DATA_POOL_H
#include"prereq.h"
#include<iostream>
#include<map>
#include<vector>
#include"mesh_object.h"
#include<memory>
class DATACOLLE_CLASS DataPool
{
protected:
	static std::map<int, std::shared_ptr<CMeshObject>> mesh_object_pool_;//pool of CMeshObject
	static int mesh_object_max_id_;// current max id in mesh_object_pool_
public:
	static int AddMeshObject(CMeshObject *mesh);// add CMeshObject to mesh_object_pool_
	static void Init();//initialization
	static std::map<int, std::shared_ptr<CMeshObject>>& GetMeshObjectPool();
};
#endif
