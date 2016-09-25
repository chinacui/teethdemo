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
	static std::map<int, std::shared_ptr<CMeshObject>> mesh_object_pool_;
	static int mesh_object_max_id_;
public:
	static int AddMeshObject(CMeshObject *mesh);
	static void Init();
	static std::map<int, std::shared_ptr<CMeshObject>>& GetMeshObjectPool();
};
#endif
