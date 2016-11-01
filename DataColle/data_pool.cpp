// Copyright 2016_9 by ChenNenglun
#include"data_pool.h"

std::map<int, std::shared_ptr<CMeshObject>> DataPool::mesh_object_pool_;
int DataPool::mesh_object_max_id_;

std::map<int, std::shared_ptr<CCurveObject>> DataPool::curve_object_pool_;
int DataPool::curve_object_max_id_;
std::map<int, std::shared_ptr<CMeshObject>>& DataPool::GetMeshObjectPool()
{
	return mesh_object_pool_;
}
CMeshObject * DataPool::GetMeshObject(int id)
{
	auto iter=mesh_object_pool_.find(id);
	if (iter == mesh_object_pool_.end())
	{
		return NULL;
	}
	else
	{
		return iter->second.get();
	}
}
int DataPool::AddMeshObject(CMeshObject *mesh)
{

	if (mesh->GetId() != -1 && (mesh_object_pool_.find(mesh->GetId()) == mesh_object_pool_.end()))
	{
		if (mesh_object_max_id_ < mesh->GetId())
		{
			mesh_object_max_id_ = mesh->GetId() + 1;
		}
		
		mesh_object_pool_.insert(std::pair<int, std::shared_ptr<CMeshObject>>(mesh->GetId(), std::shared_ptr<CMeshObject>(mesh)));
		return mesh->GetId();
	}
	else
	{
		mesh->SetId(mesh_object_max_id_);
		mesh_object_pool_.insert(std::pair<int, std::shared_ptr<CMeshObject>>(mesh_object_max_id_, std::shared_ptr<CMeshObject>(mesh)));
		mesh_object_max_id_++;
		return mesh_object_max_id_ - 1;
	}
}
void DataPool::Init()
{
	mesh_object_max_id_ = 0;
	curve_object_max_id_ = 0;
}

bool DataPool::DeleteCurveObject(int id)
{
	if (curve_object_pool_.find(id) != curve_object_pool_.end())
	{
		curve_object_pool_.erase(id);
		return true;
	}
	return false;
}
int DataPool::AddCurveObject(CCurveObject *curve)
{
	if (curve->GetId() != -1 && (curve_object_pool_.find(curve->GetId()) == curve_object_pool_.end()))
	{
		if (curve_object_max_id_ < curve->GetId())
		{
			curve_object_max_id_ = curve->GetId() + 1;
		}

		curve_object_pool_.insert(std::pair<int, std::shared_ptr<CCurveObject>>(curve->GetId(), std::shared_ptr<CCurveObject>(curve)));
		return curve->GetId();
	}
	else
	{
		curve->SetId(curve_object_max_id_);
		curve_object_pool_.insert(std::pair<int, std::shared_ptr<CCurveObject>>(curve_object_max_id_, std::shared_ptr<CCurveObject>(curve)));
		curve_object_max_id_++;
		return curve_object_max_id_ - 1;
	}
}
CCurveObject *DataPool::GetCurveObject(int id)
{
	auto iter = curve_object_pool_.find(id);
	if (iter == curve_object_pool_.end())
	{
		return NULL;
	}
	else
	{
		return iter->second.get();
	}
}
std::map<int, std::shared_ptr<CCurveObject>>& DataPool::GetCurveObjectPool()
{
	return curve_object_pool_;
}