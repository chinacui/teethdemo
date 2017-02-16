// Copyright 2016_9 by ChenNenglun
#include"data_pool.h"

std::map<int, std::shared_ptr<CMeshObject>> DataPool::mesh_object_pool_;
int DataPool::mesh_object_max_id_;

std::map<int, std::shared_ptr<CCurveObject>> DataPool::curve_object_pool_;
int DataPool::curve_object_max_id_;

std::map<int, std::shared_ptr<CVolumeDataObject>> DataPool::volume_data_object_pool_;
int DataPool::volume_data_object_max_id_;
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
bool DataPool::DeleteMeshObject(int id)
{
	if (mesh_object_pool_.find(id) != mesh_object_pool_.end())
	{
		mesh_object_pool_.erase(id);
		return true;
	}
	return false;
}
void DataPool::DeleteAllCurveObjects()
{
	curve_object_pool_.clear();
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




int DataPool::AddVolumeDataObject(CVolumeDataObject *volume_data_obj)
{
	if (volume_data_obj->GetId() != -1 && (volume_data_object_pool_.find(volume_data_obj->GetId()) == volume_data_object_pool_.end()))
	{
		if (volume_data_object_max_id_ < volume_data_obj->GetId())
		{
			volume_data_object_max_id_ = volume_data_obj->GetId() + 1;
		}

		volume_data_object_pool_.insert(std::pair<int, std::shared_ptr<CVolumeDataObject>>(volume_data_obj->GetId(), std::shared_ptr<CVolumeDataObject>(volume_data_obj)));
		return volume_data_obj->GetId();
	}
	else
	{
		volume_data_obj->SetId(volume_data_object_max_id_);
		volume_data_object_pool_.insert(std::pair<int, std::shared_ptr<CVolumeDataObject>>(volume_data_object_max_id_, std::shared_ptr<CVolumeDataObject>(volume_data_obj)));
		volume_data_object_max_id_++;
		return volume_data_object_max_id_ - 1;
	}
}
bool DataPool::DeleteVolumeDataObject(int id)
{
	if (volume_data_object_pool_.find(id) != volume_data_object_pool_.end())
	{
		volume_data_object_pool_.erase(id);
		return true;
	}
	return false;
}
void DataPool::DeleteAllVolumeDataObjects()
{
	volume_data_object_pool_.clear();
}
CVolumeDataObject *DataPool::GetVolumeDataObject(int id)
{
	auto iter = volume_data_object_pool_.find(id);
	if (iter == volume_data_object_pool_.end())
	{
		return NULL;
	}
	else
	{
		return iter->second.get();
	}
}
std::map<int, std::shared_ptr<CVolumeDataObject>>& DataPool::GetVolumeDataObjectPool()
{
	return volume_data_object_pool_;
}