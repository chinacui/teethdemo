#include"data_pool.h"

std::map<int, std::shared_ptr<CMeshObject>> DataPool::mesh_object_pool_;
int DataPool::mesh_object_max_id_;
std::map<int, std::shared_ptr<CMeshObject>>& DataPool::GetMeshObjectPool()
{
	return mesh_object_pool_;
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
}