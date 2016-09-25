#include"scene.h"
#include<fstream>
CScene::CScene()
{
	scene_mesh_.clear();
	

}
CScene::~CScene()
{

}

void CScene::Render( CCamera camera)
{
	UpdateScene();
	for (auto iter = scene_mesh_.begin(); iter != scene_mesh_.end(); iter++)
	{
		iter->second->Render(camera);
	}
}
void CScene::UpdateScene()
{
	auto mesh_obj_pool=DataPool::GetMeshObjectPool();
	for (auto iter = mesh_obj_pool.begin(); iter != mesh_obj_pool.end(); iter++)
	{
		//find new mesh
		if (scene_mesh_.find(iter->first) == scene_mesh_.end())
		{
			CSceneMeshObject *p_scene_mesh_obj=new CSceneMeshObject((std::weak_ptr<CMeshObject>)iter->second);
			scene_mesh_[iter->first] = p_scene_mesh_obj;
		}
	}
	for (auto iter = scene_mesh_.begin(); iter != scene_mesh_.end(); iter++)
	{
		//need to be removed from scene
		if (mesh_obj_pool.find(iter->first) == mesh_obj_pool.end())
		{
			scene_mesh_.erase(iter->first);
		}
	}
}