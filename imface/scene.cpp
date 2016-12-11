// Copyright 2016_9 by ChenNenglun
#include"scene.h"
#include<fstream>

CScene::CScene()
{
	scene_mesh_.clear();
	scene_curve_.clear();

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
	for (auto iter = scene_curve_.begin(); iter != scene_curve_.end(); iter++)
	{
		iter->second->Render(camera);
	}

	for (auto iter = scene_volumedata_.begin(); iter != scene_volumedata_.end(); iter++)
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
	std::vector<int>mto_be_removed;
	for (auto iter = scene_mesh_.begin(); iter != scene_mesh_.end(); ++iter)
	{
		//need to be removed from scene
		if (mesh_obj_pool.find(iter->first) == mesh_obj_pool.end())
		{
			mto_be_removed.push_back(iter->first);
		}
	}
	for (int i = 0; i < mto_be_removed.size(); i++)
		scene_mesh_.erase(mto_be_removed[i]);


	auto curve_obj_pool = DataPool::GetCurveObjectPool();
	for (auto iter = curve_obj_pool.begin(); iter != curve_obj_pool.end(); iter++)
	{
		//find new mesh
		if (scene_curve_.find(iter->first) == scene_curve_.end())
		{
			CSceneCurveObject *p_scene_curve_obj = new CSceneCurveObject((std::weak_ptr<CCurveObject>)iter->second);
			scene_curve_[iter->first] = p_scene_curve_obj;
		}
	}
	std::vector<int>cto_be_removed;
	for (auto iter = scene_curve_.begin(); iter != scene_curve_.end(); ++iter)
	{
		//need to be removed from scene
		if (curve_obj_pool.find(iter->first) == curve_obj_pool.end())
		{
			cto_be_removed.push_back(iter->first);
		}
	}
	for (int i = 0; i < cto_be_removed.size(); i++)
	{
		scene_curve_.erase(cto_be_removed[i]);
	}



	auto volume_data_obj_pool = DataPool::GetVolumeDataObjectPool();
	for (auto iter = volume_data_obj_pool.begin(); iter != volume_data_obj_pool.end(); iter++)
	{
		
		if (scene_volumedata_.find(iter->first) == scene_volumedata_.end())
		{
			CSceneVolumeObject *p_scene_volumedata_obj = new CSceneVolumeObject((std::weak_ptr<CVolumeDataObject>)iter->second);
			scene_volumedata_[iter->first] = p_scene_volumedata_obj;
		}
	}
	std::vector<int>vto_be_removed;
	for (auto iter = scene_volumedata_.begin(); iter != scene_volumedata_.end(); ++iter)
	{
		//need to be removed from scene
		if (scene_volumedata_.find(iter->first) == scene_volumedata_.end())
		{
			vto_be_removed.push_back(iter->first);
		}
	}
	for (int i = 0; i < vto_be_removed.size(); i++)
	{
		scene_volumedata_.erase(vto_be_removed[i]);
	}
}