// Copyright 2016_9 by ChenNenglun
#include"mesh_object.h"
#include<iostream>
#include"Polyhedron_type.h"
#include"aabb_type.h"
int CMeshObject::GetId()
{
	return mesh_id_;
}
void CMeshObject::SetId(int id)
{
	mesh_id_ = id;
}
bool CMeshObject::IsChanged()
{
	return is_changed_;
}

void CMeshObject::SetAttrChanged(bool is_attrchanged)
{
	is_changed_ = is_attrchanged;

}
void CMeshObject::SetChanged(bool is_changed)
{
	is_changed_ = is_changed;
	if (is_changed)
	{
		if (aabb_tree_ != NULL)
		{
			delete aabb_tree_;
			aabb_tree_ = NULL;
		}
		if (geodesic_model_ != NULL)
		{
			delete geodesic_model_;
			geodesic_model_ = NULL;
		}
		tot_surface_area_ = -1;
		mesh_.update_normals();
	}
	
	if (vtag_.size() != mesh_.n_vertices())
	{
		vtag_.clear();

		for (auto viter = mesh_.vertices_begin(); viter != mesh_.vertices_end(); viter++)
		{
			vtag_[viter->idx()] = -1;
		}

	}
	
	
}
void CMeshObject::CopyFrom(CMeshObject& b)
{
	mesh_ = b.mesh_;
	is_visiable_ = b.is_visiable_;
	vtag_ = b.vtag_;
	is_shinning_ = b.is_shinning_;
	SetChanged();
}
CMeshObject::CMeshObject(CMeshObject &b)
{
	mesh_ = b.mesh_;
	texture_id_ = b.TextureId();
	is_visiable_ = b.is_visiable_;
	use_texture_ = false;
	vtag_ = b.vtag_;
	is_shinning_ = b.is_shinning_;
	SetChanged();
}
CMeshObject::CMeshObject():CBaseObject()
{

	mesh_id_ = -1;
	texture_id_ = -1;
	use_texture_ = false;
	is_changed_ = true;
	tot_surface_area_ = -1;
	is_visiable_ = true;
	vtag_.clear();
	is_shinning_ = false;
	for (auto viter = mesh_.vertices_begin(); viter != mesh_.vertices_end(); viter++)
	{
		vtag_[viter->idx()] = -1;
	}
}
void CMeshObject::SetMeshColor(OpenMesh::Vec3d color)
{
	for (auto viter = mesh_.vertices_begin(); viter != mesh_.vertices_end(); viter++)
	{
		mesh_.set_color(viter, color);
	}
}
void CMeshObject::NormalizeUV()
{
	OpenMesh::Vec2f mmin(std::numeric_limits<float>::max(), std::numeric_limits<float>::max());
	OpenMesh::Vec2f mmax(std::numeric_limits<float>::min(), std::numeric_limits<float>::min());
	for (auto hiter = mesh_.halfedges_begin(); hiter != mesh_.halfedges_end(); hiter++)
	{
		OpenMesh::Vec2f uv=mesh_.data(hiter).GetUV();
		for (int i = 0; i < 2; i++)
		{
			if (mmin[i] > uv[i])
				mmin[i] = uv[i];
			if (mmax[i] < uv[i])
				mmax[i] = uv[i];
		}
	}
	float len0 = mmax[0] - mmin[0];
	float len1 = mmax[1] - mmin[1];
	for (auto hiter = mesh_.halfedges_begin(); hiter != mesh_.halfedges_end(); hiter++)
	{
		OpenMesh::Vec2f uv = mesh_.data(hiter).GetUV();
		if (len0 != 0)
			uv[0] = (uv[0] - mmin[0]) / len0;
		if (len1 != 0)
			uv[1] = (uv[1] - mmin[1]) / len1;
		mesh_.data(hiter).SetUV(uv);
	}
}
double CMeshObject::GetTotSurfaceArea()
{
	if (tot_surface_area_ == -1)
	{
		tot_surface_area_ = 0;
		for (auto fiter = mesh_.faces_begin(); fiter != mesh_.faces_end(); fiter++)
		{
			std::vector<OpenMesh::Vec3d>fvs;
			for (auto fviter = mesh_.fv_begin(fiter); fviter != mesh_.fv_end(fiter); fviter++)
			{
				fvs.push_back(mesh_.point(fviter));
			}
			OpenMesh::Vec3d va = fvs[1] - fvs[0];
			OpenMesh::Vec3d vb = fvs[2] - fvs[0];
			tot_surface_area_+=OpenMesh::cross(va, vb).length()/2;
		}
	}
	return tot_surface_area_;
}
CMeshObject::~CMeshObject()
{
	if (aabb_tree_ != NULL)
		delete aabb_tree_;
	if (geodesic_model_ != NULL)
		delete geodesic_model_;
}