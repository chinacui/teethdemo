// Copyright 2016_9 by ChenNenglun
#include"mesh_object.h"
#include<iostream>



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
void CMeshObject::SetChanged(bool is_changed)
{
	is_changed_ = is_changed;
	

	
	
}
void CMeshObject::CopyFrom(CMeshObject& b)
{
	mesh_ = b.mesh_;
	SetChanged();
}
CMeshObject::CMeshObject(CMeshObject &b)
{
	mesh_ = b.mesh_;
}
CMeshObject::CMeshObject()
{

	mesh_id_ = -1;
	is_changed_ = true;
}
CMeshObject::~CMeshObject()
{

}