#include"mesh_object.h"
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
CMeshObject::CMeshObject(CMeshObject &b)
{
	vertexs_ = b.vertexs_;
	faces_ = b.faces_;
	vertex_colors_ = b.vertex_colors_;
}
CMeshObject::CMeshObject()
{
	mesh_id_ = -1;
	is_changed_ = true;
}
CMeshObject::~CMeshObject()
{

}