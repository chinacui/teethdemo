#include"cgal_arap_deform.h"
#include "../DataColle/Polyhedron_type.h"
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/properties_Polyhedron_3.h>
CgalArapDeform::CgalArapDeform(Polyhedron& mesh):mesh_(mesh)
{

	set_halfedgeds_items_id(mesh_);
	p_deform_mesh_ = new Surface_mesh_deformation(mesh_);
}

CgalArapDeform::~CgalArapDeform()
{
	delete(p_deform_mesh_);
}

bool CgalArapDeform::SetRoiVertices(const std::vector<Polyhedron::Vertex_handle>& vertices_handle_vec)
{
	
	Polyhedron::Vertex_iterator vb = mesh_.vertices_begin();
	Polyhedron::Vertex_iterator ve = mesh_.vertices_end();
	for (auto itr = vertices_handle_vec.begin(); itr != vertices_handle_vec.end(); itr++)
	{
		p_deform_mesh_->insert_roi_vertex(*itr);
	}
	return true;
}

bool CgalArapDeform::SetRoiVertices()
{

	for (auto itr = mesh_.vertices_begin(); itr != mesh_.vertices_end(); itr++)
	{
		p_deform_mesh_->insert_roi_vertex(itr);
	}
	return true;
}

bool CgalArapDeform::Deform()
{
	p_deform_mesh_->deform();
	return true;
}

bool CgalArapDeform::SetDeformMap(std::map<Polyhedron::Vertex_handle, Kernel::Point_3>& deform_handle_map)
{
	for (auto itr = deform_handle_map.begin(); itr != deform_handle_map.end(); itr++)
	{
		p_deform_mesh_->insert_control_vertex(itr->first);
	}
	p_deform_mesh_->preprocess();
	for (auto itr = deform_handle_map.begin(); itr != deform_handle_map.end(); itr++)
	{
		p_deform_mesh_->set_target_position(itr->first, itr->second);
	}
	return true;
}