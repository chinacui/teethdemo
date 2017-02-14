#ifndef CGAL_ARAP_DEFORM_H
#define CGAL_ARAP_DEFORM_H
#include"prereq.h"
#include <CGAL/Surface_mesh_deformation.h>
#include "../DataColle/Polyhedron_type.h"
class CgalArapDeform
{
public:
	//CARAPDeformation() : CMeshDeformation() {}
	CgalArapDeform(Polyhedron& mesh);
	virtual ~CgalArapDeform();
	bool SetRoiVertices();//Set all of vertices on polyhedron as ROI
	bool SetRoiVertices(const std::vector<Polyhedron::Vertex_handle>& vertices_handle_vec);//Set specified vertices as ROI
	bool SetDeformMap(std::map<Polyhedron::Vertex_handle, Kernel::Point_3>& deform_handle_map);


	bool Deform();
private:
	typedef CGAL::Surface_mesh_deformation<Polyhedron> Surface_mesh_deformation;
	Surface_mesh_deformation* p_deform_mesh_;
	Polyhedron& mesh_;
};

#endif
