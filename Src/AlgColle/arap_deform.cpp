#include"arap_deform.h"

#include "cgal_arap_deform.h"
#include<map>
//CARAPDeformation::CARAPDeformation(COpenMeshT& mesh):mesh_(mesh)
//{
//	
//	std::vector<COpenMeshT::FaceHandle>id_fh_map;
//	CConverter::ConvertFromOpenMeshToIGL(mesh, v_, f_);
//	arap_data_.energy = igl::ARAP_ENERGY_TYPE_SPOKES;
//	arap_data_.max_iter = 100;
//
//	
//}
//CARAPDeformation::~CARAPDeformation()
//{
//
//}
////bool CARAPDeformation::SetRoiVertices();//Set all of vertices on polyhedron as ROI
////{
////
////}
////bool CARAPDeformation::SetRoiVertices(const std::vector<COpenMeshT::VertexHandle>& vertices_handle_vec);//Set specified vertices as ROI
////{
////
////}
//bool CARAPDeformation::SetDeformHandle(std::vector<COpenMeshT::VertexHandle>& deform_handles)
//{
//	deform_handles_ = deform_handles;
//	handle_id_.resize(deform_handles.size());
//	for (int i = 0; i < deform_handles.size(); i++)
//	{
//		handle_id_(i)=deform_handles[i].idx();
//	}
//	igl::arap_precomputation(v_, f_, v_.cols(), handle_id_, arap_data_);
//	handle_targt_.resize(handle_id_.rows(),v_.cols());
//	return true;
//}
//bool CARAPDeformation::Deform(std::map<COpenMeshT::VertexHandle, OpenMesh::Vec3d>&handle_map)
//{
//	std::cerr << "start arap deform" << std::endl;
//	for(int i=0;i<handle_id_.rows();i++)
//	{
//		int id = handle_id_(i);
//		COpenMeshT::VertexHandle vh= mesh_.vertex_handle(id);
//		OpenMesh::Vec3d p=handle_map[vh];
//		handle_targt_(i, 0) = p[0];
//		handle_targt_(i, 1) = p[1];
//		handle_targt_(i, 2) = p[2];
//		
//	}
//	
//	if (!igl::arap_solve(handle_targt_, arap_data_, v_))
//	{
//		std::cerr << "arap solve failed" << std::endl;
//	}
//	else
//	{
//		std::cerr << "arap solve" << std::endl;
//	}
//	for (auto viter = mesh_.vertices_begin(); viter != mesh_.vertices_end(); viter++)
//	{
//		mesh_.point(viter)=OpenMesh::Vec3d(v_(viter->idx(), 0), v_(viter->idx(), 1), v_(viter->idx(), 2));
//	}
//	std::cerr << "end arap" << std::endl;
//	return true;
//}


CARAPDeformation::CARAPDeformation(COpenMeshT& mesh):mesh_(mesh)
{
	std::map< COpenMeshT::VertexHandle, Polyhedron::Vertex_handle> vh_map;
	std::map<COpenMeshT::FaceHandle, Polyhedron::Facet_handle> fh_map;

	CConverter::ConvertFromOpenMeshToCGAL(mesh, poly_, vh_map, fh_map);
	set_halfedgeds_items_id(poly_);
	cgal_id_to_om_vh_map_.resize(mesh.n_vertices());
	om_id_to_cgal_vh_map_.resize(mesh.n_vertices());
	for(auto viter=mesh.vertices_begin();viter!=mesh.vertices_end();viter++)
	{
		cgal_id_to_om_vh_map_[vh_map[viter]->id()] = viter;
		om_id_to_cgal_vh_map_[viter->idx()] = vh_map[viter];
	}
	cgal_arap_deform_ = new CgalArapDeform(poly_);
	cgal_arap_deform_->SetRoiVertices();
}
CARAPDeformation::~CARAPDeformation()
{
	if (cgal_arap_deform_ == NULL)
		delete cgal_arap_deform_;
}

bool CARAPDeformation::SetDeformMap(std::map<COpenMeshT::VertexHandle, OpenMesh::Vec3d>& deform_handle_map)
{
	std::map<Polyhedron::Vertex_handle, Kernel::Point_3>cgal_deform_handle_map;
	for (auto iter = deform_handle_map.begin(); iter != deform_handle_map.end(); iter++)
	{
		cgal_deform_handle_map[om_id_to_cgal_vh_map_[iter->first.idx()]] =Point(iter->second[0], iter->second[1], iter->second[2]);
	}
	
	cgal_arap_deform_->SetDeformMap(cgal_deform_handle_map);
	return true;
}


bool CARAPDeformation::Deform()
{
	cgal_arap_deform_->Deform();
	for (auto viter = poly_.vertices_begin(); viter != poly_.vertices_end(); viter++)
	{
		Point p = viter->point();
		OpenMesh::Vec3d om_p(p[0], p[1], p[2]);
		COpenMeshT::VertexHandle vh=cgal_id_to_om_vh_map_[viter->id()];
		mesh_.point(vh) = om_p;
	}
	return true;
}




