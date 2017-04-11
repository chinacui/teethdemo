#include "cgal_igl_converter.h"
#include <CGAL/Polyhedron_incremental_builder_3.h>

bool CConverter::ConvertFromCGALToIGL(Polyhedron& mesh, Eigen::MatrixXd& v, Eigen::MatrixXi& f)
{
	set_halfedgeds_items_id(mesh);
	v.resize(mesh.size_of_vertices(), 3);
	for (auto vitr = mesh.vertices_begin(); vitr != mesh.vertices_end(); vitr++)
	{
		v(vitr->id(), 0) = vitr->point().x();
		v(vitr->id(), 1) = vitr->point().y();
		v(vitr->id(), 2) = vitr->point().z();
	}
	f.resize(mesh.size_of_facets(), 3);
	for (auto fitr = mesh.facets_begin(); fitr != mesh.facets_end(); fitr++)
	{
		assert(fitr->facet_degree() == 3);//triangle mesh
		auto vitr = fitr->facet_begin();
		for (size_t i = 0; i < 3; i++, vitr++)
		{
			f(fitr->id(), i) = vitr->vertex()->id();
		}
	}
	return true;
}

bool CConverter::ConvertFromOpenMeshToIGL(COpenMeshT &mesh, Eigen::MatrixXd& v, Eigen::MatrixXi& f, std::vector<COpenMeshT::FaceHandle>*id_fh_map,std::vector<COpenMeshT::VertexHandle>*id_vh_map )
{
	COpenMeshT::FaceIter fiter = mesh.faces_begin();
	f.resize(mesh.n_faces(), 3);
	v.resize(mesh.n_vertices(), 3);
	if (id_vh_map != NULL)
	{
		id_vh_map->resize(mesh.n_vertices());
	}
	for (int i = 0; i < mesh.n_vertices(); i++)
	{
		auto vpos=mesh.points()[i];
		for (int j = 0; j < 3; j++)
		{
			v(i, j) = vpos[j];
		}
		if (id_vh_map != NULL)
		{
			(*id_vh_map)[i] = mesh.vertex_handle(i);
		}
	}
	if (id_fh_map != NULL)
	{
		id_fh_map->resize(mesh.n_faces());
	}
	for (int i = 0; i < mesh.n_faces(); i++, fiter++)
	{
		COpenMeshT::FaceHalfedgeIter fhi = mesh.fh_begin(*fiter), fhend = mesh.fh_end(*fiter);
		if (id_fh_map != NULL)
		{
			(*id_fh_map)[i] = fiter;
		}
		int j = 0;
		while (fhi != fhend)
		{
			int vid = mesh.to_vertex_handle(*fhi).idx();
			
			f(i, j) = vid;
			j++;
			fhi++;
		}



		
	}
	return true;
}

bool CConverter::ConvertFromIGLToOpenMesh(Eigen::MatrixXd& v, Eigen::MatrixXi& f, COpenMeshT& openmesh)
{

	openmesh.clear();
	std::vector<COpenMeshT::VertexHandle>vh_map;
	vh_map.resize(v.rows());
	for (int i = 0; i < v.rows(); i++)
	{
		OpenMesh::Vec3d vertex(v(i, 0), v(i, 1), v(i, 2));
		vh_map[i] = openmesh.add_vertex(vertex);
	}
	for (int i = 0; i < f.rows(); i++)
	{
		std::vector<OpenMesh::VertexHandle> vhandles;
		for (int j = 0; j < f.cols(); j++)
		{
			vhandles.push_back(vh_map[f(i, j)]);
		}
		openmesh.add_face(vhandles);
	}
	return true;
}

bool CConverter::ConvertFromIGLToCGAL(Eigen::MatrixXd& v, Eigen::MatrixXi& f, Polyhedron &mesh, std::vector<Face_handle>&fid_fh_map, std::vector<Vertex_handle>&vid_vhs_map)
{
	typedef Polyhedron::HalfedgeDS             HalfedgeDS;
	typedef CGAL::Polyhedron_incremental_builder_3<HalfedgeDS> HDSIncreBuilder;
	HDSIncreBuilder poly_builder(mesh.hds(), true);
	poly_builder.begin_surface(v.rows(), f.rows());
	vid_vhs_map.resize(v.rows());
	fid_fh_map.resize(f.rows());
	for (int i = 0; i < v.rows(); i++)
	{
		Point p(v(i, 0), v(i, 1), v(i, 2));
		auto vh=poly_builder.add_vertex(p);
		vid_vhs_map[i] = vh;
	}
	for (int i = 0; i < f.rows(); i++)
	{
		fid_fh_map[i]=poly_builder.begin_facet();
		for (int j = 0; j < 3; j++)
		{
			poly_builder.add_vertex_to_facet(f(i, j));
		}
		poly_builder.end_facet();
	}
	poly_builder.end_surface();
	set_halfedgeds_items_id(mesh);
	return true;

}
bool CConverter::ConvertFromCGALRoiToIGL(Polyhedron& mesh, std::vector<Face_handle>&fhs, Eigen::MatrixXd& v, Eigen::MatrixXi&f)
{

	set_halfedgeds_items_id(mesh);
	std::vector<int>vertices_mask(mesh.size_of_vertices(), -1);
	int num_v = 0;
	int num_f = fhs.size();
	for (int i = 0; i < fhs.size(); i++)
	{
		Face_handle fh = fhs[i];
		int count = 0;
		for (auto hiter = fh->facet_begin(); count < fh->facet_degree(); count++, hiter++)
		{
			if (vertices_mask[hiter->vertex()->id()] == -1)
			{
				vertices_mask[hiter->vertex()->id()] = num_v;
				num_v++;
			}
		}
	}


	v.resize(num_v, 3);
	f.resize(num_f, 3);
	for (int i = 0; i < fhs.size(); i++)
	{
		Face_handle fh = fhs[i];
		int count = 0;
		for(auto hiter = fh->facet_begin(); count < fh->facet_degree(); count++, hiter++)
		{
			int old_vid = hiter->vertex()->id();
			if (vertices_mask[old_vid] == -1 || vertices_mask[old_vid] >= num_v)
				std::cerr << "vertices_mask[old_vid] == -1||vertices_mask[old_vid]>=num_v" << std::endl;
				
			f(i, count) = vertices_mask[old_vid];
			v(vertices_mask[old_vid], 0) = hiter->vertex()->point().x();
			v(vertices_mask[old_vid], 1) = hiter->vertex()->point().y();
			v(vertices_mask[old_vid], 2) = hiter->vertex()->point().z();
		}
	}
	
	return true;
}


bool CConverter::ConvertFromCGALToOpenMesh(Polyhedron& mesh, COpenMeshT& openmesh, std::map<Polyhedron::Vertex_handle, COpenMeshT::VertexHandle>& vh_map, std::map<Polyhedron::Facet_handle, COpenMeshT::FaceHandle>& fh_map)
{
	openmesh.clear();
	vh_map.clear();
	fh_map.clear();
	openmesh.request_halfedge_texcoords2D();
	openmesh.request_vertex_texcoords2D();

	for (Polyhedron::Vertex_iterator viter = mesh.vertices_begin(); viter != mesh.vertices_end(); viter++)
	{
		COMTraits::Point p(viter->point()[0], viter->point()[1], viter->point()[2]);
		COpenMeshT::VertexHandle vhandle = openmesh.add_vertex(p);
		vh_map[viter] = vhandle;
	}
	for (Polyhedron::Facet_iterator fiter = mesh.facets_begin(); fiter != mesh.facets_end(); fiter++)
	{
		std::vector<OpenMesh::VertexHandle> vhandles;
		std::map<OpenMesh::VertexHandle, COMTraits::TexCoord2D> uvs;
		std::map<OpenMesh::VertexHandle, COMTraits::TexCoord2D> uvspares;
		Polyhedron::Halfedge_around_facet_circulator heiter(fiter->facet_begin()), heend(heiter);
		CGAL_For_all(heiter, heend)
		{
			vhandles.push_back(vh_map[heiter->vertex()]);
			uvs[vhandles.back()] = COMTraits::TexCoord2D(heiter->uv_.x(), heiter->uv_.y());
			uvspares[vhandles.back()] = COMTraits::TexCoord2D(heiter->uv_spare_.x(), heiter->uv_spare_.y());
		}
		OpenMesh::FaceHandle om_fhandle = openmesh.add_face(vhandles);
		fh_map[fiter] = om_fhandle;
		COpenMeshT::FaceHalfedgeIter fhi = openmesh.fh_begin(om_fhandle), fhend = openmesh.fh_end(om_fhandle);
		while (fhi != fhend)
		{
		
			openmesh.data(*fhi).SetUV(uvs[openmesh.to_vertex_handle(*fhi)]);
			fhi++;
		}
	}
	return true;
}

bool CConverter::ConvertFromOpenMeshToCGAL(COpenMeshT& openmesh, Polyhedron &mesh, std::map<COpenMeshT::VertexHandle, Polyhedron::Vertex_handle>& vh_map, std::map<COpenMeshT::FaceHandle, Polyhedron::Facet_handle>& fh_map)
{
	mesh.clear();
	vh_map.clear();
	fh_map.clear();
	CGAL::Polyhedron_incremental_builder_3<Polyhedron::HalfedgeDS> meshbuilder(mesh.hds(), false);
	meshbuilder.begin_surface(openmesh.n_vertices(), openmesh.n_faces());
	for (int i = 0; i < openmesh.n_vertices(); i++)
	{
		Kernel::Point_3 p(openmesh.points()[i][0], openmesh.points()[i][1], openmesh.points()[i][2]);
		vh_map[openmesh.vertex_handle(i)] = meshbuilder.add_vertex(p);
	}
	COpenMeshT::FaceIter fiter = openmesh.faces_begin();
	for (int i = 0; i < openmesh.n_faces(); i++, fiter++)
	{
		COpenMeshT::FaceHalfedgeIter fhi = openmesh.fh_begin(*fiter), fhend = openmesh.fh_end(*fiter);
		std::map<Polyhedron::Vertex_handle, Kernel::Point_2> vertexuv_map;
		std::map<Polyhedron::Vertex_handle, Kernel::Point_2> vertexuvspare_map;
		auto poly_fhandle = meshbuilder.begin_facet();
		while (fhi != fhend)
		{
			meshbuilder.add_vertex_to_facet(openmesh.to_vertex_handle(*fhi).idx());
			COMTraits::TexCoord2D om_uv = openmesh.data(*fhi).GetUV();
			COMTraits::TexCoord2D om_uvspare;
			vertexuv_map[vh_map[openmesh.to_vertex_handle(*fhi)]] = Kernel::Point_2(om_uv[0], om_uv[1]);
			vertexuvspare_map[vh_map[openmesh.to_vertex_handle(*fhi)]] = Kernel::Point_2(om_uvspare[0], om_uvspare[1]);
			fhi++;
		}
		meshbuilder.end_facet();
		fh_map[*fiter] = poly_fhandle;
		Polyhedron::Halfedge_around_facet_circulator heiter(poly_fhandle->facet_begin()), heend(heiter);
		CGAL_For_all(heiter, heend)
		{
			heiter->uv_ = vertexuv_map[heiter->vertex()];
			heiter->uv_spare_ = vertexuvspare_map[heiter->vertex()];
		}
	}
	meshbuilder.end_surface();
	set_halfedgeds_items_id(mesh);
	return true;
}

