#include "geodesic_type.h"
void CGeodesicModel::Update()
{
	/*this->m_Verts.resize(mesh_->n_vertices());
	this->m_Faces.resize(mesh_->n_faces());
	for (auto viter = mesh_->vertices_begin(); viter != mesh_->vertices_end(); viter++)
	{
		m_Verts[viter->idx()] = CPoint3D(mesh_->point(viter)[0], mesh_->point(viter)[1], mesh_->point(viter)[2]);
	}
	for (auto fiter = mesh_->faces_begin(); fiter != mesh_->faces_end(); fiter++)
	{
		std::vector<int> fs;
		fs.clear();
		COpenMeshT::FaceHalfedgeIter fhi = mesh_->fh_begin(*fiter), fhend = mesh_->fh_end(*fiter);

		while (fhi != fhend)
		{
			int vid = mesh_->to_vertex_handle(*fhi).idx();
			fs.push_back(vid);

			fhi++;
		}
		if (fs.size() != 3)
		{
			std::cerr << "error : num of face vertexs must be 3" << std::endl;
			return;
		}
		m_Faces[fiter->idx()] = CBaseModel::CFace(fs[0], fs[1], fs[2]);
	}*/

	/////
	
}
CGeodesicModel::CGeodesicModel(COpenMeshT *mesh)
{
	mesh_ = mesh;
	
}