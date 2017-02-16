#include"teeth_template_object.h"

CTeethTemplateObject::CTeethTemplateObject()
{
	updir_ = OpenMesh::Vec3d(0, 1, 0);
	front_dir_ = OpenMesh::Vec3d(0, 0, 1);
	crown_fhs_.clear();
	type_ = "";
	crown_.clear();
}
CTeethTemplateObject::CTeethTemplateObject(CTeethTemplateObject & template_object):CMeshObject(template_object)
{

	updir_ = template_object.updir_;
	front_dir_ = template_object.front_dir_;
	crown_fhs_.clear();

	for (int i = 0; i < template_object.crown_fhs_.size(); i++)
	{
		crown_fhs_.push_back(mesh_.face_handle(template_object.crown_fhs_[i].idx()));
	}
	ComputeCrownMesh();

}
void CTeethTemplateObject::SetChanged(bool flag )
{
	CMeshObject::SetChanged(flag);

}


COpenMeshT& CTeethTemplateObject::GetCrownMesh()
{
	if (crown_.n_vertices() == 0)
	{
		ComputeCrownMesh();
	}
	return crown_;
}
std::map<COpenMeshT::VertexHandle, COpenMeshT::VertexHandle> &CTeethTemplateObject::GetCrown2TempVhsMap()
{
	if (crown_.n_vertices() == 0)
	{
		ComputeCrownMesh();
	}
	return crown2tempvhs_map_;
}
void CTeethTemplateObject::ComputeCrownMesh()
{
	std::cerr << "compute crown mesh" << std::endl;
	crown2tempvhs_map_.clear();
	std::map<COpenMeshT::VertexHandle, COpenMeshT::VertexHandle>vhs_old2new_map;

	for (int i = 0; i<crown_fhs_.size(); i++)
	{
		COpenMeshT::FaceHandle fh = crown_fhs_[i];
		for (auto viter = mesh_.fv_begin(fh); viter != mesh_.fv_end(fh); viter++)
		{
			if (vhs_old2new_map.find(viter) == vhs_old2new_map.end())
			{
				vhs_old2new_map[viter] = viter;
			}

		}
	}
	crown_.clear();
	for (auto iter = vhs_old2new_map.begin(); iter != vhs_old2new_map.end(); iter++)
	{
		OpenMesh::Vec3d p = this->TransformPointByLocalMatrix(mesh_.point(iter->first));
		COpenMeshT::VertexHandle new_vh = crown_.add_vertex(p);
		iter->second= new_vh;

		crown2tempvhs_map_[new_vh] = iter->first;
	}
	for (int i = 0; i < crown_fhs_.size(); i++)
	{
		COpenMeshT::FaceHandle fh = crown_fhs_[i];
		std::vector<COpenMeshT::VertexHandle>vhs;
		for (auto viter = mesh_.fv_begin(fh); viter != mesh_.fv_end(fh); viter++)
		{
			vhs.push_back(vhs_old2new_map[viter]);
		}
		crown_.add_face(vhs);
	}

	
	
}