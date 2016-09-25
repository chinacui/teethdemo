// Copyright 2016_9 by ChenNenglun
#include"data_io.h"
#include<vector>
int CDataIO::ReadPlyVertexCallback(p_ply_argument argument)
{
	void *ptemp = NULL;				// pointer to the PlyToVoxel object, store the data
	std::vector<Eigen::Vector3d> *pdata;	// casted pointer; same as temp_p
	long int index;

	/* pass the value from the input file to the output file */
	ply_get_argument_user_data(argument, &ptemp, &index);
	pdata = static_cast< std::vector<Eigen::Vector3d>* > (ptemp);

	// save the data to the list
	if (index == 0)
	{
		pdata->push_back(Eigen::Vector3d());
	}

	int vidx = pdata->size() - 1;
	(*pdata)[vidx][index] = ply_get_argument_value(argument);


	return 1;
}
int CDataIO::ReadPlyRgbCallback(p_ply_argument argument)
{
	void *ptemp = NULL;				// pointer to the PlyToVoxel object, store the data
	std::vector<Eigen::Vector3i> *pdata;	// casted pointer; same as temp_p
	long int index;

	/* pass the value from the input file to the output file */
	ply_get_argument_user_data(argument, &ptemp, &index);
	pdata = static_cast< std::vector<Eigen::Vector3i>* > (ptemp);

	// save the data to the list
	if (index == 0)
	{
		pdata->push_back(Eigen::Vector3i());
	}

	int cidx = pdata->size() - 1;
	(*pdata)[cidx][index] = ply_get_argument_value(argument);


	return 1;
}

int CDataIO::ReadPlyFaceCallback(p_ply_argument argument)
{
	void *ptemp = NULL;					// pointer to the PlyToVoxel object, store the data
	std::vector<Eigen::Vector3i> *pdata;		// casted pointer; same as temp_p
	long int idata;

	/* pass the value from the input file to the output file */
	ply_get_argument_user_data(argument, &ptemp, &idata);
	pdata = static_cast< std::vector<Eigen::Vector3i>* > (ptemp);

	// save the data to the list
	long length, value_index;
	ply_get_argument_property(argument, NULL, &length, &value_index);

	if (value_index == -1)
	{
		pdata->push_back(Eigen::Vector3i());
	}
	else
	{
		int tidx = pdata->size() - 1;
		(*pdata)[tidx][value_index] = ply_get_argument_value(argument);
	}


	return 1;
}

bool CDataIO::ReadPly(std::string fname, CMeshObject & res_mesh_obj)
{
	p_ply ply = ply_open(fname.c_str(), NULL, 0, NULL);
	if (!ply)
	{
		return false;
	}
	if (!ply_read_header(ply))
	{
		return false;
	}

	long nvertices, ntriangles;
	// vertex
	std::vector<Eigen::Vector3d> vlist;
	nvertices = ply_set_read_cb(ply, "vertex", "x", ReadPlyVertexCallback, &vlist, 0);
	ply_set_read_cb(ply, "vertex", "y", ReadPlyVertexCallback, &vlist, 1);
	ply_set_read_cb(ply, "vertex", "z", ReadPlyVertexCallback, &vlist, 2);
	
	// rgb color
	std::vector<Eigen::Vector3i> vcolorlist;
	ply_set_read_cb(ply, "vertex", "red", ReadPlyRgbCallback, &vcolorlist, 0);
	ply_set_read_cb(ply, "vertex", "green", ReadPlyRgbCallback, &vcolorlist, 1);
	ply_set_read_cb(ply, "vertex", "blue", ReadPlyRgbCallback, &vcolorlist, 2);
	

	// face
	std::vector<Eigen::Vector3i>facelist;
	ntriangles = ply_set_read_cb(ply, "face", "vertex_indices", ReadPlyFaceCallback, &facelist, 0);
	
	


	// read mesh info
	if (!ply_read(ply))
	{
		return false;
	}
	ply_close(ply);
	res_mesh_obj.vertexs_.resize(vlist.size(), 3);
	for (int i = 0; i < vlist.size(); i++)
	{
		res_mesh_obj.vertexs_(i, 0) = vlist[i](0);
		res_mesh_obj.vertexs_(i, 1) = vlist[i](1);
		res_mesh_obj.vertexs_(i, 2) = vlist[i](2);
	}
	res_mesh_obj.vertex_colors_.resize(vcolorlist.size(), 3);
	for (int i = 0; i < vcolorlist.size(); i++)
	{
		res_mesh_obj.vertex_colors_(i, 0) = vcolorlist[i](0)/255.0;
		res_mesh_obj.vertex_colors_(i, 1) = vcolorlist[i](1)/255.0;
		res_mesh_obj.vertex_colors_(i, 2) = vcolorlist[i](2)/255.0;
	}
	res_mesh_obj.faces_.resize(facelist.size(), 3);
	for (int i = 0; i < facelist.size(); i++)
	{
		res_mesh_obj.faces_(i, 0) = facelist[i][0];
		res_mesh_obj.faces_(i, 1) = facelist[i][1];
		res_mesh_obj.faces_(i, 2) = facelist[i][2];
	}

	return true;
}


