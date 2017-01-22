// Copyright 2016_9 by ChenNenglun
#ifndef DATA_IO_H
#define DATA_IO_H
#include"prereq.h"
#include<iostream>
#include<string>
#include"mesh_object.h"
#include"volume_data_object.h"
class DATACOLLE_CLASS CDataIO
{
protected:

public:
	
	static bool LoadCurveFromObj(std::string fname, std::vector<OpenMesh::Vec3d> &curve);
	static bool ReadMesh(std::string fname, CMeshObject & res_mesh_obj, OpenMesh::IO::Options io_options= OpenMesh::IO::Options::Default);
	static bool WriteMesh(std::string fname, CMeshObject & res_mesh_obj);
	static ItkVolumeDataType::Pointer ReadVolumeDataFromDICOMSeries(std::string dirname);//return null if failed
	static bool ReadVolumeDataObjFromDICOMSeries(std::string dirname, CVolumeDataObject& volume_data_obj);
};
#endif
