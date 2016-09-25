// Copyright 2016_9 by ChenNenglun
#ifndef DATA_IO_H
#define DATA_IO_H
#include"prereq.h"
#include<iostream>
#include<string>
#include"mesh_object.h"
#include"rply.h"
#include"rplyfile.h"
class DATACOLLE_CLASS CDataIO
{
protected:
	static int ReadPlyVertexCallback(p_ply_argument argument);
	static int ReadPlyRgbCallback(p_ply_argument argument);
	static int ReadPlyFaceCallback(p_ply_argument argument);
public:
	static bool ReadPly(std::string fname, CMeshObject & res_mesh_obj);
	
};
#endif
