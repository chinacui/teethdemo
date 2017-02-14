#ifndef CPANORAMIC_ALG_H
#define CPANORAMIC_ALG_H
#include"../DataColle/custom_openmesh_type.h"
#include"../DataColle/mesh_object.h"
#include"../DataColle/data_pool.h"
#include<opencv2\opencv.hpp>
class CPanoramicAlg
{
public:
	static bool ComputePanoramicImageOfToothMesh(COpenMeshT&tooth_mesh, cv::Mat &res_img);
};
#endif