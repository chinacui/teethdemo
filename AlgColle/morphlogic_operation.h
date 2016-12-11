#ifndef CMORPHLOGIC_OPERATION_H
#define CMORPHLOGIC_OPERATION_H
#include"prereq.h"
#include<Eigen/Dense>
#include<iostream>
#include<vector>
#include"../DataColle/mesh_object.h"
#include"prereq.h"
class ALGCOLLE_CLASS CMorphlogicOperation
{
public:
	enum CMOVertexTag{Complex,Disk,Center,NonFeature};
	static bool IsComplexVertex(COpenMeshT & mesh, std::vector<bool>&label, COpenMeshT::VertexHandle vh);
	static bool IsCenterVertex(COpenMeshT & mesh, std::vector<bool>&label, COpenMeshT::VertexHandle vh);
	static int ComputeDegree(COpenMeshT &mesh, std::vector<bool>&label,COpenMeshT::VertexHandle vh);
	static void TagAllVertexs(COpenMeshT &mesh, std::vector<bool>&labels, std::vector<CMOVertexTag>&tags, std::vector<int>&region_ids);//tags:complex,disk,center,nonfeature  region_ids:id of non feature(or the boundary of the feature region) region start from 0, -1 mean inner feature region
public:
	static void Dilate(COpenMeshT &mesh, std::vector<bool>&labels);
	static void Erode(COpenMeshT &mesh, std::vector<bool>&labels);

	static void FloodFill(COpenMeshT &mesh, std::vector<bool>&label, bool roi_label, COpenMeshT::VertexHandle start_vh, int target_tag, std::vector<int>&res_region_tag);
	static void Skeletonization(COpenMeshT& mesh, std::vector<bool>&labels);
};
#endif