#ifndef CPANORAMIC_SIMULATION_H
#define CPANORAMIC_SIMULATION_H


#include"prereq.h"
#include"../DataColle/mesh_object.h"
#include<opencv2\opencv.hpp>
#include<map>
#include<vector>
class TEETHROOTRECOALG_CLASS CPanoramicProjectorBase
{
public:
	virtual OpenMesh::Vec3d Project(OpenMesh::Vec3d p)=0;
	virtual void ProjectMesh(CMeshObject* meshobj, std::vector<OpenMesh::Vec3d> &res_pts, std::vector<COpenMeshT::VertexHandle > &res_orig_vhs= std::vector<COpenMeshT::VertexHandle >());
	virtual void ProjectMeshes(std::vector<CMeshObject*>mesh_objs, std::vector<std::vector<OpenMesh::Vec3d>>&res_pts, std::vector<std::vector<COpenMeshT::VertexHandle>>&res_orig_vhs = std::vector<std::vector<COpenMeshT::VertexHandle>>()) ;
	virtual OpenMesh::Vec2d ComputeProjParam( OpenMesh::Vec3d p)=0;
	void ComputeProjParams(std::vector<OpenMesh::Vec3d>&pts, std::vector<OpenMesh::Vec2d>&res_params);
	virtual void SetParams(std::vector<double>&params)=0;
	virtual void GetParams(std::vector<double>&params) = 0;
	
};
class TEETHROOTRECOALG_CLASS CCurveSurfaceProjector :public CPanoramicProjectorBase
{
protected:
	OpenMesh::Vec3d updir_;
	std::vector<OpenMesh::Vec3d>curve_;
	OpenMesh::Vec3d curve_center_;
	std::vector<OpenMesh::Vec3d>proj_curve_;
public:
	CCurveSurfaceProjector() {};
	void SetParams(OpenMesh::Vec3d updir, std::vector<OpenMesh::Vec3d>&curve);
	void SetParams(std::vector<double>&params);
	 void GetParams(std::vector<double>&params);
	void SetChanged();
	OpenMesh::Vec3d Project(OpenMesh::Vec3d p);
	OpenMesh::Vec2d ComputeProjParam(OpenMesh::Vec3d p);
	
	

};


/*CCircularSurfaceProjector -- get virtual x-ray image of crown using circular projection method*/
class TEETHROOTRECOALG_CLASS CCircularSurfaceProjector :public CPanoramicProjectorBase
{
protected:
	OpenMesh::Vec3d center_,updir_;
	OpenMesh::Vec3d start_dir_;
	double radius_=-1;
public:
	CCircularSurfaceProjector() {

	}
	/*
	*Summary: calculate the virtual x-ray projection center, up direction and radius parameters.
	*Parameters:--
	*returns:--
	*/
	OpenMesh::Vec3d GetCenter()
	{
		return center_;
	}
	OpenMesh::Vec3d GetUpDir()
	{
		return updir_;
	}
	double GetRadius()
	{
		return radius_;
	}
	void SetParams(std::vector<double>&params);
	void GetParams(std::vector<double>&params);
	void SetParams(OpenMesh::Vec3d updir, OpenMesh::Vec3d center, OpenMesh::Vec3d start_dir,  double radius);
	OpenMesh::Vec3d Project(OpenMesh::Vec3d p);
	OpenMesh::Vec2d ComputeProjParam(OpenMesh::Vec3d p);
	

};

class TEETHROOTRECOALG_CLASS CPanoramicSimulation
{
	
public:
	
	static void GeneratePanoramicImageByDentalArch(std::vector<CMeshObject *>crown_objs, std::vector<OpenMesh::Vec3d>&frame, int img_width, cv::Mat &res_panoramic);//to be finished
	static void GeneratePanoramicPointsByCircle(std::vector<CMeshObject *>crown_objs, OpenMesh::Vec3d updir,OpenMesh::Vec3d center, double radius,std::vector<std::vector<OpenMesh::Vec2d>>&res_pano_params,std::vector<std::vector<COpenMeshT::VertexHandle>>&res_corres_vhs= std::vector<std::vector<COpenMeshT::VertexHandle>>());
	static CCircularSurfaceProjector ConstructCircularProjector(std::vector<CMeshObject *>crown_objs, OpenMesh::Vec3d updir, OpenMesh::Vec3d center, double radius);
	static void GeneratePanoramicWholeBoundPointsByCircle(std::vector<CMeshObject *>crown_objs, OpenMesh::Vec3d updir, OpenMesh::Vec3d center, double radius, std::vector<std::vector<OpenMesh::Vec2d>>&res_bound_pano_params, std::vector<std::vector<COpenMeshT::VertexHandle>>&res_corres_vhs = std::vector<std::vector<COpenMeshT::VertexHandle>>());
	static void GeneratePanoramicBoundPointsByCircle(std::vector<CMeshObject *>crown_objs, OpenMesh::Vec3d updir, OpenMesh::Vec3d center, double radius, std::vector<std::vector<OpenMesh::Vec2d>>&res_bound_pano_params, std::vector<std::vector<COpenMeshT::VertexHandle>>&res_corres_vhs = std::vector<std::vector<COpenMeshT::VertexHandle>>());
};



#endif