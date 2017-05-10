#pragma once
#ifndef TEETH_PROJECTION
#define TEETH_PROJECTION

#include <iostream>
#include <Eigen\Dense>
#include <vector>
#include <map>
#include "..\..\Src\DataColle\mesh_object.h"
#include "..\..\Src\TeethRootRecoAlg\prereq.h"

class  TEETHROOTRECOALG_CLASS teethProjection
{

public:
	teethProjection();
	static teethProjection* GetInstance();
	static void DeleteInstance();
	~teethProjection();
public:
	void test();
	double computeCurveValue(double y);
	double computeDerivateCurve(double y);
	double computeNeareastPointOfCurve(double y, OpenMesh::Vec3d & tooth_point);
	double computeLenghOfCurve(double min_i_y, double max_i_y, double y);
	bool setTeethMesh(std::map<int, COpenMeshT*>&sep_meshes);
	bool setCurveCeoffs(std::vector<double> &ceoffs);
	bool curveDiscretization();
	bool establishMeshCurveMap(double min_i_y_, double max_i_y_);
	bool computeProjectionImage(double min_i_y, double max_i_y);
	bool computeProjectionImageContour();
	bool calculateTeethTop(std::map<int, std::map<int, std::map<int, double>>>& map_toothId_map_meshId_tipLen_);
	bool sortTeethTop();
	std::map<int, std::vector<OpenMesh::Vec2d>> getProjectionImage();
	std::map<int, OpenMesh::Vec3d> getTest();
	std::vector<std::map<int, int>> gettest1();
private:
	static teethProjection* instance_;
	//The mesh of each tooth
	std::map<int, COpenMeshT*>sep_meshes_;
	//The ceofficient of the curve function
	std::vector<double> curve_ceoffs_;
	//the map that connect each tooth mesh with the curve
	std::map<int, std::vector<OpenMesh::Vec2d>> mesh_curve_map_;
	//the projection image
	std::map<int, std::vector<OpenMesh::Vec2d>> projection_teeth_image_;
	//the contour of projection image
	std::map<int, std::vector<OpenMesh::Vec2d>> projection_teeth_image_contour_;
	//the map between the 3D vertices and 2D virsual point
	std::vector<std::map<int, int>> map_3dvertice_2dpoints_;
	std::map<int, std::vector<int>> judge_image_contour_;
	std::vector<std::vector<int>> map_vids_;
	std::map<int, std::map<int, double>> map_toothId_map_meshId_tipLen_;
	std::map<int, std::map<int, double>> map_toothId_map_verticeHandle_tipLen_;
	std::map<int,OpenMesh::Vec3d> map_toothId_TopPoint_;
};

#endif

