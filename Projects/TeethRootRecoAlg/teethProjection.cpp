#include "teethProjection.h"
#include "..\..\Src\AlgColle\geo_base_alg.h"
#include "..\..\Src\AlgColle\geo_alg.h"
#include <math.h>

teethProjection* teethProjection::instance_ = nullptr;

teethProjection::teethProjection()
{
}


teethProjection::~teethProjection()
{
}

teethProjection* teethProjection::GetInstance() {
	if (nullptr == instance_) {
		instance_ = new teethProjection();
	}
	return instance_;
}

void teethProjection::DeleteInstance() {
	if (nullptr != instance_) {
		delete instance_;
		instance_ = nullptr;
	}
}

bool teethProjection::setTeethMesh(std::map<int, COpenMeshT*>&sep_meshes)
{
	this->sep_meshes_ = sep_meshes;
	return true;
}

bool teethProjection::setCurveCeoffs(std::vector<double> &ceoffs)
{
	this->curve_ceoffs_ = ceoffs;
	return true;
}

bool teethProjection::curveDiscretization()
{
	if (this->sep_meshes_.size() == NULL) return false;
	double min_y_ = 1e10; //the min value of y of the curve
	double max_y_ = -1e10;//the max ...
	std::vector<double> min_i_y_; //the min value of y of each tooth
	std::vector<double> max_i_y_; //the max ...
	for (auto iter = this->sep_meshes_.begin(); iter != this->sep_meshes_.end(); iter++)
	{
		double min_i_y_temp_ = 1e10;
		double max_i_y_temp_ = -1e10;
		COpenMeshT mesh_i = *iter->second;
		for (auto viter = mesh_i.vertices_begin(); viter != mesh_i.vertices_end(); viter++)
		{
			if (mesh_i.point(viter)[1] < min_y_) min_y_ = mesh_i.point(viter)[1];
			if (mesh_i.point(viter)[1] > max_y_) max_y_ = mesh_i.point(viter)[1];
			
			if (mesh_i.point(viter)[1] < min_i_y_temp_) min_i_y_temp_ = mesh_i.point(viter)[1];
			if (mesh_i.point(viter)[1] > max_i_y_temp_) max_i_y_temp_ = mesh_i.point(viter)[1];
		}
		min_i_y_.push_back(min_i_y_temp_);
		max_i_y_.push_back(max_i_y_temp_);
	}
	std::cerr << "21" << std::endl;
	this->establishMeshCurveMap(min_y_, max_y_);
	std::cerr << "22" << std::endl;
	this->computeProjectionImage(min_y_, max_y_);
	std::cerr << "23" << std::endl;
	return true;
}

bool teethProjection::establishMeshCurveMap(double min_i_y_, double max_i_y_)
{
	double limit = 0.0001;
	double judge, min_i_y_temp, max_i_y_temp;
	for (auto iter = this->sep_meshes_.begin(); iter != this->sep_meshes_.end(); iter++)
	{
		COpenMeshT mesh_tooth_ = *iter->second;
		for (auto i = mesh_tooth_.vertices_begin(); i != mesh_tooth_.vertices_end(); i++)
		{
			min_i_y_temp = min_i_y_;
			max_i_y_temp = max_i_y_;
			std::cerr << "30" << std::endl;
			judge = this->computeNeareastPointOfCurve(min_i_y_temp, mesh_tooth_.point(i)) * this->computeNeareastPointOfCurve(max_i_y_temp, mesh_tooth_.point(i));
			std::cerr << "31" << std::endl;
			if (judge > 0) std::cerr << "error" << std::endl;
			else
			{
				while (max_i_y_temp - min_i_y_temp > limit)
				{
					if (this->computeNeareastPointOfCurve((min_i_y_temp + max_i_y_temp) / 2, mesh_tooth_.point(i)) * this->computeNeareastPointOfCurve(max_i_y_temp, mesh_tooth_.point(i)) < 0)
					{
						min_i_y_temp = (min_i_y_temp + max_i_y_temp) / 2;
					}
					else
					{
						max_i_y_temp = (min_i_y_temp + max_i_y_temp) / 2;
					}
				}
			}
			OpenMesh::Vec2d temp;
			temp[0] = min_i_y_temp;
			temp[1] = this->computeCurveValue(min_i_y_temp);
			this->mesh_curve_map_[iter->first].push_back(temp);
		}
	}

	return true;
}

double teethProjection::computeDerivateCurve(double y)
{
	double curve_derivate_x = 0;
	double curve_y = 1;
	for (auto i = 1; i < this->curve_ceoffs_.size(); i++)
	{
		curve_derivate_x += this->curve_ceoffs_[i] * curve_y * i;
		curve_y = curve_y * y;
	}
	return curve_derivate_x;
}

double teethProjection::computeCurveValue(double y)
{
	double curve_x = 0;
	double curve_y = 1;
	for (auto i = 0; i < this->curve_ceoffs_.size(); i++)
	{
		curve_x += this->curve_ceoffs_[i] * curve_y;
		curve_y = curve_y*y;
	}
	return curve_x;
}

double teethProjection::computeLenghOfCurve(double min_i_y, double max_i_y, double y)
{
	OpenMesh::Vec2d p1, p2;
	p1[0] = y;
	p1[1] = this->computeCurveValue(y);
	double len_curve = 0;
	for (double i = y; i < max_i_y;)
	{
		p2[0] = i;
		p2[1] = this->computeCurveValue(i);
		len_curve = len_curve + sqrt(pow(p2[1] - p1[1], 2) + pow(p2[0] - p1[0], 2));
		p1 = p2;
		i = i + 0.01;
	}
	return len_curve;
}

double teethProjection::computeNeareastPointOfCurve(double y, OpenMesh::Vec3d & tooth_point)
{
	return (y - tooth_point[1]) + this->computeDerivateCurve(y)*(this->computeCurveValue(y) - tooth_point[0]);
}

bool teethProjection::computeProjectionImage(double min_i_y, double max_i_y)
{
	for (auto iter = this->sep_meshes_.begin(); iter != this->sep_meshes_.end(); iter++)
	{
		COpenMeshT mesh_i = *iter->second;
		OpenMesh::Vec2d temp;
		for (auto viter = mesh_i.vertices_begin(); viter != mesh_i.vertices_end(); viter++)
		{
			temp[0] = this->computeLenghOfCurve(min_i_y, max_i_y, this->mesh_curve_map_[iter->first][viter->idx()][0]);
			temp[1] = mesh_i.point(viter)[2];
			this->projection_teeth_image_[iter->first].push_back(temp);
		}
	}
	return true;
}
std::map<int, std::vector<OpenMesh::Vec2d>> teethProjection::getProjectionImage()
{
	this->computeProjectionImageContour();
	return this->projection_teeth_image_contour_;
}
std::map<int, OpenMesh::Vec3d> teethProjection::getTest()
{
	return map_toothId_TopPoint_;
}

std::vector<std::map<int, int>> teethProjection::gettest1()
{
	return this->map_3dvertice_2dpoints_;
}


bool teethProjection::computeProjectionImageContour()
{
	std::vector<std::vector<int>> vids;
	if (this->projection_teeth_image_contour_.size() != 0) this->projection_teeth_image_contour_.clear();
	if (this->map_vids_.size() != 0) this->map_vids_.clear();
	for (auto iter = this->projection_teeth_image_.begin(); iter != this->projection_teeth_image_.end(); iter++)
	{
		std::vector<OpenMesh::Vec2d> tooth_contour_ = iter->second;
		CGeoAlg::AlphaShape2d(tooth_contour_, 0.5, vids);
		this->map_vids_.push_back(vids[0]);
		for (auto i = 0; i < vids[0].size(); i++)
		{
			projection_teeth_image_contour_[iter->first].push_back(tooth_contour_[vids[0][i]]);
			std::map<int, int> temp;
			temp[iter->first] = vids[0][i];
			this->map_3dvertice_2dpoints_.push_back(temp);
		}
	}
	return true;
}

bool teethProjection::calculateTeethTop(std::map<int, std::map<int, std::map<int, double>>>& map_toothId_map_meshId_tipLen_)
{
	int count = 0;
	for (auto ii = map_toothId_map_meshId_tipLen_.begin(); ii != map_toothId_map_meshId_tipLen_.end(); ii++)
	{
		for (auto iter = ii->second.begin(); iter != ii->second.end(); iter++) {
			for (auto iiter = iter->second.begin(); iiter != iter->second.end(); iiter++)
			{
				COpenMeshT mesh_i = *this->sep_meshes_[iter->first];
				for (auto viter = mesh_i.vertices_begin(); viter != mesh_i.vertices_end(); viter++) {
					if (viter->idx() == iiter->first)
					{
						this->map_toothId_TopPoint_[count] = mesh_i.point(viter);
						this->map_toothId_TopPoint_[count][2] = this->map_toothId_TopPoint_[count][2] + iiter->second;
						std::cerr << mesh_i.point(viter)[0] << " " << mesh_i.point(viter)[1] << " " << mesh_i.point(viter)[2] << " " << iiter->second << std::endl;
					}
				}
				count++;
			}
		}
	}
	//this->sortTeethTop();
	return true;
}

bool teethProjection::sortTeethTop()
{
	if (this->map_toothId_TopPoint_.size() == 0) return false;
	std::vector<OpenMesh::Vec3d> sort_teethtop;
	OpenMesh::Vec3d temp;
	for (auto iter = this->map_toothId_TopPoint_.begin(); iter != this->map_toothId_TopPoint_.end(); iter++)
	{
		sort_teethtop.push_back(iter->second);
	}
	for (auto i = 1; i < sort_teethtop.size(); i++)
	{
		for (auto j = sort_teethtop.size() - 1; j >= i; j--)
		{
			if (sort_teethtop[j][1] > sort_teethtop[j-1][1])
			{
				temp = sort_teethtop[j - 1];
				sort_teethtop[j - 1] = sort_teethtop[j];
				sort_teethtop[j] = temp;
			}
		}
	}
	for (auto iter = this->map_toothId_TopPoint_.begin(); iter != this->map_toothId_TopPoint_.end(); iter++)
	{
		this->map_toothId_TopPoint_[iter->first] = sort_teethtop[iter->first];
	}
	return true;
}

void teethProjection::test()
{
	this->curveDiscretization();
	std::cerr << "test" << std::endl;
}