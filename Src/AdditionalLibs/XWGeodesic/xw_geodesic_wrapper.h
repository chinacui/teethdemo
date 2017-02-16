#ifndef CXW_GEODESIC_WRAPPER_H
#define CXW_GEODESIC_WRAPPER_H
#include<iostream>
#include<vector>
#include"BaseModel.h"

class CRichModel;
class CGeoFacePoint
{
public:
	int fid_;
	double ls_[3];
	double pos_[3];
	int vids_[3];
};
class ADDITIONALLIBS_CLASS CXWGeodesic
{
protected:
	CRichModel* model_=NULL;

public:
	CXWGeodesic() {  }
	void SetModel(std::vector<CPoint3D>& vertexs, std::vector<CBaseModel::CFace>&faces);
	void GeodesicDis(int svid, std::vector<int>&tvids, std::vector<double>&dis);
	void GeodesicDis(std::vector<int>&svids, std::vector<double>&dis);
	void GeodesicPath(int svid, int tvid, std::vector<CGeoFacePoint>&path);
	void GeodesicPath(int svid, std::vector<int>&tvids, std::vector<CGeoFacePoint>&path);
	void GeodesicPath(std::vector<int>&svids, std::vector<int>&tvids, std::vector<std::vector<CGeoFacePoint>>&path);
	~CXWGeodesic();

};
#endif