#ifndef CCURVE_OBJECT_H
#define CCURVE_OBJECT_H
#include"prereq.h"
#include"base_object.h"
#include"custom_openmesh_type.h"
class DATACOLLE_CLASS CCurveObject :public CBaseObject
{
protected:
	std::vector<OpenMesh::Vec3d>curve_;
	int curve_id_;
	bool is_changed_;
	OpenMesh::Vec3d color_;
public:
	int GetId();//get id of curve object
	void SetId(int id);//set id of curve object
	bool IsChanged();
	void SetChanged(bool is_changed = true); 
	std::vector<OpenMesh::Vec3d>&GetCurve();
	OpenMesh::Vec3d GetColor();
	void SetColor(OpenMesh::Vec3d color);

	CCurveObject();
	CCurveObject(CCurveObject&b);
	~CCurveObject();

};
#endif