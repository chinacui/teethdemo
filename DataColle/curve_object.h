#ifndef CCURVE_OBJECT_H
#define CCURVE_OBJECT_H
#include"prereq.h"
#include"base_object.h"
#include"custom_openmesh_type.h"
class DATACOLLE_CLASS CCurveObject :public CBaseObject
{
public:
	enum CurveType { Default, Dots };
protected:
	std::vector<OpenMesh::Vec3d>curve_;
	int curve_id_;
	bool is_changed_;
	OpenMesh::Vec3d color_;
	CurveType render_type_;
public:
	int GetId();//get id of curve object
	void SetId(int id);//set id of curve object
	bool IsChanged();
	void SetChanged(bool is_changed = true); 
	std::vector<OpenMesh::Vec3d>&GetCurve();
	void SetCurve(std::vector<OpenMesh::Vec3d>& curve);
	OpenMesh::Vec3d GetColor();
	void SetColor(OpenMesh::Vec3d color);
	CurveType& RendereType() { return render_type_; }
	CCurveObject();
	CCurveObject(CCurveObject&b);
	~CCurveObject();

};
#endif