#include"curve_object.h"
int CCurveObject::GetId()
{
	return curve_id_;
}
void CCurveObject::SetId(int id)
{
	curve_id_ = id;
}
bool CCurveObject::IsChanged()
{
	return is_changed_;
}
void CCurveObject::SetChanged(bool is_changed)
{
	is_changed_ = is_changed;
}
void CCurveObject::SetCurve(std::vector<OpenMesh::Vec3d>&curve)
{
	curve_ = curve;
}
std::vector<OpenMesh::Vec3d>& CCurveObject::GetCurve()
{
	return curve_;
}

CCurveObject::CCurveObject()
{
	curve_id_ = -1;
	color_ = OpenMesh::Vec3d(0, 0, 0);
	render_type_ = Default;
	is_changed_ = true;
}
CCurveObject::CCurveObject(CCurveObject&b)
{
	curve_id_ = b.curve_id_;
	curve_ = b.curve_;
	is_changed_ = true;
	color_ = b.color_;
	render_type_ = b.render_type_;
}
OpenMesh::Vec3d CCurveObject::GetColor()
{
	return color_;
}
void CCurveObject::SetColor(OpenMesh::Vec3d color)
{
	color_ = color;
}
CCurveObject::~CCurveObject()
{

}