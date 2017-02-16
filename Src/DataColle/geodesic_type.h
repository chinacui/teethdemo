#ifndef CGEODESIC_TYPE_H
#define CGEODESIC_TYPE_H
#include"../AdditionalLibs/XWGeodesic/xw_geodesic_wrapper.h"
#include"custom_openmesh_type.h"
class CGeodesicModel /*:public CRichModel*/
{
protected:
	COpenMeshT *mesh_=NULL;
public:
	CGeodesicModel(COpenMeshT *mesh);
	void Update();
	~CGeodesicModel() {}
};
#endif