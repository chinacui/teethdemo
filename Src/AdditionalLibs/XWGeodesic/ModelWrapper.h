#ifndef MODEL_WRAPPER_H
#define MODEL_WRAPPER_H
#include"BaseModel.h"
class ModelWrapper:public CBaseModel
{
public:
	void setMesh(std::vector<CPoint3D>& vertexs, std::vector<CFace>&faces);
};
#endif