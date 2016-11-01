#include"ModelWrapper.h"

void ModelWrapper::setMesh(std::vector<CPoint3D>& vertexs, std::vector<CFace>&faces)
{
	m_Verts = vertexs;
	m_Faces = faces;

}