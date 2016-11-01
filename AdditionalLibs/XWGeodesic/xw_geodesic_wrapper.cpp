#include"xw_geodesic_wrapper.h"
#include"Xin_Wang.h"
void CXWGeodesic::SetModel(std::vector<CPoint3D>& vertexs, std::vector<CBaseModel::CFace>&faces)
{
	if (model_ != NULL)
		delete model_;
	model_ = new CRichModel();
	model_->SetModel(vertexs, faces);
}
void CXWGeodesic::GeodesicPath(int svid, int tvid, std::vector<CGeoFacePoint>&path)
{
	CXin_Wang alg(*model_, svid);
	std::cerr << "start exe" << std::endl;
	alg.Execute();
	std::cerr << "end exe" << std::endl;
	std::vector<EdgePoint> epath=alg.BacktraceShortestPath(tvid);
	path.resize(epath.size());
	for (int i = 0; i < epath.size(); i++)
	{
		int eid = epath[i].index;
		auto edge = model_->Edge(eid);
		int fid = edge.indexOfFrontFace;
		auto face=model_->Face(fid);
		for (int j = 0; j < 3; j++)
		{
			if (face.verts[j] == edge.indexOfOppositeVert)
			{
				path[i].fid_ = fid;
				path[i].ls_[j] = 0;
				path[i].ls_[(j + 1) % 3] = 1 - epath[i].proportion;
				path[i].ls_[(j + 2) % 3] = epath[i].proportion;
				break;
			}
		}
	}
}
CXWGeodesic::~CXWGeodesic()
{
	if (model_ != NULL)
		delete model_;
}