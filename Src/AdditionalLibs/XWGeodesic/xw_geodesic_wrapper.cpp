#include"xw_geodesic_wrapper.h"
#include"Xin_Wang.h"
#include<numeric>
void CXWGeodesic::SetModel(std::vector<CPoint3D>& vertexs, std::vector<CBaseModel::CFace>&faces)
{
	if (model_ != NULL)
		delete model_;
	model_ = new CRichModel();
	model_->SetModel(vertexs, faces);
	//model_->SaveObjFile("test.obj");
}
void CXWGeodesic::GeodesicDis(std::vector<int>&svids, std::vector<double>&dis)
{
	std::set<int>set_tvids;
	for (int i = 0; i < svids.size(); i++)
		set_tvids.insert(svids[i]);
	CXin_Wang alg(*model_, set_tvids);
	alg.Execute();
	const std::vector<double>&dis_field = alg.GetDistanceField();
	dis.resize(dis_field.size());
	for (int i = 0; i < dis.size(); i++)
	{
		dis[i] = dis_field[i];
	}
}
void CXWGeodesic::GeodesicDis(int svid, std::vector<int>&tvids, std::vector<double>&dis)
{
	std::set<int>srcs, dsts;
	srcs.insert(svid);
	for (int i = 0; i < tvids.size(); i++)
	{
		dsts.insert(tvids[i]);
	}
	CXin_Wang alg(*model_, srcs,dsts);
	//CXin_Wang(const CRichModel& model, const set<int>& sources, const set<int>& destinations);
	alg.Execute();
	const std::vector<double>&dis_field=alg.GetDistanceField();
	dis.resize(tvids.size(), -1);
	for (int i = 0; i < tvids.size(); i++)
	{
		dis[i] = dis_field[tvids[i]];
	}

}
void CXWGeodesic::GeodesicPath(std::vector<int>&svids, std::vector<int>&tvids, std::vector<std::vector<CGeoFacePoint>>&paths)
{
	std::set<int>set_tvids;
	for (int i = 0; i < tvids.size(); i++)
		set_tvids.insert(tvids[i]);
	CXin_Wang alg(*model_, set_tvids);
	alg.Execute();
	paths.resize(svids.size());
	const std::vector<double>&dis_field = alg.GetDistanceField();
	for (int i = 0; i < svids.size(); i++)
	{
		std::vector<EdgePoint> epath = alg.BacktraceShortestPath(svids[i]);
		std::reverse(epath.begin(), epath.end());
		paths[i].resize(epath.size());
		for (int j = 0; j < epath.size(); j++)
		{
			int eid = epath[j].index;
			auto edge = model_->Edge(eid);
			if (epath[j].isVertex)
			{
				eid = model_->Neigh(epath[j].index)[0].first;
				edge = model_->Edge(eid);
			}
			int fid = edge.indexOfFrontFace;
			auto point3d = epath[j].Get3DPoint(*model_);
			auto face = model_->Face(fid);
			//point3d = model_->m_Verts[face.verts[0]];
			paths[i][j].pos_[0] = point3d.x;
			paths[i][j].pos_[1] = point3d.y;
			paths[i][j].pos_[2] = point3d.z;

			paths[i][j].vids_[0] = face.verts[0];
			paths[i][j].vids_[1] = face.verts[1];
			paths[i][j].vids_[2] = face.verts[2];

			for (int k = 0; k < 3; k++)
			{
				if (face.verts[k] == edge.indexOfOppositeVert)
				{
					paths[i][j].fid_ = fid;
					paths[i][j].ls_[k] = 0;
					paths[i][j].ls_[(k + 1) % 3] = 1 - epath[j].proportion;
					paths[i][j].ls_[(k + 2) % 3] = epath[j].proportion;


					break;
				}
			}
		}
	}
	

	
}
void CXWGeodesic::GeodesicPath(int svid, std::vector<int>&tvids, std::vector<CGeoFacePoint>&path)
{

	//CXin_Wang alg(*model_, srcs, dsts);

	CXin_Wang alg(*model_, svid);
	alg.Execute();
	const std::vector<double>&dis_field = alg.GetDistanceField();
	double min_dis = 999999999;
	int mi;
	for (int i = 0; i < tvids.size(); i++)
	{
		if (min_dis > dis_field[tvids[i]])
		{
			min_dis = dis_field[tvids[i]];
			mi = i;
		}
	}
	//std::cerr << "geo mi " << mi<<" " <<tvids[mi]<<" dis:"<<min_dis<< std::endl;
	std::vector<EdgePoint> epath = alg.BacktraceShortestPath(tvids[mi]);
	path.resize(epath.size());
	for (int i = 0; i < epath.size(); i++)
	{
		int eid = epath[i].index;
		auto edge = model_->Edge(eid);
		if (epath[i].isVertex)
		{
			 eid = model_->Neigh(epath[i].index)[0].first;
			edge = model_->Edge(eid);
		}
		int fid = edge.indexOfFrontFace;
		auto point3d = epath[i].Get3DPoint(*model_);
		auto face = model_->Face(fid);
		//point3d = model_->m_Verts[face.verts[0]];
		path[i].pos_[0] = point3d.x;
		path[i].pos_[1] = point3d.y;
		path[i].pos_[2] = point3d.z;
	
		path[i].vids_[0] = face.verts[0];
		path[i].vids_[1] = face.verts[1];
		path[i].vids_[2] = face.verts[2];

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
	//std::cerr << path.front().vids_[0] << " " << path.front().vids_[1] << " " << path.front().vids_[2] << std::endl;
	//std::cerr << path.back().vids_[0] << " " << path.back().vids_[1] << " " << path.back().vids_[2] << std::endl;
}
void CXWGeodesic::GeodesicPath(int svid, int tvid, std::vector<CGeoFacePoint>&path)
{
	
	CXin_Wang alg(*model_, svid);
	//std::cerr << "start exe" << std::endl;
	alg.Execute();
	//std::cerr << "end exe" << std::endl;
	std::vector<EdgePoint> epath=alg.BacktraceShortestPath(tvid);
	path.resize(epath.size());
	for (int i = 0; i < epath.size(); i++)
	{
		int eid = epath[i].index;
		auto edge = model_->Edge(eid);
		if (epath[i].isVertex)
		{
		
			int eid = model_->Neigh(epath[i].index)[0].first;
			//std::cerr << "is vert " << epath[i].index<<" "<< nvid << std::endl;
			 //eid = model_->GetEdgeIndexFromTwoVertices(epath[i].index, nvid);
			edge = model_->Edge(eid);
			std::cerr << "is vert " <<eid<< std::endl;
		}
		int fid = edge.indexOfFrontFace;
		auto face=model_->Face(fid);
		path[i].fid_ = fid;
		auto point3d = epath[i].Get3DPoint(*model_);
		point3d = model_->m_Verts[face.verts[0]];
		path[i].vids_[0]=face.verts[0];
		path[i].vids_[1] = face.verts[1];
		path[i].vids_[2] = face.verts[2];
	
		path[i].pos_[0] = point3d.x;
		path[i].pos_[1] = point3d.y;
		path[i].pos_[2] = point3d.z;
		for (int j = 0; j < 3; j++)
		{
			if (face.verts[j] == edge.indexOfOppositeVert)
			{
				
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