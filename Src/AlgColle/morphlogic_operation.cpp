#include"morphlogic_operation.h"
#include<queue>
void CMorphlogicOperation::Dilate(COpenMeshT &mesh, std::vector<bool>&labels)
{
	std::vector<bool>in_labels = labels;
	for (int i = 0; i < labels.size(); i++)
	{
		labels[i] = 0;
	}
	for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); viter++)
	{
		int vid = viter->idx();
		if (in_labels[vid] == true)
		{
			labels[vid] = true;
			continue;
		}
	
		for (auto nviter = mesh.vv_begin(viter); nviter != mesh.vv_end(viter); nviter++)
		{
			int nvid = nviter->idx();
			if (in_labels[nvid] == true)
			{
				labels[vid] = true;
				break;
			}
		}
	}
}

void CMorphlogicOperation::Erode(COpenMeshT &mesh, std::vector<bool>&labels)
{
std::vector<bool>in_labels = labels;
for (int i = 0; i < labels.size(); i++)
{
	labels[i] = true;
}
for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); viter++)
{
	int vid = viter->idx();
	if (in_labels[vid] == false)
	{
		labels[vid] = false;
		continue;
	}

	for (auto nviter = mesh.vv_begin(viter); nviter != mesh.vv_end(viter); nviter++)
	{
		int nvid = nviter->idx();
		if (in_labels[nvid] == false)
		{
			labels[vid] = false;
			break;
		}
	}
}
}
bool CMorphlogicOperation::IsCenterVertex(COpenMeshT & mesh, std::vector<bool>&label, COpenMeshT::VertexHandle vh)
{
	if (!label[vh.idx()])
		return false;
	for (auto vviter = mesh.vv_begin(vh); vviter != mesh.vv_end(vh); vviter++)
	{
		if (label[vviter->idx()] == false)
			return false;
	}
	return true;
}
int CMorphlogicOperation::ComputeDegree(COpenMeshT &mesh, std::vector<bool>&label, COpenMeshT::VertexHandle vh)
{
	int count = 0;
	for (auto vviter = mesh.vv_begin(vh); vviter != mesh.vv_end(vh); vviter++)
	{
		if (label[vviter->idx()] == true)
			count++;
	}
	return count;
}
bool CMorphlogicOperation::IsComplexVertex(COpenMeshT & mesh, std::vector<bool>&label, COpenMeshT::VertexHandle vh)
{
	if (!label[vh.idx()])
		return false;
	auto prevviter = mesh.vv_ccwbegin(vh);
	auto vviter = prevviter;
	vviter++;
	int count = 0;
	for (; vviter != mesh.vv_ccwend(vh); vviter++)
	{
		if (label[vviter->idx()] != label[prevviter->idx()])
			count++;
		prevviter = vviter;
	}
	vviter = mesh.vv_ccwbegin(vh);
	if (label[prevviter->idx()] != label[vviter->idx()])
		count++;
	if (count >= 4)
		return true;
	else
		return false;
}
void CMorphlogicOperation::Skeletonization(COpenMeshT& mesh, std::vector<bool>&labels)
{
	std::vector<bool>in_labels = labels;
	for (int i = 0; i < labels.size(); i++)
	{
		labels[i] = true;
	}
	std::vector<int>region_ids;
	std::vector<CMOVertexTag>region_tags;
	TagAllVertexs(mesh, in_labels, region_tags, region_ids);
	for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); viter++)
	{
		int vid = viter->idx();
		if (in_labels[vid] == false)
		{
			labels[vid] = false;

		}
		/*else if (mesh.is_boundary(viter))
		{
			labels[vid] = true;
		}*/
		else if (IsComplexVertex(mesh, in_labels, viter))
		{
			labels[vid] = true;

		}
		else if (region_ids[vid] == CMOVertexTag::Center)
		{
			labels[vid] = true;
		}
		else if (region_tags[vid] == CMOVertexTag::Disk)
		{
			/*if (ComputeDegree(mesh, in_labels, viter)==1)
			{
				labels[vid] = true;
			}
			else
			{*/
				//labels[vid] = false;
				//auto nviter = mesh.vv_begin(viter);
				//for (; nviter != mesh.vv_end(viter); nviter++)
				//{
				//	int nvid = nviter->idx();
				//	if (region_tags[nvid] == CMOVertexTag::Disk&&region_ids[nvid] != region_ids[vid])
				//	{
				//		//labels[vid] = true;
				//		break;
				//	}
				//}
				//if (nviter != mesh.vv_end(viter))
				//{
					in_labels[vid] = false;
					labels[vid] = false;
					for (auto nviter = mesh.vv_begin(viter); nviter != mesh.vv_end(viter); nviter++)
					{
						int nvid = nviter->idx();
						if (IsComplexVertex(mesh, in_labels, nviter))
						{
							region_ids[nviter->idx()] = -1;
							region_tags[nviter->idx()] = CMOVertexTag::Complex;
						}

					}
				//}
			//}
			
		}
		
		
	}
}
void CMorphlogicOperation::FloodFill(COpenMeshT &mesh, std::vector<bool>&label, bool roi_label, COpenMeshT::VertexHandle start_vh, int target_tag, std::vector<int>&res_region_tag)
{
	if (label[start_vh.idx()] != roi_label)
		return;
	if (mesh.n_vertices() != res_region_tag.size())
		res_region_tag.resize(mesh.n_vertices());
	std::vector<bool>mmark(mesh.n_vertices(), 0);
	std::queue<COpenMeshT::VertexHandle> Q;
	Q.push(start_vh);
	res_region_tag[start_vh.idx()] = target_tag;
	mmark[start_vh.idx()] = true;
	while (!Q.empty())
	{
		auto pviter = Q.front();
		Q.pop();
		for (auto vviter = mesh.vv_begin(pviter); vviter != mesh.vv_end(pviter); vviter++)
		{
			if (label[vviter->idx()] == roi_label&&mmark[vviter->idx()]==false)
			{
				Q.push(vviter);
				mmark[vviter->idx()] = true;
				res_region_tag[vviter->idx()] = target_tag;
			}
		}
	}
}
void CMorphlogicOperation::TagAllVertexs(COpenMeshT &mesh, std::vector<bool>&labels, std::vector<CMOVertexTag>&tags, std::vector<int>&region_ids)
{
	region_ids.resize(mesh.n_vertices(), -1);
	tags.resize(mesh.n_vertices(), CMOVertexTag::NonFeature);
	for (int i = 0; i < mesh.n_vertices(); i++)
	{
		tags[i] = CMOVertexTag::NonFeature;
		region_ids[i] = -1;
	}
		
	
	int t = 0;
	for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); viter++)
	{
		if (region_ids[viter->idx()] == -1&&labels[viter->idx()]==false)
		{
			FloodFill(mesh, labels, false, viter, t, region_ids);
			t++;
		}
		
	}
	for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); viter++)
	{
		if (labels[viter->idx()] == true)
		{
			if (IsComplexVertex(mesh, labels, viter))
			{
				tags[viter->idx()] = CMOVertexTag::Complex;
			}
			else if (IsCenterVertex(mesh, labels, viter))
			{
				tags[viter->idx()] = CMOVertexTag::Center;
			}
			else
			{
				tags[viter->idx()] = CMOVertexTag::Disk;
			
				for (auto vviter = mesh.vv_begin(viter); vviter != mesh.vv_end(viter); vviter++)
				{
					if (region_ids[vviter->idx()]!=-1&&tags[vviter->idx()]==CMOVertexTag::NonFeature)
					{
						region_ids[viter->idx()] = region_ids[vviter->idx()];
						break;
					}
				}
			}
		}
	}

}