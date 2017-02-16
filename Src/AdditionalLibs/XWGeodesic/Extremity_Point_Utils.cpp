#include"Extremity_Point_Utils.h"
#include<omp.h>



ExtremityPointUtils::ExtremityPointUtils(CRichModel &model) :m_model(model)
{

}
void ExtremityPointUtils::ComputeMeanDistances(double radius, std::vector<double>& meanDis,int num_thread)
{
	meanDis.resize(m_model.GetNumOfVerts());
#pragma omp parallel for num_threads(num_thread)
	for (int i = 0; i < m_model.GetNumOfVerts(); i++)
	{
		cerr << "compute mean dis at vert: " << i << endl;
		ExtremityPointUtils alg(m_model);
		meanDis[i] = alg.ComputeMeanDistanceAtSourcePoint(i, radius);
	}
}
double ExtremityPointUtils::ComputeMeanDistanceAtSourcePoint(int source, double radius)
{
	CXin_Wang alg(m_model, source,radius);
	alg.Execute();
	set<int> vIDs;
	queue<int>vQue;
	vQue.push(source);
	while (!vQue.empty())
	{
		int p = vQue.front();
		vQue.pop();
		for (int i = 0; i < m_model.Neigh(p).size(); i++)
		{

			int id = m_model.Edge(m_model.Neigh(p)[i].first).indexOfRightVert;
			if (alg.m_InfoAtVertices[id].disUptodate <= radius&&vIDs.find(id) == vIDs.end())
			{
				vIDs.insert(id);
				vQue.push(id);
			}
			
		}
	}
	
	double mean_dis = 0;
	for (auto iter = vIDs.begin(); iter != vIDs.end(); iter++)
	{
		mean_dis += alg.m_InfoAtVertices[*iter].disUptodate;
	}
	mean_dis /= vIDs.size();
	return mean_dis;
}