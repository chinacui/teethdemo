#include "Xin_Wang.h"
#include"Extremity_Point_Utils.h"

int main(int argc, char* argv[])
{
	CRichModel model("camel.obj");
	model.LoadModel();
	//example: compute the geodesic distance for a single source
	std::vector<double> meanDis;
	ExtremityPointUtils alg(model);
	alg.ComputeMeanDistances(0.2, meanDis);
	model.SaveScalarFieldObjFile(meanDis, "camel_meanDis.obj");
	getchar();
	//see more examples:
	//CXin_Wang(const CRichModel& model, int source);
	//CXin_Wang(const CRichModel& model, int source, int destination);
	//CXin_Wang(const CRichModel& model, int source, double R);
	//CXin_Wang(const CRichModel& model, const map<int, double>& sources);
	//CXin_Wang(const CRichModel& model, const map<int, double>& sources, const set<int> &destinations);	
	//CXin_Wang(const CRichModel& model, const set<int>& sources);	
	//CXin_Wang(const CRichModel& model, const set<int>& sources, double R);
	//CXin_Wang(const CRichModel& model, const set<int>& sources, const set<int>& destinations);
}

