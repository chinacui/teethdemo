#ifndef DATA_POOL_H
#define DATA_POOL_H
#include"prereq.h"
#include<iostream>
#include<map>
#include<vector>
class DATACOLLE_CLASS DataPool
{
protected:
	static std::map<int, boost::shared_ptr<PolyhedronObject>> polyhedron_object_pool_;
	static std::map<int, boost::shared_ptr<PolyhedronObject>>auxiliary_shape_pool_;
	static std::vector<boost::shared_ptr<SkeletonObject>> skeleton_list_;
	static std::map<int, boost::shared_ptr<CurveObject>> curve_list_;
	static int polyhedron_object_max_id_;
	static int curve_object_max_id_;
	static int auxiliary_shape_max_id_;
};
#endif
