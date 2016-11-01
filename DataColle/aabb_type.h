#ifndef CAABB_TYPE_H
#define CAABB_TYPE_H
#include"Polyhedron_type.h"
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Aff_transformation_3.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include"custom_openmesh_type.h"
namespace CGAL
{
	template<class T>
	class AABB_tree;

	template<class T, class VertexPointPMap, class OneFaceGraphPerTree, class CacheDatum>
	class AABB_face_graph_triangle_primitive;

	template<class T, class AABBPrimitive>
	class AABB_traits;
	template<class T>
	class Aff_transformation_3;

}
typedef CGAL::AABB_face_graph_triangle_primitive<Polyhedron, CGAL::Default, CGAL::Tag_true, CGAL::Tag_false> Primitive;
typedef CGAL::AABB_traits<Kernel, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> AABBTree;
class CAABBTree  {
public:
	AABBTree *tree_=NULL;
	Polyhedron poly;
	std::vector<COpenMeshT::VertexHandle>id_vh_map_;
	std::vector<COpenMeshT::FaceHandle>id_fh_map_;
	CAABBTree() {
		tree_ = NULL;
		id_fh_map_.clear();
		id_fh_map_.clear();
	};
	~CAABBTree() {
		if (tree_ != NULL)
		{
			delete tree_;

		}
	}
};
typedef Kernel::Ray_3 Ray;
typedef boost::optional< AABBTree::Intersection_and_primitive_id<Ray>::Type > Ray_intersection;
typedef AABBTree::Point_and_primitive_id AABB_Point_and_primitive_id;
#endif