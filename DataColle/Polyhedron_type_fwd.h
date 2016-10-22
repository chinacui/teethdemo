#ifndef POLYHEDRON_TYPE_FWD_H
#define POLYHEDRON_TYPE_FWD_H

#include <memory>  // for std::allocator
#include <utility> // for std::pair
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Polyhedron_items_with_id_3.h>
#ifdef USE_FORWARD_DECL

namespace CGAL {

  template < typename FT_ >
  struct Simple_cartesian;

  class Epick;
  class Epeck;

  class Polyhedron_items_3;

  template < class T, class I, class A>
  class HalfedgeDS_default;

  template < class PolyhedronTraits_3,
             class PolyhedronItems_3,
             template < class T, class I, class A>
             class T_HDS, 
             class Alloc
             >
  class Polyhedron_3;

  template <typename Kernel_, typename Items_, typename Mark_>
  class Nef_polyhedron_3;

  class SNC_indexed_items;

  namespace Mesh_3 {
    template <typename Patch_id>
    class Mesh_polyhedron_items;
  } // end namespace Mesh_3



} // end namespace CGAL

// kernel

typedef CGAL::Epick Kernel;
typedef CGAL::Epeck Exact_Kernel;

typedef std::pair<int, int> Patch_id;


template < class Refs>
class Halfedge_with_attributes : public CGAL::HalfedgeDS_halfedge_max_base_with_id<Refs, std::size_t>
{
public:
	CGAL::Color color_;
	Kernel::Point_2 uv_ = Kernel::Point_2(0, 0);
	Kernel::Point_2 uv_spare_ = Kernel::Point_2(0, 0);
};

class Polyhedron_items_with_attributes_3 :public  CGAL::Polyhedron_items_with_id_3{
public:
	
	template < class Refs, class Traits>
	struct Halfedge_wrapper :public CGAL::Polyhedron_items_with_id_3::Halfedge_wrapper<Refs,Traits>  {
		typedef Halfedge_with_attributes<Refs> Halfedge;
	};
};


// surface mesh
//typedef CGAL::Polyhedron_3<Kernel,
//                           CGAL::Mesh_3::Mesh_polyhedron_items<Patch_id>,
//                           CGAL::HalfedgeDS_default,
//                           std::allocator<int> > Polyhedron;
typedef CGAL::Polyhedron_3<Kernel,
						Polyhedron_items_with_attributes_3,
                           CGAL::HalfedgeDS_default,
                           std::allocator<int> > Polyhedron;
////////////////////////////////////////////////////////////////



typedef CGAL::Nef_polyhedron_3<Exact_Kernel, CGAL::SNC_indexed_items, bool> Nef_polyhedron;

typedef std::vector<Kernel::Point_3> Curve;

typedef Kernel::Point_2 Point_2;





#else // USE_FORWARD_DECL

#include "Polyhedron_type.h"

#endif // USE_FORWARD_DECL

#endif // POLYHEDRON_TYPE_FWD_H
