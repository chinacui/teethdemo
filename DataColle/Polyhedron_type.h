#ifndef POLYHEDRON_TYPE_H
#define POLYHEDRON_TYPE_H

// CGAL
// kernel
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/Bbox_3.h>
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Exact_predicates_exact_constructions_kernel Exact_Kernel;
// surface mesh
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Mesh_polyhedron_3.h>
//#include <CGAL/Polyhedron_items_with_id_3.h>
#include "Polyhedron_type_fwd.h"

// simple geometric types
typedef Kernel::FT FT;
typedef Kernel::Line_3 Line;
typedef Kernel::Point_3 Point;
typedef Kernel::Plane_3 Plane;
typedef Kernel::Sphere_3 Sphere;
typedef Kernel::Vector_3 Vector;
typedef Kernel::Triangle_3 Triangle;
typedef Kernel::Iso_cuboid_3 Iso_cuboid;
typedef Kernel::Plane_3 Plane_3;
typedef Kernel::Point_2 Point_2;
typedef Kernel::Vector_2 Vector_2;
typedef Kernel::Segment_2 Segment_2;
typedef CGAL::Bbox_3 Bbox_3;
typedef CGAL::Bbox_2 Bbox_2;
typedef Kernel::Line_2 Line_2;
typedef Kernel::Segment_2 Segment_2;
typedef Kernel::Segment_3 Segment_3;
// surface mesh, Patch_id is pair<int,int>
//typedef CGAL::Mesh_polyhedron_3<Kernel, Patch_id>::type Polyhedron;


typedef CGAL::Polyhedron_3<Kernel, Polyhedron_items_with_attributes_3> Polyhedron;


typedef Polyhedron::Vertex_handle   Vertex_handle;
typedef Polyhedron::Halfedge_handle Halfedge_handle;
typedef Polyhedron::Face_handle Face_handle;

typedef std::vector<Point> Curve;
typedef CGAL::Aff_transformation_3<Kernel> AffTrans;
/// Returns a normalized (unit) vector.
/// @param v input vector.
/// @return v normalized. 
inline Vector normalize(const Vector& v) {
	double len = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
	return Vector(v[0] / len, v[1] / len, v[2] / len);
}

#endif // POLYHEDRON_TYPE_H
