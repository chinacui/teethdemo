#ifndef FACEMOD_ALGCOLLE_VORONOI_DIAGRAM_H
#define FACEMOD_ALGCOLLE_VORONOI_DIAGRAM_H
#include "prereq.h"


#include <boost/polygon/voronoi.hpp>
#include"../DataColle/mesh_object.h"
using boost::polygon::voronoi_builder;
using boost::polygon::voronoi_diagram;

// pre typedef for using voronoi diagram in boost
struct PointVD {
	int a;
	int b;
	PointVD(int x, int y) : a(x), b(y) {}
};

struct SegmentVD {
	PointVD p0;
	PointVD p1;
	SegmentVD(int x1, int y1, int x2, int y2) : p0(x1, y1), p1(x2, y2) {}
};

namespace boost 
{
	namespace polygon 
	{

		template <>
		struct geometry_concept<PointVD> 
		{
			typedef point_concept type;
		};

		template <>
		struct point_traits<PointVD> 
		{
			typedef int coordinate_type;

			static inline coordinate_type get(
				const PointVD& point, orientation_2d orient) {
				return (orient == HORIZONTAL) ? point.a : point.b;
			}
		};

		template <>
		struct geometry_concept<SegmentVD> 
		{
			typedef segment_concept type;
		};

		template <>
		struct segment_traits<SegmentVD> 
		{
			typedef int coordinate_type;
			typedef PointVD point_type;

			static inline point_type get(const SegmentVD& segment, direction_1d dir) 
			{
				return dir.to_int() ? segment.p1 : segment.p0;
			}
		};
	}  // polygon
}  // boost
// pre typedef for using voronoi diagram in boost



class ALGCOLLE_CLASS CVoronoiDiagramWrapper
{
public:
	CVoronoiDiagramWrapper();
	~CVoronoiDiagramWrapper();

	CVoronoiDiagramWrapper(std::vector<OpenMesh::Vec2d> closed_boundary_curve);

	bool constructVoronoiDiagramSegments();
	bool constructVoronoiDiagramPoints();

	void MainGetSkeletonBranches(
		std::vector<std::vector<std::pair<OpenMesh::Vec2d, double> > > &branches_vertices,
		std::vector<std::vector<int> > &branches_associated_v_sites_index);

	void getColoredSkeletalBranches(std::vector<std::pair<OpenMesh::Vec2d, int> > &color_branches);

	void getSkeletalBranches(
		std::vector<std::vector<std::pair<OpenMesh::Vec2d, double> > > &branches_vector_point_2,
		std::vector<std::vector<int> > &branches_associated_v_sites_index);

	void getSkeletalBranches(std::vector<std::vector<std::pair<OpenMesh::Vec2d, double> > > &branches_vector_point_2);

	void getMaxLengthBranch(double &accum_length,
		std::vector<std::vector<const boost::polygon::voronoi_edge<double>*> > &accum_edge_collector);


private:
	std::vector<SegmentVD> vd_segments_;
	std::vector<PointVD> vd_sites_;

	std::vector<OpenMesh::Vec2d> bc_curve_;
	voronoi_diagram<double> vd_output_;
	std::vector<OpenMesh::Vec2d> vd_vert_list_;

	std::map<const boost::polygon::voronoi_edge<double>*, bool> vde_removal_map_;
	std::map<const boost::polygon::voronoi_edge<double>*, int> skeletal_branches_;
	std::map<const boost::polygon::voronoi_edge<double>*, std::pair<int, int> > skeletal_branches_id_;
	std::map<const boost::polygon::voronoi_vertex<double>*, double> map_ptvertex_radius_;

	//std::vector<std::vector<int> >  branch_adj_map_;


private:
	bool computeVDVRemovalMap(bool is_removal);
	bool computeVDVertList(std::vector<OpenMesh::Vec2d> &edge_skeleton);
	bool computeRemovalVertList(std::vector<OpenMesh::Vec2d> &removal_edge_skeleton);
	bool computeVertexRadius(std::vector<OpenMesh::Vec2d> &vdv_vcoord_list, std::vector<double> &vdv_radius_list);
	// acylic graph ONLY!!!
	bool computeSkeletalBranch(int &id_branch, int seq_in_single_branch, 
		const boost::polygon::voronoi_edge<double>* pt_edge_seed, double ratio_threshold = 0.5);

	// new implementation for cylic graph as well
	void computeSkeletalBranchNew(
		int &id_branch, int seq_in_single_branch,
		std::set<const boost::polygon::voronoi_vertex<double>*> &vertex_visit_set,
		const boost::polygon::voronoi_edge<double>* pt_edge_seed, double ratio_threshold = 0.5);

	void computeMaxLengthBranch(
		const boost::polygon::voronoi_edge<double>* pt_edge_seed,
		double &tmp_length,
		double &accum_length,
		std::vector<const boost::polygon::voronoi_edge<double>*> &tmp_edge_collector,
		std::vector<const boost::polygon::voronoi_edge<double>*> &accum_edge_collector);

	bool computeAllVertexRadius(std::vector<OpenMesh::Vec2d> &vdv_vcoord_list, std::vector<double> &vdv_radius_list);

	// 3 inputs/outputs: skeleton and radius function
	bool computeMedialAxisTransform(std::vector<OpenMesh::Vec2d> &edge_skeleton,
		std::vector<OpenMesh::Vec2d> &vdv_vcoord_list,
		std::vector<double> &vdv_radius_list,
		double ratio_threshold = 0.6);

	// 3 inputs/outputs: skeleton, radius function and removed edges
	bool computeMedialAxisTransform(
		std::vector<OpenMesh::Vec2d> &edge_skeleton,
		std::vector<OpenMesh::Vec2d> &vdv_vcoord_list,
		std::vector<double> &vdv_radius_list,
		std::vector<OpenMesh::Vec2d> &removal_edge_skeleton,
		double ratio_threshold = 0.6);

	//void computeTree();


};
#endif // FACEMOD_ALGCOLLE_VORONOI_DIAGRAM_H
