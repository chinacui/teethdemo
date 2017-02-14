#include "compute_voronoi_diagram.h"
#include "curve_base_alg.h"
#include"geo_base_alg.h"
CVoronoiDiagramWrapper::CVoronoiDiagramWrapper()
{
}


CVoronoiDiagramWrapper::~CVoronoiDiagramWrapper()
{
}

CVoronoiDiagramWrapper::CVoronoiDiagramWrapper(std::vector<OpenMesh::Vec2d> bc_curve)
{
	bc_curve_ = bc_curve;
	
	// convert input bc_curve to SegmentVD structure
	for (int i = 0; i < bc_curve_.size() - 1; i++)
	{
		// scale to resist the double2int loss
		PointVD p1(int(bc_curve_[i][0] * 100), int(bc_curve_[i][1] * 100));
		PointVD p2(int(bc_curve_[i + 1][0] * 100), int(bc_curve_[i + 1][1] * 100));

		vd_segments_.push_back(SegmentVD(p1.a, p1.b, p2.a, p2.b));
	}

	// convert input bc_curve to PointVD structure
	for (int i = 0; i < bc_curve_.size(); i++)
	{
		// scale to resist the double2int loss
		PointVD p1(int(bc_curve_[i][0] * 100), int(bc_curve_[i][1] * 100));

		vd_sites_.push_back(PointVD(p1.a, p1.b));
	}
}

void CVoronoiDiagramWrapper::getMaxLengthBranch(double & accum_length, std::vector<std::vector<const boost::polygon::voronoi_edge<double>*>>& accum_edge_collector)
{
	const boost::polygon::voronoi_edge<double>* pt_edge_seed = skeletal_branches_id_.begin()->first;

	double tmp_length = 0;
	double accum_length_pos = 0;
	double accum_length_neg = 0;
	
	std::vector<const boost::polygon::voronoi_edge<double>*> tmp_edge_collector;
	std::vector<const boost::polygon::voronoi_edge<double>*> accum_edge_collector_pos;
	std::vector<const boost::polygon::voronoi_edge<double>*> accum_edge_collector_neg;
	
	//std::cout << "positive edge" << std::endl;
	computeMaxLengthBranch(pt_edge_seed,
		tmp_length, accum_length_pos,
		tmp_edge_collector, accum_edge_collector_pos);

	std::cout << "final accum_edge size " << accum_edge_collector_pos.size() << std::endl;
	//std::cout << "final accum_length " << accum_length_pos << std::endl;

	//std::cout << "negative edge" << std::endl;
	computeMaxLengthBranch(pt_edge_seed->twin(),
		tmp_length, accum_length_neg,
		tmp_edge_collector, accum_edge_collector_neg);
	std::cout << "final accum_edge size " << accum_edge_collector_neg.size() << std::endl;
	//std::cout << "final accum_length " << accum_length_neg << std::endl;

	accum_edge_collector.push_back(accum_edge_collector_pos);
	std::vector<const boost::polygon::voronoi_edge<double>*> vectmp;
	for (int i = 0; i < accum_edge_collector_neg.size(); i++)
	{
		vectmp.push_back(accum_edge_collector_neg[i]);
	}
	accum_edge_collector.push_back(vectmp);
	std::cout << "final accum_length " << accum_length_neg + accum_length_pos << std::endl;

	return;
}

void CVoronoiDiagramWrapper::getSkeletalBranches(
	std::vector<std::vector<std::pair<OpenMesh::Vec2d, double> > > &branches_vector_point_2,
	std::vector<std::vector<int> > &branches_associated_v_sites_index)
{
	// compute the number of branches
	std::map<const boost::polygon::voronoi_edge<double>*, std::pair<int, int>>::const_iterator itmap = skeletal_branches_id_.begin();
	int max_id_branch, min_id_branch;
	min_id_branch = itmap->second.first;
	max_id_branch = itmap->second.first;
	//itmap++;
	for (; itmap != skeletal_branches_id_.end(); ++itmap)
	{
		int tmp_id_branch = itmap->second.first;
		if (tmp_id_branch < min_id_branch)
		{
			min_id_branch = tmp_id_branch;
		}
		else if (tmp_id_branch > max_id_branch)
		{
			max_id_branch = tmp_id_branch;
		}
	}

	// initialize branches
	std::vector<std::map<int, const boost::polygon::voronoi_edge<double>*> > branches_map;
	for (int i = min_id_branch; i < max_id_branch; i++)
	{
		std::map<int, const boost::polygon::voronoi_edge<double>*> single_branch;
		branches_map.push_back(single_branch);
	}

	// construct the data of branches as a map for intermediate data transition
	itmap = skeletal_branches_id_.begin();
	int id_branch = 0;
	for (; itmap != skeletal_branches_id_.end(); ++itmap)
	{
		int tmp_id_branch = itmap->second.first > 0 ? itmap->second.first - min_id_branch - 1 : itmap->second.first - min_id_branch;
		int tmp_seq_in_single_branch = itmap->second.second;
		if (branches_map[tmp_id_branch].find(tmp_seq_in_single_branch) != branches_map[tmp_id_branch].end())
		{
			printf(" Oops! some thing wrong... (%d, %d)\n", tmp_id_branch, tmp_seq_in_single_branch);
		}
		branches_map[tmp_id_branch][tmp_seq_in_single_branch] = itmap->first;
	}

	// sort each branch and output as vertices (OpenMesh::Vec2d)
	//std::vector<std::vector<OpenMesh::Vec2d> > branches_vector_point_2;
	//std::vector<std::vector<int> > branches_associated_v_sites_index;

	for (int i = 0; i < branches_map.size(); i++)
	{
		std::vector<std::pair<OpenMesh::Vec2d, double> > tmp_single_branch_vector;
		std::vector<int> tmp_b2sid; // branch-associated v_sites id

		std::map<int, const boost::polygon::voronoi_edge<double>*> tmp_single_branch_map = branches_map[i];
		// pruning
		if (tmp_single_branch_map.size() == 1)
		{
			double r00 = map_ptvertex_radius_[tmp_single_branch_map.begin()->second->vertex0()];
			double r01 = map_ptvertex_radius_[tmp_single_branch_map.begin()->second->vertex1()];
			double ratio_pruning = r00 < r01 ? (r00 / r01) : (r01 / r00);
			std::cout << "ratio_pruning = " << ratio_pruning << std::endl;
			if (ratio_pruning > 0.3)
			{
				// branches
				OpenMesh::Vec2d p0 = OpenMesh::Vec2d(tmp_single_branch_map.begin()->second->vertex0()->x() / 100.0, tmp_single_branch_map[0]->vertex0()->y() / 100.0);
				tmp_single_branch_vector.push_back(std::pair<OpenMesh::Vec2d, double>(p0, r00));
				OpenMesh::Vec2d p1 = OpenMesh::Vec2d(tmp_single_branch_map.begin()->second->vertex1()->x() / 100.0, tmp_single_branch_map[0]->vertex1()->y() / 100.0);
				tmp_single_branch_vector.push_back(std::pair<OpenMesh::Vec2d, double>(p1, r01));
				branches_vector_point_2.push_back(tmp_single_branch_vector);

				// branch-associated v_site id
				int siteid = tmp_single_branch_map.begin()->second->cell()->source_index(); // branch-associated sites
				tmp_b2sid.push_back(siteid);
				siteid = tmp_single_branch_map.begin()->second->twin()->cell()->source_index();
				tmp_b2sid.push_back(siteid);
				branches_associated_v_sites_index.push_back(tmp_b2sid);
			}
		}
		else
		{
			int j = 0;
			for (; j < tmp_single_branch_map.size(); j++)
			{
				// branches
				OpenMesh::Vec2d p = OpenMesh::Vec2d(tmp_single_branch_map[j]->vertex0()->x() / 100.0, tmp_single_branch_map[j]->vertex0()->y() / 100.0);
				double r = map_ptvertex_radius_[tmp_single_branch_map[j]->vertex0()];
				tmp_single_branch_vector.push_back(std::pair<OpenMesh::Vec2d, double>(p, r));

				// branch-associated sites
				int siteid = tmp_single_branch_map[j]->cell()->source_index();
				tmp_b2sid.push_back(siteid);
				siteid = tmp_single_branch_map[j]->twin()->cell()->source_index();
				tmp_b2sid.push_back(siteid);
			}
			OpenMesh::Vec2d p = OpenMesh::Vec2d(tmp_single_branch_map[j - 1]->vertex1()->x() / 100.0, tmp_single_branch_map[j - 1]->vertex1()->y() / 100.0);
			double r = map_ptvertex_radius_[tmp_single_branch_map[j - 1]->vertex1()];
			tmp_single_branch_vector.push_back(std::pair<OpenMesh::Vec2d, double>(p, r));

			// branch edges
			branches_vector_point_2.push_back(tmp_single_branch_vector);

			// branch-associated v_site id
			branches_associated_v_sites_index.push_back(tmp_b2sid);
		}
	}
}

void CVoronoiDiagramWrapper::getSkeletalBranches(std::vector<std::vector<std::pair<OpenMesh::Vec2d, double> > > &branches_vector_point_2)
{
	// compute the number of branches
	std::map<const boost::polygon::voronoi_edge<double>*, std::pair<int, int>>::const_iterator itmap = skeletal_branches_id_.begin();
	int max_id_branch, min_id_branch;
	min_id_branch = itmap->second.first;
	max_id_branch = itmap->second.first;
	//itmap++;
	for (; itmap != skeletal_branches_id_.end(); ++itmap)
	{
		int tmp_id_branch = itmap->second.first;
		if (tmp_id_branch < min_id_branch)
		{
			min_id_branch = tmp_id_branch;
		}
		else if (tmp_id_branch > max_id_branch)
		{
			max_id_branch = tmp_id_branch;
		}
	}

	// initialize branches
	std::vector<std::map<int, const boost::polygon::voronoi_edge<double>*> > branches_map;
	for (int i = min_id_branch; i < max_id_branch; i++)
	{
		std::map<int, const boost::polygon::voronoi_edge<double>*> single_branch;
		branches_map.push_back(single_branch);
	}

	// construct the data of branches as a map for intermediate data transition
	itmap = skeletal_branches_id_.begin();
	int id_branch = 0;
	for (; itmap != skeletal_branches_id_.end(); ++itmap)
	{
		int tmp_id_branch = itmap->second.first > 0 ? itmap->second.first - min_id_branch - 1 : itmap->second.first - min_id_branch;
		int tmp_seq_in_single_branch = itmap->second.second;
		if (branches_map[tmp_id_branch].find(tmp_seq_in_single_branch) != branches_map[tmp_id_branch].end())
		{
			printf(" Oops! some thing wrong... (%d, %d)\n", tmp_id_branch, tmp_seq_in_single_branch);
		}
		branches_map[tmp_id_branch][tmp_seq_in_single_branch] = itmap->first;
	}

	// sort each branch and output as vertices (OpenMesh::Vec2d)
	//std::vector<std::vector<OpenMesh::Vec2d> > branches_vector_point_2;

	for (int i = 0; i < branches_map.size(); i++)
	{
		std::vector<std::pair<OpenMesh::Vec2d, double> > tmp_single_branch_vector;
		std::vector<int> tmp_b2sid; // branch-associated v_sites id
		std::map<int, const boost::polygon::voronoi_edge<double>*> tmp_single_branch_map = branches_map[i];

		// pruning
		if (tmp_single_branch_map.size() == 0)
		{
			std::cout << "tmp_single_branch_map.size() = " << tmp_single_branch_map.size() << std::endl;
			continue;
		}
		if (tmp_single_branch_map.size() == 1)
		{
			std::cout << "tmp_single_branch_map.size() = " << tmp_single_branch_map.size() << std::endl;
			double r00 = map_ptvertex_radius_[tmp_single_branch_map[0]->vertex0()];
			double r01 = map_ptvertex_radius_[tmp_single_branch_map[0]->vertex1()];
			double ratio_pruning = r00 < r01 ? (r00 / r01) : (r01 / r00);
			std::cout << "ratio_pruning = " << ratio_pruning << std::endl;
			if (ratio_pruning > 0.3)
			{
				OpenMesh::Vec2d p0 = OpenMesh::Vec2d(tmp_single_branch_map[0]->vertex0()->x() / 100.0, tmp_single_branch_map[0]->vertex0()->y() / 100.0);
				tmp_single_branch_vector.push_back(std::pair<OpenMesh::Vec2d, double>(p0, r00));

				OpenMesh::Vec2d p1 = OpenMesh::Vec2d(tmp_single_branch_map[0]->vertex1()->x() / 100.0, tmp_single_branch_map[0]->vertex1()->y() / 100.0);
				tmp_single_branch_vector.push_back(std::pair<OpenMesh::Vec2d, double>(p1, r01));
				
				branches_vector_point_2.push_back(tmp_single_branch_vector);
			}
		}
		else
		{
			int j = 0;
			//std::cout << "j = " << j << std::endl;
			std::cout << "tmp_single_branch_map.size() = " << tmp_single_branch_map.size() << std::endl;
			for (; j < tmp_single_branch_map.size(); j++)
			{
				// branches
				OpenMesh::Vec2d p = OpenMesh::Vec2d(tmp_single_branch_map[j]->vertex0()->x() / 100.0, tmp_single_branch_map[j]->vertex0()->y() / 100.0);
				double r = map_ptvertex_radius_[tmp_single_branch_map[j]->vertex0()];
				tmp_single_branch_vector.push_back(std::pair<OpenMesh::Vec2d, double>(p, r));
			}
			//std::cout << "j = " << j << std::endl;
			OpenMesh::Vec2d p = OpenMesh::Vec2d(tmp_single_branch_map[j - 1]->vertex1()->x() / 100.0, tmp_single_branch_map[j - 1]->vertex1()->y() / 100.0);
			double r = map_ptvertex_radius_[tmp_single_branch_map[j - 1]->vertex1()];
			tmp_single_branch_vector.push_back(std::pair<OpenMesh::Vec2d, double>(p, r));

			// branch edges
			branches_vector_point_2.push_back(tmp_single_branch_vector);
		}
	}

	//for (int i = 0; i < branches_map.size(); i++)
	//{
	//	std::vector<std::pair<OpenMesh::Vec2d, double> > tmp_single_branch_vector;
	//	std::map<int, const boost::polygon::voronoi_edge<double>*> tmp_single_branch_map = branches_map[i];
	//	int j = 0;
	//	for (; j < tmp_single_branch_map.size(); j++)
	//	{
	//		// branches
	//		OpenMesh::Vec2d p = OpenMesh::Vec2d(tmp_single_branch_map[j]->vertex0()->x() / 100.0, tmp_single_branch_map[j]->vertex0()->y() / 100.0);
	//		double r = map_ptvertex_radius_[tmp_single_branch_map[j]->vertex0()];
	//		tmp_single_branch_vector.push_back(std::pair<OpenMesh::Vec2d, double>(p, r));
	//	}
	//	OpenMesh::Vec2d p = OpenMesh::Vec2d(tmp_single_branch_map[j - 1]->vertex1()->x() / 100.0, tmp_single_branch_map[j - 1]->vertex1()->y() / 100.0);
	//	double r = map_ptvertex_radius_[tmp_single_branch_map[j - 1]->vertex1()];
	//	tmp_single_branch_vector.push_back(std::pair<OpenMesh::Vec2d, double>(p, r));

	//	// branch edges
	//	branches_vector_point_2.push_back(tmp_single_branch_vector);
	//}
}

void CVoronoiDiagramWrapper::getColoredSkeletalBranches(std::vector<std::pair<OpenMesh::Vec2d, int> > &color_branches)
{
	std::map<const boost::polygon::voronoi_edge<double>*, int>::const_iterator itmap = skeletal_branches_.begin();
	for (; itmap != skeletal_branches_.end(); ++itmap)
	{
		std::pair<OpenMesh::Vec2d, int> tmp_pair_1
		= std::pair<OpenMesh::Vec2d, int>(OpenMesh::Vec2d(itmap->first->vertex0()->x() / 100.0, itmap->first->vertex0()->y() / 100.0), itmap->second);
		color_branches.push_back(tmp_pair_1);
		std::pair<OpenMesh::Vec2d, int> tmp_pair_2
			= std::pair<OpenMesh::Vec2d, int>(OpenMesh::Vec2d(itmap->first->vertex1()->x() / 100.0, itmap->first->vertex1()->y() / 100.0), itmap->second);
		color_branches.push_back(tmp_pair_2);
	}
}

bool CVoronoiDiagramWrapper::constructVoronoiDiagramSegments()
{
	// construct the vd
	boost::polygon::construct_voronoi(vd_segments_.begin(), vd_segments_.end(), &vd_output_);
	if (vd_output_.num_edges() > 0)
	{
		return true;
	}
	return false;
}

bool CVoronoiDiagramWrapper::constructVoronoiDiagramPoints()
{
	// construct the vd
	boost::polygon::construct_voronoi(vd_sites_.begin(), vd_sites_.end(), &vd_output_);
	if (vd_output_.num_edges() > 0)
	{
		return true;
	}
	return false;
}

void CVoronoiDiagramWrapper::MainGetSkeletonBranches(
	std::vector<std::vector<std::pair<OpenMesh::Vec2d, double> > > &branches_vertices,
	std::vector<std::vector<int> > &branches_associated_v_sites_index)
{
	double ratio_threshold = 0.01;
	std::vector<OpenMesh::Vec2d> vdv_vcoord_list;
	std::vector<double> vdv_radius_list;
	std::vector<OpenMesh::Vec2d> edge_skeleton;
	std::vector<OpenMesh::Vec2d> removal_edge_skeleton;

	computeMedialAxisTransform(edge_skeleton, vdv_vcoord_list, vdv_radius_list, removal_edge_skeleton, ratio_threshold);
	getSkeletalBranches(branches_vertices, branches_associated_v_sites_index);



}

bool CVoronoiDiagramWrapper::computeMedialAxisTransform(
	std::vector<OpenMesh::Vec2d> &edge_skeleton,
	std::vector<OpenMesh::Vec2d> &vdv_vcoord_list,
	std::vector<double> &vdv_radius_list,
	double ratio_threshold)
{
	computeVDVRemovalMap(true);
	computeVDVertList(edge_skeleton);
	computeVertexRadius(vdv_vcoord_list, vdv_radius_list);
	return true;
};

bool CVoronoiDiagramWrapper::computeMedialAxisTransform(
	std::vector<OpenMesh::Vec2d> &edge_skeleton,
	std::vector<OpenMesh::Vec2d> &vdv_vcoord_list,
	std::vector<double> &vdv_radius_list,
	std::vector<OpenMesh::Vec2d> &removal_edge_skeleton,
	double ratio_threshold)
{
	computeVDVRemovalMap(true);
	computeVDVertList(edge_skeleton);
	computeVertexRadius(vdv_vcoord_list, vdv_radius_list);
	computeRemovalVertList(removal_edge_skeleton);

	const boost::polygon::voronoi_edge<double>* pt_edge_seed;
	std::map<const boost::polygon::voronoi_edge<double>*, bool>::const_iterator itmap = vde_removal_map_.begin();
	for (; itmap != vde_removal_map_.end(); ++itmap)
	{
		if (itmap->second)
		{
			pt_edge_seed = itmap->first;
			break;
		}
	}

	std::set<const boost::polygon::voronoi_vertex<double>*> vertex_visit_set;
	int starting_id_p = 1;
	int starting_id_n = -1;
	computeSkeletalBranchNew(starting_id_p, 0, vertex_visit_set, pt_edge_seed->twin(), ratio_threshold);
	computeSkeletalBranchNew(starting_id_n, 0, vertex_visit_set, pt_edge_seed, ratio_threshold);

	return true;
};

bool CVoronoiDiagramWrapper::computeVDVRemovalMap(bool is_removal)
{
	voronoi_diagram<double>::const_edge_iterator ite = vd_output_.edges().begin();

	double dist_threshold = 0.03;
	
	for (; ite != vd_output_.edges().end(); ++ite)
	{
		if (ite->is_finite())
		{
			const boost::polygon::voronoi_edge<double>* itet = ite->twin()->twin();
			const boost::polygon::voronoi_edge<double>* itet_twin = ite->twin();

			const boost::polygon::voronoi_vertex<double> *itv0 = ite->vertex0();
			const boost::polygon::voronoi_vertex<double> *itv1 = ite->vertex1();

			double p0x = itv0->x() / 100.0;
			double p0y = itv0->y() / 100.0;

			double p1x = itv1->x() / 100.0;
			double p1y = itv1->y() / 100.0;

			// outside the closed curve region
			vde_removal_map_[itet] = true;
			vde_removal_map_[itet_twin] = true;

			if (is_removal)
			{
				// remove secondary edges ONLY for segment inputs
				if (ite->is_secondary())
				{
					vde_removal_map_[itet] = false;
					vde_removal_map_[itet_twin] = false;
				}

				// remove edges that intersect with the boundary curve
				if (vde_removal_map_[itet])
				{
					bool b_InRegion = CCurveBaseAlg::IsInRegion(bc_curve_, OpenMesh::Vec2d(p0x, p0y)) && CCurveBaseAlg::IsInRegion(bc_curve_, OpenMesh::Vec2d(p1x, p1y));
					vde_removal_map_[itet] = b_InRegion;
					vde_removal_map_[itet_twin] = b_InRegion;
				}

				// remove edges whose vertices is close to boundary

			}
			
		}	
	}
	
	return true;
}


bool CVoronoiDiagramWrapper::computeVDVertList(std::vector<OpenMesh::Vec2d> &vd_vert_list_new)
{
	// update vd_vert_list
	std::map<const boost::polygon::voronoi_edge<double>*, bool>::iterator itmap = vde_removal_map_.begin();
	for (; itmap != vde_removal_map_.end(); ++itmap)
	{
		if (itmap->second)
		{
			const boost::polygon::voronoi_edge<double>* ite = itmap->first;

			double p0x = ite->vertex0()->x()/100.0;
			double p0y = ite->vertex0()->y()/100.0;

			double p1x = ite->vertex1()->x()/100.0;
			double p1y = ite->vertex1()->y()/100.0;

			vd_vert_list_new.push_back(OpenMesh::Vec2d(p0x, p0y));
			vd_vert_list_new.push_back(OpenMesh::Vec2d(p1x, p1y));
		}
	}

	return true;
}


bool CVoronoiDiagramWrapper::computeRemovalVertList(std::vector<OpenMesh::Vec2d> &vd_vert_list_removal)
{
	// update vd_vert_list
	std::map<const boost::polygon::voronoi_edge<double>*, bool>::iterator itmap = vde_removal_map_.begin();
	for (; itmap != vde_removal_map_.end(); ++itmap)
	{
		if (!itmap->second)
		{
			const boost::polygon::voronoi_edge<double>* ite = itmap->first;
			if (ite->is_finite())
			{
				double p0x = ite->vertex0()->x() / 100.0;
				double p0y = ite->vertex0()->y() / 100.0;

				double p1x = ite->vertex1()->x() / 100.0;
				double p1y = ite->vertex1()->y() / 100.0;

				vd_vert_list_removal.push_back(OpenMesh::Vec2d(p0x, p0y));
				vd_vert_list_removal.push_back(OpenMesh::Vec2d(p1x, p1y));
			}

		}
	}

	return true;
}


bool CVoronoiDiagramWrapper::computeVertexRadius(std::vector<OpenMesh::Vec2d> &vdv_vcoord_list, std::vector<double> &vdv_radius_list)
{
	if (vde_removal_map_.empty())
	{
		return false;
	}

	// iterate the edges
	voronoi_diagram<double>::const_edge_iterator ite = vd_output_.edges().begin();
	for (; ite != vd_output_.edges().end(); ++ite)
	{
		const boost::polygon::voronoi_edge<double>* pt_edge = ite->twin()->twin();

		if (!vde_removal_map_[pt_edge])
		{
			continue;
		}

		const boost::polygon::voronoi_vertex<double>* pt_vert_0 = pt_edge->vertex0();
		const boost::polygon::voronoi_vertex<double>* pt_vert_1 = pt_edge->vertex1();

		OpenMesh::Vec2d p0 = OpenMesh::Vec2d(pt_vert_0->x() / 100.0, pt_vert_0->y() / 100.0);
		OpenMesh::Vec2d p1 = OpenMesh::Vec2d(pt_vert_1->x() / 100.0, pt_vert_1->y() / 100.0);

		double min_dist = std::numeric_limits<double>::max();


		// transverse the starting point
		min_dist = std::numeric_limits<double>::max();
		std::vector<const boost::polygon::voronoi_edge<double>* > search_edges;
		search_edges.push_back(pt_edge);
		search_edges.push_back(pt_edge->prev());
		search_edges.push_back(pt_edge->twin());
		search_edges.push_back(pt_edge->twin()->next());
		for (int i = 0; i < 4; i++)
		{
			int qid = search_edges[i]->cell()->source_index();
			double dist_prev_tmp;
			double dist_next_tmp;
			if (qid == 0)
			{
				dist_prev_tmp = CGeoBaseAlg::Min_Dist_Point2_to_Line2(p0, bc_curve_[qid], bc_curve_.back());
				dist_next_tmp = CGeoBaseAlg::Min_Dist_Point2_to_Line2(p0, bc_curve_[qid], bc_curve_[qid + 1]);
			}
			else if (qid == bc_curve_.size() - 1)
			{
				dist_prev_tmp = CGeoBaseAlg::Min_Dist_Point2_to_Line2(p0, bc_curve_[qid], bc_curve_[qid - 1]);
				dist_next_tmp = CGeoBaseAlg::Min_Dist_Point2_to_Line2(p0, bc_curve_[qid], bc_curve_.front());
			}
			else
			{
				dist_prev_tmp = CGeoBaseAlg::Min_Dist_Point2_to_Line2(p0, bc_curve_[qid], bc_curve_[qid - 1]);
				dist_next_tmp = CGeoBaseAlg::Min_Dist_Point2_to_Line2(p0, bc_curve_[qid], bc_curve_[qid + 1]);
			}
			if (dist_prev_tmp < min_dist)
			{
				min_dist = dist_prev_tmp;
			}
			if (dist_next_tmp < min_dist)
			{
				min_dist = dist_next_tmp;
			}
		}
		// add to maps
		if (map_ptvertex_radius_.end() != map_ptvertex_radius_.find(pt_edge->vertex0()))
		{
			if (min_dist < map_ptvertex_radius_.find(pt_edge->vertex0())->second)
			{
				map_ptvertex_radius_[pt_edge->vertex0()] = min_dist;
			}
		}
		else
		{
			map_ptvertex_radius_[pt_edge->vertex0()] = min_dist;
		}
		// vdv_vcoord_list.push_back(p0);
		// vdv_radius_list.push_back(min_dist);

		// transverse the end point
		min_dist = std::numeric_limits<double>::max();
		search_edges.clear();
		search_edges.push_back(pt_edge);
		search_edges.push_back(pt_edge->next());
		search_edges.push_back(pt_edge->twin());
		search_edges.push_back(pt_edge->twin()->prev());
		for (int i = 0; i < 4; i++)
		{
			int qid = search_edges[i]->cell()->source_index();
			double dist_prev_tmp;
			double dist_next_tmp;
			if (qid == 0)
			{
				dist_prev_tmp = CGeoBaseAlg::Min_Dist_Point2_to_Line2(p1, bc_curve_[qid], bc_curve_.back());
				dist_next_tmp = CGeoBaseAlg::Min_Dist_Point2_to_Line2(p1, bc_curve_[qid], bc_curve_[qid + 1]);
			}
			else if (qid == bc_curve_.size() - 1)
			{
				dist_prev_tmp = CGeoBaseAlg::Min_Dist_Point2_to_Line2(p1, bc_curve_[qid], bc_curve_[qid - 1]);
				dist_next_tmp = CGeoBaseAlg::Min_Dist_Point2_to_Line2(p1, bc_curve_[qid], bc_curve_.front());
			}
			else
			{
				dist_prev_tmp = CGeoBaseAlg::Min_Dist_Point2_to_Line2(p1, bc_curve_[qid], bc_curve_[qid - 1]);
				dist_next_tmp = CGeoBaseAlg::Min_Dist_Point2_to_Line2(p1, bc_curve_[qid], bc_curve_[qid + 1]);
			}
			if (dist_prev_tmp < min_dist)
			{
				min_dist = dist_prev_tmp;
			}
			if (dist_next_tmp < min_dist)
			{
				min_dist = dist_next_tmp;
			}
		}
		// add to maps
		if (map_ptvertex_radius_.end() != map_ptvertex_radius_.find(pt_edge->vertex1()))
		{
			if (min_dist < map_ptvertex_radius_.find(pt_edge->vertex1())->second)
			{
				map_ptvertex_radius_[pt_edge->vertex1()] = min_dist;
			}
		}
		else
		{
			map_ptvertex_radius_[pt_edge->vertex1()] = min_dist;
		}

	}
	std::map<const boost::polygon::voronoi_vertex<double>*, double>::const_iterator itmap = map_ptvertex_radius_.begin();
	for (; itmap != map_ptvertex_radius_.end(); ++itmap)
	{
		vdv_vcoord_list.push_back(OpenMesh::Vec2d(itmap->first->x()/100.0, itmap->first->y()/100.0));
		vdv_radius_list.push_back(itmap->second);
	}
	return true;
}

bool CVoronoiDiagramWrapper::computeAllVertexRadius(std::vector<OpenMesh::Vec2d> &vdv_vcoord_list, std::vector<double> &vdv_radius_list)
{

	// iterate the vertices
	voronoi_diagram<double>::const_vertex_iterator itv = vd_output_.vertices().begin();
	for (; itv != vd_output_.vertices().end(); ++itv)
	{
		OpenMesh::Vec2d p = OpenMesh::Vec2d(itv->x()/100.0, itv->y()/100.0);
		double min_dist = std::numeric_limits<double>::max();
		const boost::polygon::voronoi_edge<double>* pt_edge = itv->incident_edge();
		// transverse the edges around itv to find the minimal radius
		std::vector<const boost::polygon::voronoi_edge<double>* > search_edges;
		search_edges.push_back(pt_edge);
		search_edges.push_back(pt_edge->prev());
		search_edges.push_back(pt_edge->twin());
		search_edges.push_back(pt_edge->twin()->next());
		for (int i = 0; i < 4; i++)
		{
			int qid = search_edges[i]->cell()->source_index();
			double dist_prev_tmp;
			double dist_next_tmp;
			if (qid == 0)
			{
				dist_prev_tmp = CGeoBaseAlg::Min_Dist_Point2_to_Line2(p, bc_curve_[qid], bc_curve_.back());
				dist_next_tmp = CGeoBaseAlg::Min_Dist_Point2_to_Line2(p, bc_curve_[qid], bc_curve_[qid + 1]);
			}
			else if (qid == bc_curve_.size() - 1)
			{
				dist_prev_tmp = CGeoBaseAlg::Min_Dist_Point2_to_Line2(p, bc_curve_[qid], bc_curve_[qid - 1]);
				dist_next_tmp = CGeoBaseAlg::Min_Dist_Point2_to_Line2(p, bc_curve_[qid], bc_curve_.front());
			}
			else
			{
				dist_prev_tmp = CGeoBaseAlg::Min_Dist_Point2_to_Line2(p, bc_curve_[qid], bc_curve_[qid - 1]);
				dist_next_tmp = CGeoBaseAlg::Min_Dist_Point2_to_Line2(p, bc_curve_[qid], bc_curve_[qid + 1]);
			}
			if (dist_prev_tmp < min_dist)
			{
				min_dist = dist_prev_tmp;
			}
			if (dist_next_tmp < min_dist)
			{
				min_dist = dist_next_tmp;
			}
		}
		// add to maps
		map_ptvertex_radius_[pt_edge->vertex0()] = min_dist;
		vdv_vcoord_list.push_back(p);
		vdv_radius_list.push_back(min_dist);
	}
	return true;
}

void CVoronoiDiagramWrapper::computeSkeletalBranchNew(
	int &id_branch, int seq_in_single_branch,
	std::set<const boost::polygon::voronoi_vertex<double>*> &vertex_visit_set,
	const boost::polygon::voronoi_edge<double>* pt_edge_seed, double ratio_threshold)
{
	// using a set that records the visit of nodes to guide the graph transverse
	if (vertex_visit_set.find(pt_edge_seed->vertex1()) == vertex_visit_set.end())
	{
		vertex_visit_set.insert(pt_edge_seed->vertex1());
	}
	else
	{
		// update skeletal_branches_
		skeletal_branches_[pt_edge_seed] = id_branch;
		skeletal_branches_id_[pt_edge_seed] = std::pair<int, int>(id_branch, seq_in_single_branch);
		return;
	}

	// storage the to-go edges
	std::vector<const boost::polygon::voronoi_edge<double>*> togo_edges;

	const boost::polygon::voronoi_edge<double>* pt_edge_seed_next = pt_edge_seed->next();
	if (pt_edge_seed->twin() != pt_edge_seed_next) // not the same edge with reversed dir
	{
		const boost::polygon::voronoi_edge<double>* pt_edge_seed_cir = pt_edge_seed_next;
		do
		{
			pt_edge_seed_cir = pt_edge_seed_cir->rot_next();
			if (vde_removal_map_[pt_edge_seed_cir] && pt_edge_seed_cir != pt_edge_seed->twin())
			{
				togo_edges.push_back(pt_edge_seed_cir);
			}
		} while (pt_edge_seed_next != pt_edge_seed_cir);
	}
	// update skeletal_branches_
	skeletal_branches_[pt_edge_seed] = id_branch;
	skeletal_branches_id_[pt_edge_seed] = std::pair<int, int>(id_branch, seq_in_single_branch);

	// go to random next edge around current vertex seed
	if (togo_edges.size() > 1)
	{
		for (int i = 0; i < togo_edges.size(); i++)
		{
			if (id_branch > 0)
			{
				id_branch++;
			}
			else if (id_branch < 0)
			{
				id_branch--;
			}
			seq_in_single_branch = 0;
			computeSkeletalBranchNew(id_branch, seq_in_single_branch, vertex_visit_set, togo_edges[i], ratio_threshold);
		}
	}
	else if (!togo_edges.empty())
	{
		//double radius_0, radius_1;
		//radius_0 = map_ptvertex_radius_[togo_edges[0]->vertex0()];
		//radius_1 = map_ptvertex_radius_[togo_edges[0]->vertex1()];
		//double ratio = radius_0 < radius_1 ? radius_0 / radius_1 : radius_1 / radius_0;
		////std::cout << "ratio: " << ratio << std::endl;
		//if (ratio < ratio_threshold)
		//{
		//	if (id_branch > 0)
		//	{
		//		id_branch++;
		//	}
		//	else if (id_branch < 0)
		//	{
		//		id_branch--;
		//	}
		//	seq_in_single_branch = 0;
		//	computeSkeletalBranchNew(id_branch, seq_in_single_branch, vertex_visit_set, togo_edges[0], ratio_threshold);
		//}
		//else
		//{
			seq_in_single_branch++;
			computeSkeletalBranchNew(id_branch, seq_in_single_branch, vertex_visit_set, togo_edges[0], ratio_threshold);
		//}
	}

}


// recursive for acylic graph
bool CVoronoiDiagramWrapper::computeSkeletalBranch(int &id_branch, int seq_in_single_branch,
	const boost::polygon::voronoi_edge<double>* pt_edge_seed, double ratio_threshold)
{
	//std::cout << "id_branch = " << id_branch << std::endl;
	// transverse the edges
	std::vector<const boost::polygon::voronoi_edge<double>*> togo_edges;

	// circulate the next edge to check if it is the only valid edge next to the seed edge; 
	// if so, the end point of the seed edge is a normal point,
	// if not, the end point of the seed edge is a joint point. 
	const boost::polygon::voronoi_edge<double>* pt_edge_seed_next = pt_edge_seed->next();

	if (pt_edge_seed->twin() != pt_edge_seed_next) // not the same edge
	{
		const boost::polygon::voronoi_edge<double>* pt_edge_seed_cir = pt_edge_seed_next;
		int cnt = 0;
		do 
		{
			pt_edge_seed_cir = pt_edge_seed_cir->rot_next();
			if (vde_removal_map_[pt_edge_seed_cir] && pt_edge_seed_cir != pt_edge_seed->twin())
			{
				togo_edges.push_back(pt_edge_seed_cir);
				cnt++;
			}
		} while (pt_edge_seed_next != pt_edge_seed_cir || cnt == 5);
	}
	// update skeletal_branches_
	skeletal_branches_[pt_edge_seed] = id_branch;
	skeletal_branches_id_[pt_edge_seed] = std::pair<int, int>(id_branch, seq_in_single_branch);

	// go to random next edge around current vertex seed
	if (togo_edges.size() > 1)
	{
		for (int i = 0; i < togo_edges.size(); i++)
		{
			if (id_branch > 0)
			{
				id_branch++;
			}
			else if (id_branch < 0)
			{
				id_branch--;
			}
			seq_in_single_branch = 0;
			computeSkeletalBranch(id_branch, seq_in_single_branch, togo_edges[i], ratio_threshold);
		}
	}
	else if (!togo_edges.empty())
	{
		double radius_0, radius_1;
		radius_0 = map_ptvertex_radius_[togo_edges[0]->vertex0()];
		radius_1 = map_ptvertex_radius_[togo_edges[0]->vertex1()];
		double ratio = radius_0 < radius_1 ? radius_0 / radius_1 : radius_1 / radius_0;
		//std::cout << "ratio: " << ratio << std::endl;
		if (ratio < ratio_threshold)
		{
			if (id_branch > 0)
			{
				id_branch++;
			}
			else if (id_branch < 0)
			{
				id_branch--;
			}
			seq_in_single_branch = 0;
			computeSkeletalBranch(id_branch, seq_in_single_branch, togo_edges[0], ratio_threshold);
		}
		else
		{
			seq_in_single_branch++;
			computeSkeletalBranch(id_branch, seq_in_single_branch, togo_edges[0], ratio_threshold);
		}
	}
	
	return true;
}

void CVoronoiDiagramWrapper::computeMaxLengthBranch(
	const boost::polygon::voronoi_edge<double>* pt_edge_seed,
	double &tmp_length,
	double &accum_length,
	std::vector<const boost::polygon::voronoi_edge<double>*> &tmp_edge_collector,
	std::vector<const boost::polygon::voronoi_edge<double>*> &accum_edge_collector)
{
	OpenMesh::Vec2d p0(pt_edge_seed->vertex0()->x() / 100.0, pt_edge_seed->vertex0()->y() / 100.0);
	OpenMesh::Vec2d p1(pt_edge_seed->vertex1()->x() / 100.0, pt_edge_seed->vertex1()->y() / 100.0);
	double length_eh_seed = (p0-p1).length();
	
	// update tmp_length
	tmp_length += length_eh_seed;
	tmp_edge_collector.push_back(pt_edge_seed);

	// storage the to-go edges
	std::vector<const boost::polygon::voronoi_edge<double>*> togo_edges;
	const boost::polygon::voronoi_edge<double>* pt_edge_seed_next = pt_edge_seed->next();
	if (pt_edge_seed->twin() != pt_edge_seed_next || pt_edge_seed->vertex1() != pt_edge_seed_next->vertex1()) // not the same edge with reversed dir
	{
		const boost::polygon::voronoi_edge<double>* pt_edge_seed_cir = pt_edge_seed_next;
		//togo_edges.push_back(pt_edge_seed_cir);
		do
		{
			pt_edge_seed_cir = pt_edge_seed_cir->rot_next();
			// removal map true means not removed
			if (vde_removal_map_[pt_edge_seed_cir] && pt_edge_seed_cir != pt_edge_seed->twin())
			{
				togo_edges.push_back(pt_edge_seed_cir);
			}
		} while (pt_edge_seed_next != pt_edge_seed_cir);
	}

	
	if (togo_edges.empty()) // degree-1 vertex; reach to an extremity
	{
		/*std::cout << "accum_length = " << accum_length << std::endl;
		std::cout << "tmp_length = " << tmp_length << std::endl;*/
		if (accum_length < tmp_length)
		{
			accum_length = tmp_length;
			accum_edge_collector.clear();
			accum_edge_collector = tmp_edge_collector;
			tmp_length = 0;
		}
		return; // go back to previous bifurcation
	}
	else if (togo_edges.size() == 1) // degree-2 vertex; no bifurcation; 
	{
		// go into next recursion
		computeMaxLengthBranch(togo_edges[0], tmp_length, accum_length, tmp_edge_collector, accum_edge_collector);
	}
	else if (togo_edges.size() > 1) // case 2: bifurcation
	{
		// storage some temp information for recursive use
		double tmp_tmp_length = tmp_length;
		std::vector<const boost::polygon::voronoi_edge<double>*> tmp_tmp_edge_collector = tmp_edge_collector;

		// go into next recursion
		for (int i = 0; i < togo_edges.size(); i++)
		{
			// instant storage
			tmp_length = tmp_tmp_length;
			tmp_edge_collector.clear();
			tmp_edge_collector = tmp_tmp_edge_collector;
			
			// compute
			computeMaxLengthBranch(
				togo_edges[i], 
				tmp_length, accum_length,
				tmp_edge_collector, accum_edge_collector);
			//

		}
			
	}



}