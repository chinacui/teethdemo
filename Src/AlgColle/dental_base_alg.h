
#ifndef CDENTAL_BASE_ALG_H
#define CDENTAL_BASE_ALG_H
#include"../DataColle/mesh_object.h"
#include"../DataColle/geo_primit.h"
#include"prereq.h"
#include"curve_base_alg.h"
#include"geo_base_alg.h"
#include"geo_alg.h"
class ALGCOLLE_CLASS CDentalBaseAlg
{
protected:
	class CVertexValue
	{
	public:

		COpenMeshT::VertexHandle vh_;
		double value_;

		CVertexValue(COpenMeshT::VertexHandle vh,double value)
		{
			vh_ = vh;
			value_ = value;
		}
		bool operator<(CVertexValue& b)
		{
			if (value_ > b.value_)
				return true;
			else
				return false;
		}
	};

	class CUnionFind 
	{
	public:
		int num_;
		std::vector<int>roots_;
		CUnionFind(int num):num_(num)
		{
			roots_.resize(num);
			for (int i = 0; i < num; i++)
			{
				roots_[i] = i;
			}
		}
		int Find(int id)
		{
			if (roots_[id] == id)
				return id;
			else
			{
				roots_[id]=Find(roots_[id]);
				return roots_[id];
			}
		}
		void Union(int a, int b)
		{
			int ra = Find(a);
			roots_[b] = ra;
		}
		bool IsRoot(int id)
		{
			if (roots_[id] == id)
				return true;
			else
				return false;
		}
	};
	static void ComputeVertexPenaltyWeight(COpenMeshT&mesh,std::vector<COpenMeshT::VertexHandle>&feature_points, Eigen::VectorXd &mean_curvature,Eigen::VectorXd& res_weight);
	static double ComputePenaltyValue(COpenMeshT&mesh, CPlane plane, std::vector<COpenMeshT::VertexHandle>&feature_points, Eigen::VectorXd &mean_curvature, Eigen::VectorXd& res_weight);

public:
	class CCuttingPath
	{
		double len_=-1;
	public:
		COpenMeshT::VertexHandle start_vh_, end_vh_;
		std::vector<COpenMeshT::FaceHandle> path_fhs_;
		std::vector<OpenMesh::Vec3d> path_barycoords_;
		COpenMeshT* mesh_;
		CCuttingPath CCuttingPath::operator=(CCuttingPath &b)
		{
			mesh_ = b.mesh_;
			start_vh_ = b.start_vh_;
			end_vh_ = b.end_vh_;
			path_fhs_ = b.path_fhs_;
			path_barycoords_ = b.path_barycoords_;
			return *this;
		}
		CCuttingPath()
		{

		}
		CCuttingPath(const CCuttingPath&b)
		{
			mesh_ = b.mesh_;
			start_vh_ = b.start_vh_;
			end_vh_ = b.end_vh_;
			path_fhs_ = b.path_fhs_;
			path_barycoords_ = b.path_barycoords_;
			len_ = b.len_;
		}
		CCuttingPath(COpenMeshT&mesh,COpenMeshT::VertexHandle start_vh, COpenMeshT::VertexHandle end_vh, std::vector<COpenMeshT::FaceHandle>&path_fhs, std::vector<OpenMesh::Vec3d>& path_barycoords)
		{
				mesh_ = &mesh;
			   start_vh_ = start_vh;
			   end_vh_ = end_vh;
			   path_fhs_ = path_fhs;
			   path_barycoords_ = path_barycoords;
			   OpenMesh::Vec3d prev = CGeoBaseAlg::ComputePointFromBaryCoord(mesh, path_fhs_[0], path_barycoords_[0]);
			   for (int i = 1; i < path_fhs_.size(); i++)
			   {
				   OpenMesh::Vec3d pv = CGeoBaseAlg::ComputePointFromBaryCoord(mesh, path_fhs_[i], path_barycoords_[i]);
				   len_ += (pv - prev).length();

			   }
		 };
		void GetPathPoints(std::vector<OpenMesh::Vec3d>&res_path_points)
		{
			res_path_points.clear();
			for (int i = 0; i < path_fhs_.size(); i++)
			{
				OpenMesh::Vec3d pv = CGeoBaseAlg::ComputePointFromBaryCoord(*mesh_, path_fhs_[i], path_barycoords_[i]);
				res_path_points.push_back(pv);
			}
		}
		double GetStraightLength()
		{
			OpenMesh::Vec3d sv = mesh_->point(start_vh_);
			OpenMesh::Vec3d ev = mesh_->point(end_vh_);
			return (sv - ev).length();
		}
		int GetSize()
		{
			return path_fhs_.size();
		}
		double GetLength()
		{
			if (len_ == -1)
			{
				len_ = 0;
				OpenMesh::Vec3d prev = CGeoBaseAlg::ComputePointFromBaryCoord(*mesh_, path_fhs_[0], path_barycoords_[0]);

				for (int i = 1; i < path_fhs_.size(); i++)
				{
					OpenMesh::Vec3d pv = CGeoBaseAlg::ComputePointFromBaryCoord(*mesh_, path_fhs_[i], path_barycoords_[i]);
					len_ += (pv - prev).length();

				}
			}
			return len_;
		}

	};

	static void ComputePCAFrameFromHighCurvaturePoints(COpenMeshT&mesh, double threshold, OpenMesh::Vec3d&mean,std::vector<OpenMesh::Vec3d>&res_frame);
	static void ComputeTeethFeaturePoints(COpenMeshT &mesh, std::vector<COpenMeshT::VertexHandle>&res_vhs);
	static void ComputeTeethFeaturePointsUsingSmoothedMesh(COpenMeshT &mesh, std::vector<COpenMeshT::VertexHandle>&res_vhs);
	static void ComputeGingivaVhs(COpenMeshT &mesh, std::vector<COpenMeshT::VertexHandle>&res_vhs);
	static void ComputeHoleVertMean(COpenMeshT&mesh, OpenMesh::Vec3d &res_mean);
	static void ComputeExtremePointsOfClosedCurve(std::vector<OpenMesh::Vec2d>&curve, std::vector<int>&res_ids);
	static void DetectRootOfTeethSilhouette(std::vector<OpenMesh::Vec2d>&curve,std::vector<int>&root_pids);
	static void ComputeTwoSideBoundsOfToothMesh(COpenMeshT&mesh, std::vector<std::vector<COpenMeshT::VertexHandle>>&res_inside_bounds, std::vector<std::vector<COpenMeshT::VertexHandle>>&res_outside_bounds);
	static void ComputeBoundCuttingPointsOfToothMesh(CMeshObject&mesh, std::vector<std::vector<COpenMeshT::VertexHandle>>&res_inside_vhs, std::vector<std::vector<COpenMeshT::VertexHandle>>&res_outside_vhs);
	static void MergeCuttingPointsByDis(COpenMeshT&mesh, std::vector<std::vector<COpenMeshT::VertexHandle>>&inside_vhs, std::vector<std::vector<COpenMeshT::VertexHandle>>&outside_vhs,double threshold);
	static void ComputeCuttingPath(CMeshObject&meshobj, std::vector<std::vector<COpenMeshT::VertexHandle>>&inside_cutting_vhs, std::vector<std::vector<COpenMeshT::VertexHandle>>&outside_cutting_vhs,std::vector<CCuttingPath>&res_cuttingpath);
	static int TagToothByCuttingPath(COpenMeshT&mesh, std::vector<CCuttingPath>&cutting_pathes,std::vector<int>&res_tags);
	static Eigen::Matrix4d  ComputeTeethLocalCoordinateFromFaPointAndLongAxis(OpenMesh::Vec3d fa_point, OpenMesh::Vec3d mean, OpenMesh::Vec3d longaxis);

	static void PCABasedOrientationCorrection(COpenMeshT& mesh);
	
};
#endif