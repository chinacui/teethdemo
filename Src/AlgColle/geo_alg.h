#ifndef CGEO_ALG_H
#define CGEO_ALG_H
#include"prereq.h"
#include<Eigen/Dense>
#include<iostream>
#include<vector>
#include<map>
#include"../DataColle/geo_primit.h"
#include"../DataColle/mesh_object.h"
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>


class ALGCOLLE_CLASS CGeoAlg
{

protected:
	
public:
	//reomve self intersection face
	static bool SelfIntersectionRemoval(COpenMeshT& mesh);

	//fill holes of the mesh
	static void FillHoles(Eigen::MatrixXd& vertexs, Eigen::MatrixXi& faces);

	//fill hole of the mesh, the largest hole can be remained by option
	static void FillHoles(COpenMeshT &mesh, bool remain_largest=false);

	static void SeparateMeshByVhsLoop(COpenMeshT &mesh, std::vector<COpenMeshT::VertexHandle>&vhs, std::vector< COpenMeshT*>&res_meshes,std::vector<std::map<COpenMeshT::VertexHandle,COpenMeshT::VertexHandle>>&new2orig_vhsmap);
	static void FindEdgePointsPath(CMeshObject &mesh_obj, std::pair<COpenMeshT::FaceHandle, OpenMesh::Vec3d>fpa, std::pair<COpenMeshT::FaceHandle, OpenMesh::Vec3d>fpb, std::vector<std::pair<COpenMeshT::EdgeHandle, OpenMesh::Vec3d>>&res_eps, std::vector<bool>&is_vertex,int recursive_num);
	static void SmoothMesh(COpenMeshT &mesh, std::vector<COpenMeshT::VertexHandle>&roi_vhs,int num);
	static void FindPointsPath(CMeshObject &mesh_obj, std::vector<std::pair<COpenMeshT::FaceHandle, OpenMesh::Vec3d>>face_points_path, std::vector<OpenMesh::Vec3d>&res_points, std::map<int, COpenMeshT::FaceHandle>&res_fhs_map, std::map<int, COpenMeshT::EdgeHandle>&res_ehs_map, std::map<int,COpenMeshT::VertexHandle>&res_vhs_map,std::vector<int>&orig_res_map, bool is_closed = false);
	static void ComputeClosestPoint(OpenMesh::Vec3d p, std::pair<OpenMesh::Vec3d,OpenMesh::Vec3d>seg, OpenMesh::Vec3d &res_p, double &res_dis);
	static bool ComputeClosestPointOfTwoSegments(std::pair<OpenMesh::Vec3d, OpenMesh::Vec3d>sega, std::pair<OpenMesh::Vec3d, OpenMesh::Vec3d>segb, OpenMesh::Vec3d&res_pa, OpenMesh::Vec3d &res_pb);//is they are paralle return false
	static void GetSilhouetteEdges(COpenMeshT&mesh,OpenMesh::Vec3d view_dir, std::vector<std::pair<COpenMeshT::VertexHandle,COpenMeshT::VertexHandle>>&res_edges);
	static void PointSetPCA3D(std::vector<OpenMesh::Vec3d> &pts, OpenMesh::Vec3d&res_mean, std::vector<OpenMesh::Vec3d> &res_eigen_vects, std::vector<double>& res_eigen_values);//row order,sorted from largest to smallest
	static void CutByPlane(COpenMeshT &mesh, CPlane plane, COpenMeshT &res_mesh,std::map<COpenMeshT::VertexHandle,COpenMeshT::VertexHandle>&vid_orig);//approximate return positive side of plane
	static void SeparateDisconnectedParts(COpenMeshT &mesh, std::vector<COpenMeshT*>&res_meshes, std::vector<std::map<COpenMeshT::VertexHandle, COpenMeshT::VertexHandle>>&vid_orig);
	static void GetFeatureGroupsByConnectivity(COpenMeshT&mesh, std::vector<bool>&is_feature,std::vector<std::vector<COpenMeshT::VertexHandle>>&feature_groups);
	static void ComputeMeshPCA(COpenMeshT&mesh, OpenMesh::Vec3d&mean, std::vector<OpenMesh::Vec3d>&res_frame);
	static void SeparateMeshByVertexTag(COpenMeshT &mesh_obj, std::vector<int>&v_tag,std::map<int,COpenMeshT*>&res_meshes, std::map<int,std::map<COpenMeshT::VertexHandle, COpenMeshT::VertexHandle>>&vid_new_2_orig);//v_tag should be continues, start from 0
	static void SeparateMeshByFaceTag(COpenMeshT &mesh, std::vector<int>&f_tag, std::map<int, COpenMeshT*>&res_meshes, std::map<int, std::map<COpenMeshT::VertexHandle, COpenMeshT::VertexHandle>>&vid_orig);//f_tag should be continues, start from 0
	static bool RayMeshIntersection(OpenMesh::Vec3d  source, OpenMesh::Vec3d dir, CMeshObject &mesh_obj, COpenMeshT::FaceHandle & res_fh, OpenMesh::Vec3d &res_bary_coord);
	static bool RayMeshIntersection(OpenMesh::Vec3d  source, OpenMesh::Vec3d dir, CMeshObject &mesh_obj, COpenMeshT::VertexHandle & res_vh);
	static int PickMesh(OpenMesh::Vec3d  source, OpenMesh::Vec3d dir, std::map<int, std::shared_ptr<CMeshObject>>&data_pool,bool is_visiable);
	static void ComputeClosestVertex(OpenMesh::Vec3d source, CMeshObject &mesh_obj, OpenMesh::VertexHandle &res_vh);
	static void SeparateVhsByConnectivity(std::vector<COpenMeshT::VertexHandle>&vhs, COpenMeshT &mesh, std::vector<std::vector<COpenMeshT::VertexHandle>>&res_vhs);
	static double ComputeClosestPoint(OpenMesh::Vec3d source, CMeshObject &mesh_obj, OpenMesh::FaceHandle &res_fh, OpenMesh::Vec3d &res_p);//res_p is in world space
	static bool ComputeGeodesicPath(CMeshObject &mesh_obj,int svid,int tvid, std::vector<COpenMeshT::FaceHandle>&res_fhs, std::vector<OpenMesh::Vec3d>&res_bary_coords);
	static bool ComputeGeodesicPath(CMeshObject &mesh_obj, int svid, int tvid, std::vector<OpenMesh::Vec3d>&path);
	static bool ComputeGeodesicPath(CMeshObject &mesh_obj, COpenMeshT::VertexHandle svh, std::vector<COpenMeshT::VertexHandle>& tvhs, std::vector<OpenMesh::Vec3d>&path);
	static bool ComputeGeodesicPath(CMeshObject &mesh_obj, COpenMeshT::VertexHandle svh, std::vector<COpenMeshT::VertexHandle>& tvhs, std::vector<COpenMeshT::FaceHandle>&res_fhs, std::vector<OpenMesh::Vec3d>&res_bary_coords);
	static bool ComputeGeodesicPath(CMeshObject &mesh_obj, std::vector<COpenMeshT::VertexHandle>&svhs, std::vector<COpenMeshT::VertexHandle>&tvhs, std::vector<std::vector<COpenMeshT::FaceHandle>>&res_fhs, std::vector<std::vector<OpenMesh::Vec3d>>&res_bary_coords);
	static bool ComputeShortestPathAlongEdge(COpenMeshT &mesh, COpenMeshT::VertexHandle svh, COpenMeshT::VertexHandle tvh, std::vector<COpenMeshT::VertexHandle>&path);
	static void ComputeGeodesicDis(CMeshObject &mesh_obj, COpenMeshT::VertexHandle svh, std::vector<COpenMeshT::VertexHandle>&dst_vhs, std::vector<double>&res_dis);
	static void ComputeGeodesicDis(CMeshObject &mesh_obj, std::vector<COpenMeshT::VertexHandle>&svhs, std::vector<double>&res_dis);
	static void InitGeodesicModel(CMeshObject &mesh_obj);
	static void ComputeGradientOfScalarField(COpenMeshT &mesh, Eigen::VectorXd &scalars, Eigen::MatrixXd & res_grad);
	static void SimplifyMesh(OpenMesh::TriMesh_ArrayKernelT<COMTraits> &mesh,int edgenum);
	static void ExtractNRing(COpenMeshT &mesh, COpenMeshT::VertexHandle vh,int ringnum, std::vector<COpenMeshT::VertexHandle>&res_vhs);//except vh
	static void ExtractNeiByDis(COpenMeshT &mesh, COpenMeshT::VertexHandle vh, double ring_dis, std::vector<COpenMeshT::VertexHandle>&res_vhs);
	static void ComputeLocalExtremum(COpenMeshT &mesh, Eigen::VectorXd &scalars, int neighbor_num,std::vector<COpenMeshT::VertexHandle>&res_vhs);
	static void RefineMeshBoundary(COpenMeshT &mesh);
	static double ComputeAverageEdgeLength(COpenMeshT &mesh);
	static void Remeshing(COpenMeshT &mesh);
	static void ComputeScaleAndTransformFrom2SetOfRegisteredPoints(std::vector<OpenMesh::Vec2d>&src_points, std::vector<OpenMesh::Vec2d>&tgt_points, double &res_scale, OpenMesh::Vec2d &res_trans);
	static void ComputeScaleAndTransMatformFrom2SetOfRegisteredPoints(std::vector<OpenMesh::Vec2d>&src_points, std::vector<OpenMesh::Vec2d>&tgt_points, Eigen::Matrix3d &res_mat);
	static void AlphaShape2d(std::vector<OpenMesh::Vec2d>&pts,double alpha,std::vector<std::vector<int>>&bound_vids);
	static void SampleCircle(OpenMesh::Vec3d updir, OpenMesh::Vec3d center,double radius, double degree_step, std::vector<OpenMesh::Vec3d>&res_pts);
	static void ComputeConvexHull(std::vector<OpenMesh::Vec2d>&pts, std::vector<int>&res_convex_ids);
	static void Remeshing(COpenMeshT &mesh, std::vector<COpenMeshT::VertexHandle>&roi);
	static void SplitFaceByPoints(COpenMeshT &mesh, COpenMeshT::FaceHandle fh, std::vector<OpenMesh::Vec3d> &pts,std::vector<OpenMesh::VertexHandle>&new_vhs=std::vector<OpenMesh::VertexHandle>());
	static void SplitFaceByVhs(COpenMeshT &mesh, COpenMeshT::FaceHandle fh, std::vector<COpenMeshT::VertexHandle>&vhs, std::vector<std::pair<COpenMeshT::VertexHandle, COpenMeshT::VertexHandle>>&cons,std::map<COpenMeshT::VertexHandle,COpenMeshT::EdgeHandle>&vh_eh_map);//if a vh is on an edge point , it should be added to vh_eh_map,trick: is cons <0  it indicates the orig vertex on triangle
	static void Triangulate(std::vector<Eigen::Vector2d>&boundv, std::vector<Eigen::Vector2d>&inner_pts, Eigen::MatrixXd &res_v, Eigen::MatrixXi &res_f,std::vector<int>&bound_vids=std::vector<int>(),std::vector<int>&inner_vids = std::vector<int>(),double max_edge_len=-1);
	static void DelaunayTriangulation(std::vector<Eigen::Vector2d>&points, Eigen::MatrixXd &res_v, Eigen::MatrixXi &res_f, std::vector<int>&orig_vidmap = std::vector<int>());//res_v: 2d vertice
	static void DelaunayTriangulation3d(std::vector<Eigen::Vector3d>&points, Eigen::MatrixXd &res_v, Eigen::MatrixXi &res_f, std::vector<int>&orig_vidmap = std::vector<int>());//res_v:3d vertice
	static void ConstrainedTriangulation(std::vector<Eigen::Vector2d>&points,std::vector<std::pair<int,int>>&cons ,Eigen::MatrixXd &res_v, Eigen::MatrixXi &res_f, std::vector<int>&orig_vidmap = std::vector<int>());//res_v: 2d vertice)
	static void ConstrainedTriangulation3d(std::vector<Eigen::Vector3d>&points, std::vector<std::pair<int, int>>&cons, Eigen::MatrixXd &res_v, Eigen::MatrixXi &res_f, std::vector<int>&res_orig_vidmap);
	static void ComputeLocalExtremum(COpenMeshT &mesh, std::vector<double> &scalars, int neighbor_num, std::vector<COpenMeshT::VertexHandle>&res_vhs);
	static void ComputeBaryCoord(COpenMeshT&mesh,COpenMeshT::FaceHandle  fh, OpenMesh::Vec3d p, double &w0, double &w1, double &w2);
	static void ComputeBaryCoord(OpenMesh::Vec3d a, OpenMesh::Vec3d b, OpenMesh::Vec3d c, OpenMesh::Vec3d p, double &res_w0,  double &res_w1, double &res_w2);
	static void ComputeBaryCoord(OpenMesh::Vec2d a, OpenMesh::Vec2d b, OpenMesh::Vec2d c, OpenMesh::Vec2d p, double &res_w0,  double &res_w1, double &res_w2);
	static void LaplacianSmooth(COpenMeshT &mesh, int epoch, double step);
	static void GetOrderedRegionBound(COpenMeshT&mesh, std::vector<int>&region_tags,COpenMeshT::VertexHandle start_vh, std::vector<COpenMeshT::VertexHandle>&res_bound_vhs);
};
#endif