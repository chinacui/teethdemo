#include"dental_template_fitting.h"
#include"../AdditionalLibs/obb/cpqp_obb_wrapper.h"
#include"../DataColle/cgal_igl_converter.h"
#include"../AlgColle/geo_alg.h"
#include"../AlgColle/geo_base_alg.h"
#include "../AdditionalLibs/Remesher/UniformRemesherT.hh"
#include "../AdditionalLibs/Remesher/AdaptiveRemesherT.hh"
#include "../AdditionalLibs/Remesher/MeshSelectionT.hh"
#include"../DataColle/curve_object.h"
#include"../DataColle/data_pool.h"
#include<igl\writeOBJ.h>
#include"../AlgColle/curve_base_alg.h"
#include"../DataColle/data_io.h"
#include"../DataColle/teeth_template_object.h"
#include"../AlgColle/numerical_base_alg.h"
void CDentalTemplateFitting::RefineFittingTemplate(COpenMeshT &mesh_template, COpenMeshT &mesh_target)
{
	Eigen::MatrixXd temp_v, tgt_v;
	Eigen::MatrixXi temp_f, tgt_f;
	std::vector<COpenMeshT::FaceHandle>temp_idfh_map, tgt_idfh_map;
	CConverter::ConvertFromOpenMeshToIGL(mesh_template, temp_v, temp_f, &temp_idfh_map);
	CConverter::ConvertFromOpenMeshToIGL(mesh_target, tgt_v, tgt_f, &tgt_idfh_map);

	std::cerr << "compute obb" << std::endl;
	CPqpObb tmp_obb(temp_v, temp_f);
	CPqpObb tgt_obb(tgt_v, tgt_f);
	std::cerr << "compute obb finish" << std::endl;
	std::map<OpenMesh::FaceHandle, std::vector<OpenMesh::Vec3d>>temp_face_add_points;

	for (int i = 0; i < temp_f.rows(); i++)
	{
		Eigen::Vector3d center(0, 0, 0);
		for (int j = 0; j < 3; j++)
		{
			int vid = temp_f(i, j);
			center = center + temp_v.row(vid).transpose();
		}
		center =center/3;
		temp_face_add_points[temp_idfh_map[i]] = std::vector<OpenMesh::Vec3d>();
		temp_face_add_points[temp_idfh_map[i]].push_back(OpenMesh::Vec3d(center.x(), center.y(), center.z()));
	}
	
	std::set<OpenMesh::VertexHandle>new_vhs;
	for (auto iter = temp_face_add_points.begin(); iter != temp_face_add_points.end(); iter++)
	{
		std::vector<OpenMesh::VertexHandle>vhs;
		//	CGeoAlg::SplitFaceByPoints(mesh_template, iter->first, iter->second, vhs);
		std::vector<COpenMeshT::VertexHandle>fvhs;
		COpenMeshT::FaceHandle fh = iter->first;
		for (auto viter = mesh_template.fv_begin(fh); viter != mesh_template.fv_end(fh); viter++)
		{
			fvhs.push_back(viter);
		}
		mesh_template.delete_face(fh);
		vhs.push_back(mesh_template.add_vertex(iter->second[0]));

		mesh_template.add_face(vhs[0], fvhs[0], fvhs[1]);
		mesh_template.add_face(vhs[0], fvhs[2], fvhs[0]);
		mesh_template.add_face(vhs[0], fvhs[1], fvhs[2]);

		for (int i = 0; i < vhs.size(); i++)
		{
			new_vhs.insert(vhs[i]);
		}
	}
	/*for (auto iter = new_vhs.begin(); iter != new_vhs.end(); iter++)
	{
		OpenMesh::Vec3d p = mesh_template.point(*iter);
		Eigen::Vector3d ep(p[0], p[1], p[2]);
		CPqpObb::QueryResult res = tgt_obb.QueryClosestPoint(ep);
		Eigen::Vector3d tp = res.closest_pnt_;
		p = OpenMesh::Vec3d(tp.x(), tp.y(), tp.z());
		mesh_template.point(*iter) = p;
	}*/
	mesh_template.request_face_normals();
	mesh_template.request_vertex_normals();
	mesh_template.garbage_collection();
	
	std::cerr << "RefineFittingTemplate finish" << std::endl;

}
void CDentalTemplateFitting::ComputeCrownFrontDirs(std::vector<CMeshObject*>&crowns, std::vector<OpenMesh::Vec3d>&res_front_dirs)
{
	res_front_dirs.clear();
	std::vector<Eigen::Matrix4d>trans_mats;
	CDentalTemplateFitting::ComputeStretchCrowns2LineMatrix(crowns, trans_mats);
	
	for (int i = 0; i < trans_mats.size(); i++)
	{


		OpenMesh::Vec3d center_before = crowns[i]->TransformPointByLocalMatrix(CGeoBaseAlg::ComputeMeshCenter(crowns[i]->GetMesh()));
		OpenMesh::Vec3d center_after = CGeoBaseAlg::Transform(trans_mats[i], center_before);
		OpenMesh::Vec3d dir(0, 0, 1);
		Eigen::Matrix4d inv_mat = trans_mats[i].inverse();
		dir = CGeoBaseAlg::Transform(inv_mat, center_after + dir) - center_before;

		res_front_dirs.push_back(dir);
		/*CCurveObject *co = new CCurveObject();

		co->GetCurve().push_back(center_before);
		co->GetCurve().push_back(center_before + dir * 10);
		co->SetChanged();
		DataPool::AddCurveObject(co);*/

		//crowns[i]->GetMatrix() = (trans_mats[i] * crowns[i]->GetMatrix());
		//crowns[i]->SetChanged();
	}
}
void CDentalTemplateFitting::ComputeStretchCrowns2LineMatrix(std::vector<CMeshObject*>&crowns, std::vector<Eigen::Matrix4d>&res_mats)
{
	std::vector<OpenMesh::Vec2d>center_points;
	std::vector<OpenMesh::Vec3d>center_points_3d;
	double mean_y = 0;
	double mean_z = 0;
	double min_x = std::numeric_limits<double>::max();
	double max_x = std::numeric_limits<double>::min();
	//std::vector<OpenMesh::Vec3d>curve0;
	center_points.resize(crowns.size());
	center_points_3d.resize(crowns.size());
	for (int i=0;i<crowns.size();i++)
	{
		OpenMesh::Vec3d center = crowns[i]->TransformPointByLocalMatrix(CGeoBaseAlg::ComputeMeshCenter(crowns[i]->GetMesh()));
		center_points[i] = OpenMesh::Vec2d(center[0], center[2]);
		center_points_3d[i] = center;
		mean_y += center[1];
		mean_z += center[2];
		if (min_x > center[0])
			min_x = center[0];
		if (max_x < center[0])
			max_x = center[0];
		//curve0.push_back(center);
	}
	/*CCurveObject *curve_obj0 = new CCurveObject();
	curve_obj0->SetCurve(curve0);
	curve_obj0->SetColor(OpenMesh::Vec3d(0, 0, 1));
	curve_obj0->SetChanged();
	DataPool::AddCurveObject(curve_obj0);*/
	std::cerr << "max x " << max_x << " min x" << min_x << std::endl;
	mean_y /= crowns.size();
	mean_z /= crowns.size();
	std::vector<double>coeffs;
	std::vector<OpenMesh::Vec2d>vec_center_points;
	vec_center_points.resize(center_points.size());
	for (int i=0;i<vec_center_points.size();i++)
	{
		vec_center_points[i] = center_points[i];
	}
	CCurveBaseAlg::PolynomialFitting(vec_center_points, 4, coeffs);
	std::vector<OpenMesh::Vec2d>curve;
	std::vector<OpenMesh::Vec3d>curve3d;
	for (double xi = min_x - 0.02; xi <= max_x + 0.02; xi = xi + 0.01)
	{
		//std::cerr << xi << std::endl;
		double zi = 0;
		double nxt = 1;
		for (int di = 0; di < coeffs.size(); di++)
		{
			zi += coeffs[di] * nxt;
			nxt *= xi;
		}
		curve.push_back(OpenMesh::Vec2d(xi, zi));
		curve3d.push_back(OpenMesh::Vec3d(xi, 0.15, zi));
	}
	/*	CCurveObject *curveobj = new CCurveObject();
	curveobj->SetCurve(curve3d);
	curveobj->SetColor(OpenMesh::Vec3d(0, 1, 0));
	curveobj->SetChanged();
	std::cerr << "curve size " << curveobj->GetCurve().size() << std::endl;
	DataPool::AddCurveObject(curveobj);*/
	double curve_len = CCurveBaseAlg::ComputeLenOfCurve(curve);

	std::vector<double>curve_lens(curve.size(), 0);
	curve_lens[0] = 0;
	for (int i = 0; i < curve.size(); i++)
	{
		int prei = i - 1 >= 0 ? i - 1 : i;
		curve_lens[i] = curve_lens[prei] + ((curve[i] - curve[prei]).length());
	}
	res_mats.resize(crowns.size());
	for (int i=0;i<center_points.size();i++)
	{
	
		
		int closest_pid;
		CCurveBaseAlg::ComputeClosestPoint(curve, center_points[i], closest_pid);
		OpenMesh::Vec2d new_centerpoint = OpenMesh::Vec2d(curve_lens[closest_pid] - curve_len / 2.0, mean_z);
		OpenMesh::Vec3d trans(new_centerpoint[0] - center_points[i][0], 0, new_centerpoint[1] - center_points[i][1]);
		
		OpenMesh::Vec2d normal_dir = CCurveBaseAlg::ComputeNormalOfPolynomial(coeffs, curve[closest_pid][0]);
		OpenMesh::Vec3d rot_axis = OpenMesh::cross(OpenMesh::Vec3d(normal_dir[0], 0, normal_dir[1]), OpenMesh::Vec3d(0, 0, 1));
		double rot_angle = std::acos(OpenMesh::dot(normal_dir, OpenMesh::Vec2d(0, 1)));
		Eigen::Vector3d rot_axis_eg;

		Eigen::Matrix4d rot_mat = CGeoBaseAlg::ComputeRotMat(rot_axis, rot_angle, OpenMesh::Vec3d(new_centerpoint[0], mean_y, new_centerpoint[1]));
		Eigen::Matrix4d trans_mat;
		trans_mat.setIdentity();
		for (int j = 0; j < 3; j++)
		{
			trans_mat(j, 3) = trans[j];
		}
	

		res_mats[i]= rot_mat*trans_mat;
	}
}
void CDentalTemplateFitting::OrderCrowns(std::vector<CMeshObject*>&crowns, OpenMesh::Vec3d crown_updir, bool from_left2right)
{
	std::vector<Eigen::Matrix4d>mats;
	ComputeStretchCrowns2LineMatrix(crowns, mats);
	std::vector<OpenMesh::Vec3d>centers;
	for (int i = 0; i < crowns.size(); i++)
	{
		OpenMesh::Vec3d p = crowns[i]->TransformPointByLocalMatrix(CGeoBaseAlg::ComputeMeshCenter(crowns[i]->GetMesh()));
		centers.push_back(CGeoBaseAlg::Transform(mats[i], p));
	}
	std::vector<double>center_x;
	for (int i = 0; i < centers.size(); i++)
	{
		center_x.push_back(centers[i][0]);
	}
	std::vector<int>order;
	CNumericalBaseAlg::ComputeOrder(center_x, order);

	std::vector<CMeshObject*>tmp_crowns = crowns;
	for (int i = 0; i < crowns.size(); i++)
	{
		crowns[order[i]] = tmp_crowns[i];
	}
	if (!from_left2right)
	{
		std::reverse(crowns.begin(), crowns.end());
	}
	


}
void CDentalTemplateFitting::ComputeCrownsType(std::vector<CMeshObject*>crowns, std::map<CMeshObject*, int>&fixed_id, std::map<CMeshObject *, std::string>&res_crown_type)
{
	std::map<int,int>left_teeth,right_teeth;
	left_teeth[25] = 29; left_teeth[21] = 25;
	left_teeth[17] = 21; left_teeth[13] = 17;
	left_teeth[9] = 13; left_teeth[5] = 9;
	left_teeth[1] = 5; left_teeth[2] = 1;
	left_teeth[6] = 2; left_teeth[10] = 6;
	left_teeth[14] = 10; left_teeth[18] = 14;
	left_teeth[22] = 18; left_teeth[26] = 22;
	left_teeth[30] = 26; 

	left_teeth[27] = 31; left_teeth[23] = 27;
	left_teeth[19] = 23; left_teeth[15] = 19;
	left_teeth[11] = 15; left_teeth[7] = 11;
	left_teeth[3] = 7; left_teeth[4] = 3;
	left_teeth[8] = 4; left_teeth[12] = 8;
	left_teeth[16] = 12; left_teeth[20] = 16;
	left_teeth[24] = 20; left_teeth[28] = 24;
	left_teeth[32] = 28; 

	for (auto iter = left_teeth.begin(); iter != left_teeth.end(); iter++)
	{
		right_teeth[iter->second] = iter->first;
	}
	std::map<int, std::string>crown_id2type_map;
	crown_id2type_map[5] = crown_id2type_map[1] = crown_id2type_map[2] = crown_id2type_map[6] = "incisor";
	crown_id2type_map[7] = crown_id2type_map[3] = crown_id2type_map[4] = crown_id2type_map[8] = "incisor";
	crown_id2type_map[9] = crown_id2type_map[10]=  crown_id2type_map[11]=  crown_id2type_map[12] = "canine";
	crown_id2type_map[17] = crown_id2type_map[13] = crown_id2type_map[14] = crown_id2type_map[18] = "premolar";
	crown_id2type_map[19] = crown_id2type_map[15] = crown_id2type_map[16] = crown_id2type_map[20] = "premolar";
	crown_id2type_map[25] = crown_id2type_map[21] = crown_id2type_map[26] = crown_id2type_map[22] = "molar";
	crown_id2type_map[23] = crown_id2type_map[27] = crown_id2type_map[24] = crown_id2type_map[28] = "molar";
	crown_id2type_map[31] = crown_id2type_map[32] = crown_id2type_map[29] = crown_id2type_map[30] = "molar";
	CDentalTemplateFitting::OrderCrowns(crowns, OpenMesh::Vec3d(0, 1, 0), true);
	std::cerr << std::endl;
	for (int i = 0; i < crowns.size(); i++)
	{
		OpenMesh::Vec3d c = CGeoBaseAlg::ComputeMeshCenter(crowns[i]->GetMesh());
		std::cerr << c << std::endl;
	}
	std::vector<int>crown_ids;//id of teeth in dentical
	crown_ids.resize(crowns.size(),-1);
	for (int i = 0; i < crown_ids.size(); i++)
	{
		if (fixed_id.find(crowns[i]) != fixed_id.end())
		{
			crown_ids[i] = fixed_id[crowns[i]];
		}
	}
	for (int i = 0; i < crown_ids.size(); i++)
	{
	
		if (crown_ids[i] != -1)
		{
			for (int j = i - 1; j >= 0; j--)
			{
				if (crown_ids[j] != -1)
				{
					break;
				}
				crown_ids[j] = left_teeth[crown_ids[j + 1]];
			}
			for (int j = i + 1; j<crown_ids.size(); j++)
			{
				if (crown_ids[j] != -1)
				{
					break;
				}
				crown_ids[j] = right_teeth[crown_ids[j - 1]];
			}
		}
	}
	std::cerr << std::endl;
	for (int i = 0; i < crown_ids.size(); i++)
	{
		std::cerr << "id " << crown_ids[i] << std::endl;
	}

	res_crown_type.clear();
	std::cerr << std::endl;
	for (int i = 0; i < crowns.size(); i++)
	{
		res_crown_type[crowns[i]] = crown_id2type_map[crown_ids[i]];
		std::cerr <<"id "<< crowns[i]->GetId() << " type " << res_crown_type[crowns[i]] << std::endl;
	}
}
void CDentalTemplateFitting::LoadTemplates(std::string dirname, std::map<std::string,  CTeethTemplateObject*>&res_templates)
{
	std::string index_fname = dirname + "//index.txt";
	std::ifstream in(index_fname);
	char buf[256];
	res_templates.clear();
	while (in.getline(buf, sizeof buf))
	{
		std::istringstream line(buf);
		std::string type;
		line >> type;
		std::string fname;
		line >> fname;
		fname = dirname + "//" + fname;
		CDataIO dataio;
		CTeethTemplateObject *temp_obj = new CTeethTemplateObject();
		dataio.ReadMesh(fname, *temp_obj, OpenMesh::IO::Options::VertexColor);
		COpenMeshT&mesh = temp_obj->GetMesh();
		std::vector<COpenMeshT::FaceHandle>roifhs;
		for (auto fiter = mesh.faces_begin(); fiter != mesh.faces_end(); fiter++)
		{
			bool flag = true;
			for (auto viter = mesh.fv_begin(fiter); viter != mesh.fv_end(fiter); viter++)
			{
				OpenMesh::Vec3d c = mesh.color(viter);
				if (c[0] > 0.1 || c[1] > 0.1&&c[2] > 0.1)
				{
					flag = false;
					break;
				}
			}
			if (flag)
			{
				roifhs.push_back(fiter);
			}
		}
		temp_obj->CrownFhs() = roifhs;
		res_templates[type] = temp_obj;
		temp_obj->SetMeshColor(OpenMesh::Vec3d(0.8, 0.8, 0.8));
		temp_obj->SetAttrChanged();
	}
}
void CDentalTemplateFitting::MergeBoundaryRegisteredMeshes(COpenMeshT &mesh_a, Eigen::Matrix4d mat_a, COpenMeshT&mesh_b, Eigen::Matrix4d mat_b, std::map<COpenMeshT::VertexHandle, COpenMeshT::VertexHandle>&bound_map_a2b, COpenMeshT&res_mesh, std::vector<COpenMeshT::VertexHandle>&res_seam_vhs)
{

	res_mesh.clear();
	res_seam_vhs.clear();
	std::vector<COpenMeshT::VertexHandle>bound_a;
	CGeoBaseAlg::GetLargestOrderedBoundary(mesh_a, bound_a);
	std::map<COpenMeshT::VertexHandle, COpenMeshT::VertexHandle>b2new_vhsmap, a2new_vhsmap;
	for (auto viter = mesh_b.vertices_begin(); viter != mesh_b.vertices_end(); viter++)
	{
		COpenMeshT::VertexHandle new_vh=res_mesh.add_vertex(CGeoBaseAlg::Transform(mat_b,mesh_b.point(viter)));
		b2new_vhsmap[viter] = new_vh;
	}
	for (auto viter = mesh_a.vertices_begin(); viter != mesh_a.vertices_end(); viter++)
	{
		if (bound_map_a2b.find(viter) == bound_map_a2b.end())
		{
			COpenMeshT::VertexHandle new_vh = res_mesh.add_vertex(CGeoBaseAlg::Transform(mat_a, mesh_a.point(viter)));
			a2new_vhsmap[viter] = new_vh;
		}
		else
		{
			a2new_vhsmap[viter] = b2new_vhsmap[bound_map_a2b[viter]];
		}
	
	}
	for (int i = 0; i < bound_a.size(); i++)
	{
		res_seam_vhs.push_back(a2new_vhsmap[bound_a[i]]);
	}
	for (auto fiter = mesh_a.faces_begin(); fiter != mesh_a.faces_end(); fiter++)
	{
		std::vector<COpenMeshT::VertexHandle>vhs;
		for (auto viter = mesh_a.fv_begin(fiter); viter != mesh_a.fv_end(fiter); viter++)
		{
			vhs.push_back(a2new_vhsmap[viter]);
		}
		res_mesh.add_face(vhs);
	}
	for (auto fiter = mesh_b.faces_begin(); fiter != mesh_b.faces_end(); fiter++)
	{
		std::vector<COpenMeshT::VertexHandle>vhs;
		for (auto viter = mesh_b.fv_begin(fiter); viter != mesh_b.fv_end(fiter); viter++)
		{
			vhs.push_back(b2new_vhsmap[viter]);
		}
		res_mesh.add_face(vhs);
	}
	res_mesh.garbage_collection();
}

void CDentalTemplateFitting::AlignTemplate2Crown(CMeshObject &crown, OpenMesh::Vec3d front_dir, OpenMesh::Vec3d up_dir, CTeethTemplateObject& teeth_template)
{
	teeth_template.ApplyTransform();
	COpenMeshT &temp_crown= teeth_template.GetCrownMesh();
	COpenMeshT crown_mesh = crown.GetMesh();
	OpenMesh::Vec3d temp_center=CGeoBaseAlg::ComputeMeshCenter(temp_crown);
	OpenMesh::Vec3d crown_center = crown.TransformPointByLocalMatrix(CGeoBaseAlg::ComputeMeshCenter(crown.GetMesh()));

	
	OpenMesh::Vec3d crown_axisx, crown_axisy, crown_axisz;
	OpenMesh::Vec3d temp_axisx, temp_axisy, temp_axisz;
	crown_axisx = OpenMesh::cross(up_dir, front_dir);
	crown_axisy = up_dir;
	crown_axisz = front_dir;
	temp_axisy = teeth_template.UpDir();
	temp_axisz = teeth_template.FrontDir();
	temp_axisx = OpenMesh::cross(temp_axisy, temp_axisz);
	crown_axisx.normalize(), crown_axisy.normalize(), crown_axisz.normalize();
	temp_axisx.normalize(), temp_axisy.normalize(), temp_axisz.normalize();

	double crown_width, crown_w_max = std::numeric_limits<double>::min(), crown_w_min = std::numeric_limits<double>::max();
	for (auto viter = crown_mesh.vertices_begin(); viter != crown_mesh.vertices_end(); viter++)
	{
		OpenMesh::Vec3d p = crown.TransformPointByLocalMatrix(crown_mesh.point(viter));
		p = p - crown_center;
		double tmp_w=OpenMesh::dot(p, crown_axisx);
		if (crown_w_min > tmp_w)
		{
			crown_w_min = tmp_w;
		}
		if (crown_w_max < tmp_w)
		{
			crown_w_max = tmp_w;
		}
	}
	crown_width = crown_w_max-crown_w_min;

	double temp_width, temp_w_max = std::numeric_limits<double>::min(), temp_w_min = std::numeric_limits<double>::max();
	for (auto viter = temp_crown.vertices_begin(); viter != temp_crown.vertices_end(); viter++)
	{
		OpenMesh::Vec3d p = temp_crown.point(viter);
		p = p - temp_center;
		double tmp_w = OpenMesh::dot(p, temp_axisx);
		if (temp_w_min > tmp_w)
		{
			temp_w_min = tmp_w;
		}
		if (temp_w_max < tmp_w)
		{
			temp_w_max = tmp_w;
		}
	}
	temp_width = temp_w_max - temp_w_min;
	double scale = crown_width / temp_width;
	std::cerr << "crown width " << crown_width << " temp_width " << temp_width << std::endl;

	Eigen::Matrix4d scale_mat =CGeoBaseAlg::ComputeScaleMat(scale, temp_center);

	std::vector<OpenMesh::Vec3d>crown_frame, temp_frame;
	crown_frame.push_back(crown_axisx);
	crown_frame.push_back(crown_axisy);
	crown_frame.push_back(crown_axisz);
	temp_frame.push_back(temp_axisx);
	temp_frame.push_back(temp_axisy);
	temp_frame.push_back(temp_axisz);
	Eigen::Matrix4d transform_mat = CGeoBaseAlg::ComputeFrameTransMatrix( temp_center, temp_frame, crown_center, crown_frame);
	
	
	Eigen::Matrix4d &teeth_temp_mat = teeth_template.GetMatrix();
	teeth_temp_mat = transform_mat*scale_mat*teeth_temp_mat;


	//std::pair<OpenMesh::Vec3d, OpenMesh::Vec3d>temp_bound, crown_bound;
	//CGeoBaseAlg::ComputeAABB(temp_crown, temp_bound.first, temp_bound.second);
	//CGeoBaseAlg::ComputeAABB(crown.GetMesh(), crown_bound.first, crown_bound.second);
	//
	//crown_bound.first = crown.TransformPointByLocalMatrix(crown_bound.first);
	//crown_bound.second = crown.TransformPointByLocalMatrix(crown_bound.second);
	//

	//std::cerr << "temp bound " << temp_bound.first << " " << temp_bound.second << std::endl;
	//std::cerr << "crown bound " << crown_bound.first << " " << crown_bound.second << std::endl;

	//teeth_template.Transform(-crown_center);
	//Eigen::Matrix4d& temp_mat = teeth_template.GetMatrix();
	//OpenMesh::Vec3d vec_scale = (crown_bound.second - crown_bound.first) / (temp_bound.second-temp_bound.first);
	//double scale = (std::abs(vec_scale[0])+std::abs(vec_scale[0]))/2.0;
	//std::cerr << "vec scale " << vec_scale << std::endl;
	//std::cerr << "scale " << scale << std::endl;
	//Eigen::Matrix4d scale_mat;
	//scale_mat.setIdentity();

	//scale_mat(0, 0) = scale;
	//scale_mat(1, 1) = scale;
	//scale_mat(2, 2) = scale;
	//temp_mat = scale_mat*temp_mat;
	//teeth_template.Transform(crown_center);
	////std::cerr << teeth_template.GetMatrix() << std::endl;
	teeth_template.ApplyTransform();
	teeth_template.SetChanged();

}
void CDentalTemplateFitting::ReplaceTemplateCrown(CMeshObject &teeth_template_obj, CMeshObject &mesh_crown_obj)
{
	std::cerr << "ReplaceTemplateCrown" << std::endl;
	std::vector<COpenMeshT::VertexHandle>crown_bound;
	CGeoBaseAlg::GetLargestOrderedBoundary(mesh_crown_obj.GetMesh(), crown_bound);
	
	COpenMeshT &temp_mesh = teeth_template_obj.GetMesh();
	COpenMeshT &crown_mesh = mesh_crown_obj.GetMesh();
	OpenMesh::Vec3d pre_temp_mesh_color = temp_mesh.color(temp_mesh.vertices_begin());
	std::vector<std::pair<COpenMeshT::FaceHandle, OpenMesh::Vec3d>>face_points;

	for (int i = 0; i < crown_bound.size(); i++)
	{
		OpenMesh::Vec3d p = crown_mesh.point(crown_bound[i]);
		p = mesh_crown_obj.TransformPointByLocalMatrix(p);
		COpenMeshT::FaceHandle res_fh;
		OpenMesh::Vec3d res_p;
		double dis=CGeoAlg::ComputeClosestPoint(p, teeth_template_obj,res_fh,res_p);
		
		face_points.push_back(std::pair<COpenMeshT::FaceHandle, OpenMesh::Vec3d>(res_fh, res_p));
	

	}
	std::vector<OpenMesh::Vec3d>res_ps;
	std::vector<bool>is_vertex;
	
	std::map<int, COpenMeshT::FaceHandle>res_fhs_map;
	std::map<int, COpenMeshT::EdgeHandle>res_ehs_map;
	std::map<int, COpenMeshT::VertexHandle>res_vhs_map;
	std::vector<int>face_points2res_ps_map;
	CGeoAlg::FindPointsPath(teeth_template_obj, face_points, res_ps, res_fhs_map, res_ehs_map, res_vhs_map, face_points2res_ps_map,true);


	CCurveObject *c_obj = new CCurveObject();
	std::vector<OpenMesh::Vec3d>curve;
	for (int i = 0; i < res_ps.size(); i++)
	{
		curve.push_back(res_ps[i]);
	}
	//curve.push_back(curve[0]);
	c_obj->SetColor(OpenMesh::Vec3d(0, 1, 0));
	c_obj->SetCurve(curve);
	c_obj->SetChanged();
	std::cerr << "curve num " << c_obj->GetCurve().size() << std::endl;
	DataPool::AddCurveObject(c_obj);
//	return;
	std::map<COpenMeshT::FaceHandle, std::vector<COpenMeshT::VertexHandle>>temp_fv_map;
	std::map < COpenMeshT::FaceHandle, std::vector<std::pair<COpenMeshT::VertexHandle, COpenMeshT::VertexHandle>>>temp_fcons_map;
	std::map<COpenMeshT::VertexHandle, COpenMeshT::EdgeHandle>vh_eh_map;
	std::vector<COpenMeshT::VertexHandle>vhs_loop;
	std::set<COpenMeshT::VertexHandle>dup_set;
	std::set<COpenMeshT::FaceHandle> pre_fhs,first_fhs;
	COpenMeshT::VertexHandle pre_vh,first_vh;
	for (int i = 0; i < res_ps.size(); i++)
	{
		auto vhiter = res_vhs_map.find(i);
		auto ehiter = res_ehs_map.find(i);
		auto fhiter = res_fhs_map.find(i);
		if (vhiter != res_vhs_map.end())
		{
			
			for (auto vfiter = temp_mesh.vf_begin(vhiter->second);vfiter!= temp_mesh.vf_end(vhiter->second);vfiter++)
			{
				if (pre_fhs.find(vfiter) != pre_fhs.end())
				{
					if (temp_fcons_map.find(vfiter) == temp_fcons_map.end())
					{
						temp_fcons_map[vfiter] = std::vector<std::pair<COpenMeshT::VertexHandle, COpenMeshT::VertexHandle>>();
					}
					temp_fcons_map[vfiter].push_back(std::make_pair(vhiter->second, pre_vh));
				}
				if (i == res_ps.size() - 1)
				{
					if (first_fhs.find(vfiter) != first_fhs.end())
					{
						if (temp_fcons_map.find(vfiter) == temp_fcons_map.end())
						{
							temp_fcons_map[vfiter] = std::vector<std::pair<COpenMeshT::VertexHandle, COpenMeshT::VertexHandle>>();
						}
						temp_fcons_map[vfiter].push_back(std::make_pair(vhiter->second, first_vh));
					}
				}
				
			}

			pre_fhs.clear();
			for (auto vfiter = temp_mesh.vf_begin(vhiter->second); vfiter != temp_mesh.vf_end(vhiter->second); vfiter++)
			{
				pre_fhs.insert(vfiter);
			}
			pre_vh = vhiter->second;

			vhs_loop.push_back(vhiter->second);
			if (dup_set.find(vhiter->second) == dup_set.end())
			{
				dup_set.insert(vhiter->second);
			}
			else
			{
				std::cerr << "error vh already exists in dup_set" << std::endl;
			}
		}
		else if (ehiter != res_ehs_map.end())
		{
			COpenMeshT::EdgeHandle eh = ehiter->second;
			COpenMeshT::HalfedgeHandle he0 = temp_mesh.halfedge_handle(eh, 0), he1 = temp_mesh.halfedge_handle(eh, 1);
			COpenMeshT::FaceHandle fh0=temp_mesh.face_handle(he0);
			COpenMeshT::FaceHandle fh1=temp_mesh.face_handle(he1);
		
			COpenMeshT::VertexHandle new_vh = temp_mesh.add_vertex(res_ps[i]);
			vhs_loop.push_back(new_vh);
			if (dup_set.find(new_vh) == dup_set.end())
			{
				dup_set.insert(new_vh);
			}
			else
			{
				std::cerr << "error evh already exists in dup_set" << std::endl;
			}
			vh_eh_map[new_vh] = eh;
			if (temp_fv_map.find(fh0) == temp_fv_map.end())
			{
				temp_fv_map[fh0] = std::vector<COpenMeshT::VertexHandle>();
				temp_fcons_map[fh0] = std::vector<std::pair<COpenMeshT::VertexHandle, COpenMeshT::VertexHandle>>();
			}
			if (temp_fv_map.find(fh1) == temp_fv_map.end())
			{
				temp_fv_map[fh1] = std::vector<COpenMeshT::VertexHandle>();
				temp_fcons_map[fh1] = std::vector<std::pair<COpenMeshT::VertexHandle, COpenMeshT::VertexHandle>>();
			}
		
			if (pre_fhs.find(fh0)!=pre_fhs.end())
			{
				
				temp_fcons_map[fh0].push_back(std::make_pair(pre_vh, new_vh));
			}
			
			if (pre_fhs.find(fh1) != pre_fhs.end())
			{
				temp_fcons_map[fh1].push_back(std::make_pair(pre_vh, new_vh));
			}

			if (i == res_ps.size() - 1)
			{

				if (first_fhs.find(fh0) != first_fhs.end())
				{
					if (temp_fcons_map.find(fh0) == temp_fcons_map.end())
					{
						temp_fcons_map[fh0] = std::vector<std::pair<COpenMeshT::VertexHandle, COpenMeshT::VertexHandle>>();
					}
					temp_fcons_map[fh0].push_back(std::make_pair(first_vh, new_vh));
				}

				if (first_fhs.find(fh1) != first_fhs.end())
				{
					if (temp_fcons_map.find(fh1) == temp_fcons_map.end())
					{
						temp_fcons_map[fh1] = std::vector<std::pair<COpenMeshT::VertexHandle, COpenMeshT::VertexHandle>>();
					}
					temp_fcons_map[fh1].push_back(std::make_pair(first_vh, new_vh));
				}
			}

			temp_fv_map[fh0].push_back(new_vh);
			temp_fv_map[fh1].push_back(new_vh);
			pre_fhs.clear();
			pre_fhs.insert(fh0);
			pre_fhs.insert(fh1);
			pre_vh = new_vh;
			
		}
		else if(fhiter!=res_fhs_map.end())
		{
			COpenMeshT::FaceHandle fh = fhiter->second;
			if (temp_fv_map.find(fh) == temp_fv_map.end())
			{
				temp_fv_map[fh] = std::vector<COpenMeshT::VertexHandle>();
				temp_fcons_map[fh] = std::vector<std::pair<COpenMeshT::VertexHandle, COpenMeshT::VertexHandle>>();
			}
			COpenMeshT::VertexHandle new_vh = temp_mesh.add_vertex(res_ps[i]);
			vhs_loop.push_back(new_vh);
			if (dup_set.find(new_vh) == dup_set.end())
			{
				dup_set.insert(new_vh);
			}
			else
			{
				std::cerr << "error fvh already exists in dup_set" << std::endl;
			}
			
			if (pre_fhs.find(fh)!= pre_fhs.end())
			{
				temp_fcons_map[fh].push_back(std::make_pair(pre_vh, new_vh));
			}
			if (i == res_ps.size() - 1)
			{

				if (first_fhs.find(fh) != first_fhs.end())
				{
					if (temp_fcons_map.find(fh) == temp_fcons_map.end())
					{
						temp_fcons_map[fh]= std::vector<std::pair<COpenMeshT::VertexHandle, COpenMeshT::VertexHandle>>();
					}
					temp_fcons_map[fh].push_back(std::make_pair(first_vh, new_vh));
				}
			}
			
			temp_fv_map[fh].push_back(new_vh);
			pre_fhs.clear();
			pre_fhs.insert(fh);
			pre_vh = new_vh;
		}
		else
		{
			std::cerr << "size res_fhs_map " << res_fhs_map.size() << std::endl;
			std::cerr << "size res_ehs_map " << res_ehs_map.size() << std::endl;
			std::cerr << "error in the result of FindPointsPath" << std::endl;
		}
		if (i == 0)
		{
			first_fhs = pre_fhs;
			first_vh = pre_vh;
		}
	
	}
	for (auto iter = temp_fv_map.begin(); iter != temp_fv_map.end(); iter++)
	{
	
		CGeoAlg::SplitFaceByVhs(temp_mesh, iter->first, iter->second, temp_fcons_map[iter->first], vh_eh_map);
	}

	temp_mesh.garbage_collection();
//	teeth_template_obj.SetChanged();
	//return;

	std::vector<COpenMeshT*>sep_meshes;

	std::cerr << "separate mesh" << std::endl;
	std::vector<std::map<COpenMeshT::VertexHandle, COpenMeshT::VertexHandle>>temp_new2orig_vhsmap,temp_orig2new_vhsmap;
	std::vector<std::vector<COpenMeshT::VertexHandle>>crown_bound2new_vhsmaps;
	CGeoAlg::SeparateMeshByVhsLoop(temp_mesh, vhs_loop, sep_meshes, temp_new2orig_vhsmap);
	temp_orig2new_vhsmap.resize(temp_new2orig_vhsmap.size());
	crown_bound2new_vhsmaps.resize(temp_new2orig_vhsmap.size());
	for (int i = 0; i < temp_new2orig_vhsmap.size(); i++)
	{
		temp_orig2new_vhsmap[i] = std::map<COpenMeshT::VertexHandle, COpenMeshT::VertexHandle>();
		for (auto iter = temp_new2orig_vhsmap[i].begin(); iter != temp_new2orig_vhsmap[i].end(); iter++)
		{
			temp_orig2new_vhsmap[i][iter->second] = iter->first;
			if (sep_meshes[i]->vertex_handle(iter->first.idx()).idx() == -1)
			{
				std::cerr << "error sep_meshes[i]->vertex_handle() == -1,first second " << iter->first << " " << iter->second << std::endl;
			}
		//	
		}
	}

	for (int i = 0; i < crown_bound2new_vhsmaps.size(); i++)
	{
		//std::cerr << "temp_mesh vsize " << temp_mesh.n_vertices() << " temp_orig2new_vhsmap size " << temp_orig2new_vhsmap[i].size() << std::endl;
		
		crown_bound2new_vhsmaps[i].clear();
		for (int j = 0; j < face_points2res_ps_map.size(); j++)
		{
			if (vhs_loop.size()<face_points2res_ps_map[j])
			{
				std::cerr << "error vhs_loop.size()<face_points2res_ps_map[j]" << std::endl;
			}
			
			crown_bound2new_vhsmaps[i].push_back(temp_orig2new_vhsmap[i][vhs_loop[face_points2res_ps_map[j]]]);
			if (crown_bound2new_vhsmaps[i].back().idx() == -1)
			{
				std::cerr << "error crown_bound2new_vhsmaps[i].back().idx==-1 vhs_loop[face_points2res_ps_map[j]]:" << vhs_loop[face_points2res_ps_map[j]].idx() << std::endl;
				

			/*	CCurveObject *c_obj = new CCurveObject();
				std::vector<OpenMesh::Vec3d>curve;
				
				curve.push_back(temp_mesh.point(vhs_loop[face_points2res_ps_map[j]]));
				curve.push_back(OpenMesh::Vec3d(0, 0, 0));
				c_obj->SetColor(OpenMesh::Vec3d(0, 1, 1));
				c_obj->SetCurve(curve);
				c_obj->SetChanged();
				DataPool::AddCurveObject(c_obj);*/
			}
		}
	}
	CCurveObject *co_loop = new CCurveObject();
	//std::cerr << "detecting interval start" << std::endl;
	for (int i = 0; i < vhs_loop.size(); i++)
	{
		COpenMeshT::EdgeHandle eh;
		if (!CGeoBaseAlg::GetCommonEdge(temp_mesh, vhs_loop[i], vhs_loop[(i + 1) % vhs_loop.size()], eh))
		{
			if (temp_mesh.point(vhs_loop[i]) != res_ps[i])
			{
				std::cerr << "error temp_mesh.point(vhs_loop[i]) != res_ps[i], i: " << i << std::endl;
			}
			//std::cerr << i << std::endl;
			//std::cerr << "res ps i:" << res_ps[i] << std::endl;
			CCurveObject *co = new CCurveObject();
			co->GetCurve().push_back(OpenMesh::Vec3d(0, 0, 0));
			co->GetCurve().push_back(temp_mesh.point(vhs_loop[i]));
			co->SetColor(OpenMesh::Vec3d(1, 0, 0));
			DataPool::AddCurveObject(co);

			CCurveObject *co2 = new CCurveObject();
			co2->GetCurve().push_back(OpenMesh::Vec3d(0, 0, 0));
			co2->GetCurve().push_back(temp_mesh.point(vhs_loop[(i + 1) % vhs_loop.size()]));
			co2->SetColor(OpenMesh::Vec3d(1, 0, 0));
			DataPool::AddCurveObject(co2);
			std::cerr << temp_mesh.point(vhs_loop[i]) << std::endl;
			std::cerr << temp_mesh.point(vhs_loop[(i + 1) % vhs_loop.size()]) << std::endl;
		
		}
		co_loop->GetCurve().push_back(temp_mesh.point(vhs_loop[i]));
	}
	//std::cerr << "detecting interval end" << std::endl;
	co_loop->SetColor(OpenMesh::Vec3d(0, 0, 1));
	co_loop->SetChanged();
	DataPool::AddCurveObject(co_loop);
	//return;
	OpenMesh::Vec3d temp_center=CGeoBaseAlg::ComputeMeshCenter(crown_mesh);
	
	double max_dis = std::numeric_limits<double>::min();
	//std::cerr << "sep meshes " << sep_meshes.size() << std::endl;
	int mi;
	for (int i = 0; i < sep_meshes.size(); i++)
	{
		OpenMesh::Vec3d sep_mesh_center=CGeoBaseAlg::ComputeMeshCenter(*sep_meshes[i]);
		OpenMesh::Vec3d tmp_dir = sep_mesh_center - temp_center;
		double dis = std::sqrt(tmp_dir.sqrnorm());
		if (dis > max_dis)
		{
			mi = i;
			max_dis = dis;
		}
		//std::cerr << "sep mesh v num " << sep_meshes[i]->n_vertices() << std::endl;
	}
	//std::cerr << "root temp part mi " << mi << std::endl;
	teeth_template_obj.GetMatrix().setIdentity();
	teeth_template_obj.SetMesh(*sep_meshes[mi]);
	teeth_template_obj.SetMeshColor(pre_temp_mesh_color);
	teeth_template_obj.SetChanged();
	

	COpenMeshT &temp_root_mesh = teeth_template_obj.GetMesh();
	//crown_bound
	std::vector<COpenMeshT::VertexHandle>crown_bound2temp_vhsmap;
	
	for (int i = 0; i < crown_bound2new_vhsmaps[mi].size(); i++)
	{
		//std::cerr << "crown_bound2temp_vhsmap " << crown_bound2temp_vhsmap[i].idx() << std::endl;
		crown_bound2temp_vhsmap.push_back( temp_root_mesh.vertex_handle(crown_bound2new_vhsmaps[mi][i].idx()));
		if (crown_bound2temp_vhsmap.back().idx() == -1)
		{
			std::cerr << "crown_bound2temp_vhsmap.back().idx() == -1,crown_bound2new_vhsmaps[mi][i].idx: " << crown_bound2new_vhsmaps[mi][i].idx() << std::endl;
		}
	}
	
	//srand((unsigned)time(NULL));
	//for (int i = 0; i < crown_bound2temp_vhsmap.size(); i++)
	//{
	//	double r = rand() % 1000 / 1000.0;
	//	double g = rand() % 1000 / 1000.0;
	//	double b = rand() % 1000 / 1000.0;
	//	crown_mesh.set_color(crown_bound[i], OpenMesh::Vec3d(r, g, b));
	//	/*std::cerr << "bbb" << std::endl;
	//	std::cerr << crown_bound2temp_vhsmap[i].idx() << std::endl;
	//	std::cerr << "aaa" << std::endl;*/
	//	temp_root_mesh.set_color(crown_bound2temp_vhsmap[i], OpenMesh::Vec3d(r, g, b));
	//}

	std::set<COpenMeshT::VertexHandle>mapped_temp_bound_set;
	for (int i = 0; i < crown_bound2temp_vhsmap.size(); i++)
	{
		mapped_temp_bound_set.insert(crown_bound2temp_vhsmap[i]);
	}
	std::vector<COpenMeshT::VertexHandle>temp_bound;
	CGeoBaseAlg::GetLargestOrderedBoundary(temp_root_mesh, temp_bound);
	std::vector<OpenMesh::Vec3d>temp_bound_points, crown_bound_points;
	for (int i = 0; i < temp_bound.size(); i++)
	{
		temp_bound_points.push_back(teeth_template_obj.TransformPointByLocalMatrix(temp_root_mesh.point(temp_bound[i])));
	}
	for (int i = 0; i < crown_bound.size(); i++)
	{
		crown_bound_points.push_back(mesh_crown_obj.TransformPointByLocalMatrix(crown_mesh.point(crown_bound[i])));
	}
	std::vector<std::pair<int, OpenMesh::Vec3d>>temp_bound_anchor, crown_bound_anchor;
	for (int i = 0; i < crown_bound2temp_vhsmap.size(); i++)
	{
		crown_bound_anchor.push_back(std::make_pair(i, crown_bound_points[i]));
		COpenMeshT::VertexHandle tb_vh = crown_bound2temp_vhsmap[i];
		bool flag = true;
		for (int j = 0; j < temp_bound.size(); j++)
		{
			if (tb_vh == temp_bound[j])
			{
				temp_bound_anchor.push_back(std::make_pair(j, temp_bound_points[j]));
				flag = false;
				break;
			}
		}
		if (flag)
		{
			std::cerr << "error in iter "<<i<<"tb_vh:"<<tb_vh.idx()<<" temp_bound_anchor.size£º "<< temp_bound_anchor.size()<<"  crown_bound_anchor.size£º "<< crown_bound_anchor.size() << std::endl;
			
		}
	}
	std::vector<std::pair<int, double>>temp2crown_correspond;
	CCurveBaseAlg::ComputeMatchingWithAnchorsFixed(temp_bound_points, temp_bound_anchor,crown_bound_points,crown_bound_anchor, temp2crown_correspond);

	std::map<COpenMeshT::VertexHandle, COpenMeshT::EdgeHandle>crown_new_vh_eh_map;
	std::map<COpenMeshT::FaceHandle, std::vector<COpenMeshT::VertexHandle>>crown_new_fh_vh_map;
	std::map<COpenMeshT::VertexHandle, COpenMeshT::VertexHandle>temp_bound2crown_bound_vhs_map;
	for (int i = 0; i < temp2crown_correspond.size(); i++)
	{
		
		if (mapped_temp_bound_set.find(temp_bound[i]) != mapped_temp_bound_set.end())
		{

			temp_bound2crown_bound_vhs_map[temp_bound[i]] = crown_bound[temp2crown_correspond[i].first];
			//std::cerr << "continue" << std::endl;
			continue;
		}
	//	std::cerr << temp2crown_correspond[i].first << " " << temp2crown_correspond[i].second << std::endl;
		COpenMeshT::VertexHandle vh0 = crown_bound[temp2crown_correspond[i].first];
		COpenMeshT::VertexHandle vh1 = crown_bound[(temp2crown_correspond[i].first + 1) % crown_bound.size()];
		COpenMeshT::HalfedgeHandle hh;
		if (CGeoBaseAlg::GetHalfedgeHandle(crown_mesh, vh0, vh1, hh))
		{
			if (crown_mesh.is_boundary(hh))
			{
				hh = crown_mesh.opposite_halfedge_handle(hh);
			}
		}
		else
		{
			std::cerr << "cannot find halfedge of vh0 and vh1" << std::endl;
		}
		COpenMeshT::FaceHandle fh = crown_mesh.face_handle(hh);

		OpenMesh::Vec3d p0 = crown_mesh.point(vh0);
		OpenMesh::Vec3d p1 = crown_mesh.point(vh1);
		OpenMesh::Vec3d p = (p1 - p0)*temp2crown_correspond[i].second + p0;
		COpenMeshT::VertexHandle new_vh=crown_mesh.add_vertex(p);
		temp_bound2crown_bound_vhs_map[temp_bound[i]] = new_vh;
		crown_new_vh_eh_map[new_vh] = crown_mesh.edge_handle(hh);
		if (crown_new_fh_vh_map.find(fh) == crown_new_fh_vh_map.end())
		{
			crown_new_fh_vh_map[fh] = std::vector<COpenMeshT::VertexHandle>();
		}
		crown_new_fh_vh_map[fh].push_back(new_vh);
	}
	for (auto iter = crown_new_fh_vh_map.begin(); iter != crown_new_fh_vh_map.end(); iter++)
	{
	
		CGeoAlg::SplitFaceByVhs(crown_mesh, iter->first, iter->second, std::vector<std::pair<COpenMeshT::VertexHandle, COpenMeshT::VertexHandle>>(), crown_new_vh_eh_map);
	}
	
	crown_mesh.garbage_collection();


	COpenMeshT res_mesh;
	std::vector<COpenMeshT::VertexHandle>seam_vhs;
	MergeBoundaryRegisteredMeshes(temp_mesh, teeth_template_obj.GetMatrix(), crown_mesh, mesh_crown_obj.GetMatrix(), temp_bound2crown_bound_vhs_map, res_mesh, seam_vhs);
	std::vector<COpenMeshT::VertexHandle>seam_nei_vhs;
	CGeoBaseAlg::GetNeighborVhs(res_mesh, seam_vhs, 4, seam_nei_vhs);
	CGeoAlg::SmoothMesh(res_mesh, seam_nei_vhs, 5);
	CGeoBaseAlg::GetNeighborVhs(res_mesh, seam_vhs, 5, seam_nei_vhs);
	CGeoAlg::Remeshing(res_mesh, seam_nei_vhs);
	CGeoAlg::SmoothMesh(res_mesh, seam_nei_vhs, 3);
	mesh_crown_obj.SetMesh(res_mesh);
	mesh_crown_obj.SetChanged();
	//temp_mesh = teeth_template_obj.GetMesh();
	mesh_crown_obj.SetMeshColor(pre_temp_mesh_color);
	teeth_template_obj.IsVisiable() = false;
	for (int i = 0; i < sep_meshes.size(); i++)
	{
		delete sep_meshes[i];
	}
	
}
void CDentalTemplateFitting::RefineCrownBoundaryAfterFitting(CMeshObject &mesh_crown_obj, CMeshObject &teeth_template_obj, std::map<COpenMeshT::VertexHandle, COpenMeshT::FaceHandle>&temp_crown_map)
{
	std::cerr << "RefineCrownBoundaryAfterFitting" << std::endl;
	COpenMeshT& mesh_crown = mesh_crown_obj.GetMesh();
	Eigen::Matrix4d mat = mesh_crown_obj.GetMatrix();
	std::vector<COpenMeshT::VertexHandle>to_del;
	for (auto viter = mesh_crown.vertices_begin(); viter != mesh_crown.vertices_end(); viter++)
	{
		OpenMesh::Vec3d p = mesh_crown_obj.TransformPointByLocalMatrix(mesh_crown.point(viter));
		COpenMeshT::FaceHandle res_fh,res_inv_fh;
		OpenMesh::Vec3d res_p,res_inv_p;
		double dis=CGeoAlg::ComputeClosestPoint(p, teeth_template_obj, res_fh, res_p);
		double inv_dis = CGeoAlg::ComputeClosestPoint(res_p, mesh_crown_obj, res_inv_fh, res_inv_p);
		if (dis < inv_dis / 4.0*3.0 || inv_dis < dis / 4.0*3.0)
		{
			to_del.push_back(viter);
		}

	}
	std::vector<std::vector<COpenMeshT::VertexHandle>>sep_vhs;
	CGeoAlg::SeparateVhsByConnectivity(to_del, mesh_crown, sep_vhs);
	to_del.clear();
	for (int i = 0; i < sep_vhs.size(); i++)
	{
		bool flag = false;
		for (int j = 0; j < sep_vhs[i].size(); j++)
		{
			if (mesh_crown.is_boundary(sep_vhs[i][j]))
			{
				flag = true;
				break;
			}
		}
		if (flag)
		{
			for (int j = 0; j < sep_vhs[i].size(); j++)
			{
				to_del.push_back(sep_vhs[i][j]);
			}

		}
	}
	for (auto iter = to_del.begin(); iter != to_del.end(); iter++)
	{
	
		mesh_crown.delete_vertex(*iter); 

	}
	mesh_crown.garbage_collection();
	std::vector<COpenMeshT*>res_meshes;
	std::vector<std::map<COpenMeshT::VertexHandle, COpenMeshT::VertexHandle>>vid_orig;
	CGeoAlg::SeparateDisconnectedParts(mesh_crown, res_meshes, vid_orig);
	int max_size = -1;
	int mi=0;
	for (int i = 0; i < res_meshes.size(); i++)
	{
		if (max_size < res_meshes[i]->n_vertices())
		{
			max_size = res_meshes[i]->n_vertices();
			mi = i;
		}
	}

	for (int i = 0; i < res_meshes.size(); i++)
	{
		if (i != mi)
		{
			for (auto viter = res_meshes[i]->vertices_begin(); viter != res_meshes[i]->vertices_end(); viter++)
			{
				mesh_crown.delete_vertex(vid_orig[i][viter]);
			}
		}
	}
	CGeoBaseAlg::RemoveNonManifold(mesh_crown);
	mesh_crown.garbage_collection();
	mesh_crown_obj.SetChanged();

}
//void CDentalTemplateFitting::RefineFittingTemplate(COpenMeshT &mesh_template, COpenMeshT &mesh_target)
//{
//	Eigen::MatrixXd temp_v, tgt_v;
//	Eigen::MatrixXi temp_f, tgt_f;
//	std::vector<COpenMeshT::FaceHandle>temp_idfh_map, tgt_idfh_map;
//	CConverter::ConvertFromOpenMeshToIGL(mesh_template, temp_v, temp_f, &temp_idfh_map);
//	CConverter::ConvertFromOpenMeshToIGL(mesh_target, tgt_v, tgt_f, &tgt_idfh_map);
//
//	std::cerr << "compute obb" << std::endl;
//	CPqpObb tmp_obb(temp_v, temp_f);
//	CPqpObb tgt_obb(tgt_v, tgt_f);
//	std::cerr << "compute obb finish" << std::endl;
//	std::map<OpenMesh::FaceHandle,std::vector<OpenMesh::Vec3d>>temp_face_add_points;
//	for (int i = 0; i < tgt_v.rows(); i++)
//	{
//		Eigen::Vector3d p = tgt_v.row(i);
//		CPqpObb::QueryResult res=tmp_obb.QueryClosestPoint(p);
//		if (temp_face_add_points.find(temp_idfh_map[res.fid_]) == temp_face_add_points.end())
//		{
//			temp_face_add_points[temp_idfh_map[res.fid_]] = std::vector<OpenMesh::Vec3d>();
//		}
//		Eigen::Vector3d cp = res.closest_pnt_;
//		temp_face_add_points[temp_idfh_map[res.fid_]].push_back(OpenMesh::Vec3d(cp.x(),cp.y(),cp.z()));
//	}
//	std::set<OpenMesh::VertexHandle>new_vhs;
//	for (auto iter = temp_face_add_points.begin(); iter != temp_face_add_points.end(); iter++)
//	{
//		std::vector<OpenMesh::VertexHandle>vhs;
//	//	CGeoAlg::SplitFaceByPoints(mesh_template, iter->first, iter->second, vhs);
//		std::vector<COpenMeshT::VertexHandle>fvhs;
//		COpenMeshT::FaceHandle fh = iter->first;
//		for (auto viter = mesh_template.fv_begin(fh); viter != mesh_template.fv_end(fh); viter++)
//		{
//			fvhs.push_back(viter);
//		}
//		mesh_template.delete_face(fh);
//		vhs.push_back(mesh_template.add_vertex(iter->second[0]));
//
//		mesh_template.add_face(vhs[0], fvhs[0], fvhs[1]);
//		mesh_template.add_face(vhs[0], fvhs[2], fvhs[0]);
//		mesh_template.add_face(vhs[0], fvhs[1], fvhs[2]);
//
//		for (int i = 0; i < vhs.size(); i++)
//		{
//			new_vhs.insert(vhs[i]);
//		}
//	}
//	for (auto iter = new_vhs.begin(); iter != new_vhs.end(); iter++)
//	{
//		OpenMesh::Vec3d p = mesh_template.point(*iter);
//		Eigen::Vector3d ep(p[0], p[1], p[2]);
//		CPqpObb::QueryResult res = tgt_obb.QueryClosestPoint(ep);
//		Eigen::Vector3d tp = res.closest_pnt_;
//		p = OpenMesh::Vec3d(tp.x(), tp.y(), tp.z());
//		mesh_template.point(*iter) = p;
//	}
//	mesh_template.garbage_collection();
//
//	std::cerr << "RefineFittingTemplate finish" << std::endl;
//
//}