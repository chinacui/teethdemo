#include <single_Teeth_Projection_Action.h>
#include <igl\readDMAT.h>
#include <igl\writeDMAT.h>
#include <algorithm> 
#include "qfiledialog.h"
#include "../DataColle/mesh_object.h"
#include "../DataColle/data_io.h"
#include "..\DataColle\data_pool.h"
#include "../AlgColle/dental_base_alg.h"
#include "../AlgColle/geo_alg.h"
#include "..\..\Src\AlgColle\curve_base_alg.h"


#include "ui_context.h"
#include "qtextstream.h"
#include "../AlgColle/geo_base_alg.h"
#include<boost\filesystem.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include<igl/writeOFF.h>
#include<igl/readOFF.h>


#include"../AlgColle/non_rigid_icp.h"
#include"../AlgColle/arap_deform.h"

#include <igl/colon.h>
#include <igl/harmonic.h>

#include "c:\Users\admin\Desktop\workspace\Src\DataColle\cgal_igl_converter.h"

#include <igl/grad.h>
#include <igl/massmatrix.h>
#include <igl/barycenter.h>
#include <igl/per_vertex_normals.h>

CSingleTeethProjectionAction* CSingleTeethProjectionAction::instance_ = nullptr;
int CSingleTeethProjectionAction::id_1_num_ = 3;
int CSingleTeethProjectionAction::id_2_num_ = 2;
int CSingleTeethProjectionAction::id_3_num_ = 1;
int CSingleTeethProjectionAction::id_4_num_ = 4;
int CSingleTeethProjectionAction::id_5_num_ = 1;
int CSingleTeethProjectionAction::id_6_num_ = 2;
int CSingleTeethProjectionAction::id_7_num_ = 3;

void CSingleTeethProjectionAction::Init()
{
}

CSingleTeethProjectionAction* CSingleTeethProjectionAction::GetInstance() {
	if (nullptr == instance_) {
		instance_ = new CSingleTeethProjectionAction();
	}
	return instance_;
}

void CSingleTeethProjectionAction::DeleteInstance() {
	if (nullptr != instance_) {
		delete instance_;
		instance_ = nullptr;
	}
}

CSingleTeethProjectionAction::CSingleTeethProjectionAction()
{
	type_ = CSingleTeethProjection;
}
void CSingleTeethProjectionAction::templateIdNum(int id_1_, int id_2_, int id_3_, int id_4_, int id_5_, int id_6_, int id_7_)
{
	id_1_num_ = id_1_;
	id_2_num_ = id_2_;
	id_3_num_ = id_3_;
	id_4_num_ = id_4_;
	id_5_num_ = id_5_;
	id_6_num_ = id_6_;
	id_7_num_ = id_7_;
}
void CSingleTeethProjectionAction::KeyPressEvent(QKeyEvent *e)
{
	//ui.setupUi(this);
	switch (e->key())
	{
	case Qt::Key_2:
	{
		QString path = QFileDialog::getOpenFileName(NULL, "load dental mesh", ".", "Files(*.dmat *.obj *.stl *.off )");
		if (path.length() == 0)
		{
			std::cerr << "unable to load the mesh!\n" << std::endl;
			break;
		}
		CMeshObject *meshobj = new CMeshObject();
		COpenMeshT &mesh = meshobj->GetMesh();
		Eigen::VectorXd scalars;
		std::map<int, int>&tags = meshobj->GetVertexTags();
		std::vector<int> gingiva_seg_mark_;
		std::cerr << path.toStdString() << std::endl;
		if (!CDataIO::ReadMesh(path.toStdString(), *meshobj))
		{
			std::cerr << "unable to load mesh\n" << std::endl;
		}
		else
		{
			path[path.length() - 3] = 'd';
			path[path.length() - 2] = 'm';
			path[path.length() - 1] = 'a';
			path = path + 't';
			if (igl::readDMAT(path.toStdString(), scalars))
			{
				std::map<int, int>&tags = meshobj->GetVertexTags();
				tags.clear();
				for (int i = 0; i < scalars.size(); i++)
				{
					tags[i] = scalars(i);
					if (tags[i] == -1) { gingiva_seg_mark_.push_back(-1); tag_mesh.push_back(tags[i]); }
					if (tags[i] != -1) { gingiva_seg_mark_.push_back(1); tag_mesh.push_back(tags[i]); }

				}
			}
			/*segment the teeth into individual and calculate the center of each tooth and the tooth tag*/
			std::cerr << "1" << std::endl;
			segTeethToIndividual(mesh, tag_mesh);
			for (auto iter = sep_meshes.begin(); iter != sep_meshes.end(); iter++)
			{
				COpenMeshT mesh_i = *iter->second;
				int count = 0;
				OpenMesh::Vec3d center_temp = OpenMesh::Vec3d(0, 0, 0);
				for (auto viter = mesh_i.vertices_begin(); viter != mesh_i.vertices_end(); viter++)
				{
					center_temp = center_temp + mesh_i.point(viter);
					count++;
				}
				teethCenter.push_back(center_temp / count);
				tagTeeth.push_back(iter->first);
			}
			std::cerr << "11" << std::endl;
			CCurveBaseAlg * a;
			std::vector<OpenMesh::Vec2d> curve_teeth_center_;
			OpenMesh::Vec2d teeth_center_2d;
			for (auto i = 0; i < teethCenter.size(); i++)
			{
				teeth_center_2d[0] = teethCenter[i][0];
				teeth_center_2d[1] = teethCenter[i][1];
				curve_teeth_center_.push_back(teeth_center_2d);
			}
			std::vector<double> ceoffs;
			a->PolynomialFitting(curve_teeth_center_, 3, ceoffs);
			ceoffs[0] = ceoffs[0] + 7;
			teeth_projection_ = teethProjection::GetInstance();
			teeth_projection_->setTeethMesh(sep_meshes);
			teeth_projection_->setCurveCeoffs(ceoffs);
			std::cerr << "111111" << std::endl;
			this->teeth_projection_->test();
              		
			std::cerr << "2" << std::endl;
			//test
			std::vector<OpenMesh::Vec3d> t;
			std::vector<OpenMesh::Vec2d> tt;
			OpenMesh::Vec3d t1;
			OpenMesh::Vec2d t0;
			for (auto i = this->projection_teeth_image_.begin(); i != this->projection_teeth_image_.end(); i++)
			{
				tt = i->second;
				for (auto j = 0; j != tt.size(); j++)
				{
					t1[0] = tt[j][0];
					t1[1] = tt[j][1];
					t1[2] = 0;
					t.push_back(t1);
				}

			}



			std::cerr << "3" << std::endl;
			std::map<int, COpenMeshT*>sep_meshes;
			std::map<int, std::map<COpenMeshT::VertexHandle, COpenMeshT::VertexHandle>> vid_org_map_;
			CGeoAlg::SeparateMeshByVertexTag(mesh, gingiva_seg_mark_, sep_meshes, vid_org_map_);
			CMeshObject* p_seg_tooth_mesh = new CMeshObject();
			for (auto iter = sep_meshes.begin(); iter != sep_meshes.end(); iter++)
			{
				int tid = iter->first;
				if (tid == 1)
				{

					//	CGeoBaseAlg::RemoveNonManifold(*iter->second);
					p_seg_tooth_mesh->GetMesh() = *iter->second;


					delete sep_meshes[tid];

					p_seg_tooth_mesh->SetMeshColor(OpenMesh::Vec3d(0.8, 0.8, 0.8));
					p_seg_tooth_mesh->SetChanged();
					p_seg_tooth_mesh->SetMeshColor(OpenMesh::Vec3d(0.8, 0.8, 0.8));

					//seg_tooth_mesh_id_ = DataPool::AddMeshObject(p_seg_tooth_mesh);
				}


			}

			std::cerr << "4" << std::endl;
			/*show the teeth curve on the screen*/
			std::vector<OpenMesh::Vec3d> poly_teeth_cruve;
			QFile fp("C:/polyfit.txt");
			if (fp.open(fp.ReadOnly))
			{
				int count = 0;
				OpenMesh::Vec3d point;
				while (!fp.atEnd())
				{
					if (count % 2 == 0)
					{
						QString lineString = QString(fp.readLine()).trimmed();
						point[0] = lineString.toDouble();
					}
					if (count % 2 == 1)
					{
						QString lineString = QString(fp.readLine()).trimmed();
						point[1] = lineString.toDouble();
						point[2] = 0;
						poly_teeth_cruve.push_back(point);
					}
					count++;
				}

			}



			std::vector<int>::iterator it;
			std::cerr << "5" << std::endl;
			/*project the teeth of the teethcurve*/
			bubblesort(tagTeeth.size());
			int teeth_select_id;
			for (teeth_select_id = 0; teeth_select_id < tagTeeth.size(); teeth_select_id++)
			{
				calProjectDir(teeth_select_id, teethCenter, poly_teeth_cruve);
				//calProjectPlane(teeth_select_id, teethCenter, tagTeeth, mesh, tags);
				//CGeoAlg::AlphaShape2d(projection_Plane_front, 0.5, vids);
				//calToothEdge(judge_volume_edge_front, projection_Plane_front, teeth_select_id, vids);
				//vids.clear();
				//CGeoAlg::AlphaShape2d(projection_Plane_back, 0.5, vids);
				//calToothEdge(judge_volume_edge_back, projection_Plane_back, teeth_select_id, vids);
				//vids.clear();
				//judge_volume_edge_front.clear();
				//judge_volume_edge_back.clear();
				//projection_Plane_front.clear();
				//projection_Plane_back.clear();
				//edge_Point.clear();
			}

			std::cerr << "6" << std::endl;
			/*std::vector<OpenMesh::Vec3d> temp;
			temp.push_back(teethCenter[0]);*/
			CCurveObject *center = new CCurveObject();
			OpenMesh::Vec3d color;
			color[0] = 200.0;
			color[1] = 200.0;
			color[2] = 0;
			center->SetCurve(t);
			center->RendereType() = CCurveObject::CurveType::Dots;
			center->SetColor(color);
			//co->SetChanged();
			DataPool::AddCurveObject(center);



			/*p_seg_tooth_mesh->SetChanged();
			int id = DataPool::AddMeshObject(p_seg_tooth_mesh);
			p_seg_tooth_mesh->RestoreCurrentVPos();
			CUIContext::SetSelectedMeshObjectId(id);
			p_seg_tooth_mesh->UseTexture() = false;
			p_seg_tooth_mesh->SetAttrChanged();*/
		}
		break;
	}
	case Qt::Key_3: {
		std::cerr << "the begin" << std::endl;
		for (int i = 0; i < 16; i++)
		{
			if (i != 3 && i != 8 && i != 9)
			{
				getCrownRootPoint(i);
			}
		}
		break;
	}

	case Qt::Key_4: {
		int num_teeth;
		this->edge_crown_root_top_point = this->teeth_projection_->getTest();
		//test
		/*for (auto i = this->edge_crown_root_top_point.begin(); i != this->edge_crown_root_top_point.end(); i++)
		{
			OpenMesh::Vec3d t = teethCenter[i->first] - i->second;
			OpenMesh::Vec3d t1;
			double t2 = OpenMesh::dot(t, project_dir[i->first]);
			t1 = i->second + t2*project_dir[i->first];
			i->second = t1;
		}*/
		//test



		QString path = QFileDialog::getOpenFileName(NULL, "load dental mesh", ".", "Files(*.dmat *.obj *.stl *.off )");
		if (path.length() == 0)
		{
			std::cerr << "unable to load the mesh!\n" << std::endl;
			break;
		}
		std::cerr << "path " << path.toStdString() << std::endl;
		for (auto i = 1; i < 5; i++)
		{
			if (i == 1) path[path.length() - 5] = '1';
			if (i == 2) path[path.length() - 5] = '2';
			if (i == 3) path[path.length() - 5] = '3';
			if (i == 4) path[path.length() - 5] = '4';
			std::cerr << path.toStdString() << std::endl;
			LoadRawData(path.toStdString(), i);
		}

		for (int i = 0; i < 16; i++)
		{
			//if (i != 3 && i != 8 && i != 9)
			{
				crownRegistration(i);
			}
		}
		break;
	}
	}
}

void CSingleTeethProjectionAction::segTeethToIndividual(COpenMeshT &mesh, std::vector<int>&tag_mesh)
{
	std::map<int, std::map<COpenMeshT::VertexHandle, COpenMeshT::VertexHandle>> vid_org_map_;
	CGeoAlg::SeparateMeshByVertexTag(mesh, tag_mesh, sep_meshes, vid_org_map_);
	for (auto iter = sep_meshes.begin(); iter != sep_meshes.end(); iter++)
	{
		if (vid_org_map_[iter->first].size() < 50 || vid_org_map_[iter->first].size() > 5000)
		{
			//delete sep_meshes[iter->first];
			sep_meshes.erase(iter->first);
		}
	}
	vid_org_map_.clear();
}
void CSingleTeethProjectionAction::calProjectDir(int teeth_select_id, std::vector<OpenMesh::Vec3d> teethCenter, std::vector<OpenMesh::Vec3d> poly_teeth_cruve)
{
	for (auto i = 0; i < poly_teeth_cruve.size(); i++)
	{
		if (teethCenter[teeth_select_id][1] - poly_teeth_cruve[i][1] <= 0.1)
		{
			project_tangent_temp[0] = poly_teeth_cruve[i + 1][0] - poly_teeth_cruve[i][0];
			project_tangent_temp[1] = poly_teeth_cruve[i + 1][1] - poly_teeth_cruve[i][1];
			project_tangent_temp[2] = poly_teeth_cruve[i + 1][2] - poly_teeth_cruve[i][2];
			break;
		}
	}
	project_tangent_temp = project_tangent_temp.normalized();
	if (project_tangent_temp[1] > 0)
	{
		project_tangent_temp = -1 * project_tangent_temp;
	}
	project_tangent.push_back(project_tangent_temp);

	project_dir_temp[0] = 1;
	project_dir_temp[1] = -1 * project_tangent_temp[0] * project_dir_temp[0] / project_tangent_temp[1];
	project_dir_temp[2] = 0;
	project_dir_temp = project_dir_temp.normalized();
	project_dir.push_back(project_dir_temp);

	project_z_temp = OpenMesh::cross(project_dir_temp, project_tangent_temp);
	project_z_temp = project_z_temp.normalized();
	if (project_z_temp[2] > 0)
	{
		project_z_temp = -1 * project_z_temp;
	}
	project_z.push_back(project_z_temp);
}
void CSingleTeethProjectionAction::calProjectPlane(int teeth_select_id, std::vector<OpenMesh::Vec3d> teethCenter, std::vector<int> tagTeeth, COpenMeshT mesh, std::map<int, int> tags)
{
	int count = 0;
	std::map<int, OpenMesh::Vec3d>single_Teeth;
	OpenMesh::Vec3d origin = teethCenter[teeth_select_id];
	bool flag = FALSE;
	for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); viter++)
	{
		if (tags[viter->idx()] == tagTeeth[teeth_select_id])
		{
			single_Teeth[count] = mesh.point(viter) - origin;
			for (auto it = mesh.vv_begin(viter); it != mesh.vv_end(viter); ++it)
			{
				//auto v = it.handle();
				if (tags[it->idx()] != tagTeeth[teeth_select_id])
				{
					judge_volume_edge[teeth_select_id].push_back(TRUE);
					flag = TRUE;
					break;
				}
			}
			if (flag == FALSE)
			{
				judge_volume_edge[teeth_select_id].push_back(FALSE);
			}
			flag = FALSE;
			count++;
		}
	}
	for (auto viter = 0; viter < single_Teeth.size(); viter++)
	{
		OpenMesh::Vec2d point;
		if (OpenMesh::dot(single_Teeth[viter], project_dir[teeth_select_id]) > 0)
		{
			point[0] = OpenMesh::dot(single_Teeth[viter], project_tangent[teeth_select_id]);
			point[1] = OpenMesh::dot(single_Teeth[viter], project_z[teeth_select_id]);
			projection_Plane_front.push_back(point);
			judge_volume_edge_front[teeth_select_id].push_back(judge_volume_edge[teeth_select_id][viter]);
		}
		if (OpenMesh::dot(single_Teeth[viter], project_dir[teeth_select_id]) <= 0)
		{
			point[0] = OpenMesh::dot(single_Teeth[viter], project_tangent[teeth_select_id]);
			point[1] = OpenMesh::dot(single_Teeth[viter], project_z[teeth_select_id]);
			projection_Plane_back.push_back(point);
			judge_volume_edge_back[teeth_select_id].push_back(judge_volume_edge[teeth_select_id][viter]);
		}
	}
}
void CSingleTeethProjectionAction::calToothEdge(std::map<int, std::vector<bool>>&judge_volume_edge_temp, std::vector<OpenMesh::Vec2d> projection_Plane, int teeth_select_id, std::vector<std::vector<int>>vids)
{
	QString str;
	str = "C:/Projection_edge" + QString::number(teeth_select_id, 10) + ".txt";
	QByteArray ba = str.toLatin1();
	char *path = ba.data();
	FILE *fp;
	if ((fp = fopen(path, "a+")) != NULL)  //判断文件是否已经被打开
		for (auto i = 0; i < vids[0].size(); i++)
		{
			if (judge_volume_edge_temp[teeth_select_id][vids[0][i]] == TRUE)
			{
				continue;
			}
			OpenMesh::Vec3d p;
			p[0] = projection_Plane[vids[0][i]][0];
			p[1] = projection_Plane[vids[0][i]][1];
			p[2] = 0.;
			fprintf(fp, "%f", p[0]);
			fprintf(fp, "\n");
			fprintf(fp, "%f", p[1]);
			fprintf(fp, "\n");
			edge_Point.push_back(p);
		}
	fclose(fp);
}
void CSingleTeethProjectionAction::getCrownRootPoint(int ii)
{
	std::vector<OpenMesh::Vec3d> edge_crown_root_Point;
	QString path = "C:/project_edge/imageToModelPoints0.txt";
	char s[10];
	itoa(ii, s, 10);
	if (ii<10)
	{
		path[path.length() - 5] = s[0];
	}
	else
	{
		path[path.length() - 5] = s[0];
		path[path.length() - 4] = s[1];
		path[path.length() - 3] = '.';
		path[path.length() - 2] = 't';
		path[path.length() - 1] = 'x';
		path = path + 't';
	}
	QFile fp(path);
	if (fp.open(fp.ReadOnly))
	{
		int count = 0;
		OpenMesh::Vec3d point;
		while (!fp.atEnd())
		{
			if (count % 2 == 0)
			{
				QString lineString = QString(fp.readLine()).trimmed();
				point[0] = lineString.toDouble();
			}
			if (count % 2 == 1)
			{
				QString lineString = QString(fp.readLine()).trimmed();
				point[1] = lineString.toDouble();
				point[2] = 0;
				edge_crown_root_Point.push_back(point);
			}
			count++;
		}

	}

	edge_crown_root_top_point[ii] = edge_crown_root_Point[0];
	for (auto i = 0; i < edge_crown_root_Point.size(); ++i)
	{
		edge_crown_root_Point[i] = edge_crown_root_Point[i][0] * project_tangent[ii] + edge_crown_root_Point[i][1] * project_z[ii];
		edge_crown_root_Point[i] = edge_crown_root_Point[i] + teethCenter[ii];
		id_edge_crown_root_Point[ii].push_back(edge_crown_root_Point[i]);
		for (auto j = 0; j < 16; j++)
		{
			if (edge_crown_root_Point[i][2] > edge_crown_root_top_point[ii][2]) edge_crown_root_top_point[ii] = edge_crown_root_Point[i];
		}
	}

/*	CCurveObject *co = new CCurveObject();
	OpenMesh::Vec3d color;
	color[0] = 255.0;
	color[1] = 0.0;
	color[2] = 0;
	co->SetCurve(edge_crown_root_Point);
	co->RendereType() = CCurveObject::CurveType::Dots;
	co->SetColor(color);
	//co->SetChanged();
	DataPool::AddCurveObject(co);*/
}
void CSingleTeethProjectionAction::bubblesort(int count)
{
	OpenMesh::Vec3d iTemp;
	int tagsTemp;
	for (int i = 1; i < count; i++)
	{
		for (int j = count-1; j>=i; j--)
		{
			if (teethCenter[j][1]>teethCenter[j-1][1])
			{
				iTemp = teethCenter[j - 1];
				teethCenter[j - 1] = teethCenter[j];
				teethCenter[j] = iTemp;
				tagsTemp = tagTeeth[j - 1];
				tagTeeth[j - 1] = tagTeeth[j];
				tagTeeth[j] = tagsTemp;
			}
		}
	}
}
bool CSingleTeethProjectionAction::LoadRawData(std::string dir_name, int id_kind)
{
	std::pair<Eigen::MatrixXd, Eigen::MatrixXi>mean_shapes;

	boost::filesystem::path path(dir_name);
	if (!boost::filesystem::exists(path))
	{
		return false;
	}
	/*	Eigen::MatrixXd tmp_N;*/
	igl::readOFF(dir_name, mean_shapes.first, mean_shapes.second);
	std::cerr << "load down" << std::endl;

	int num_template = 16;
	if (id_kind == 1)
	{
		for (auto i = 0; i < id_1_num_; i++)
		{
			mesh_template = new CMeshObject();
			template_map[i] = mesh_template;
			CConverter::ConvertFromIGLToOpenMesh(mean_shapes.first, mean_shapes.second, template_map[i]->GetMesh());
		}
		int temp_value = id_1_num_ + id_2_num_ + id_3_num_ + id_4_num_ + id_5_num_ + id_6_num_;
		for (auto i = temp_value; i < temp_value + id_7_num_; i++)
		{
			mesh_template = new CMeshObject();
			template_map[i] = mesh_template;
			CConverter::ConvertFromIGLToOpenMesh(mean_shapes.first, mean_shapes.second, template_map[i]->GetMesh());
		}
	}

	if (id_kind == 2)
	{
		for (auto i = id_1_num_; i < id_1_num_ + id_2_num_; i++)
		{
			mesh_template = new CMeshObject();
			template_map[i] = mesh_template;
			CConverter::ConvertFromIGLToOpenMesh(mean_shapes.first, mean_shapes.second, template_map[i]->GetMesh());
		}
		int temp_value = id_1_num_ + id_2_num_ + id_3_num_ + id_4_num_ + id_5_num_;
		for (auto i = temp_value; i < temp_value + id_6_num_; i++)
		{
			mesh_template = new CMeshObject();
			template_map[i] = mesh_template;
			CConverter::ConvertFromIGLToOpenMesh(mean_shapes.first, mean_shapes.second, template_map[i]->GetMesh());
		}
	}

	if (id_kind == 3)
	{
		for (auto i = id_1_num_ + id_2_num_; i < id_1_num_ + id_2_num_ + id_3_num_; i++)
		{
			mesh_template = new CMeshObject();
			template_map[i] = mesh_template;
			CConverter::ConvertFromIGLToOpenMesh(mean_shapes.first, mean_shapes.second, template_map[i]->GetMesh());
		}
		int temp_value = id_1_num_ + id_2_num_ + id_3_num_ + id_4_num_;
		for (auto i = temp_value; i < temp_value + id_5_num_; i++)
		{
			mesh_template = new CMeshObject();
			template_map[i] = mesh_template;
			CConverter::ConvertFromIGLToOpenMesh(mean_shapes.first, mean_shapes.second, template_map[i]->GetMesh());
		}
	}

	if (id_kind == 4)
	{
		for (auto i = id_1_num_ + id_2_num_ + id_3_num_; i < id_1_num_ + id_2_num_ + id_3_num_ + id_4_num_; i++)
		{
			std::cerr << i << std::endl;
			mesh_template = new CMeshObject();
			template_map[i] = mesh_template;
			CConverter::ConvertFromIGLToOpenMesh(mean_shapes.first, mean_shapes.second, template_map[i]->GetMesh());
		}
	}


	return true;

}
void CSingleTeethProjectionAction::crownRegistration(int id)
{
	std::cerr << id << std::endl;
	mesh_template = template_map[id];
	/*separate template crown from tooth*/
	std::map<int, COpenMeshT*>temp_crown_root_mesh;
	std::map<int, std::map<COpenMeshT::VertexHandle, COpenMeshT::VertexHandle>> temp_crown_root_map;
	if (temp_crown_root_mesh.size() != 0) temp_crown_root_mesh.clear();
	if (temp_crown_root_map.size() != 0) temp_crown_root_map.clear();
	segTemplateCrownRoot(temp_crown_root_mesh, temp_crown_root_map);


	CMeshObject *mesh_crown = new CMeshObject();
	CMeshObject *temp_mesh_crown = new CMeshObject();
	CARAPDeformation* crown_temp_arap_;
	CNonRigidICP* non_rigid_icp_;
	mesh_crown->GetMesh() = *sep_meshes[tagTeeth[id]];
	temp_mesh_crown->GetMesh() = *temp_crown_root_mesh[1];

	tempPositionInit((temp_mesh_crown->GetMesh()), (mesh_crown->GetMesh()), id);//Initlizae the postion of template according to the project axis
	if (temp_mesh_crown != NULL && mesh_crown != NULL)
	{
		std::map<COpenMeshT::VertexHandle, OpenMesh::Vec3d>deform_map;
		if (crown_temp_arap_ == NULL) delete crown_temp_arap_;
		if (non_rigid_icp_ == NULL) delete non_rigid_icp_;
		non_rigid_icp_ = new CNonRigidICP(&temp_mesh_crown->GetMesh(), &(mesh_crown->GetMesh()));
		crown_temp_arap_ = new CARAPDeformation(mesh_template->GetMesh());                                //note:errors in time
		for (auto i = temp_mesh_crown->GetMesh().vertices_begin(); i != temp_mesh_crown->GetMesh().vertices_end(); i++)
		{
			deform_map[temp_crown_root_map[1][i]] = temp_mesh_crown->GetMesh().point(i);
		}
		crown_temp_arap_->SetDeformMap(deform_map);
		crown_temp_arap_->Deform();
		for (auto viter = mesh_template->GetMesh().vertices_begin(); viter != mesh_template->GetMesh().vertices_end(); viter++)
		{
			mesh_template->GetMesh().point(viter) = mesh_template->GetMesh().point(viter) + OpenMesh::Vec3d(0, 0, 1);
		}
		TempToMeshDisplacement(mesh_template->GetMesh(), mesh_crown->GetMesh(), id);


		int teeth_id = id;
		topPointHarmonic(temp_mesh_crown->GetMesh(), temp_crown_root_map, id);
		meshSmooth(mesh_template->GetMesh());
		handlePointBoundary(mesh_template->GetMesh(), teeth_id);
		//boundPointHarmonic(mesh_template->GetMesh());
	}

	double i, j, k;
	i = rand() / double(RAND_MAX);
	j = rand() / double(RAND_MAX);
	k = rand() / double(RAND_MAX);
	mesh_template->SetMeshColor(OpenMesh::Vec3d(i, j, k));
	DataPool::AddMeshObject(mesh_template);
	mesh_template->SetAttrChanged();

}
void CSingleTeethProjectionAction::segTemplateCrownRoot(std::map<int, COpenMeshT*>&temp_crown_root_mesh, std::map<int, std::map<COpenMeshT::VertexHandle, COpenMeshT::VertexHandle>>&temp_crown_root_map)
{
	std::vector<int>tag_temp_crown;
	for (auto i = mesh_template->GetMesh().vertices_begin(); i != mesh_template->GetMesh().vertices_end(); i++)
	{
		if (mesh_template->GetMesh().point(i)[1] > 0)
			tag_temp_crown.push_back(1);
		else
			tag_temp_crown.push_back(-1);
	}
	CGeoAlg::SeparateMeshByVertexTag(mesh_template->GetMesh(), tag_temp_crown, temp_crown_root_mesh, temp_crown_root_map);
}

void CSingleTeethProjectionAction::tempPositionInit(COpenMeshT &template_mesh_crown, COpenMeshT& tooth_mesh, int id)
{
	OpenMesh::Vec3d temp_project_z = OpenMesh::Vec3d(0.0, 1, 0);
	temp_project_z = temp_project_z.normalized();
	OpenMesh::Vec3d temp_project_tangent = OpenMesh::Vec3d(-1.0, 0.0, -1.0);
	Eigen::Vector3d u0, v0;
	Eigen::Vector3d u2, v2, test;

	test[0] = temp_project_z[0];
	test[1] = temp_project_z[1];
	test[2] = temp_project_z[2];

	u0[0] = temp_project_z[0];
	u0[1] = temp_project_z[1];
	u0[2] = temp_project_z[2];

	v0[0] = temp_project_tangent[0];
	v0[1] = temp_project_tangent[1];
	v0[2] = temp_project_tangent[2];

	u2[0] = project_z[id][0];
	u2[1] = project_z[id][1];
	u2[2] = project_z[id][2];

	v2[0] = project_tangent[id][0];
	v2[1] = project_tangent[id][1];
	v2[2] = project_tangent[id][2];



	Eigen::Quaterniond q2 = Eigen::Quaterniond::FromTwoVectors(u0, u2);
	Eigen::Vector3d v1 = q2._transformVector(v0);
	Eigen::Quaterniond q1 = Eigen::Quaterniond::FromTwoVectors(v1, v2);
	Eigen::Quaterniond q = q2*q1;
	//q = q.normalized();

	Eigen::Vector3d temp_center = Eigen::Vector3d(0.0, 0.0, 0.0);

	for (auto i = template_mesh_crown.vertices_begin(); i != template_mesh_crown.vertices_end(); i++)
	{
		Eigen::Vector3d temporary;
		temporary[0] = template_mesh_crown.point(i)[0];
		temporary[1] = template_mesh_crown.point(i)[1];
		temporary[2] = template_mesh_crown.point(i)[2];
		temporary = q2._transformVector(temporary);
		temporary = q1._transformVector(temporary);
		template_mesh_crown.point(i)[0] = temporary[0];
		template_mesh_crown.point(i)[1] = temporary[1];
		template_mesh_crown.point(i)[2] = temporary[2];
	}

	int count = 0;
	for (auto i = template_mesh_crown.vertices_begin(); i != template_mesh_crown.vertices_end(); i++)
	{
		temp_center[0] = temp_center[0] + template_mesh_crown.point(i)[0];
		temp_center[1] = temp_center[1] + template_mesh_crown.point(i)[1];
		temp_center[2] = temp_center[2] + template_mesh_crown.point(i)[2];
		count++;
	}
	temp_center = temp_center / count;

	for (auto i = mesh_template->GetMesh().vertices_begin(); i != mesh_template->GetMesh().vertices_end(); i++)
	{
		Eigen::Vector3d temporary;
		temporary[0] = mesh_template->GetMesh().point(i)[0];
		temporary[1] = mesh_template->GetMesh().point(i)[1];
		temporary[2] = mesh_template->GetMesh().point(i)[2];
		temporary = q2._transformVector(temporary);
		temporary = q1._transformVector(temporary);
		mesh_template->GetMesh().point(i)[0] = temporary[0] + teethCenter[id][0] - temp_center[0];
		mesh_template->GetMesh().point(i)[1] = temporary[1] + teethCenter[id][1] - temp_center[1];
		mesh_template->GetMesh().point(i)[2] = temporary[2] + teethCenter[id][2] - temp_center[2];
	}

	for (auto i = template_mesh_crown.vertices_begin(); i != template_mesh_crown.vertices_end(); i++)
	{
		template_mesh_crown.point(i)[0] = template_mesh_crown.point(i)[0] + teethCenter[id][0] - temp_center[0];
		template_mesh_crown.point(i)[1] = template_mesh_crown.point(i)[1] + teethCenter[id][1] - temp_center[1];
		template_mesh_crown.point(i)[2] = template_mesh_crown.point(i)[2] + teethCenter[id][2] - temp_center[2];
	}
}
void CSingleTeethProjectionAction::topPointHarmonic(COpenMeshT &template_mesh_crown, std::map<int, std::map<COpenMeshT::VertexHandle, COpenMeshT::VertexHandle>>&temp_crown_root_map, int id)
{
	std::pair<Eigen::MatrixXd, Eigen::MatrixXi>mean_shapes;
	Eigen::MatrixXd V, U, V_bc, U_bc;
	Eigen::VectorXd Z;
	Eigen::MatrixXi F;
	Eigen::VectorXi b;
	COpenMeshT &mesh_test = mesh_template->GetMesh();
	double maxZ = -0.1;
	auto topPoint = mesh_test.vertices_begin();
	for (auto i = mesh_test.vertices_begin(); i != mesh_test.vertices_end(); i++)
	{
		if (maxZ < mesh_test.point(i)[2])
		{
			maxZ = mesh_test.point(i)[2];
			topPoint = i;
		}
	}
	/*
	double scale;
	Eigen::RowVector3d top_displacement;
	scale = (project_dir[id][0] * (mesh_test.point(topPoint)[0] - edge_crown_root_top_point[id][0]) +
		project_dir[id][1] * (mesh_test.point(topPoint)[1] - edge_crown_root_top_point[id][1]) +
		project_dir[id][2] * (mesh_test.point(topPoint)[2] - edge_crown_root_top_point[id][2])) /
		(pow(project_dir[id][0], 2) + pow(project_dir[id][1], 2) + pow(project_dir[id][2], 2));
	top_displacement[0] = scale*project_dir[id][0] + edge_crown_root_top_point[id][0] - mesh_test.point(topPoint)[0];
	top_displacement[1] = scale*project_dir[id][1] + edge_crown_root_top_point[id][1] - mesh_test.point(topPoint)[1];
	top_displacement[2] = scale*project_dir[id][2] + edge_crown_root_top_point[id][2] - mesh_test.point(topPoint)[2];
	*/

	Eigen::RowVector3d top_displacement;
	CConverter::ConvertFromOpenMeshToIGL(mesh_test, mean_shapes.first, mean_shapes.second);
	V = mean_shapes.first;
	F = mean_shapes.second;
	U = V;
	Eigen::VectorXi S(V.rows());
	for (auto i = 0; i < V.rows(); i++)
	{
		S(i) = -1;
	}

	for (auto iter = this->template_crown_displacement_.begin(); iter != this->template_crown_displacement_.end(); iter++){
		S(iter->first) = 1;
	}


	igl::colon<int>(0, V.rows() - 1, b);
	b.conservativeResize(std::stable_partition(b.data(), b.data() + b.size(),
		[&S](int i)->bool {return S(i) >= 0; }) - b.data());

	// Boundary conditions directly on deformed positions
	U_bc.resize(b.size(), V.cols());
	V_bc.resize(b.size(), V.cols());
	for (int bi = 0; bi < b.size(); bi++)
	{
		V_bc.row(bi) = V.row(b(bi));
		switch (S(b(bi)))
		{
		case 0:
			// Don't move handle 0
			U_bc.row(bi) = V.row(b(bi));
			break;
		case 1:
			// move handle 1 down
			top_displacement[0] = this->template_crown_displacement_[b(bi)][0];
			top_displacement[1] = this->template_crown_displacement_[b(bi)][1];
			top_displacement[2] = this->template_crown_displacement_[b(bi)][2];
			U_bc.row(bi) = V.row(b(bi)) + top_displacement;

			break;
		case 2:
		default:
			// move other handles forward
			U_bc.row(bi) = V.row(b(bi)) + Eigen::RowVector3d(0, 0, -25);
			break;
		}
	}
	Eigen::MatrixXd C(F.rows(), 3);
	Eigen::RowVector3d purple(80.0 / 255.0, 64.0 / 255.0, 255.0 / 255.0);
	Eigen::RowVector3d gold(255.0 / 255.0, 228.0 / 255.0, 58.0 / 255.0);
	for (int f = 0; f < F.rows(); f++)
	{
		if (S(F(f, 0)) >= 0 && S(F(f, 1)) >= 0 && S(F(f, 2)) >= 0)
		{
			C.row(f) = purple;
		}
		else
		{
			C.row(f) = gold;
		}
	}
	const Eigen::MatrixXd U_bc_anim = V_bc + (U_bc - V_bc);
	Eigen::MatrixXd D;
	Eigen::MatrixXd D_bc = U_bc_anim - V_bc;
	igl::harmonic(V, F, b, D_bc, 2, D);
	U = V + D;
	CConverter::ConvertFromIGLToOpenMesh(U, F, mesh_test);
}
bool CSingleTeethProjectionAction::TempToMeshDisplacement(COpenMeshT &template_mesh, COpenMeshT& tooth_mesh, int id)
{
	if (this->template_crown_displacement_.size() != NULL)  this->template_crown_displacement_.clear();
	OpenMesh::Vec3d template_crown_vector_;
	OpenMesh::Vec3d template_minilen_;
	OpenMesh::Vec3d temp_one, temp_two;
	double template_crown_len_;
	double template_crown_;
	tooth_mesh.request_vertex_normals();
	tooth_mesh.update_normals();
	template_mesh.request_vertex_normals();
	template_mesh.update_normals();
	for (auto viter_crown = tooth_mesh.vertices_begin(); viter_crown != tooth_mesh.vertices_end(); viter_crown++){
		int viter_crown_id;
		template_crown_len_ = 1e10;
		for (auto viter_template = template_mesh.vertices_begin(); viter_template != template_mesh.vertices_end(); viter_template++){
			template_crown_vector_ = tooth_mesh.point(viter_crown) - template_mesh.point(viter_template);
			temp_one[0] = template_mesh.normal(viter_template.handle()).data()[0];
			temp_one[1] = template_mesh.normal(viter_template.handle()).data()[1];
			temp_one[2] = template_mesh.normal(viter_template.handle()).data()[2];
	        
			temp_two[0] = tooth_mesh.normal(viter_crown.handle()).data()[0];
			temp_two[1] = tooth_mesh.normal(viter_crown.handle()).data()[1];
			temp_two[2] = tooth_mesh.normal(viter_crown.handle()).data()[2];
			
			template_crown_ = template_crown_vector_.length();
			if (template_crown_len_ > template_crown_){
				viter_crown_id = viter_template->idx();
				template_minilen_ = template_crown_vector_;
				template_crown_len_ = template_crown_;
			}
		}
		template_crown_displacement_[viter_crown_id] = template_minilen_;
	}

	return true;
}
bool CSingleTeethProjectionAction::meshSmooth(COpenMeshT &template_mesh)
{
	std::pair<Eigen::MatrixXd, Eigen::MatrixXi>mean_shapes;
	COpenMeshT &mesh_test = mesh_template->GetMesh();
	Eigen::MatrixXd V, U;
	Eigen::MatrixXi F;
	Eigen::SparseMatrix<double> L;
	CConverter::ConvertFromOpenMeshToIGL(mesh_test, mean_shapes.first, mean_shapes.second);
	V = mean_shapes.first;
	F = mean_shapes.second;
	igl::cotmatrix(V, F, L);
	// Alternative construction of same Laplacian
	Eigen::SparseMatrix<double> G, K;
	// Gradient/Divergence
	igl::grad(V, F, G);
	// Diagonal per-triangle "mass matrix"
	Eigen::VectorXd dblA;
	igl::doublearea(V, F, dblA);
	// Place areas along diagonal #dim times
	const auto & T = 1.*(dblA.replicate(3, 1)*0.5).asDiagonal();
	// Laplacian K built as discrete divergence of gradient or equivalently
	// discrete Dirichelet energy Hessian
	K = -G.transpose() * T * G;
	std::cerr << "|K-L|: " << (K - L).norm() << std::endl;
	// Use original normals as pseudo-colors
	Eigen::MatrixXd N;
	igl::per_vertex_normals(V, F, N);
	Eigen::MatrixXd C = N.rowwise().normalized().array()*0.5 + 0.5;
	// Initialize smoothing with base mesh
	U = V;

	for (int i = 0; i < 60; i++)
	{
		// switch case for interation
		// Recompute just mass matrix on each step
		Eigen::SparseMatrix<double> M;
		igl::massmatrix(U, F, igl::MASSMATRIX_TYPE_BARYCENTRIC, M);
		// Solve (M-delta*L) U = M*U
		const auto & S = (M - 0.001*L);
		Eigen::SimplicialLLT<Eigen::SparseMatrix<double > > solver(S);
		assert(solver.info() == Eigen::Success);
		U = solver.solve(M*U).eval();
		//CConverter::ConvertFromIGLToOpenMesh(U, F, mesh_test);
		// Compute centroid and subtract (also important for numerics)
		/*Eigen::VectorXd dblAA;
		igl::doublearea(U, F, dblAA);
		double area = 0.5*dblAA.sum();
		Eigen::MatrixXd BC;
		igl::barycenter(U, F, BC);
		Eigen::RowVector3d centroid(0, 0, 0);
		for (int i = 0; i<BC.rows(); i++)
		{
			centroid += 0.5*dblAA(i) / area*BC.row(i);
		}
		U.rowwise() -= centroid;
		// Normalize to unit surface area (important for numerics)
		U.array() /= sqrt(area);
		//U.array() = U.array()*sqrt(area);
		//U.rowwise() += centroid;*/
	}
	
	CConverter::ConvertFromIGLToOpenMesh(U, F, mesh_test);

	return true;
}
bool CSingleTeethProjectionAction::handlePointBoundary(COpenMeshT &template_mesh, int teeth_id)
{
	std::pair<Eigen::MatrixXd, Eigen::MatrixXi>mean_shapes;
	Eigen::MatrixXd V, U, V_bc, U_bc;
	Eigen::VectorXd Z;
	Eigen::MatrixXi F;
	Eigen::VectorXi b;
	COpenMeshT &mesh_test = mesh_template->GetMesh();
	double maxZ = -0.1;
	auto topPoint = mesh_test.vertices_begin();
	for (auto i = mesh_test.vertices_begin(); i != mesh_test.vertices_end(); i++)
	{
		if (maxZ < mesh_test.point(i)[2])
		{
			maxZ = mesh_test.point(i)[2];
			topPoint = i;
		}
	}

	for (auto i = this->edge_crown_root_top_point.begin(); i != this->edge_crown_root_top_point.end(); i++)
	{
		OpenMesh::Vec3d t = mesh_test.point(topPoint) - i->second;
		OpenMesh::Vec3d t1;
		double t2 = OpenMesh::dot(t, project_dir[i->first]);
		t1 = i->second + t2*project_dir[i->first];
		i->second = t1;
	}


	CConverter::ConvertFromOpenMeshToIGL(mesh_test, mean_shapes.first, mean_shapes.second);
	V = mean_shapes.first;
	F = mean_shapes.second;
	U = V;
	Eigen::VectorXi S(V.rows());
	for (auto i = 0; i < V.rows(); i++)
	{
		S(i) = -1;
	}
	for (auto iter = this->template_crown_displacement_.begin(); iter != this->template_crown_displacement_.end(); iter++) {
		S(iter->first) = 1;
	}
	S(topPoint->idx()) = 1;

	Eigen::RowVector3d top_displacement;
	top_displacement[0] = this->edge_crown_root_top_point[teeth_id][0] - mesh_test.point(topPoint)[0];
	top_displacement[1] = this->edge_crown_root_top_point[teeth_id][1] - mesh_test.point(topPoint)[1];
	top_displacement[2] = this->edge_crown_root_top_point[teeth_id][2] - mesh_test.point(topPoint)[2];

	igl::colon<int>(0, V.rows() - 1, b);
	b.conservativeResize(std::stable_partition(b.data(), b.data() + b.size(),
		[&S](int i)->bool {return S(i) >= 0; }) - b.data());

	// Boundary conditions directly on deformed positions
	U_bc.resize(b.size(), V.cols());
	V_bc.resize(b.size(), V.cols());
	for (int bi = 0; bi < b.size(); bi++)
	{
		V_bc.row(bi) = V.row(b(bi));
		switch (S(b(bi)))
		{
		case 0:
			// Don't move handle 0
			U_bc.row(bi) = V.row(b(bi));
			break;
		case 1:
			// move handle 1 down
			if (b(bi) == topPoint->idx())
			{
				U_bc.row(bi) = V.row(b(bi)) + top_displacement;
			}
			else
			{
				U_bc.row(bi) = V.row(b(bi));
			}

			break;
		case 2:
		default:
			// move other handles forward
			U_bc.row(bi) = V.row(b(bi)) + Eigen::RowVector3d(0, 0, -25);
			break;
		}
	}
	Eigen::MatrixXd C(F.rows(), 3);
	Eigen::RowVector3d purple(80.0 / 255.0, 64.0 / 255.0, 255.0 / 255.0);
	Eigen::RowVector3d gold(255.0 / 255.0, 228.0 / 255.0, 58.0 / 255.0);
	for (int f = 0; f < F.rows(); f++)
	{
		if (S(F(f, 0)) >= 0 && S(F(f, 1)) >= 0 && S(F(f, 2)) >= 0)
		{
			C.row(f) = purple;
		}
		else
		{
			C.row(f) = gold;
		}
	}
	const Eigen::MatrixXd U_bc_anim = V_bc + (U_bc - V_bc);
	Eigen::MatrixXd D;
	Eigen::MatrixXd D_bc = U_bc_anim - V_bc;
	igl::harmonic(V, F, b, D_bc, 2, D);
	U = V + D;
	CConverter::ConvertFromIGLToOpenMesh(U, F, mesh_test);
	return true;
}