#include"../TeethRootRecoAlg/panoramic_simulation.h"
#include"panoramic_simulation_test_action.h"
#include"action_manager.h"
#include "qfiledialog.h"
#include"../DataColle/mesh_object.h"
#include"../DataColle/data_io.h"
#include<igl/writeDMAT.h>
#include<igl/readDMAT.h>
#include"../DataColle/data_pool.h"
#include"visual_utils.h"
#include"../AlgColle/geo_alg.h"
#include"cmodelviewer.h"
#include"../AlgColle/dental_base_alg.h"
#include"../DataColle/data_io.h"
#include"ui_context.h"
#include"../AlgColle/non_rigid_icp.h"
#include"../DataColle/aux_geo_utils.h"
#include"../AlgColle/correspondence_builder.h"
#include"../AlgColle/compute_voronoi_diagram.h"
#include"../AlgColle/arap_deform.h" 
#include"../TeethRootRecoAlg/dental_template_fitting.h"
#include"../TeethRootRecoAlg/panoramic_fitting.h"
#include"texture_wrapper.h"
#include"cvmatandqimage.h"
void CPanoramicSimulationTestAction::Init()
{
	radius_ = 62;
	rot_degree_step_ = M_PI / 72.0;
	center_move_step_ = 1;
	radius_step_ = 1;

	is_selecting_pair_ = false;
	panoramic_pic_id_ = -1;
	pick_pts_on_panorama_.clear();
	pick_vhs_on_crown_.clear();

}
void CPanoramicSimulationTestAction::SavePoints()
{
	QString path = QFileDialog::getSaveFileName((QWidget*)CUIContext::GetMainWindow(), "save proj result");

	if (path.length() == 0)
	{
		std::cerr << "unable to open path\n" << std::endl;
		return;
	}

	std::ofstream ofstream(path.toStdString());
	for (auto iter= proj_coords_.begin();iter!= proj_coords_.end();iter++)
	{
		for (int i = 0; i < iter->second.size(); i++)
		{
			ofstream <<"v "<< iter->second[i].second[0] << " " << iter->second[i].second[1] << " " << 0 << std::endl;
		}
	}
	ofstream.close();
}
void CPanoramicSimulationTestAction::MousePressEvent(QMouseEvent *e)
{
	if (e->modifiers() == Qt::Modifier::CTRL || e->modifiers() == Qt::Modifier::ALT)
	{
		e->setAccepted(false);
	}
	else
	{
		if (is_selecting_pair_)
		{
			auto camera = viewer_->GetCamera();
			OpenMesh::Vec3d orig, dir;

			camera.ConvertClickToLine(e->pos(), orig, dir);
			if (pick_pts_on_panorama_.size()== pick_vhs_on_crown_.size())
			{
				CMeshObject *panoramic_mesh_obj = DataPool::GetMeshObject(panoramic_pic_id_);
				
				COpenMeshT::FaceHandle fh;
				OpenMesh::Vec3d res_bary_coord;
				if (CGeoAlg::RayMeshIntersection(orig, dir, *panoramic_mesh_obj, fh, res_bary_coord))
				{
					OpenMesh::Vec3d aabb_min, aabb_max;
					CGeoBaseAlg::ComputeAABB(*panoramic_mesh_obj, aabb_min, aabb_max);
					OpenMesh::Vec3d pick_p = CGeoBaseAlg::ComputeVPosFromBaryCoord(panoramic_mesh_obj->GetMesh(), fh, res_bary_coord);
					pick_p = panoramic_mesh_obj->TransformPointByLocalMatrix(pick_p);
					double x = (pick_p[0] /*- aabb_min[0]*/);
					double y = (pick_p[1] /*- aabb_min[1]*/);
					pick_pts_on_panorama_.push_back(OpenMesh::Vec2d(x, y));
					std::cerr << "pick panorama" << std::endl;
				}
			}
			else
			{
				int min_crown_id = -1;
				OpenMesh::VertexHandle min_vh;
				double min_dis= std::numeric_limits<double>::max();
				for (auto iter = proj_coords_.begin(); iter != proj_coords_.end(); iter++)
				{
				
					std::vector<std::pair<COpenMeshT::VertexHandle, OpenMesh::Vec2d>>&params = iter->second;
					
					
					for (int i = 0; i < params.size(); i++)
					{
						OpenMesh::Vec3d p(params[i].second[0] + param_rendering_offset_[0], params[i].second[1] + param_rendering_offset_[1], param_rendering_offset_[2]);
						OpenMesh::Vec3d diff = p - orig;
						double dis = diff.length();
						if (dis < min_dis)
						{
							min_dis = dis;
							min_crown_id = iter->first;
							min_vh = params[i].first;
						}


					}
				}
				pick_vhs_on_crown_.push_back(std::make_pair(min_crown_id,min_vh));
				std::cerr << "pick crown" << std::endl;
			}
			
			
		}
	}
}

void CPanoramicSimulationTestAction::MouseMoveEvent(QMouseEvent *e)
{
	if (e->modifiers() == Qt::Modifier::CTRL || e->modifiers() == Qt::Modifier::ALT)
	{
		e->setAccepted(false);
	}
}
void CPanoramicSimulationTestAction::MouseReleaseEvent(QMouseEvent *e)
{
	if (e->modifiers() == Qt::Modifier::CTRL || e->modifiers() == Qt::Modifier::ALT)
	{
		e->setAccepted(false);
	}
}
void CPanoramicSimulationTestAction::KeyReleaseEvent(QKeyEvent *e)
{
	switch (e->key())
	{
		case Qt::Key_Z:
		{
		if (e->isAutoRepeat())
			break;
		is_selecting_pair_ = false;
		break;
		}
	}
}
void CPanoramicSimulationTestAction::RotateFrame(std::vector<OpenMesh::Vec3d>&frame,OpenMesh::Vec3d axis, double degree)
{
	Eigen::Matrix4d rot_mat=CGeoBaseAlg::ComputeRotMat(axis, degree, OpenMesh::Vec3d(0, 0, 0));
	for (int i = 0; i < frame.size(); i++)
	{
		frame[i]=CGeoBaseAlg::Transform(rot_mat, frame[i]);

	}
}
void CPanoramicSimulationTestAction::ComputeFrameAndCenter()
{
	std::vector<OpenMesh::Vec3d>mesh_centers;
	for (auto iter = crowns_.begin(); iter != crowns_.end(); iter++)
	{
		CMeshObject *mesh_obj = iter->second;
		mesh_centers.push_back(mesh_obj->TransformPointByLocalMatrix(CGeoBaseAlg::ComputeMeshCenter(mesh_obj->GetMesh())));
	}
	OpenMesh::Vec3d frame_mean;

	std::vector<double>eg_values;


	CGeoAlg::PointSetPCA3D(mesh_centers, frame_mean, frame_, eg_values);
	
	OpenMesh::Vec3d tmp = frame_[2];
	frame_[2] = frame_[1];
	frame_[1] = tmp;
	frame_[2] = -frame_[2];
	

	for (int i = 0; i < 3; i++)
	{
		frame_[i].normalize();
	}

	center_ = OpenMesh::Vec3d(0, 0, 0);
	for (int i = 0; i < mesh_centers.size(); i++)
	{
		center_ += mesh_centers[i];
	}
	center_ /= mesh_centers.size();
	center_ = center_ -frame_[2] * 20;
	for (auto iter = crowns_.begin(); iter != crowns_.end(); iter++)
	{
		CMeshObject *mesh_obj = iter->second;

	}
	
}
void CPanoramicSimulationTestAction::KeyPressEvent(QKeyEvent *e)
{
	
	switch (e->key())
	{
		case Qt::Key_Z:
		{
			if (e->isAutoRepeat())
				break;
			is_selecting_pair_ = true;
			break;
		}
		case Qt::Key_L:
		{
			QString path = QFileDialog::getOpenFileName(NULL, "load dental mesh", ".", "Files(*.obj *.stl *.off )");

			if (path.length() == 0)
			{
				std::cerr << "unable to load mesh\n" << std::endl;
				break;
			}
			CMeshObject *meshobj = new CMeshObject();
			COpenMeshT &mesh = meshobj->GetMesh();
			std::cerr << path.toStdString() << std::endl;
			if (!CDataIO::ReadMesh(path.toStdString(), *meshobj))
			{
				std::cerr << "unable to load mesh\n" << std::endl;
			}
			else
			{
				//CDentalBaseAlg::PCABasedOrientationCorrection(meshobj->GetMesh());
				//CGeoBaseAlg::NormalizeMeshSize(meshobj->GetMesh());
				path[path.length() - 3] = 'd';
				path[path.length() - 2] = 'm';
				path[path.length() - 1] = 'a';
				path = path + 't';
				Eigen::VectorXd scalars;
				std::map<int, int>&tags = meshobj->GetVertexTags();
				if (igl::readDMAT(path.toStdString(), scalars))
				{
					std::map<int, int>&tags = meshobj->GetVertexTags();
					tags.clear();
					for (int i = 0; i < scalars.size(); i++)
					{
						tags[i] = scalars(i);

					}
				}
				else
				{
					std::vector<COpenMeshT*>sep_meshes;
					std::vector<std::map<COpenMeshT::VertexHandle, COpenMeshT::VertexHandle>> vid_orig;
					CGeoAlg::SeparateDisconnectedParts(meshobj->GetMesh(), sep_meshes, vid_orig);

					for (int i = 0; i < vid_orig.size(); i++)
					{
						for (auto iter = vid_orig[i].begin(); iter != vid_orig[i].end(); iter++)
						{
							int id = iter->second.idx();
							tags[id] = i;
						}
					}

				}
				for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); viter++)
				{

					OpenMesh::Vec3d color = CVisualUtils::GetRandColorByTag(tags[viter->idx()]);
					if (tags[viter->idx()] == -1)
						mesh.set_color(viter, OpenMesh::Vec3d(0.8, 0.8, 0.8));
					else
						mesh.set_color(viter, color);
				}
				segmented_jaws_.push_back(meshobj);
				meshobj->SetChanged();
				DataPool::AddMeshObject(meshobj);

			}
			
			break;
		}
		case Qt::Key_A://adjust center
		{
			center_ -= center_move_step_*frame_[0];
			GenPanorama();
			RenderAuxilliaryShape();
			break;
		}
		case Qt::Key_D://adjust center
		{
			center_ += center_move_step_*frame_[0];
			GenPanorama();
			RenderAuxilliaryShape();
			break;
		}
		case Qt::Key_W://adjust center
		{
			center_ -= center_move_step_*frame_[2];
			GenPanorama();
			RenderAuxilliaryShape();
			break;
		}
		case Qt::Key_S://adjust center
		{
			center_ += center_move_step_*frame_[2];
			GenPanorama();
			RenderAuxilliaryShape();
			break;
		}
		case Qt::Key_Left://adjust updir
		{
			OpenMesh::Vec3d rot_axis = frame_[2];
			RotateFrame(frame_, rot_axis, rot_degree_step_);
			GenPanorama();
			RenderAuxilliaryShape();
			break;
		}
		case Qt::Key_Right://adjust updir
		{
			OpenMesh::Vec3d rot_axis = frame_[2];
			RotateFrame(frame_, rot_axis, -rot_degree_step_);
			GenPanorama();
			RenderAuxilliaryShape();
			break;
		}
		case Qt::Key_Up://adjust updir
		{

			OpenMesh::Vec3d rot_axis = frame_[0];
			RotateFrame(frame_, rot_axis, rot_degree_step_);
			GenPanorama();
			RenderAuxilliaryShape();
			break;
		}
		case Qt::Key_Down://adjust updir
		{
			OpenMesh::Vec3d rot_axis = frame_[0];
			RotateFrame(frame_, rot_axis, -rot_degree_step_);
			GenPanorama();
			RenderAuxilliaryShape();
			break;
		}
	
		case Qt::Key_Comma://adjust radius
		{
			radius_ += radius_step_;
		
			GenPanorama();
			RenderAuxilliaryShape();
			break;
		}
		case Qt::Key_Period://adjust radius
		{
			radius_ -= radius_step_;
			
			GenPanorama();
			RenderAuxilliaryShape();
			break;
		}
		case Qt::Key_V:
		{
			std::cerr << "clear picking " << std::endl;
			pick_pts_on_panorama_.clear();
			pick_vhs_on_crown_.clear();
			
			break;
		}
		case Qt::Key_O:
		{
			SavePoints();
			break;
		}
		case Qt::Key_G://separate teeth
		{
		
		
			


			crowns_.clear();
			jaw_crowns_.clear();
			for (int i = 0; i < segmented_jaws_.size(); i++)
			{
				jaw_crowns_[i] = std::vector<int>();
				segmented_jaws_[i]->IsVisiable() = false;
				std::map<int, int>tags;
				tags = segmented_jaws_[i]->GetVertexTags();
				std::vector<int>vtags(segmented_jaws_[i]->GetMesh().n_vertices());
				std::map<int, COpenMeshT*>res_meshes;
				for (int i = 0; i < vtags.size(); i++)
				{
					vtags[i] = tags[i];
				}
				CGeoAlg::SeparateMeshByVertexTag(segmented_jaws_[i]->GetMesh(), vtags, res_meshes, std::map<int, std::map<COpenMeshT::VertexHandle, COpenMeshT::VertexHandle>>());
				

				for (auto iter = res_meshes.begin(); iter != res_meshes.end(); iter++)
				{
					if (iter->first != -1&&iter->second->n_vertices()>30)
					{
						CMeshObject *tmp_mesh_obj = new CMeshObject();

						tmp_mesh_obj->GetMesh() = *res_meshes[iter->first];
						tmp_mesh_obj->SetMeshColor(OpenMesh::Vec3d(0.8, 0.8, 0.8));
						tmp_mesh_obj->IsShinning() = false;
						DataPool::AddMeshObject(tmp_mesh_obj);
						tmp_mesh_obj->SetChanged();
						crowns_[tmp_mesh_obj->GetId()] = tmp_mesh_obj;
						jaw_crowns_[i].push_back(tmp_mesh_obj->GetId());
						delete res_meshes[iter->first];
					}

				}

			}
			ComputeFrameAndCenter();
			RenderAuxilliaryShape();

			break;
		}
		case Qt::Key_Space:
		{

			e->setAccepted(false);
			break;
		}
		case Qt::Key_Y://load panoramic image
		{
			QString path = QFileDialog::getOpenFileName(NULL, "load panoramic radiograph", ".", "Files(*.png *.jpg *.bmp )");

			if (path.length() != 0)
			{
				QImage panorama(path);
				std::cerr << "load panoramic radiograph " << path.toStdString() << std::endl;
				int tid = CUIContext::AddTexture(&panorama);
				double w = panorama.width();
				double h = panorama.height();

				if (w > h)
				{
					h = h / w*4;
					w = 1.0 * 4;
				}
				else
				{
					w = w / h * 4;
					h = 1.0 * 4;
				}
				std::cerr << "panoramic size " << w << " " << h << std::endl;
				COpenMeshT plain_mesh;
				CAuxGeoUtils::GetPlainMeshFromPointAndAxis(OpenMesh::Vec3d(0, 0, -10), OpenMesh::Vec3d(w, 0, 0), OpenMesh::Vec3d(0, h, 0), OpenMesh::Vec3d(0, 0, 1), 2, plain_mesh);
				for (auto viter = plain_mesh.vertices_begin(); viter != plain_mesh.vertices_end(); viter++)
				{
					OpenMesh::Vec3d p = plain_mesh.point(viter);
					COpenMeshT::TexCoord2D texcoord;

					for (int i = 0; i < 2; i++)
					{
						if (p[i] > 0)
						{
							texcoord[i] = 1;
						}
						else if (p[i] < 0)
						{
							texcoord[i] = 0;
						}
					}
					for (auto hiter = plain_mesh.vih_begin(viter); hiter != plain_mesh.vih_end(viter); hiter++)
					{
						plain_mesh.data(*hiter).SetUV(texcoord);
					}

				}
				CMeshObject *plain_obj = new CMeshObject();
				//plain_obj->IsPickAble() = false;
				plain_obj->GetMesh() = plain_mesh;
				plain_obj->SetChanged();
				plain_obj->TextureId() = tid;
				plain_obj->UseTexture() = true;
				panoramic_pic_id_=DataPool::AddMeshObject(plain_obj);
			}

			break;
		}
		case Qt::Key_C:
		{
			ComputeFrameAndCenter();
			RenderAuxilliaryShape();
			if (pick_vhs_on_crown_.size() == pick_pts_on_panorama_.size() && pick_vhs_on_crown_.size() >= 3)
			{
				std::cerr << "optimizing parameters of projection" << std::endl;
				std::vector<CMeshObject*>vec_crowns;
				for (auto iter = crowns_.begin(); iter != crowns_.end(); iter++)
				{
					CMeshObject *mesh_obj = iter->second;
					vec_crowns.push_back(mesh_obj);
				}
				std::cerr << "frame updir " << frame_[1] << "¡¡center: " << center_ << " radius: " << radius_ << std::endl;
				CCircularSurfaceProjector cp = CPanoramicSimulation::ConstructCircularProjector(vec_crowns, frame_[1], center_, radius_);
				std::vector<OpenMesh::Vec3d>picked_points;
				for (int i = 0; i < pick_vhs_on_crown_.size(); i++)
				{
					int mid = pick_vhs_on_crown_[i].first;
					COpenMeshT::VertexHandle vh = pick_vhs_on_crown_[i].second;
					OpenMesh::Vec3d p= crowns_[mid]->GetMesh().point(vh);
					p = crowns_[mid]->TransformPointByLocalMatrix(p);
					picked_points.push_back(p);
					
				}
				CMeshObject *pan_obj=DataPool::GetMeshObject(panoramic_pic_id_);
				if (pan_obj != NULL)
				{
					int tid = pan_obj->TextureId();
					TextureWraper *tw = CUIContext::GetTexture(tid);
					QImage& img=tw->GetImage();
					cv::Mat mat_img=QtOcv::image2Mat(img, CV_64FC1);
					CMeshObject *panoramic_obj = DataPool::GetMeshObject(panoramic_pic_id_);
					panoramic_obj->ApplyTransform();
					OpenMesh::Vec3d res_min_p,res_max_p;
					CGeoBaseAlg::ComputeAABB(panoramic_obj->GetMesh(), res_min_p, res_max_p);

					CPanoramicFittingOptimization pfoptor(cp, picked_points, pick_pts_on_panorama_, mat_img, res_max_p[0]- res_min_p[0], res_max_p[1] - res_min_p[1]);
					pfoptor.RunFitting();
					pfoptor.GetResProjector(cp);
					radius_ = cp.GetRadius();
					center_ = cp.GetCenter();
					OpenMesh::Vec3d updir = cp.GetUpDir();
					Eigen::Matrix4d rot_mat = CGeoBaseAlg::ComputeRotMat(frame_[1], updir);
					for (int i = 0; i < frame_.size(); i++)
					{
						frame_[i] = CGeoBaseAlg::Transform(rot_mat, frame_[i]);
						frame_[i].normalize();
					}
					GenPanorama();
					
					double res_scale;
					OpenMesh::Vec2d res_trans;
					pfoptor.GetPanoScaleAndTrans(res_scale, res_trans);
					
					COpenMeshT& mesh = panoramic_obj->GetMesh();
					for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); viter++)
					{
						OpenMesh::Vec3d p = mesh.point(viter);
						p[0] = p[0] * res_scale + res_trans[0] + param_rendering_offset_[0];
						p[1] = p[1] * res_scale + res_trans[1] + param_rendering_offset_[1];
						mesh.point(viter) = p;
					}

					panoramic_obj->SetChanged();
		
				}
			
				

			}
			else
			{
				GenPanorama();
			}
	
			RenderAuxilliaryShape();

			
			break;
		}
		case Qt::Key_E://gen panoramic
		{
	
			GenPanorama();
			RenderAuxilliaryShape();
			break;
		}
	}
}
void CPanoramicSimulationTestAction::GenPanorama()
{
	//center_ = center_ - frames[2] * 0.2;
	//CPanoramicSimulation::GeneratePanoramicImageByDentalArch(vec_crowns, frames,1000,cv::Mat());
	//std::cerr << "orig center " << center_ << std::endl;
	std::vector<CMeshObject*>vec_crowns;
	for (auto iter = crowns_.begin(); iter != crowns_.end(); iter++)
	{
		CMeshObject *mesh_obj = iter->second;
		vec_crowns.push_back(mesh_obj);
	}
	std::vector<std::vector<OpenMesh::Vec2d>>proj_coords;
	std::vector<std::vector<COpenMeshT::VertexHandle>>corres_vhs;
	//CPanoramicSimulation::GeneratePanoramicBoundPointsByCircle(vec_crowns, frame_[1], center_, radius_, proj_coords,corres_vhs);
	CPanoramicSimulation::GeneratePanoramicPointsByCircle(vec_crowns, frame_[1], center_, radius_, proj_coords, corres_vhs);
	
	
	proj_coords_.clear();
	for (int i = 0; i < vec_crowns.size(); i++)
	{
		int id = vec_crowns[i]->GetId();
		proj_coords_[id] = std::vector<std::pair<COpenMeshT::VertexHandle, OpenMesh::Vec2d>>();
		for (int j = 0; j < proj_coords[i].size(); j++)
		{
			proj_coords_[id].push_back(std::make_pair(corres_vhs[i][j], proj_coords[i][j]));
		}
	}
	double mminx = std::numeric_limits<double>::max();
	double mminy = std::numeric_limits<double>::max();
	for (auto iter = proj_coords_.begin(); iter != proj_coords_.end(); iter++)
	{
		for (int j = 0; j < iter->second.size(); j++)
		{
			OpenMesh::Vec2d p = iter->second[j].second;
			if (mminx > p[0])
			{
				mminx = p[0];
			}
			if (mminy > p[1])
			{
				mminy = p[1];
			}
		}
	}
	param_rendering_offset_[0] = -mminx+20;
	param_rendering_offset_[1] = -mminy;
	param_rendering_offset_[0] =0;
	param_rendering_offset_[1] = 0;
	
}
void CPanoramicSimulationTestAction::RenderAuxilliaryShape()
{
	std::vector<OpenMesh::Vec3d>frame_colors;
	frame_colors.push_back(OpenMesh::Vec3d(1, 0, 0));
	frame_colors.push_back(OpenMesh::Vec3d(0, 1, 0));
	frame_colors.push_back(OpenMesh::Vec3d(0, 0, 1));
	DataPool::DeleteAllCurveObjects();
	srand((unsigned)time(NULL));
	for (int i = 0; i < frame_.size(); i++)
	{
		CCurveObject *co = new CCurveObject();
		co->GetCurve().push_back(center_);
		co->GetCurve().push_back(frame_[i]*50+center_);
		co->SetChanged();
	
		co->SetColor(frame_colors[i]);
		DataPool::AddCurveObject(co);
	}
	CCurveObject *co = new CCurveObject();
	std::vector<OpenMesh::Vec3d>cc;
	CGeoAlg::SampleCircle(frame_[1], center_, radius_, rot_degree_step_, cc);
	co->SetCurve(cc);
	co->SetColor(OpenMesh::Vec3d(1, 1, 0));
	co->SetChanged();
	DataPool::AddCurveObject(co);




	for (auto iter= proj_coords_.begin();iter!= proj_coords_.end();iter++)
	{
		std::vector<OpenMesh::Vec3d>p3d;
		OpenMesh::Vec3d color;
		color[0] = rand() % 1000 / 1000.0;
		color[1] = rand() % 1000 / 1000.0;
		color[2] = rand() % 1000 / 1000.0;
		for (int j = 0; j < iter->second.size(); j++)
		{
	
			p3d.push_back(OpenMesh::Vec3d(iter->second[j].second[0] + param_rendering_offset_[0], iter->second[j].second[1] + param_rendering_offset_[1], param_rendering_offset_[2]));
		}
		CCurveObject *co = new CCurveObject();
		co->SetCurve(p3d);
		co->RendereType() = CCurveObject::CurveType::Dots;
		co->SetColor(color);
		co->SetChanged();
		DataPool::AddCurveObject(co);
	}
}
CPanoramicSimulationTestAction::CPanoramicSimulationTestAction()
{
	type_ = CPanoramicSimulationTest;
}