#include"teeth_reconstruction_test_action.h"
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
void CTeethReconstructionTestAction::Init()
{
	sel_crown_id_ = -1;
	sel_temp_id_ = -1;
	is_picking_ = true;
}
bool CTeethReconstructionTestAction::IsTemplateTeeth(int id)
{
	if (template_tooth_.find(id) != template_tooth_.end())
		return true;
	else
		return false;

}

void CTeethReconstructionTestAction::MousePressEvent(QMouseEvent *e)
{
	if (e->modifiers() == Qt::Modifier::CTRL || e->modifiers() == Qt::Modifier::ALT)
	{
		e->setAccepted(false);
	}
	else if (sel_temp_id_ != -1&&false)
	{
		CMeshObject *temp_obj = DataPool::GetMeshObject(sel_temp_id_);

		auto camera = viewer_->GetCamera();
		OpenMesh::Vec3d orig, dir;

		camera.ConvertClickToLine(e->pos(), orig, dir);
		auto data_pool = DataPool::GetMeshObjectPool();
		COpenMeshT::FaceHandle res_fh;
		OpenMesh::Vec3d bary_coord;
		if (CGeoAlg::RayMeshIntersection(orig, dir, *temp_obj, res_fh, bary_coord))
		{
			OpenMesh::Vec3d res_p=CGeoBaseAlg::ComputePointFromBaryCoord(temp_obj->GetMesh(), res_fh, bary_coord);
			res_p = temp_obj->TransformPointByLocalMatrix(res_p);
			test_sel_ps.push_back(std::make_pair(res_fh, res_p));
			if (test_sel_ps.size() >= 2)
			{
				std::vector<std::pair<COpenMeshT::EdgeHandle, OpenMesh::Vec3d>>res_eps;
				std::vector<bool>is_vertex;
				CGeoAlg::FindEdgePointsPath(*temp_obj, test_sel_ps.back(), test_sel_ps[test_sel_ps.size() - 2], res_eps, is_vertex,10);
				CCurveObject *co = new CCurveObject();
				std::vector<OpenMesh::Vec3d>cc;
				cc.push_back(test_sel_ps.back().second);
				for (int i = 0; i < res_eps.size(); i++)
				{
					cc.push_back(res_eps[i].second);
				}
				cc.push_back(test_sel_ps[test_sel_ps.size() - 2].second);
				co->SetCurve(cc);
				DataPool::AddCurveObject(co);
				std::cerr << "curve size " << co->GetCurve().size() << std::endl;
			}
		
		}
		
	}
	if (is_picking_)
	{
		auto camera = viewer_->GetCamera();
		OpenMesh::Vec3d orig, dir;
	
		camera.ConvertClickToLine(e->pos(), orig, dir);
		auto data_pool=DataPool::GetMeshObjectPool();
		std::cerr<<DataPool::GetMeshObjectPool().size() << std::endl;
		int meshid=CGeoAlg::PickMesh(orig, dir, DataPool::GetMeshObjectPool(), true);
	//	std::cerr << orig << " " << dir << std::endl;
		//std::cerr << meshid << std::endl;
	
		if (meshid != -1)
		{
			CMeshObject *meshobj = DataPool::GetMeshObject(meshid);
		
			bool is_template = IsTemplateTeeth(meshid);
			for (auto iter = data_pool.begin(); iter != data_pool.end(); iter++)
			{
				if (IsTemplateTeeth(iter->first) == is_template)
				{	
					iter->second->IsShinning() = false;
				}
			
			}
			
			if (is_template)
			{
				sel_temp_id_ = meshobj->GetId();
				for (auto iter= template_tooth_.begin();iter!= template_tooth_.end();iter++)
				{
					if (iter->second->GetId() == sel_temp_id_)
					{
						sel_temp_crown_id_ = crown_of_template_tooth_[iter->second->GetId()]->GetId();
						break;
					}
				}
				std::cerr << "sel_temp_crown_id " << sel_temp_crown_id_ << std::endl;
				std::cerr << "sel template" << std::endl;
			}
			else
			{
				sel_crown_id_ = meshobj->GetId();
				std::cerr << "sel crown" << std::endl;
			}
			
			meshobj->IsShinning() = true;

		}
	}

	
}
void CTeethReconstructionTestAction::MouseMoveEvent(QMouseEvent *e)
{
	if (e->modifiers() == Qt::Modifier::CTRL || e->modifiers() == Qt::Modifier::ALT)
	{
		e->setAccepted(false);
	}
}
void CTeethReconstructionTestAction::MouseReleaseEvent(QMouseEvent *e)
{
	if (e->modifiers() == Qt::Modifier::CTRL || e->modifiers() == Qt::Modifier::ALT)
	{
		e->setAccepted(false);
	}
}
void CTeethReconstructionTestAction::KeyPressEvent(QKeyEvent *e)
{
	if (e->modifiers() == Qt::Modifier::CTRL|| e->modifiers() == Qt::Modifier::ALT)
	{
		e->setAccepted(false);
	}
	switch (e->key())
	{
		case Qt::Key_Space:
		{

			e->setAccepted(false);
			break;
		}
		
		case Qt::Key_O:
		{//save selected crown and template

			QString path = QFileDialog::getSaveFileName((QWidget*)CUIContext::GetMainWindow(), "save sel models");

			if (path.length() == 0)
			{
				std::cerr << "unable to open path\n" << std::endl;
				return;
			}
			CMeshObject *crown_obj=DataPool::GetMeshObject(sel_crown_id_);
			CMeshObject *temp_obj = DataPool::GetMeshObject(sel_temp_id_);
			std::string fname = path.toStdString();
			std::string crownfname,tempfname;
			if (fname.length() - 4 >= 0 && fname[fname.length() - 4] != '.')
			{
				crownfname = fname + "_crown.obj";
				tempfname = fname+ "_temp.obj";
			}
			else
			{
				crownfname = fname.substr(0, fname.length() - 4) + "_crown.obj";
				tempfname = fname.substr(0, fname.length() - 4) + "_temp.obj";
			}
			if (crown_obj != NULL)
			{
				CDataIO::WriteMesh(crownfname, *crown_obj);
			}
			if (temp_obj != NULL)
			{
				//CGeoBaseAlg::NormalizeMeshSize(temp_obj->GetMesh());
				CDataIO::WriteMesh(tempfname, *temp_obj);
			}
			auto curve_pool=DataPool::GetCurveObjectPool();
			
			for (auto iter = curve_pool.begin(); iter != curve_pool.end(); iter++)
			{
				std::string fcurvename ;
				std::stringstream sstream;
				sstream << path.toStdString() + "_curve"<< iter->first << ".obj";
				fcurvename = sstream.str();
	
				CDataIO::SaveCurveToObj(fcurvename, iter->second->GetCurve());
			}
			break;
		}
		case Qt::Key_V:
		{
			QString path = QFileDialog::getOpenFileName(NULL, "load dental mesh", ".", "Files(*.obj *.stl *.off )");

			if (path.length() == 0)
			{
				std::cerr << "unable to load mesh\n" << std::endl;
				break;
			}
			CMeshObject *meshobj = new CMeshObject();
			//COpenMeshT &mesh = meshobj->GetMesh();
			//std::cerr << path.toStdString() << std::endl;
			if (!CDataIO::ReadMesh(path.toStdString(), *meshobj, OpenMesh::IO::Options::VertexColor))
			{
				std::cerr << "unable to load mesh\n" << std::endl;
			}
			DataPool::AddMeshObject(meshobj);
			template_tooth_[meshobj->GetId()]=meshobj;
			COpenMeshT &mesh = meshobj->GetMesh();
			std::vector<COpenMeshT::FaceHandle>roifhs;
		
			for (auto fiter = mesh.faces_begin(); fiter != mesh.faces_end(); fiter++)
			{
				bool flag = true;
				for (auto viter = mesh.fv_begin(fiter); viter != mesh.fv_end(fiter); viter++)
				{
					OpenMesh::Vec3d c = mesh.color(viter);
					if (c[0] > 0.1||c[1] > 0.1&&c[2] > 0.1)
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
	
			COpenMeshT roimesh;
			crown_to_template_fmap_[meshobj->GetId()]=std::map<COpenMeshT::FaceHandle, COpenMeshT::FaceHandle>();
			crown_to_template_vmap_[meshobj->GetId()] = std::map<COpenMeshT::VertexHandle, COpenMeshT::VertexHandle>();
			CGeoBaseAlg::ConvertFromOpenMeshROIToOpenMesh(mesh, roifhs, roimesh, &crown_to_template_fmap_[meshobj->GetId()]);

			for (auto iter = crown_to_template_fmap_[meshobj->GetId()].begin(); iter != crown_to_template_fmap_[meshobj->GetId()].end(); iter++)
			{
				COpenMeshT::FaceHandle fh_crown = iter->first;
				COpenMeshT::FaceHandle fh_temp = iter->second;
				auto cviter = roimesh.fv_begin(fh_crown);
				auto tviter = mesh.fv_begin(fh_temp);
				for (; cviter != roimesh.fv_end(fh_crown); tviter++,cviter++)
				{

					crown_to_template_vmap_[meshobj->GetId()][cviter] = tviter;
				}
			}
		
			CMeshObject *roimesh_obj = new CMeshObject();
			roimesh_obj->GetMesh() = roimesh;
			DataPool::AddMeshObject(roimesh_obj);
			
			crown_of_template_tooth_[meshobj->GetId()]=roimesh_obj;
			roimesh_obj->SetMeshColor(OpenMesh::Vec3d(0.8, 0.8, 0.8));
			roimesh_obj->SetChanged(true);

			meshobj->SetMeshColor(OpenMesh::Vec3d(0.8, 0.8, 0.8));
			meshobj->SetChanged();
	
			break;
			
		}
		case Qt::Key_C:
		{
			QString path = QFileDialog::getOpenFileName(NULL, "load dental mesh", ".", "Files(*.obj *.stl *.off )");

			if (path.length() == 0)
			{
				std::cerr << "unable to load mesh\n" << std::endl;
				break;
			}
			CCurveObject *ccobj = new CCurveObject();
			
			CDataIO::LoadCurveFromObj(path.toStdString(), ccobj->GetCurve());
			std::vector<OpenMesh::Vec3d>&cc = ccobj->GetCurve();
			for (int i = 0; i < cc.size(); i++)
			{
				cc[i][2] = -0.9;
			}
			ccobj->SetColor(OpenMesh::Vec3d(1, 0, 0));
			ccobj->SetCurve(cc);
			ccobj->SetChanged();
			//DataPool::DeleteCurveObject(sel_curve_id_);
			sel_curve_id_=DataPool::AddCurveObject(ccobj);
			break;
		}
		//case Qt::Key_G:
		//{
		//	CMeshObject *tmp_obj = DataPool::GetMeshObject(sel_temp_id_);
		//	CMeshObject *crown_obj = DataPool::GetMeshObject(sel_crown_id_);
		//	if (tmp_obj != NULL&&crown_obj!=NULL)
		//	{
		//	
		//		OpenMesh::Vec3d crown_mean= crown_obj->TransformPointByLocalMatrix(CGeoBaseAlg::ComputeMeshCenter(crown_obj->GetMesh()));
		//		OpenMesh::Vec3d temp_mean = tmp_obj->TransformPointByLocalMatrix(CGeoBaseAlg::ComputeMeshCenter(tmp_obj->GetMesh()));
		//		std::vector<OpenMesh::Vec3d>temp_frame;
		//		tmp_obj->Transform(crown_mean - temp_mean);
		//		tmp_obj->ApplyTransform();
		//		tmp_obj->SetChanged();
		//		/*CGeoAlg::ComputeMeshPCA(tmp_obj->GetMesh(), temp_mean, temp_frame);
		//		temp_mean = tmp_obj->TransformPointByLocalMatrix(temp_mean);
		//		for (int i = 0; i < temp_frame.size(); i++)
		//		{
		//			temp_frame[i] = tmp_obj->TransformPointByLocalMatrix(temp_frame[i]);
		//		}

		//		std::vector<OpenMesh::Vec3d>points(2);
		//		points[0] = temp_mean;
		//		for (int i = 0; i < temp_frame.size(); i++)
		//		{
		//			CCurveObject *ccobj = new CCurveObject();
		//			points[1] = temp_frame[i]+temp_mean;
		//			ccobj->SetCurve(points);
		//			ccobj->SetChanged(true);
		//			ccobj->SetColor(OpenMesh::Vec3d(1, 0, 0));
		//			DataPool::AddCurveObject(ccobj);
		//		}*/
		//	}
		//	
		//	break;
		//}
		case Qt::Key_T:
		{//load teeth template
			QString path = QFileDialog::getOpenFileName(NULL, "load dental mesh", ".", "Files(*.obj *.stl *.off )");

			if (path.length() == 0)
			{
				std::cerr << "unable to load mesh\n" << std::endl;
				break;
			}
			CMeshObject *meshobj = new CMeshObject();
			//COpenMeshT &mesh = meshobj->GetMesh();
			//std::cerr << path.toStdString() << std::endl;
			if (!CDataIO::ReadMesh(path.toStdString(), *meshobj))
			{
				std::cerr << "unable to load mesh\n" << std::endl;
			}
			CGeoAlg::FillHoles(meshobj->GetMesh(), true);
			meshobj->SetMeshColor(OpenMesh::Vec3d(0.8, 0.8, 0.8));
			meshobj->SetChanged();
			DataPool::AddMeshObject(meshobj);
			orig_tooth_template_ = meshobj;
			break;
		}
		case Qt::Key_M:
		{
			CMeshObject *sel_crown = DataPool::GetMeshObject(sel_crown_id_);
			CMeshObject *sel_temp = DataPool::GetMeshObject(sel_temp_id_);
			CMeshObject *sel_temp_crown = DataPool::GetMeshObject(sel_temp_crown_id_);
			
			std::cerr<<sel_crown->GetMesh().n_vertices() << std::endl;
			std::cerr << sel_temp_crown->GetMesh().n_vertices() << std::endl;
			if (sel_crown != NULL&&sel_temp != NULL&&sel_temp_crown!=NULL)
			{

				/*OpenMesh::Vec3d crown_mean = sel_crown->TransformPointByLocalMatrix(CGeoBaseAlg::ComputeMeshCenter(sel_crown->GetMesh()));
				OpenMesh::Vec3d temp_mean = sel_temp->TransformPointByLocalMatrix(CGeoBaseAlg::ComputeMeshCenter(sel_temp->GetMesh()));
				std::vector<OpenMesh::Vec3d>temp_frame;
				sel_temp->Transform(crown_mean - temp_mean);
				sel_temp->ApplyTransform();
				sel_temp->SetChanged();*/
				if (non_rigid_icp_ != NULL)
					delete non_rigid_icp_;
				if (template_arap_ != NULL)
					delete template_arap_;
				non_rigid_icp_=new CNonRigidICP(&(sel_temp_crown->GetMesh()),&(sel_crown->GetMesh()));
				template_arap_ = new CARAPDeformation(sel_temp->GetMesh());
				/*std::vector<COpenMeshT::VertexHandle>vec_deform_handle;
				for(auto iter= crown_to_template_vmap_[sel_temp_id_].begin(); iter != crown_to_template_vmap_[sel_temp_id_].end(); iter++)
				{
					vec_deform_handle.push_back(iter->second);
				}

				template_arap_->SetDeformHandle(vec_deform_handle);*/
				sel_temp_crown->IsVisiable() = false;
				sel_temp_crown->SetChanged();
				sel_crown->SetChanged();
			}
			
			break;
		}
		case Qt::Key_U:
		{
			break;
		}
		case Qt::Key_B:
		{
			CMeshObject *sel_crown = DataPool::GetMeshObject(sel_crown_id_);
			CMeshObject *sel_temp = DataPool::GetMeshObject(sel_temp_id_);
			CMeshObject *sel_temp_crown = DataPool::GetMeshObject(sel_temp_crown_id_);
			
			static bool is_finished = false;
			if ((!is_finished)&&sel_crown != NULL&&sel_temp != NULL&&sel_temp_crown != NULL&&non_rigid_icp_!=NULL)
			{
				is_finished = non_rigid_icp_->Run(1);
			
				
		
				sel_crown->SetChanged();
				sel_temp_crown->SetChanged();


				std::map<COpenMeshT::VertexHandle, OpenMesh::Vec3d>deform_map;
				COpenMeshT& temp_crown_mesh = sel_temp_crown->GetMesh();
				for (auto iter = crown_to_template_vmap_[sel_temp_id_].begin(); iter != crown_to_template_vmap_[sel_temp_id_].end(); iter++)
				{
					deform_map[iter->second] = temp_crown_mesh.point(iter->first);
				}
			
				template_arap_->SetDeformMap(deform_map);
				template_arap_->Deform();
				
				//std::cerr << "start change" << std::endl;


				if(is_finished)
				{
					//CGeoAlg::Remeshing(sel_temp->GetMesh());
					////CDentalTemplateFitting::RefineFittingTemplate(sel_temp->GetMesh(),sel_crown->GetMesh());
					////sel_crown->SetMeshColor(OpenMesh::Vec3d(0.8, 0.8, 0.8));
					std::map<COpenMeshT::VertexHandle, COpenMeshT::FaceHandle>vh_fh_map;
					non_rigid_icp_->GetTgt2SrcMap(vh_fh_map);
					CDentalTemplateFitting::RefineCrownBoundaryAfterFitting(*sel_crown, *sel_temp, vh_fh_map);
					/*CGeoAlg::SelfIntersectionRemoval(sel_temp->GetMesh());
					CGeoAlg::FillHoles(sel_temp->GetMesh(), false);*/
					CDentalTemplateFitting::ReplaceTemplateCrown(*sel_temp, *sel_crown);
					
					std::cerr << "finish" << std::endl;
					
				}
				
				sel_temp->SetChanged();

			
				std::cerr << "end change" << std::endl;
			}
			break;
		}
		case Qt::Key_H:
		{
			for (auto iter= crowns_.begin();iter!= crowns_.end();iter++)
			{
				if (sel_crown_id_ != -1&&sel_crown_id_ != iter->second->GetId())
				{
					iter->second->IsVisiable() = !iter->second->IsVisiable();
				}
			}
			for (auto iter=template_tooth_.begin();iter!=template_tooth_.end();iter++)
			{
				if (sel_temp_id_ != -1&&sel_temp_id_ != iter->first)
				{
					iter->second->IsVisiable() = !iter->second->IsVisiable();
				}
			}
			break;
		}
		case Qt::Key_S://separate teeth
		{
			if (orig_tooth_template_ != NULL)
			{
				orig_tooth_template_->IsVisiable() = false;
			}
			if (segmented_crown_ != NULL)
			{
				segmented_crown_->IsVisiable() = false;
			}
			
			
			std::vector<COpenMeshT*>tmp;
			if (orig_tooth_template_ != NULL)
			{
				CGeoAlg::SeparateDisconnectedParts(orig_tooth_template_->GetMesh(), tmp, std::vector<std::map<COpenMeshT::VertexHandle, COpenMeshT::VertexHandle>>());
				template_tooth_.clear();
				for (int i = 0; i < tmp.size(); i++)
				{
					CMeshObject *tmp_mesh_obj = new CMeshObject();
				
					
					tmp_mesh_obj->GetMesh() = *tmp[i];
					DataPool::AddMeshObject(tmp_mesh_obj);
					tmp_mesh_obj->SetChanged();
					template_tooth_[tmp_mesh_obj->GetId()] = tmp_mesh_obj;
					delete tmp[i];
				}
				//template_tooth_[0]->IsShinning() = true;
			}
			
			
			if (segmented_crown_ != NULL)
			{
				std::map<int, int>tags;
				tags = segmented_crown_->GetVertexTags();
				std::vector<int>vtags(segmented_crown_->GetMesh().n_vertices());
				std::map<int, COpenMeshT*>res_meshes;
				for (int i = 0; i < vtags.size(); i++)
				{
					vtags[i] = tags[i];
				}
				CGeoAlg::SeparateMeshByVertexTag(segmented_crown_->GetMesh(), vtags, res_meshes, std::map<int, std::map<COpenMeshT::VertexHandle, COpenMeshT::VertexHandle>>());
				crowns_.clear();

				for (auto iter = res_meshes.begin(); iter!=res_meshes.end(); iter++)
				{
					if (iter->first != -1)
					{
						CMeshObject *tmp_mesh_obj = new CMeshObject();
					
						tmp_mesh_obj->GetMesh() = *res_meshes[iter->first];
						tmp_mesh_obj->SetMeshColor(OpenMesh::Vec3d(0.8, 0.8, 0.8));
						tmp_mesh_obj->IsShinning() = false;
						DataPool::AddMeshObject(tmp_mesh_obj);
						tmp_mesh_obj->SetChanged();
						crowns_[tmp_mesh_obj->GetId()]=tmp_mesh_obj;
						delete res_meshes[iter->first];
					}
					
				}
				
			}
			
			
			break;
		}
		case Qt::Key_A://select teeth
		{
			is_picking_ = true;
		
			break;
		}
		case Qt::Key_F://fitting polynomial
		{
			std::vector<Eigen::Matrix4d>trans_mats;
			std::vector<CMeshObject*>crown_vec;
			for (auto iter = crowns_.begin(); iter != crowns_.end(); iter++)
			{
				crown_vec.push_back(iter->second);
			}
			CDentalTemplateFitting::ComputeStretchCrowns2LineMatrix(crown_vec, trans_mats);
			int count = 0;
			for (auto iter = crowns_.begin(); iter != crowns_.end(); iter++)
			{
				std::cerr << "mat:" << std::endl << trans_mats[count] << std::endl;
				crown_vec[count]->GetMatrix() = (trans_mats[count] * crown_vec[count]->GetMatrix());
				crown_vec[count]->ApplyTransform();
				crown_vec[count]->SetChanged();
				count++;
			}
			break;
			std::map<int,OpenMesh::Vec2d>center_points;
			std::map<int,OpenMesh::Vec3d>center_points_3d;
			double mean_y = 0;
			double mean_z = 0;
			double min_x = std::numeric_limits<double>::max();
			double max_x = std::numeric_limits<double>::min();
			//std::vector<OpenMesh::Vec3d>curve0;
			
			for (auto iter=crowns_.begin();iter!=crowns_.end();iter++)
			{
				OpenMesh::Vec3d center= iter->second->TransformPointByLocalMatrix(CGeoBaseAlg::ComputeMeshCenter(iter->second->GetMesh()));
				center_points[iter->first]=OpenMesh::Vec2d(center[0],center[2]);
				center_points_3d[iter->first]=center;
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
			mean_y /= crowns_.size();
			mean_z /= crowns_.size();
			std::vector<double>coeffs;
			std::vector<OpenMesh::Vec2d>vec_center_points;
			for (auto iter = center_points.begin(); iter != center_points.end(); iter++)
			{
				vec_center_points.push_back(iter->second);
			}
			CCurveBaseAlg::PolynomialFitting(vec_center_points,4, coeffs);
			std::vector<OpenMesh::Vec2d>curve;
			std::vector<OpenMesh::Vec3d>curve3d;
			for (double xi = min_x-0.02; xi <= max_x+0.02; xi =xi+0.01)
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
			double curve_len=CCurveBaseAlg::ComputeLenOfCurve(curve);
	
			std::vector<double>curve_lens(curve.size(),0);
			curve_lens[0] = 0;
			for (int i = 0; i < curve.size(); i++)
			{
				int prei = i - 1 >= 0 ? i - 1 : i;
				curve_lens[i]= curve_lens[prei]+((curve[i]-curve[prei]).length());
			}
			for (auto iter=center_points.begin();iter!=center_points.end();iter++)
			{
				int closest_pid;
				CCurveBaseAlg::ComputeClosestPoint(curve, iter->second, closest_pid);
				OpenMesh::Vec2d new_centerpoint = OpenMesh::Vec2d(curve_lens[closest_pid] - curve_len / 2.0, mean_z);
				OpenMesh::Vec3d trans(new_centerpoint[0] - iter->second[0], 0, new_centerpoint[1] - iter->second[1]);
				crowns_[iter->first]->Transform(trans);
				OpenMesh::Vec2d normal_dir = CCurveBaseAlg::ComputeNormalOfPolynomial(coeffs, curve[closest_pid][0]);
				OpenMesh::Vec3d rot_axis=OpenMesh::cross(OpenMesh::Vec3d(normal_dir[0], 0, normal_dir[1]), OpenMesh::Vec3d(0, 0, 1));
				double rot_angle = std::acos(OpenMesh::dot(normal_dir,OpenMesh::Vec2d(0,1)));
				crowns_[iter->first]->Rotate(rot_axis, rot_angle,OpenMesh::Vec3d(new_centerpoint[0], mean_y,new_centerpoint[1]));
				crowns_[iter->first]->ApplyTransform();
				crowns_[iter->first]->SetChanged();
				std::cerr << "trans " << trans << std::endl;
			/*	CCurveObject *tmp_cc = new CCurveObject();
				tmp_cc->GetCurve().push_back(center_points_3d[i]);
				tmp_cc->GetCurve().push_back(OpenMesh::Vec3d(center_points_3d[i][0]+normal_dir[0]*0.2, center_points_3d[i][1], center_points_3d[i][2]+normal_dir[1]*0.2));
				DataPool::AddCurveObject(tmp_cc);*/
			}
			
			break;
		}
		case Qt::Key_Y:
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
					h = h / w;
					w = 1.0;
				}
				else
				{
					w = w / h;
					h = 1.0;
				}
				std::cerr << "panoramic size " << w << " " << h << std::endl;
				COpenMeshT plain_mesh;
				CAuxGeoUtils::GetPlainMeshFromPointAndAxis(OpenMesh::Vec3d(0, 0, -1), OpenMesh::Vec3d(w, 0, 0), OpenMesh::Vec3d(0, h, 0), OpenMesh::Vec3d(0, 0, 1), 2, plain_mesh);
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
				plain_obj->IsPickAble() = false;
				plain_obj->GetMesh() = plain_mesh;
				plain_obj->SetChanged();
				plain_obj->TextureId() = tid;
				plain_obj->UseTexture() = true;
				DataPool::AddMeshObject(plain_obj);
			}
			
			break;
		}
		case Qt::Key_K:
		{
			CMeshObject *mesh_obj = DataPool::GetMeshObject(sel_temp_id_);
			CCurveObject *cc_obj = DataPool::GetCurveObject(sel_curve_id_);
			if (/*mesh_obj != NULL&&*/cc_obj != NULL)
			{
				std::vector<OpenMesh::Vec2d>curve_2d,res_skel_points;
				std::vector<int>res_root_pids,convex_pids;
				for (int i = 0; i < cc_obj->GetCurve().size(); i++)
				{
					curve_2d.push_back(OpenMesh::Vec2d(cc_obj->GetCurve()[i][0], cc_obj->GetCurve()[i][1]));
				}
				
				CCurveBaseAlg::ComputLocalMinimalConcavityPoints(curve_2d,10, 10, convex_pids);

				CDentalBaseAlg::DetectRootOfTeethSilhouette(curve_2d, res_root_pids);
				CCurveObject*cc = new CCurveObject();
				for (int i = 0; i < res_root_pids.size(); i++)
				{
					if (convex_pids.size() >=res_root_pids.size())
					{
						int dis = std::numeric_limits<int>::max();
						int pid = 0;
						for (int j = 0; j < convex_pids.size(); j++)
						{
							if (dis > std::abs(convex_pids[j] - res_root_pids[i]))
							{
								dis = std::abs(convex_pids[j] - res_root_pids[i]);
								pid = convex_pids[j];
							}
						}
						
						res_root_pids[i] = pid;
					}
					
					OpenMesh::Vec3d tmp_p= cc_obj->GetCurve()[res_root_pids[i]];
					tmp_p[2] = 0.6;
					cc->GetCurve().push_back(tmp_p);
				}
				cc->SetColor(OpenMesh::Vec3d(1, 1, 0));
				cc->RendereType() = CCurveObject::Dots;
				cc->SetChanged();
				DataPool::AddCurveObject(cc);

				//CVoronoiDiagramWrapper voronoi_wrapper(curve_2d);
				//voronoi_wrapper.constructVoronoiDiagramPoints();
				//std::vector<std::vector<std::pair<OpenMesh::Vec2d, double> > > branches_vertices;
				//std::vector<std::vector<int> > branches_associated_v_sites_index; // this is the voronoi sites associated with branch vertices
				//voronoi_wrapper.MainGetSkeletonBranches(branches_vertices, branches_associated_v_sites_index);
				//for (int i = 0; i < branches_vertices.size(); i++)
				//{
				//	CCurveObject *tmp_c = new CCurveObject();
				//	for (int j = 0; j < branches_vertices[i].size(); j++)
				//	{
				//		OpenMesh::Vec2d pv=branches_vertices[i][j].first;
				//		tmp_c->GetCurve().push_back(OpenMesh::Vec3d(pv[0], pv[1], 0));
				//	}
				//	tmp_c->SetColor(OpenMesh::Vec3d(0, 1, 0));
				//	tmp_c->SetChanged();
				//	DataPool::AddCurveObject(tmp_c);
				//}
			}
			break;
		}
		case Qt::Key_R://generate template
		{
		
			QString path = QFileDialog::getExistingDirectory();

			QString dir = path+QString("\\template");
			std::cerr << dir.toStdString() << std::endl;
			QDir qdir;
			qdir.mkpath(dir);
			std::string std_dir = dir.toStdString();
			CDataIO data_io;
			int count = 0;
			for (auto iter=crowns_.begin();iter!=crowns_.end();iter++)
			{
				CMeshObject *crown_mesh = iter->second;
				OpenMesh::Vec3d mean;
				std::vector<OpenMesh::Vec3d>axis;
				CGeoAlg::ComputeMeshPCA(crown_mesh->GetMesh(), mean, axis);

				OpenMesh::Vec3d updir(0, 1, 0);
				OpenMesh::Vec3d rot_axis = OpenMesh::cross(axis[0], updir);
				double degree = CGeoBaseAlg::ComputeAngleDegreeOfVector(axis[0], updir);

		
				Eigen::Vector3d eg_rot_axis(rot_axis[0], rot_axis[1], rot_axis[2]);
				eg_rot_axis.normalize();
				Eigen::Matrix3d rot_mat;
				rot_mat = Eigen::AngleAxisd(degree, eg_rot_axis);

				for (auto viter = crown_mesh->GetMesh().vertices_begin(); viter != crown_mesh->GetMesh().vertices_end(); viter++)
				{
					OpenMesh::Vec3d p = crown_mesh->GetMesh().point(viter);
					p = p - mean;
					Eigen::Vector3d ep(p[0], p[1], p[2]);
					ep = rot_mat*(ep);
					p = OpenMesh::Vec3d(ep[0], ep[1], ep[2]);
					crown_mesh->GetMesh().set_point(viter, p);
				}
				CGeoBaseAlg::NormalizeMeshSize(crown_mesh->GetMesh());
				std::stringstream sstream;
				sstream << std_dir << "\\" << count ++<<".obj";
				
				data_io.WriteMesh(sstream.str(), *crown_mesh);
			}

			break;
		}
		case Qt::Key_J:
		{
			CMeshObject *mesh_obj = DataPool::GetMeshObject(sel_temp_id_);
			CCurveObject *cc_obj = DataPool::GetCurveObject(sel_curve_id_);
			if (mesh_obj != NULL&&cc_obj!=NULL)
			{
		
				mesh_obj->ApplyTransform();
				mesh_obj->SetChanged(true);
				COpenMeshT &mesh = mesh_obj->GetMesh();
				CCorrespondenceBuilder corr_builder;
				std::map<COpenMeshT::VertexHandle, OpenMesh::Vec3d> deform_map;
				std::vector<std::pair<COpenMeshT::VertexHandle, COpenMeshT::VertexHandle>>edges;
				CGeoAlg::GetSilhouetteEdges(mesh, OpenMesh::Vec3d(0, 0, -1), edges);
				std::vector<COpenMeshT::VertexHandle>sil_vhs;
				std::set<COpenMeshT::VertexHandle>sil_vhs_set;
				for (int i = 0; i < edges.size(); i++)
				{
					sil_vhs_set.insert(edges[i].first);
					sil_vhs_set.insert(edges[i].second);
				}
				CCurveObject * edge_dots = new CCurveObject();
				for (auto iter = sil_vhs_set.begin(); iter != sil_vhs_set.end(); iter++)
				{
					sil_vhs.push_back(*iter);
					edge_dots->GetCurve().push_back(mesh.point(sil_vhs.back()));
				}
				edge_dots->RendereType() = CCurveObject::CurveType::Dots;
				edge_dots->SetColor(OpenMesh::Vec3d(1, 0, 0));
				edge_dots->SetChanged(true);
				DataPool::AddCurveObject(edge_dots);
				std::cerr << "sil vhs " << sil_vhs.size() << std::endl;
				deform_map = corr_builder.FindCorrespondenceMap(mesh_obj, cc_obj->GetCurve());
				for (auto iter = deform_map.begin(); iter != deform_map.end(); iter++)
				{
					std::cerr << "deform map" << std::endl;
					OpenMesh::Vec3d pa = mesh.point(iter->first);
					OpenMesh::Vec3d pb = iter->second;
			
					CCurveObject *tmp_c = new CCurveObject();
					tmp_c->GetCurve().resize(2);
					tmp_c->GetCurve()[0] = pa;
					tmp_c->GetCurve()[1] = pb;
					tmp_c->SetColor(OpenMesh::Vec3d(1, 0,1));
					tmp_c->SetChanged(true);
					DataPool::AddCurveObject(tmp_c);
				}
			}
		
			break;
		}
		case Qt::Key_N:
		{
			std::map<int, std::shared_ptr<CMeshObject>>& pool=DataPool::GetMeshObjectPool();
			for (auto iter = pool.begin(); iter != pool.end(); iter++)
			{
				iter->second.get()->IsVisiable() = !iter->second.get()->IsVisiable();
			}
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
				segmented_crown_ = meshobj;
				meshobj->SetChanged();
				DataPool::AddMeshObject(meshobj);

			}
			
			
			break;
		}
		case Qt::Key_Q:
		{
			std::cerr << "switch to common action" << std::endl;
			manager_->SetCurrentActionType(CActionType::Common);
			break;
		}
	}
}
void CTeethReconstructionTestAction::KeyReleaseEvent(QKeyEvent *e)
{
	switch (e->key())
	{
		case Qt::Key_A:
		{
			is_picking_ = false;
			break;
		}
	}
}

CTeethReconstructionTestAction::CTeethReconstructionTestAction()
{
	is_picking_ = false;
	type_ = CTeethReconstructionTest;
	orig_tooth_template_ = NULL;
	segmented_crown_ = NULL;
}