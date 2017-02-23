#include"teeth_reconstruction_action.h"
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

#include"../DataColle/aux_geo_utils.h"
#include"../AlgColle/correspondence_builder.h"
#include"../AlgColle/compute_voronoi_diagram.h"
#include"../AlgColle/numerical_base_alg.h"

#include"../TeethRootRecoAlg/dental_template_fitting.h"
void CTeethReconstructionAction::Init()
{

}
void CTeethReconstructionAction::MousePressEvent(QMouseEvent *e)
{
	if (e->modifiers() == Qt::Modifier::CTRL|| e->modifiers() == Qt::Modifier::ALT)
	{
		e->setAccepted(false);
	}
	if (is_selecting_teeth_)
	{
		std::cerr << "picking crown" << std::endl;
		auto camera = viewer_->GetCamera();
		OpenMesh::Vec3d orig, dir;

		camera.ConvertClickToLine(e->pos(), orig, dir);
		
		for (int i = 0; i < crown_ids_.size(); i++)
		{
			CMeshObject *mesh_obj = DataPool::GetMeshObject(crown_ids_[i]);
			if (mesh_obj == NULL)
			{
				continue;
			}
			COpenMeshT::VertexHandle res_vh;
			if (CGeoAlg::RayMeshIntersection(orig, dir, *mesh_obj, res_vh))
			{
				if (sel_teeth_ids_.find(crown_ids_[i]) == sel_teeth_ids_.end())
				{
					sel_teeth_ids_.insert(crown_ids_[i]);
					mesh_obj->IsShinning() = true;
				}
				else
				{
					sel_teeth_ids_.erase(crown_ids_[i]);
					mesh_obj->IsShinning() = false;
				}
				
			
			}
		}
	}
	else
	{
		auto camera = viewer_->GetCamera();
		OpenMesh::Vec3d orig, dir;

		camera.ConvertClickToLine(e->pos(), orig, dir);
		int id=CGeoAlg::PickMesh(orig, dir, DataPool::GetMeshObjectPool(), true);
		if (id != -1)
		{
			if (crown2temptype_map_.find(id) != crown2temptype_map_.end())
			{
				std::cerr <<"pick id "<<id<< " type " << crown2temptype_map_[id] << std::endl;
			}
		}
	}
}
void CTeethReconstructionAction::MouseMoveEvent(QMouseEvent *e)
{
	if (e->modifiers() == Qt::Modifier::CTRL || e->modifiers() == Qt::Modifier::ALT)
	{
		e->setAccepted(false);
	}
}
void CTeethReconstructionAction::MouseReleaseEvent(QMouseEvent *e)
{
	if (e->modifiers() == Qt::Modifier::CTRL || e->modifiers() == Qt::Modifier::ALT)
	{
		e->setAccepted(false);
	}
}
void CTeethReconstructionAction::KeyPressEvent(QKeyEvent *e)
{
	if (e->modifiers() == Qt::Modifier::CTRL || e->modifiers() == Qt::Modifier::ALT)
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
		case Qt::Key_T:
		{
			QString path = QFileDialog::getExistingDirectory(NULL, "load dental mesh", ".");

			if (path.length() == 0)
			{
				std::cerr << "unable to load template meshes\n" << std::endl;
				break;
			}
			
			CDentalTemplateFitting::LoadTemplates(path.toStdString(), template_meshes_);
			
			//for (auto iter = template_meshes_.begin(); iter != template_meshes_.end(); iter++)
			//{
			//	CMeshObject *tmp_obj = new CMeshObject();
			//	iter->second->GetCrownMesh(tmp_obj->GetMesh(), std::map<COpenMeshT::VertexHandle, COpenMeshT::VertexHandle>());
			//	tmp_obj->SetMeshColor(OpenMesh::Vec3d(0.8, 0.8, 0.8));
			////	DataPool::AddMeshObject(tmp_obj);
			//}
			break;
		}
		case Qt::Key_O:
		{

			QString path = QFileDialog::getSaveFileName((QWidget*)CUIContext::GetMainWindow(), "save sel models");
		
			if (path.size() != 0)
			{
				std::string stdpath = path.toStdString();
				QDir qdir(QString::fromLocal8Bit(stdpath.data()));
				if (!qdir.exists())
				{



					if (!qdir.mkpath(QString::fromLocal8Bit(stdpath.data())))
					{
						std::cerr << "failed to create " << stdpath.data() << std::endl;
					}
				}
				
				
				for (int i = 0; i < crown_ids_.size(); i++)
				{
					std::stringstream crownss;
					crownss << stdpath << "\\crown_" << i << ".obj";
					std::stringstream tempss;
					tempss << stdpath << "\\temp_" << i << ".obj";
					CMeshObject *crown_obj = DataPool::GetMeshObject(crown_ids_[i]);
					CTeethTemplateObject  *temp_obj=crown_temp_map_[crown_ids_[i]];
				
						CDataIO dataio;
						dataio.WriteMesh(crownss.str(), *crown_obj);
						if (temp_obj != NULL)
						{
							COpenMeshT &mesh = temp_obj->GetMesh();
							auto vhs_map=temp_obj->GetCrown2TempVhsMap();
							for (auto iter = vhs_map.begin(); iter != vhs_map.end(); iter++)
							{
								mesh.set_color(iter->second, OpenMesh::Vec3d(0, 0, 0));
							
							}
							temp_obj->SetAttrChanged();
						dataio.WriteMesh(tempss.str(), *temp_obj);
						temp_obj->SetMeshColor(OpenMesh::Vec3d(0.8, 0.8, 0.8));
						}
					
				}
			}
			
			break;
		}
		case Qt::Key_S:
		{
			if (e->isAutoRepeat())
				break;
			std::cerr << "start select crown" << std::endl;
			if (is_selecting_teeth_ == false)
			{
				auto pool = DataPool::GetMeshObjectPool();
				for (auto iter = pool.begin(); iter != pool.end(); iter++)
				{
					iter->second->IsShinning() = false;
				}
				is_selecting_teeth_ = true;
				sel_teeth_ids_.clear();
					
			}
			
			break;
		}
		case Qt::Key_B:
		{
			std::vector<CMeshObject*>crowns;
			if (sel_teeth_ids_.size() != 0)
			{
				for (auto iter=sel_teeth_ids_.begin();iter!=sel_teeth_ids_.end();iter++)
				{
					CMeshObject * mesh_obj = DataPool::GetMeshObject(*iter);
					crowns.push_back(mesh_obj);
				}
			}
			else
			{
				for (int i = 0; i < crown_ids_.size(); i++)
				{
					CMeshObject * mesh_obj = DataPool::GetMeshObject(crown_ids_[i]);
					crowns.push_back(mesh_obj);
				}
			}
			

			for (int i = 0; i < crowns.size(); i++)
			{
				CMeshObject *sel_crown =crowns[i];
				int crown_id = sel_crown->GetId();
				CTeethTemplateObject *sel_temp = crown_temp_map_[crown_id];
				
				
				bool &is_finished = is_temp_matching_finished_[crown_id];
				if (!is_finished)
				{
					is_finished = non_rigid_icp_[crown_id]->Run(1);
					sel_crown->SetChanged();


					
					
					std::map<COpenMeshT::VertexHandle, COpenMeshT::VertexHandle> &crown_to_template_vmap_ =sel_temp->GetCrown2TempVhsMap();

					std::map<COpenMeshT::VertexHandle, OpenMesh::Vec3d>deform_map;
					COpenMeshT& temp_crown_mesh = sel_temp->GetCrownMesh();
					std::cerr << "crown_to_template_vmap_ num " << crown_to_template_vmap_.size() << std::endl;
					for (auto iter = crown_to_template_vmap_.begin(); iter != crown_to_template_vmap_.end(); iter++)
					{
						deform_map[iter->second] = temp_crown_mesh.point(iter->first);
						//std::cerr << "deform_map  "<< sel_temp->GetMesh().point(iter->second)<<" to " << temp_crown_mesh.point(iter->first) << std::endl;
					}

					crown_temp_arap_[crown_id]->SetDeformMap(deform_map);
					crown_temp_arap_[crown_id]->Deform();

					//std::cerr << "start change" << std::endl;


					if (is_finished)
					{
						//CGeoAlg::Remeshing(sel_temp->GetMesh());
						////CDentalTemplateFitting::RefineFittingTemplate(sel_temp->GetMesh(),sel_crown->GetMesh());
						////sel_crown->SetMeshColor(OpenMesh::Vec3d(0.8, 0.8, 0.8));
						std::map<COpenMeshT::VertexHandle, COpenMeshT::FaceHandle>vh_fh_map;
						non_rigid_icp_[crown_id]->GetTgt2SrcMap(vh_fh_map);
						CDentalTemplateFitting::RefineCrownBoundaryAfterFitting(*sel_crown, *sel_temp, vh_fh_map);
						/*CGeoAlg::SelfIntersectionRemoval(sel_temp->GetMesh());
						CGeoAlg::FillHoles(sel_temp->GetMesh(), false);*/
						CDentalTemplateFitting::ReplaceTemplateCrown(*sel_temp, *sel_crown);

						std::cerr << "finish" << std::endl;

					}

					sel_temp->SetChanged();


					//std::cerr << "end change" << std::endl;
				}
				
			}
			break;
		}
		case Qt::Key_I:
		{
			std::vector<CMeshObject*>crowns;
			for (int i = 0; i < crown_ids_.size(); i++)
			{
				CMeshObject * mesh_obj = DataPool::GetMeshObject(crown_ids_[i]);
				crowns.push_back(mesh_obj);
			}
			for (int i = 0; i < crowns.size(); i++)
			{
				CMeshObject *crown_obj = crowns[i];
				
				CTeethTemplateObject *crown_temp = crown_temp_map_[crown_obj->GetId()];

				non_rigid_icp_[crown_obj->GetId()] = new CNonRigidICP(&crown_temp->GetCrownMesh(), &(crown_obj->GetMesh()));
				crown_temp_arap_[crown_obj->GetId()] = new CARAPDeformation(crown_temp->GetMesh());
				crown_temp->SetChanged();
				is_temp_matching_finished_[crown_obj->GetId()] = false;
			}
			break;
		}
		case Qt::Key_P:
		{
			if (sel_teeth_ids_.size() == 2)
			{
				std::vector<CMeshObject*>crowns;
				for (int i = 0; i < crown_ids_.size(); i++)
				{
					CMeshObject * mesh_obj = DataPool::GetMeshObject(crown_ids_[i]);
					crowns.push_back(mesh_obj);
				}
				auto iter = sel_teeth_ids_.begin();
				CMeshObject *sel_mesh0 = DataPool::GetMeshObject(*iter);
				OpenMesh::Vec3d center0 = CGeoBaseAlg::ComputeMeshCenter(sel_mesh0->GetMesh());
				iter++;
				CMeshObject *sel_mesh1 = DataPool::GetMeshObject(*iter);
				OpenMesh::Vec3d center1 = CGeoBaseAlg::ComputeMeshCenter(sel_mesh1->GetMesh());
				std::map<CMeshObject*, int>fix_id;
				if (center0[0] < center1[0])
				{
					fix_id[sel_mesh0] = 1;
					fix_id[sel_mesh1] = 2;
				}
				else
				{
					fix_id[sel_mesh0] = 2;
					fix_id[sel_mesh1] = 1;
				}
				
				std::map<CMeshObject *, std::string>crown_type;
				CDentalTemplateFitting::ComputeCrownsType(crowns, fix_id, crown_type);
				crown2temptype_map_.clear();
				for (auto iter = crown_type.begin(); iter != crown_type.end(); iter++)
				{
					crown2temptype_map_[iter->first->GetId()] = iter->second;
				}
				for (auto iter=sel_teeth_ids_.begin();iter!=sel_teeth_ids_.end();iter++)
				{
					DataPool::GetMeshObject(*iter)->IsShinning() = false;
				}
				sel_teeth_ids_.clear();
			
				std::vector<OpenMesh::Vec3d>crown_front_dirs;
				CDentalTemplateFitting::ComputeCrownFrontDirs(crowns, crown_front_dirs);
				crown_temp_arap_.clear();
				non_rigid_icp_.clear();
				is_temp_matching_finished_.clear();
				for (int i = 0; i < crowns.size(); i++)
				{
					CMeshObject *crown_obj = crowns[i];
					crown_obj->ApplyTransform();
					crown_obj->SetChanged();
					std::cerr << "crown id " << crown_ids_[i] << std::endl;
					
					std::string temp_type = crown2temptype_map_[crown_obj->GetId()];
					CTeethTemplateObject *crown_temp = new CTeethTemplateObject(*template_meshes_[temp_type]);
					
					crown_temp_map_[crown_obj->GetId()] = crown_temp;
					CDentalTemplateFitting::AlignTemplate2Crown(*crown_obj, crown_front_dirs[i],OpenMesh::Vec3d(0,1,0), *crown_temp);
					crown_temp->ComputeCrownMesh();
					crown_temp->SetId(-1);
					//std::cerr << "crown temp size " << crown_temp->GetMesh().n_vertices() << std::endl;
					DataPool::AddMeshObject(crown_temp);
				//	crown_temp->IsVisiable()=false;

					/*CMeshObject *tmp_crown = new CMeshObject();
					tmp_crown->SetMesh(crown_temp->GetCrownMesh());
					tmp_crown->SetMeshColor(OpenMesh::Vec3d(0.8, 0.8, 0.8));
					tmp_crown->SetChanged();
					DataPool::AddMeshObject(tmp_crown);*/


					non_rigid_icp_[crown_obj->GetId()] = new CNonRigidICP(&crown_temp->GetCrownMesh(), &(crown_obj->GetMesh()));
					crown_temp_arap_[crown_obj->GetId()] = new CARAPDeformation(crown_temp->GetMesh());
					crown_temp->SetChanged();
					is_temp_matching_finished_[crown_obj->GetId()] = false;
					
			
				}
				
			}
			
			
			break;
		}
		case Qt::Key_L:
		{
			//	std::cerr << "l" << std::endl;
			QString path = QFileDialog::getOpenFileName(NULL, "load dental mesh", ".", "Stl Files(*.obj)");

			if (path.length() == 0)
			{
				std::cerr << "unable to load mesh\n" << std::endl;
				return;
			}
			CMeshObject *meshobj = new CMeshObject();
			COpenMeshT&mesh = meshobj->GetMesh();
			if (!CDataIO::ReadMesh(path.toStdString(), *meshobj))
			{
				std::cerr << "unable to load mesh\n" << std::endl;
			}
			

			meshobj->RestoreCurrentVPos();
			CDentalBaseAlg::PCABasedOrientationCorrection(meshobj->GetMesh());
			path[path.length() - 3] = 'd';
			path[path.length() - 2] = 'm';
			path[path.length() - 1] = 'a';
			path = path + 't';
			Eigen::VectorXd tags;
			igl::readDMAT(path.toStdString(), tags);
			std::vector<int>vtags;
			for (int i = 0; i < tags.rows(); i++)
			{
				vtags.push_back(tags(i));
			}
			std::map<int, COpenMeshT*>sep_meshes;
			CGeoAlg::SeparateMeshByVertexTag(mesh, vtags, sep_meshes, std::map<int, std::map<COpenMeshT::VertexHandle, COpenMeshT::VertexHandle>>());
		
			crown_ids_.clear();
			for (auto iter = sep_meshes.begin(); iter != sep_meshes.end(); iter++)
			{
				if (iter->first == -1)
					continue;
				if (iter->second->n_vertices() < 20)
					continue;
				CMeshObject* crown_obj = new CMeshObject();
				crown_obj->SetMesh(*iter->second);
				crown_obj->SetChanged();
				int id=DataPool::AddMeshObject(crown_obj);
				crown_ids_.push_back(id);
				delete iter->second;
				
			}
			std::cerr << "crown size " << crown_ids_.size() << std::endl;
			sep_meshes.clear();
			
			break;
			
		}
	}
}
void CTeethReconstructionAction::KeyReleaseEvent(QKeyEvent *e)
{
	if (e->modifiers() == Qt::Modifier::CTRL || e->modifiers() == Qt::Modifier::ALT)
	{
		e->setAccepted(false);
	}
	switch (e->key())
	{
	
		case Qt::Key_S:
		{
			if (e->isAutoRepeat())
				break;
			if (is_selecting_teeth_)
			{
				is_selecting_teeth_ = false;
			}
		
			break;
		}
	}
}
void CTeethReconstructionAction::QuitAction()
{
	std::cerr << "switch to common action" << std::endl;
	manager_->SetCurrentActionType(CActionType::Common);
	
}
CTeethReconstructionAction::CTeethReconstructionAction()
{
	type_ = CTeethReconstruction;
}