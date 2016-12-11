#include "action_manager.h"
#include<iostream>
#include"cmodelviewer.h"
#include"hot_key_action.h"
#include"edit_feature_edge_action.h"
#include"harmonic_field_segmentation_action.h"
#include"volume_data_segmentation_action.h"
#include"teeth_reconstruction_action.h"
void CActionManager::Init(CModelViewer* viewer)
{
	viewer_ = viewer;
	for (int i = 0; i < actions_.size(); i++)
	{
		if (actions_[i] != NULL)
			delete actions_[i];
		
	}
	actions_.clear();
	actions_.push_back(new CHotKeyAction());
	actions_.push_back(new CEditFeatureEdgeAction());
	actions_.push_back(new CHarmonicFieldSegmentation());
	actions_.push_back(new CVolumeDataSegmentationAction());
	actions_.push_back(new CTeethReconstruction());
	for (int i = 0; i < actions_.size(); i++)
	{
		actions_[i]->SetViewer(viewer_);
		actions_[i]->SetManager(this);
	}
	current_action_type_ = Common;
	
}
void CActionManager::MousePressEvent(QMouseEvent *e)
{
	for (int i = 0; i < actions_.size(); i++)
	{
		if (actions_[i]->GetType() == current_action_type_)
		{
			actions_[i]->MousePressEvent(e);
		}
	}
	if (e->isAccepted() == false)
	{
		for (int i = 0; i < actions_.size(); i++)
		{
			if (actions_[i]->GetType() == Common)
			{
				actions_[i]->MousePressEvent(e);
			}
		}
	}
	
}
void CActionManager::MouseMoveEvent(QMouseEvent *e)
{
	for (int i = 0; i < actions_.size(); i++)
	{
		if (actions_[i]->GetType() == current_action_type_)
		{
			actions_[i]->MouseMoveEvent(e);
		}
	}
	if (e->isAccepted() == false)
	{
		for (int i = 0; i < actions_.size(); i++)
		{
			if (actions_[i]->GetType() == Common)
			{
				actions_[i]->MouseMoveEvent(e);
			}
		}
	}
}
void CActionManager::MouseReleaseEvent(QMouseEvent *e)
{
	for (int i = 0; i < actions_.size(); i++)
	{
		if (actions_[i]->GetType() == current_action_type_)
		{
			actions_[i]->MouseReleaseEvent(e);
		}
	}
	if (e->isAccepted() == false)
	{
		for (int i = 0; i < actions_.size(); i++)
		{
			if (actions_[i]->GetType() == Common)
			{
				actions_[i]->MouseReleaseEvent(e);
			}
		}
	}
}
void CActionManager::KeyPressEvent(QKeyEvent *e)
{
	for (int i = 0; i < actions_.size(); i++)
	{
		if (actions_[i]->GetType() == current_action_type_)
		{
			actions_[i]->KeyPressEvent(e);
		}
	}
	if (e->isAccepted() == false)
	{
		for (int i = 0; i < actions_.size(); i++)
		{
			if (actions_[i]->GetType() == Common)
			{
				actions_[i]->KeyPressEvent(e);
			}
		}
	}
}
void CActionManager::SetCurrentActionType(CActionType t)
{
	current_action_type_ = t;
	for (int i = 0; i < actions_.size(); i++)
	{
		if (actions_[i]->GetType() == current_action_type_)
		{
			actions_[i]->Init();
		}
	}
}
void CActionManager::KeyReleaseEvent(QKeyEvent *e)
{
	for (int i = 0; i < actions_.size(); i++)
	{
		if (actions_[i]->GetType() == current_action_type_)
		{
			actions_[i]->KeyReleaseEvent(e);
		}
	}
	if (e->isAccepted() == false)
	{
		for (int i = 0; i < actions_.size(); i++)
		{
			if (actions_[i]->GetType() == Common)
			{
				actions_[i]->KeyReleaseEvent(e);
			}
		}
	}
}
void CActionManager::WheelEvent(QWheelEvent *e)
{
	for (int i = 0; i < actions_.size(); i++)
	{
		if (actions_[i]->GetType() == current_action_type_)
		{
			actions_[i]->WheelEvent(e);
		}
	}
	if (e->isAccepted() == false)
	{
		for (int i = 0; i < actions_.size(); i++)
		{
			if (actions_[i]->GetType() == Common)
			{
				actions_[i]->WheelEvent(e);
			}
		}
	}
}

void CActionManager::MouseDoubleClickEvent(QMouseEvent * e)
{
	for (int i = 0; i < actions_.size(); i++)
	{
		if (actions_[i]->GetType() == current_action_type_)
		{
			actions_[i]->MouseDoubleClickEvent(e);
		}
	}
	if (e->isAccepted() == false)
	{
		for (int i = 0; i < actions_.size(); i++)
		{
			if (actions_[i]->GetType() == Common)
			{
				actions_[i]->MouseDoubleClickEvent(e);
			}
		}
	}
}