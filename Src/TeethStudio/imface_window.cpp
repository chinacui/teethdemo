// Copyright 2016_9 by ChenNenglun
#include "imface_window.h"
#include"ui_context.h"
#include"../DataColle/data_io.h"
#include"../AlgColle/geo_base_alg.h"
#include "..\..\Src\TeethRootRecoAlg\ndt_registration.h"

#include "QFileDialog"
#include <QKeyEvent>
#include <Eigen\Dense>

CMainWindow::CMainWindow(QWidget *parent)
	: QMainWindow(parent)
{
	ui.setupUi(this);
	/*connect(ui.button_load_model, SIGNAL(clicked()), this, SLOT(OnClickButtonLoadData()));*/
	ui.model_viewer->SetScene(CUIContext::GetScene());
	SetComponents();
	SetConnections();
	//ui.base_cut_slilder->hide();
	//ui.region_threshold_slilder->hide();
//	this->connect(ui.base_cut_slilder, SIGNAL(valueChanged(int)), this, SLOT(AdjustBaseCuttingPlane(int)));
	//this->connect(ui.region_threshold_slilder, SIGNAL(valueChanged(int)), this, SLOT(AdjustSmallRegionThreshold(int)));
	
}

void CMainWindow::pushTemplateKind()
{
	single_teeth_projection_action_->templateIdNum(this->id_1_->text().toInt(), this->id_2_->text().toInt(), this->id_3_->text().toInt(),
		this->id_4_->text().toInt(), this->id_5_->text().toInt(), this->id_6_->text().toInt(), this->id_7_->text().toInt());
}

void CMainWindow::loadXImage()
{
	QString path = QFileDialog::getOpenFileName(NULL, "load primitive image", ".",
		"Files(*.png *.jpg *.bmp )");
	if (0 == path.size()) return;
	this->panoramic_image_ = new QImage(path);
	int image_width = this->panoramic_image_->width();
	int image_height = this->panoramic_image_->height();
	vector<double> image_data(image_width * image_height, 0.0);
	for (int w = 0; w < image_width; ++w) {
		for (int h = 0; h < image_height; ++h) {
			image_data[w * image_height + h] = QColor(this->panoramic_image_->pixel(w, h)).red();
		}
	}

	this->ndt_registration_->SetPrimitiveImage(image_data, image_width, image_height);
	loadCrownData();

	this->ndt_registration_->SetBoundaryPoints(this->projection_teeth_crown_);
	this->ndt_registration_->SetCellPartitionParameters(40, 40, 10, 10, 5);
	this->ndt_registration_->SetMaxIteratorTimes(60);
	this->ndt_registration_->ExecuteNDTRegister();
	this->projection_teeth_crown_ = this->ndt_registration_->GetBoundaryPoints();
	
	//...
	this->map_3dvertice_2dpoints_ = this->teeth_projection_->gettest1();
	teeth_root_tip_[0] = Eigen::Vector2d(96, 125);
	teeth_root_tip_[1] = Eigen::Vector2d(129, 139);
	teeth_root_tip_[2] = Eigen::Vector2d(166, 148);
	teeth_root_tip_[3] = Eigen::Vector2d(199, 152);
	teeth_root_tip_[4] = Eigen::Vector2d(229, 152);
	teeth_root_tip_[5] = Eigen::Vector2d(256, 151);
	teeth_root_tip_[6] = Eigen::Vector2d(290, 155);
	teeth_root_tip_[7] = Eigen::Vector2d(324, 156);
	teeth_root_tip_[8] = Eigen::Vector2d(361, 149);
	teeth_root_tip_[9] = Eigen::Vector2d(412, 151);
	teeth_root_tip_[10] = Eigen::Vector2d(449, 138);
	teeth_root_tip_[11] = Eigen::Vector2d(470, 143);
	teeth_root_tip_[12] = Eigen::Vector2d(500, 143);
	teeth_root_tip_[13] = Eigen::Vector2d(532, 140);
	teeth_root_tip_[14] = Eigen::Vector2d(574, 137);
	teeth_root_tip_[15] = Eigen::Vector2d(597, 137);

	for (auto iter = teeth_root_tip_.begin(); iter != teeth_root_tip_.end(); iter++)
	{
		for (auto i = 0; i < this->projection_teeth_crown_.size(); i++)
		{
			if (abs(iter->second[0] - int(projection_teeth_crown_[i][0])) < 3)
			{
				std::cerr << projection_teeth_crown_[i][0] << " " << projection_teeth_crown_[i][1] << std::endl;
				this->map_toothId_map_meshId_tipLen_[iter->first][this->map_3dvertice_2dpoints_[i].begin()->first][this->map_3dvertice_2dpoints_[i].begin()->second] =
					(this->projection_crown_maxW_ - this->projectiom_crown_minW_) / (this->image_crown_maxW_ - this->image_crown_minW_) * abs(iter->second[1] - projection_teeth_crown_[i][1]);
				std::cerr << this->map_toothId_map_meshId_tipLen_[iter->first][this->map_3dvertice_2dpoints_[i].begin()->first][this->map_3dvertice_2dpoints_[i].begin()->second] << std::endl;
				break;
			}
		}
	}
    
	//..
    this->teeth_projection_->calculateTeethTop(this->map_toothId_map_meshId_tipLen_);
	//panoramic_image_registration_
	std::cerr << "good" << std::endl;
}

int CMainWindow::test()
{
	return true;
}

void CMainWindow::SetComponents()
{
	this->id_1_ = this->ui.id_1;
	this->id_2_ = this->ui.id_2;
	this->id_3_ = this->ui.id_3;
	this->id_4_ = this->ui.id_4;
	this->id_5_ = this->ui.id_5;
	this->id_6_ = this->ui.id_6;
	this->id_7_ = this->ui.id_7;

	
	this->push_template_kind = this->ui.push_template_kind;
	this->push_load_image_ = this->ui.LoadImageW;
	this->single_teeth_projection_action_ = CSingleTeethProjectionAction::GetInstance();
	this->ndt_registration_ = NDTRegistration::GetInstance();
	this->teeth_projection_ = teethProjection::GetInstance();
}

void CMainWindow::SetConnections()
{
	connect(this->push_template_kind, SIGNAL(clicked(void)),
		this, SLOT(pushTemplateKind(void)));
	connect(this->push_load_image_, SIGNAL(clicked(void)),
		this, SLOT(loadXImage(void)));
}

void CMainWindow::loadCrownData()
{
	std::vector<Eigen::Vector2d> boundary_points_temp;
	this->image_crown_center_[0] = 0;
	this->image_crown_center_[1] = 0;
	this->projection_crown_center_[0] = 0;
	this->projection_crown_center_[1] = 0;
	this->projectiom_crown_minW_ = 1e10;
	this->projection_crown_maxW_ = -1e10;
	this->image_crown_minW_ = 1e10;
	this->image_crown_maxW_ = -1e10;
	//for (int i = 0; i < 15; i++)
	{
		//QString path = "C:/x-ray_crown/crown" + QString::number(i, 10) + ".txt";
		QString path = "C:/Projection_edgeTest.txt";
		QFile fpcrown(path);
		if (fpcrown.open(fpcrown.ReadOnly))
		{
			int count = 0;
			Eigen::Vector2d point;
			while (!fpcrown.atEnd())
			{
				if (count % 2 == 0)
				{
					QString lineString = QString(fpcrown.readLine()).trimmed();
					point[0] = lineString.toDouble();
					if (point[0] < this->image_crown_minW_) this->image_crown_minW_ = point[0];
					if (point[0] > this->image_crown_maxW_) this->image_crown_maxW_ = point[0];
				}
				if (count % 2 == 1)
				{
					QString lineString = QString(fpcrown.readLine()).trimmed();
					point[1] = lineString.toDouble();
					boundary_points_temp.push_back(point);
					image_crown_center_[0] = image_crown_center_[0] + point[0];
					image_crown_center_[1] = image_crown_center_[1] + point[1];
				}
				count++;
			}

		}
		fpcrown.close();
	}
	image_crown_center_[0] = image_crown_center_[0] / boundary_points_temp.size();
	image_crown_center_[1] = image_crown_center_[1] / boundary_points_temp.size();

	this->projection_tooth_crown_ = this->teeth_projection_->getProjectionImage();
	for (auto iter = projection_tooth_crown_.begin(); iter != projection_tooth_crown_.end(); iter++)
	{
		Eigen::Vector2d temp;
		for (auto i = 0; i < iter->second.size(); i++)
		{
			temp[0] = iter->second[i][0];
			temp[1] = -1 * iter->second[i][1];
			if (temp[0] < this->projectiom_crown_minW_) this->projectiom_crown_minW_ = temp[0];
			if (temp[0] > this->projection_crown_maxW_) this->projection_crown_maxW_ = temp[0];
			this->projection_crown_center_[0] = this->projection_crown_center_[0] + temp[0];
			this->projection_crown_center_[1] = this->projection_crown_center_[1] + temp[1];
			this->projection_teeth_crown_.push_back(temp);
		}
	}
	//compensation
	this->image_crown_maxW_ = this->image_crown_maxW_ + 6;
	this->projection_crown_center_ = this->projection_crown_center_ / this->projection_teeth_crown_.size();


	for (auto i = 0; i < projection_teeth_crown_.size(); i++)
	{
		this->projection_teeth_crown_[i][0] = (this->image_crown_center_[0] - this->projection_crown_center_[0]) + (this->image_crown_maxW_ - this->image_crown_minW_) / (this->projection_crown_maxW_ - this->projectiom_crown_minW_)*(this->projection_teeth_crown_[i][0] - this->projection_crown_center_[0]) + this->projection_crown_center_[0];
		this->projection_teeth_crown_[i][1] = (this->image_crown_center_[1] - this->projection_crown_center_[1]) + (this->image_crown_maxW_ - this->image_crown_minW_) / (this->projection_crown_maxW_ - this->projectiom_crown_minW_)*(this->projection_teeth_crown_[i][1] - this->projection_crown_center_[1]) + this->projection_crown_center_[1];
	}
}

CModelViewer* CMainWindow::GetModelViewer() { 
	return ui.model_viewer; 
}

void CMainWindow::UpdateRequest()
{
	this->ui.model_viewer->update();
	this->update();	
}
CMainWindow::~CMainWindow()
{
	this->single_teeth_projection_action_->DeleteInstance();
}
