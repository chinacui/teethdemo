#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include "QObject"
#include "QFileDialog"
#include "QProgressBar"
#include "QLineEdit"
#include "QLabel"
#include "QRgb"
#include <Eigen\Dense>
#include "panoramic_image_registration.h"
#include "ndt_registration.h"
#include "vcl_iostream.h"
#include "qtextstream.h"
#include<opencv2\opencv.hpp>
#include "../AlgColle/geo_alg.h"
#include "c:\Users\admin\Desktop\workspace\Src\AlgColle\image_base_alg.h"

PanoramicImageRegistration::PanoramicImageRegistration(QWidget *parent)
	: QWidget(parent)
{
	ui.setupUi(this);

	this->SetComponents();
	this->SetConnections();
}

PanoramicImageRegistration::~PanoramicImageRegistration()
{
	this->ndt_registration_->DeleteInstance();
}

void PanoramicImageRegistration::SetComponents() {
	this->button_load_primitive_image_ = this->ui.button_load_primitive_image;
	this->button_labeled_image_ = this->ui.button_labeled_image;
	this->button_boundaries_up_ = this->ui.button_boundaries_up;
	this->button_boundaries_down_ = this->ui.button_boundaries_down;
	this->button_boundaries_left_ = this->ui.button_boundaries_left;
	this->button_boundaries_right_ = this->ui.button_boundaries_right;
	this->button_registering_with_NDT_algorithm_ = this->ui.button_registering_with_NDT_algorithm;
	this->label_image_panel1_ = this->ui.label_image_panel1;
	this->label_image_panel2_ = this->ui.label_image_panel2;
	this->lineedit_cell_num_w_ = this->ui.lineedit_cell_num_w;
	this->lineedit_cell_num_h_ = this->ui.lineedit_cell_num_h;
	this->lineedit_cell_num_addition_w_ = this->ui.lineedit_cell_num_addition_w;
	this->lineedit_cell_num_addition_h_ = this->ui.lineedit_cell_num_addition_h;
	this->lineedit_cell_num_count_ = this->ui.lineedit_cell_num_count;
	this->lineedit_max_iterator_times_ = this->ui.lineedit_max_iterator_times;
	this->button_load_obj_points_ = this->ui.button_load_obj_points;
	this->ndt_registration_ = NDTRegistration::GetInstance();
}

void PanoramicImageRegistration::SetConnections() {
	connect(this->button_load_primitive_image_, SIGNAL(clicked(void)),
		this, SLOT(OnLoadPrimiticeImage(void)));
	connect(this->button_labeled_image_, SIGNAL(clicked(void)),
		this, SLOT(OnLoadLabeledImage(void)));
	connect(this->button_boundaries_up_, SIGNAL(clicked(void)),
		this, SLOT(OnTransformAndRotateBoundaries(void)));
	connect(this->button_boundaries_down_, SIGNAL(clicked(void)),
		this, SLOT(OnTransformAndRotateBoundaries(void)));
	connect(this->button_boundaries_left_, SIGNAL(clicked(void)),
		this, SLOT(OnTransformAndRotateBoundaries(void)));
	connect(this->button_boundaries_right_, SIGNAL(clicked(void)),
		this, SLOT(OnTransformAndRotateBoundaries(void)));
	connect(this->button_registering_with_NDT_algorithm_, SIGNAL(clicked(void)),
		this, SLOT(OnRegisterWithNDTAlgorithm(void)));
	connect(this->button_load_obj_points_, SIGNAL(clicked(void)),
		this, SLOT(OnLoadOBJPoints(void)));
}

void PanoramicImageRegistration::OnLoadPrimiticeImage(void) {
	QString path = QFileDialog::getOpenFileName(NULL, "load primitive image", ".", 
		"Files(*.png *.jpg *.bmp )");
	if (0 == path.size()) return;
	this->primitive_image_ = new QImage(path);
	int image_width = this->primitive_image_->width();
	int image_height = this->primitive_image_->height();
	vector<double> image_data(image_width * image_height, 0.0);
	for (int w = 0; w < image_width; ++w) {
		for (int h = 0; h < image_height; ++h) {
			image_data[w * image_height + h] = QColor(this->primitive_image_->pixel(w, h)).red();
		}
	}
	this->ndt_registration_->SetPrimitiveImage(image_data, image_width, image_height);
	this->label_image_panel1_->setPixmap(QPixmap::fromImage(*this->primitive_image_));
}

void PanoramicImageRegistration::OnLoadLabeledImage(void) {
	QString path = QFileDialog::getOpenFileName(NULL, "load labeled image", ".",
		"Files(*.png *.jpg *.bmp )");
	vector<Eigen::Vector2d> boundary_points;
	vector<Eigen::Vector2d> boundary_points_temp; // the boundary segmented from x-ray image
	vector<Eigen::Vector2d> boundary_points_root;
	Eigen::Vector2d boundary_points_center_temp; // the center of the boundary in x-ray image
	Eigen::Vector2d boundary_points_center; // the center of the boundary in projection plane
	if (0 == path.size()) return;
	this->labeled_image_ = new QImage(path);
	int image_width = this->labeled_image_->width();
	int image_height = this->labeled_image_->height();
	double minW = 1e30;
	double maxW = 1e-30;
	/*for (int w = 0; w < image_width; ++w) {
		for (int h = 0; h < image_height; ++h) {
			if (QColor(this->labeled_image_->pixel(w, h)).red() == 17 &&
				QColor(this->labeled_image_->pixel(w, h)).green() == 255 &&
				QColor(this->labeled_image_->pixel(w, h)).blue() == 238) {
				boundary_points_temp.push_back(Eigen::Vector2d(w, h));
				boundary_points_center_temp[0] = boundary_points_center_temp[0] + w;
				boundary_points_center_temp[1] = boundary_points_center_temp[1] + h;
				if (w < minW) minW = w;
				if (w > maxW) maxW = w;
			}
			if (QColor(this->labeled_image_->pixel(w, h)).red() == 255 &&
				QColor(this->labeled_image_->pixel(w, h)).green() == 0 &&
				QColor(this->labeled_image_->pixel(w, h)).blue() == 0)
			{
				boundary_points_root.push_back(Eigen::Vector2d(w, h));
			}
		}
	}*/
	//boundary_points_center_temp[0] = boundary_points_center_temp[0] / boundary_points_temp.size();
	//boundary_points_center_temp[1] = boundary_points_center_temp[1] / boundary_points_temp.size();

	QFile fp("C:/x-ray_root/root10.txt");
	if (fp.open(fp.ReadOnly))
	{
		int count = 0;
		Eigen::Vector2d point;
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
				boundary_points_root.push_back(point);
			}
			count++;
		}

	}
	fp.close();

	QFile fpcrown("C:/x-ray_crown/crown10.txt");
	if (fpcrown.open(fp.ReadOnly))
	{
		int count = 0;
		Eigen::Vector2d point;
		while (!fpcrown.atEnd())
		{
			if (count == 0)
			{
				boundary_points_center_temp[0] = 0;
				boundary_points_center_temp[1] = 0;
			}
			if (count % 2 == 0)
			{
				QString lineString = QString(fpcrown.readLine()).trimmed();
				point[0] = lineString.toDouble();
				if (point[0] < minW) minW = point[0];
				if (point[0] > maxW) maxW = point[0];
			}
			if (count % 2 == 1)
			{
				QString lineString = QString(fpcrown.readLine()).trimmed();
				point[1] = lineString.toDouble();
				boundary_points_temp.push_back(point);
				boundary_points_center_temp[0] = boundary_points_center_temp[0] + point[0];
				boundary_points_center_temp[1] = boundary_points_center_temp[1] + point[1];
			}
			count++;
		}

	}
	fpcrown.close();
	boundary_points_center_temp[0] = boundary_points_center_temp[0] / boundary_points_temp.size();
	boundary_points_center_temp[1] = boundary_points_center_temp[1] / boundary_points_temp.size();

	cv::Mat binary_root_crown = cv::Mat(image_height, image_width, CV_8UC1);
	unsigned char background_color = 0;
	unsigned char foreground_color = 255;
	binary_root_crown.setTo(background_color);
	int count = 0;
	cv::Point point_one, point_two;

	boundary_points_temp.insert(boundary_points_temp.end(), boundary_points_root.begin(), boundary_points_root.end());
	std::vector<OpenMesh::Vec2d> root_crown;
	std::vector<std::vector<int>>vids;
	for (auto i = 0; i < boundary_points_temp.size(); i++)
	{
		OpenMesh::Vec2d p;
		p[0] = boundary_points_temp[i][0];
		p[1] = boundary_points_temp[i][1];
		root_crown.push_back(p);
	}
	CGeoAlg::AlphaShape2d(root_crown, 400, vids);
	for (auto i = 0; i < vids[0].size() - 1; i++)
	{
		cv::Point2d a, b;
		a = cv::Point(root_crown[vids[0][i]][0], root_crown[vids[0][i]][1]);
		b = cv::Point(root_crown[vids[0][i + 1]][0], root_crown[vids[0][i + 1]][1]);
		std::cerr << a << " " << b << std::endl;
		cv::line(binary_root_crown, a, b, cv::Scalar(foreground_color), 0);
		if (i == vids[0].size() - 2)
		{
			a = cv::Point(root_crown[vids[0][i + 1]][0], root_crown[vids[0][i + 1]][1]);
			b = cv::Point(root_crown[vids[0][0]][0], root_crown[vids[0][0]][1]);
			cv::line(binary_root_crown, a, b, cv::Scalar(foreground_color), 0);
		}
	}
	RegionGrow(binary_root_crown, boundary_points_center_temp);
	CImageBaseAlg::MorphSkeleton(binary_root_crown);
	for (auto i = 0; i < binary_root_crown.rows; i++)
	{
		for (auto j = 0; j < binary_root_crown.cols; j++)
		{
			if (binary_root_crown.at<uchar>(i, j) == foreground_color)
			{
				OpenMesh::Vec2d p;
				p[0] = i;
				p[1] = j;
				line_skeleton.push_back(p);
				std::cerr << i << " " << j << std::endl;
				std::cerr << count++ << std::endl;
			}
		}
	}

	QFile fpproject_edge("C:/project_edge/Projection_edge11.txt");
	if (fpproject_edge.open(fpproject_edge.ReadOnly))
	{
		int count = 0;
		Eigen::Vector2d point;
		while (!fpproject_edge.atEnd())
		{
			if (count == 0)
			{
				boundary_points_center[0] = 0;
				boundary_points_center[1] = 0;
			}
			if (count%2 == 0)
			{
				QString lineString = QString(fpproject_edge.readLine()).trimmed();
				point[0] = lineString.toDouble();
			}
			if (count%2 == 1)
			{
				QString lineString = QString(fpproject_edge.readLine()).trimmed();
				point[1] = lineString.toDouble();
				boundary_points.push_back(point);
				boundary_points_center[0] = boundary_points_center[0] + point[0];
				boundary_points_center[1] = boundary_points_center[1] + point[1];
			}
			count++;
		}
 
	}
	fpproject_edge.close();
	boundary_points_center[0] = boundary_points_center[0] / boundary_points.size();
	boundary_points_center[1] = boundary_points_center[1] / boundary_points.size();
	/*set the projection boundary to be stardard*/
	double minW_P = 1e30;
	double maxW_P = -1e30;
	for (auto i = 0; i < boundary_points.size(); i++)
	{
		if (minW_P > boundary_points[i][0]) minW_P = boundary_points[i][0];
		if (maxW_P < boundary_points[i][0]) maxW_P = boundary_points[i][0];
	}
	for (auto i = 0; i < boundary_points.size(); i++)
	{
		boundary_points[i][0] = (boundary_points_center_temp[0] - boundary_points_center[0]) + (maxW - minW) / (maxW_P - minW_P) * (boundary_points[i][0] - boundary_points_center[0]) + boundary_points_center[0];
		boundary_points[i][1] = (boundary_points_center_temp[1] - boundary_points_center[1]) + (maxW - minW) / (maxW_P - minW_P) * (boundary_points[i][1] - boundary_points_center[1]) + boundary_points_center[1];
	}
	this->ndt_registration_->SetBoundaryCenterPoints(boundary_points_center_temp);
	this->ndt_registration_->SetBoundaryPoints(boundary_points);
	this->ndt_registration_->setBoundaryRootPoint(boundary_points_root);

	b_c = boundary_points_center_temp;
	b_c_P = boundary_points_center;
	scalse = (maxW_P - minW_P) / (maxW - minW);

	this->label_image_panel2_->setPixmap(QPixmap::fromImage(*this->labeled_image_));
	this->ShowBoundaryPointsInImage();
}

void PanoramicImageRegistration::OnLoadOBJPoints(void) {
	QString path = QFileDialog::getOpenFileName(NULL, "load obj points", ".",
		"Files(*.obj)");
	std::ifstream in(path.toLocal8Bit());
	double minw = 1e30, maxw = -1e30;
	double minh = 1e30, maxh = -1e30;
	double valuew, valueh;
	vector<Eigen::Vector2d> boundary_points;
	while (in >> valuew >> valueh) {
		valuew *= 3;
		valueh *= 3;
		valueh = -valueh;
		minw = std::min(minw, valuew);
		maxw = std::max(maxw, valuew);
		minh = std::min(minh, valueh);
		maxh = std::max(maxh, valueh);
		boundary_points.push_back(Eigen::Vector2d(valuew, valueh));
		std::cerr << valuew << valueh << std::endl;
	}
	int points_count = (int)boundary_points.size();
	int image_width = this->ndt_registration_->GetImageWidth();
	int image_height = this->ndt_registration_->GetImageHeight();
	QImage image(*primitive_image_);
	for (int i = 0; i < points_count; ++i) {
		boundary_points[i](0) -= minw;
		boundary_points[i](1) -= minh;
		int idx = (int)boundary_points[i](0);
		int idy = (int)boundary_points[i](1);
		if (idx >= 0 && idx < image_width && idy >= 0 && idy < image_height) {
			image.setPixel(idx, idy, QColor(255, 0, 0).rgb());
		}
	}
	this->ndt_registration_->SetBoundaryPoints(boundary_points);
	this->label_image_panel2_->setPixmap(QPixmap::fromImage(image));
}

void PanoramicImageRegistration::OnTransformAndRotateBoundaries(void) {
	QPushButton* sender = qobject_cast<QPushButton*>(this->sender());
	if (sender == this->button_boundaries_up_) {
		this->ndt_registration_->TransformAndRotateBoundaryPoints(1, 3);
	}
	else if (sender == this->button_boundaries_down_) {
		this->ndt_registration_->TransformAndRotateBoundaryPoints(2, 3);
	}
	else if (sender == this->button_boundaries_left_) {
		this->ndt_registration_->TransformAndRotateBoundaryPoints(3, 3);
	}
	else if (sender == this->button_boundaries_right_) {
		this->ndt_registration_->TransformAndRotateBoundaryPoints(4, 3);
	}
	this->ShowBoundaryPointsInImage();
}

void PanoramicImageRegistration::OnRegisterWithNDTAlgorithm(void) {
	int cell_init_num_w = this->lineedit_cell_num_w_->text().toInt();
	int cell_init_num_h = this->lineedit_cell_num_h_->text().toInt();
	int cell_num_addition_w = this->lineedit_cell_num_addition_w_->text().toInt();
	int cell_num_addition_h = this->lineedit_cell_num_addition_h_->text().toInt();
	int cell_num_count = this->lineedit_cell_num_count_->text().toInt();
	int max_iterator_times = this->lineedit_max_iterator_times_->text().toInt();
	this->ndt_registration_->SetCellPartitionParameters(cell_init_num_w,
		cell_init_num_h, cell_num_addition_w, cell_num_addition_h, cell_num_count);
	this->ndt_registration_->SetMaxIteratorTimes(max_iterator_times);
	this->ndt_registration_->ExecuteNDTRegister();
	for (int i = 0; i < this->ndt_registration_->boundary_points_image_to_model.size(); i++)
	{
		this->ndt_registration_->boundary_points_image_to_model[i](0) = (b_c_P(0) - b_c(0)) + b_c(0) + scalse * (this->ndt_registration_->boundary_points_image_to_model[i](0) - b_c(0));
		this->ndt_registration_->boundary_points_image_to_model[i](1) = (b_c_P(1) - b_c(1)) + b_c(1) + scalse * (this->ndt_registration_->boundary_points_image_to_model[i](1) - b_c(1));
	}
	FILE *fp;
	if ((fp = fopen("C:/imageToModelPoints.txt", "w+")) != NULL)  //判断文件是否已经被打开
	{
		for (int i = 0; i < this->ndt_registration_->boundary_points_image_to_model.size(); i++)
		{
			fprintf(fp, "%f", this->ndt_registration_->boundary_points_image_to_model[i](0));
			fprintf(fp, "\n");
			fprintf(fp, "%f", this->ndt_registration_->boundary_points_image_to_model[i](1));
			fprintf(fp, "\n"); 
		}
		this->ndt_registration_->boundary_points_image_to_model.clear();
	}
	fclose(fp);
	this->ShowBoundaryPointsInImage();
}

void PanoramicImageRegistration::ShowBoundaryPointsInImage() {
	QImage image(*primitive_image_);
	vector<Eigen::Vector2d>& boundary_points = this->ndt_registration_->GetBoundaryPoints();
	int image_width = this->ndt_registration_->GetImageWidth();
	int image_height = this->ndt_registration_->GetImageHeight();
	int boundary_points_count = boundary_points.size();
	/*vector<double> gradient;
	gradient = this->ndt_registration_->getImageGradient();
	if (!gradient.empty()) {
		for (int w = 0; w < image_width; ++w) {
			for (int h = 0; h < image_height; ++h) {
				std::cerr << w << " "<<h << std::endl;
				gradient[w * image_height + h] = gradient[w * image_height + h] * 10;
				if (gradient[w * image_height + h] > 255) gradient[w * image_height + h] = 255;
				image.setPixel(w, h, QColor(int(gradient[w * image_height + h]), int(gradient[w * image_height + h]), int(gradient[w * image_height + h])).rgb());
			}
		}*/
	for (int i = 0; i < line_skeleton.size(); i++)
	{
		int idx = (int)line_skeleton[i][1];
		int idy = (int)line_skeleton[i][0];
		if (idx >= 0 && idx < image_width && idy >= 0 && idy < image_height) {
			image.setPixel(idx, idy, QColor(255, 0, 0).rgb());
		}
	}
		/*for (int p = 0; p < boundary_points_count; ++p) {
			int idx = (int)boundary_points[p](0);
			int idy = (int)boundary_points[p](1);
			if (idx >= 0 && idx < image_width && idy >= 0 && idy < image_height) {
				image.setPixel(idx, idy, QColor(255, 0, 0).rgb());
			}
		}*/
		this->label_image_panel1_->setPixmap(QPixmap::fromImage(image));
	//}
}

void PanoramicImageRegistration::RegionGrow(cv::Mat &src, Eigen::Vector2d pt)
{
	Eigen::Vector2d ptGrowing;                      //待生长点位置  
	unsigned char foreground_color = 255, background_color = 0;
	int nGrowLable = 0;                             //标记是否生长过  
	int nSrcValue = 0;                              //生长起点灰度值  
	int nCurValue = 0;                              //当前生长点灰度值  
	cv::Mat matDst = cv::Mat::zeros(src.size(), CV_8UC1);   //创建一个空白区域，填充为黑色  
													//生长方向顺序数据  
	int DIR[4][2] = { { 0,-1 },{ 1,0 },{ 0,1 },{ -1,0 } };
	std::vector<Eigen::Vector2d> vcGrowPt;                     //生长点栈 
	pt[0] = int(pt[0]);
	pt[1] = int(pt[1]);
	vcGrowPt.push_back(pt);                         //将生长点压入栈中  
	src.at<uchar>(pt[1], pt[0]) = foreground_color;               //标记生长点  
	nSrcValue = src.at<uchar>(pt[1], pt[0]);           //记录生长点的灰度值  

	while (!vcGrowPt.empty())                       //生长栈不为空则生长  
	{
		pt = vcGrowPt.back();                       //取出一个生长点  
		vcGrowPt.pop_back();

		//分别对八个方向上的点进行生长  
		for (int i = 0; i<4; ++i)
		{
			ptGrowing[0] = pt[0] + DIR[i][0];
			ptGrowing[1] = pt[1] + DIR[i][1];
			if (ptGrowing[0] == 466 && ptGrowing[1] == 150)
			{
				std::cerr<<"a"<<std::endl;
			}
			//检查是否是边缘点  
			if (ptGrowing[0] < 0 || ptGrowing[1] < 0 || ptGrowing[0] >(src.cols - 1) || (ptGrowing[1] > src.rows - 1))
				continue;

			nGrowLable = src.at<uchar>(ptGrowing[1], ptGrowing[0]);      //当前待生长点的灰度值  

			if (nGrowLable == background_color)                    //如果标记点还没有被生长  
			{
				nCurValue = src.at<uchar>(ptGrowing[1], ptGrowing[0]);
				if (nSrcValue - nCurValue == foreground_color)                 //在阈值范围内则生长  
				{
					src.at<uchar>(ptGrowing[1], ptGrowing[0]) = foreground_color;     //标记为白色  
					vcGrowPt.push_back(ptGrowing);                  //将下一个生长点压入栈中  
				}
			}
		}
	}
}