#ifndef CIMAGE_BASE_ALG_H
#define CIMAGE_BASE_ALG_H
#include"prereq.h"
#include"../DataColle/custom_openmesh_type.h"
#include<opencv2/opencv.hpp>
class ALGCOLLE_CLASS CImageBaseAlg
{
public:
	static void Curve2dToGrayImage(std::vector<OpenMesh::Vec2d>&curve,int dim_lens0,int dim_len01,unsigned char background_color, unsigned char foreground_color, cv::Mat &res_gray_image);//dim_lens0 :width,dim_lens1 :height
	static void MorphSkeleton(cv::Mat& img);//morph skeleton of the uchar binary image with black background 
	static void ShortestDis(cv::Mat &img, cv::Point src, std::vector<cv::Point> dst, std::vector<double>&dis);//img must CV_8UC1, and color of dst,src as well as path must be the same
	static void GetExtremePoints(cv::Mat&img, uchar fore_color, std::vector<cv::Point>&res_points);
};
#endif