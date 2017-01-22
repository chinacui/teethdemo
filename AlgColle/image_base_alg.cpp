#include"image_base_alg.h"
#include"curve_base_alg.h"
#include<queue>
void CImageBaseAlg::Curve2dToGrayImage(std::vector<OpenMesh::Vec2d>&curve, int dim_len0, int dim_len1, unsigned char background_color, unsigned char foreground_color, cv::Mat &res_gray_image)
{
	std::cerr << "curve size" << curve.size() << std::endl;
	std::cerr << "Curve2dToGrayImage" << std::endl;
	res_gray_image=cv::Mat(dim_len1, dim_len0,CV_8UC1);
	res_gray_image.setTo(background_color);
	OpenMesh::Vec2d bbox_min, bbox_max;
	CCurveBaseAlg::ComputeBoundingBoxOf2dCurve(curve, bbox_min, bbox_max);
	double len[2];
	len[0] = bbox_max[0] - bbox_min[0];
	len[1] = bbox_max[1] - bbox_min[1];
	for (int i = 0; i < curve.size(); i++)
	{
		int j = (i + 1) % curve.size();
		cv::Point2d a((curve[i][0] - bbox_min[0]) / len[0]* dim_len0, (curve[i][1] - bbox_min[1]) / len[1]* dim_len1);
		cv::Point2d b((curve[j][0] - bbox_min[0]) / len[0]* dim_len0, (curve[j][1] - bbox_min[1]) / len[1]* dim_len1);
		//std::cerr << a << " " << b << std::endl;
		cv::line(res_gray_image, a, b,cv::Scalar(foreground_color),3);
	}
	
}
void CImageBaseAlg::GetExtremePoints(cv::Mat&img, uchar fore_color, std::vector<cv::Point>&res_points)
{
	int dir[8][2] = { 1,0,1,1,0,1,-1,1,-1,0,-1,-1,0,-1,1,-1 };
	res_points.clear();
	for (int i = 0; i < img.rows; i++)
	{
		for (int j = 0; j < img.cols; j++)
		{
			if (img.at<uchar>(i, j) == fore_color)
			{
				int count = 0;
				for (int t = 0; t < 8; t++)
				{
					cv::Point nextp(j + dir[t][1], i + dir[t][0]);
					if (img.at<uchar>(nextp) == fore_color)
						count++;
				}
				if (count == 1)
					res_points.push_back(cv::Point(j, i));
			}
			
		}
	}
}
void CImageBaseAlg::ShortestDis(cv::Mat &img, cv::Point src, std::vector<cv::Point> dst, std::vector<double>&dis)
{
	std::queue<cv::Point>Q;
	Q.push(src);
	cv::Mat tmp_dis(img.size(), CV_64FC1);
	tmp_dis.setTo(-1);
	tmp_dis.at<double>(src) = 0;
	cv::Mat mmark(img.size(), CV_8UC1);
	mmark.setTo(0);
	mmark.at<uchar>(src) = 1;
	int dir[8][2] = { 1,0,1,1,0,1,-1,1,-1,0,-1,-1,0,-1,1,-1 };
	double sqrt2 = std::sqrt(2);
	double step_dis[8] = { 1,sqrt2 ,1,sqrt2 ,1,sqrt2 ,1,sqrt2 };
	while (!Q.empty())
	{
		auto p = Q.front();
		mmark.at<uchar>(p) = 0;
		Q.pop();
		for (int i = 0; i < 8; i++)
		{
			cv::Point nextp = cv::Point(p.x + dir[i][0], p.y + dir[i][1]);
			if (img.at<uchar>(nextp) != img.at<uchar>(src))
				continue;
			if (tmp_dis.at<double>(nextp) == -1 || tmp_dis.at<double>(p) + step_dis[i] < tmp_dis.at<double>(nextp))
			{
				tmp_dis.at<double>(nextp) = tmp_dis.at<double>(p) + step_dis[i];
				if (mmark.at<uchar>(nextp) == 0)
				{
					Q.push(nextp);
					mmark.at<uchar>(nextp) = 1;
				}
			
			}
		}
	}
	dis.resize(dst.size());
	for (int i = 0; i < dst.size(); i++)
	{
		dis[i] = tmp_dis.at<double>(dst[i]);
	}

}
void CImageBaseAlg::MorphSkeleton(cv::Mat& img)
{
	cv::Mat dst;
	int iterations = 150;
	const int height = img.rows - 1;
	const int width = img.cols - 1;

	//拷贝一个数组给另一个数组
	if (img.data != dst.data)
	{
		img.copyTo(dst);
	}


	int n = 0, i = 0, j = 0;
	cv::Mat tmpImg;
	uchar *pU, *pC, *pD;
	bool isFinished = false;

	for (n = 0; n < iterations; n++)
	{
		dst.copyTo(tmpImg);
		isFinished = false;   //一次 先行后列扫描 开始
							  //扫描过程一 开始
		for (i = 1; i < height; i++)
		{
			pU = tmpImg.ptr<uchar>(i - 1);
			pC = tmpImg.ptr<uchar>(i);
			pD = tmpImg.ptr<uchar>(i + 1);
			for (int j = 1; j < width; j++)
			{
				if (pC[j] > 0)
				{
					int ap = 0;
					int p2 = (pU[j] > 0);
					int p3 = (pU[j + 1] > 0);
					if (p2 == 0 && p3 == 1)
					{
						ap++;
					}
					int p4 = (pC[j + 1] > 0);
					if (p3 == 0 && p4 == 1)
					{
						ap++;
					}
					int p5 = (pD[j + 1] > 0);
					if (p4 == 0 && p5 == 1)
					{
						ap++;
					}
					int p6 = (pD[j] > 0);
					if (p5 == 0 && p6 == 1)
					{
						ap++;
					}
					int p7 = (pD[j - 1] > 0);
					if (p6 == 0 && p7 == 1)
					{
						ap++;
					}
					int p8 = (pC[j - 1] > 0);
					if (p7 == 0 && p8 == 1)
					{
						ap++;
					}
					int p9 = (pU[j - 1] > 0);
					if (p8 == 0 && p9 == 1)
					{
						ap++;
					}
					if (p9 == 0 && p2 == 1)
					{
						ap++;
					}
					if ((p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9) > 1 && (p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9) < 7)
					{
						if (ap == 1)
						{
							if ((p2*p4*p6 == 0) && (p4*p6*p8 == 0))
							{
								dst.ptr<uchar>(i)[j] = 0;
								isFinished = true;
							}

							//   if((p2*p4*p8==0)&&(p2*p6*p8==0))
							//    {                           
							//         dst.ptr<uchar>(i)[j]=0;
							//         isFinished =TRUE;                            
							//    }

						}
					}
				}

			} //扫描过程一 结束


			dst.copyTo(tmpImg);
			//扫描过程二 开始
			for (i = 1; i < height; i++)  //一次 先行后列扫描 开始
			{
				pU = tmpImg.ptr<uchar>(i - 1);
				pC = tmpImg.ptr<uchar>(i);
				pD = tmpImg.ptr<uchar>(i + 1);
				for (int j = 1; j < width; j++)
				{
					if (pC[j] > 0)
					{
						int ap = 0;
						int p2 = (pU[j] > 0);
						int p3 = (pU[j + 1] > 0);
						if (p2 == 0 && p3 == 1)
						{
							ap++;
						}
						int p4 = (pC[j + 1] > 0);
						if (p3 == 0 && p4 == 1)
						{
							ap++;
						}
						int p5 = (pD[j + 1] > 0);
						if (p4 == 0 && p5 == 1)
						{
							ap++;
						}
						int p6 = (pD[j] > 0);
						if (p5 == 0 && p6 == 1)
						{
							ap++;
						}
						int p7 = (pD[j - 1] > 0);
						if (p6 == 0 && p7 == 1)
						{
							ap++;
						}
						int p8 = (pC[j - 1] > 0);
						if (p7 == 0 && p8 == 1)
						{
							ap++;
						}
						int p9 = (pU[j - 1] > 0);
						if (p8 == 0 && p9 == 1)
						{
							ap++;
						}
						if (p9 == 0 && p2 == 1)
						{
							ap++;
						}
						if ((p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9) > 1 && (p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9) < 7)
						{
							if (ap == 1)
							{
								//   if((p2*p4*p6==0)&&(p4*p6*p8==0))
								//   {                           
								//         dst.ptr<uchar>(i)[j]=0;
								//         isFinished =TRUE;                            
								//    }

								if ((p2*p4*p8 == 0) && (p2*p6*p8 == 0))
								{
									dst.ptr<uchar>(i)[j] = 0;
									isFinished = true;
								}

							}
						}
					}

				}

			} //一次 先行后列扫描完成          
			  //如果在扫描过程中没有删除点，则提前退出
			if (isFinished == false)
			{
				break;
			}
		}

	}
	dst.copyTo(img);
}