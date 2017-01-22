// Copyright 2016_9 by ChenNenglun
#include"data_io.h"
#include<vector>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include<igl/writeOBJ.h>
#include"cgal_igl_converter.h"
#include "itkGDCMImageIO.h"
#include "itkGDCMSeriesFileNames.h"
#include "itkImageSeriesReader.h"
#include "itkImageFileWriter.h"
bool CDataIO::ReadMesh(std::string fname, CMeshObject & res_mesh_obj, OpenMesh::IO::Options io_options)
{
	if (io_options != OpenMesh::IO::Options::Default)
	{
		if (!OpenMesh::IO::read_mesh(res_mesh_obj.GetMesh(), fname, io_options))
		{
			std::cerr << "read error\n";
			return false;
		}
	}
	else
	{
		if (!OpenMesh::IO::read_mesh(res_mesh_obj.GetMesh(), fname))
		{
			std::cerr << "read error\n";
			return false;
		}
	}
	
	COpenMeshT &mesh = res_mesh_obj.GetMesh();
	if (io_options== OpenMesh::IO::Options::Default)
	{
		for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); viter++)
		{
			mesh.set_color(viter, OpenMesh::Vec3d(0.8, 0.8, 0.8));
		}
	}
	
	res_mesh_obj.SetChanged();
	return true;
}
bool CDataIO::LoadCurveFromObj(std::string fname, std::vector<OpenMesh::Vec3d> &curve)
{
	std::ifstream in(fname);
	char buf[256];
	curve.clear();
	while (in.getline(buf, sizeof buf))
	{
		std::istringstream line(buf);
		std::string word;
		line >> word;
		if (word == "v")
		{
			double x, y, z;
			line >> x;
			line >> y;
			line >> z;
			OpenMesh::Vec3d pt(x, y, z);
			curve.push_back(pt);
		}
	}
	return true;
}
bool CDataIO::ReadVolumeDataObjFromDICOMSeries(std::string dirname, CVolumeDataObject& volume_data_obj)
{
	auto pvdata=ReadVolumeDataFromDICOMSeries(dirname);
	auto region = pvdata->GetLargestPossibleRegion();
	/*if (pvdata->)
	{*/
		volume_data_obj.SetVolumeData(pvdata);
		volume_data_obj.SetChanged();
		return true;
	//}
	////else
	//{
	//	return false;
	//}
	
}
ItkVolumeDataType::Pointer CDataIO::ReadVolumeDataFromDICOMSeries(std::string dirname)
{
	typedef itk::ImageSeriesReader< ItkVolumeDataType >        ReaderType;
	typedef itk::GDCMSeriesFileNames NamesGeneratorType;
	NamesGeneratorType::Pointer nameGenerator = NamesGeneratorType::New();

	nameGenerator->SetUseSeriesDetails(true);
	nameGenerator->AddSeriesRestriction("0008|0021");
	nameGenerator->SetGlobalWarningDisplay(false);
	nameGenerator->SetDirectory(dirname);

	try
	{
		typedef std::vector< std::string >    SeriesIdContainer;
		const SeriesIdContainer & seriesUID = nameGenerator->GetSeriesUIDs();
		SeriesIdContainer::const_iterator seriesItr = seriesUID.begin();
		SeriesIdContainer::const_iterator seriesEnd = seriesUID.end();

		if (seriesItr != seriesEnd)
		{
			std::cout << "The directory: ";
			std::cout << dirname << std::endl;
			std::cout << "Contains the following DICOM Series: ";
			std::cout << std::endl;
		}
		else
		{
			std::cout << "No DICOMs in: " << dirname << std::endl;
			return NULL;
		}

		while (seriesItr != seriesEnd)
		{
			std::cout << seriesItr->c_str() << std::endl;
			++seriesItr;
		}

		seriesItr = seriesUID.begin();
		while (seriesItr != seriesUID.end())
		{
			std::string seriesIdentifier;


			seriesIdentifier = seriesItr->c_str();
			seriesItr++;

			std::cout << "\nReading: ";
			std::cout << seriesIdentifier << std::endl;
			typedef std::vector< std::string >   FileNamesContainer;
			FileNamesContainer fileNames;
			fileNames = nameGenerator->GetFileNames(seriesIdentifier);

			 ReaderType::Pointer reader = ReaderType::New();
			typedef itk::GDCMImageIO       ImageIOType;
			ImageIOType::Pointer dicomIO = ImageIOType::New();
			reader->SetImageIO(dicomIO);
			reader->SetFileNames(fileNames);
			reader->UpdateLargestPossibleRegion();
			reader->Update();
			/*reader->SetReleaseDataFlag(false);
			reader->ReleaseDataBeforeUpdateFlagOff();
			reader->ReleaseDataFlagOff();
			*/
			ItkVolumeDataType::Pointer volume_data =reader->GetOutput();
			volume_data->Update();
			//for (int i = 0; i < 100; i++)
			//{
			//	for (int j = 0; j < 100; j++)
			//	{
			//		for (int k = 0; k < 100; k++)
			//		{
			//			ItkVolumeDataType::IndexType pixel_idx;
			//			pixel_idx[0] = j;
			//			pixel_idx[1] = k;
			//			pixel_idx[2] = i;
			//			auto pixv1 = volume_data->GetPixel(pixel_idx);
			//			//volume_data->SetPixel(pixel_idx,89);
			//			auto pixv=volume_data->GetPixel(pixel_idx);
			//			if (pixv != -1000&&pixv!=-999)
			//			{
			//				int b = 100;
			//			}
			//			int a = 1;
			//		}
			//	}
			//}
			//

			/*volume_data->SetReleaseDataFlag(false);
			volume_data->SetGlobalReleaseDataFlag(false);
			volume_data->ReleaseDataFlagOff();*/
			
			/*auto region = volume_data->GetLargestPossibleRegion();
			auto size = region.GetSize();
			std::cerr << size << std::endl;*/
			return volume_data;
		
		}
	
	}
	catch (itk::ExceptionObject &ex)
	{
		std::cout << ex << std::endl;
		return NULL;
	}

}
bool CDataIO::WriteMesh(std::string fname, CMeshObject & res_mesh_obj)
{
	Eigen::MatrixXd V;
	Eigen::MatrixXi F;
	CConverter::ConvertFromOpenMeshToIGL(res_mesh_obj.GetMesh(), V, F);
	igl::writeOBJ(fname, V, F);
	/*if (!OpenMesh::IO::write_mesh(res_mesh_obj.GetMesh(), fname))
	{
		std::cerr << "write error\n";
		return false;
	}*/
	return true;
}


