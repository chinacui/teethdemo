#include"volume_data_object.h"
int CVolumeDataObject::GetId()
{
	return id_;
}
ItkVolumeDataType::Pointer CVolumeDataObject::GetVolumeData()
{
	return volume_data_;
}
void CVolumeDataObject::SetVolumeData(ItkVolumeDataType::Pointer pdata)
{
	volume_data_ = pdata;
}

OpenMesh::Vec3d CVolumeDataObject::GetDataSize()
{
	auto region = volume_data_->GetLargestPossibleRegion();
	auto size = region.GetSize();
	return OpenMesh::Vec3d(size[0], size[1], size[2]);
}
void CVolumeDataObject::GenTestData()
{
	data_height_ = 400;
	data_depth_ = 400;
	data_width_ = 400;
	vdata_.resize(data_height_*data_width_*data_depth_ ,OpenMesh::Vec4d(0.2,0.2,0.2,0.01));
	for (int d = 100; d < 200; d++)
	{
		for (int w = 100; w < 200; w++)
		{
			for (int h = 100; h < 200; h++)
			{
				vdata_[(d*data_height_*data_width_ + w*data_height_ + h)] = OpenMesh::Vec4d(1,0,0,0.9);

			}
		}
	}
}
CVolumeDataObject::CVolumeDataObject()
{
	id_ = -1;
	is_changed_ = true;
	//GenTestData();
	//volume_data_ = NULL;
}
void CVolumeDataObject::SetChanged(bool is_changed)
{
	is_changed_ = is_changed;
}
CVolumeDataObject::~CVolumeDataObject()
{
	
	//volume_data_->ReleaseData();

}
void CVolumeDataObject::SetId(int id)
{
	id_ = id;
}