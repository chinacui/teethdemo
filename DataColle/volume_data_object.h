#ifndef CVOLUME_DATA_OBJECT_H
#define CVOLUME_DATA_OBJECT_H
#include"prereq.h"
#include"base_object.h"
#include<vector>
#include<iostream>
#include"custom_openmesh_type.h"
#include"custom_itk_type.h"

class DATACOLLE_CLASS CVolumeDataObject:public CBaseObject
{
protected:
	int id_;
	bool is_changed_;
	
	ItkVolumeDataType::Pointer volume_data_;

	//test
	void GenTestData();
public:
	int data_width_, data_height_, data_depth_;
	std::vector<OpenMesh::Vec4d>vdata_;
	CVolumeDataObject();
	~CVolumeDataObject();

	void SetChanged(bool is_changed=true);
	bool IsChanged() { return is_changed_; }
	int GetId();
	void SetId(int id);

	OpenMesh::Vec3d GetDataSize();
	
	ItkVolumeDataType::Pointer GetVolumeData();
	void SetVolumeData(ItkVolumeDataType::Pointer pdata);


	
};
#endif