#ifndef CRIGID_ICP_H
#define CRIGID_ICP_H
#include"prereq.h"
#include <Eigen/Sparse>
#include<ANN/ANN.h>
#include"../DataColle/mesh_object.h"
class  ALGCOLLE_CLASS CRigidIcp
{
protected:
	COpenMeshT* src_mesh_, *tgt_mesh_;
public:
	CRigidIcp(COpenMeshT* src_mesh);
	void FitTarget(COpenMeshT *target_mesh);
};

#endif