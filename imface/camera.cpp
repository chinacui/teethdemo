#include"camera.h"

#include"QGLViewer/quaternion.h"
#include"qquaternion.h"
#include"qvector3d.h"
CCamera::CCamera() :Camera()
{

}


CCamera::CCamera(const CCamera& c) :Camera(c)
{

	
}

CCamera & CCamera::operator=(const CCamera &b)
{

	/*((Camera)(*this)) = (Camera)b;
	this->ortho_size_ = b.ortho_size_;
	return *this;
	*/

	*((Camera*)this) = (Camera)b;

	return *this;
}
