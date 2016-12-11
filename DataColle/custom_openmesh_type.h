#ifndef CUSTOM_OPENMESH_TYPE_H
#define CUSTOM_OPENMESH_TYPE_H
#include "prereq.h"
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Mesh/Traits.hh>
struct COMTraits : public OpenMesh::DefaultTraits
{
	// store barycenter of neighbors in this member
	typedef OpenMesh::Vec3d Point; // use double-values points
	typedef OpenMesh::Vec3d Color; 
	typedef OpenMesh::Vec3d  Normal;
	HalfedgeTraits
	{
	private:
		TexCoord2D  uv_;
		TexCoord3D uv3d_;
		Color color_;
	public:
		HalfedgeT() : uv_(TexCoord2D(0.0f, 0.0f)), color_(Color(0,0, 0)) { }
		const TexCoord2D& GetUV() const { return uv_; }
		const TexCoord3D& GetUV3D() const { return uv3d_; }
		const Color& GetColor() const { return color_; }
		void SetUV(const TexCoord2D& _p) { uv_ = _p; }
		void SetUV3D(const TexCoord3D& _p) { uv3d_ = _p; }
		void SetColor(const TexCoord2D& c) { color_ = c; }
	};
	
	FaceTraits
	{
	private:
		
	public:
		FaceT() {}
		
	};
	
};

class COpenMeshT :public OpenMesh::TriMesh_ArrayKernelT<COMTraits>
{
public:
	COpenMeshT()
	{
		request_vertex_colors();
		request_face_status();
		request_edge_status();
		request_vertex_status();
		request_face_normals();
		request_vertex_normals();
		request_halfedge_texcoords3D();
		request_halfedge_texcoords2D();
		
	}
};

// ---------------------------------------------------------------------------
#endif