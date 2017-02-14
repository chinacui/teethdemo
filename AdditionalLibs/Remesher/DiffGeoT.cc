/*===========================================================================*\
*                                                                            *
*                              OpenFlipper                                   *
 *           Copyright (c) 2001-2015, RWTH-Aachen University                 *
 *           Department of Computer Graphics and Multimedia                  *
 *                          All rights reserved.                             *
 *                            www.openflipper.org                            *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * This file is part of OpenFlipper.                                         *
 *---------------------------------------------------------------------------*
 *                                                                           *
 * Redistribution and use in source and binary forms, with or without        *
 * modification, are permitted provided that the following conditions        *
 * are met:                                                                  *
 *                                                                           *
 * 1. Redistributions of source code must retain the above copyright notice, *
 *    this list of conditions and the following disclaimer.                  *
 *                                                                           *
 * 2. Redistributions in binary form must reproduce the above copyright      *
 *    notice, this list of conditions and the following disclaimer in the    *
 *    documentation and/or other materials provided with the distribution.   *
 *                                                                           *
 * 3. Neither the name of the copyright holder nor the names of its          *
 *    contributors may be used to endorse or promote products derived from   *
 *    this software without specific prior written permission.               *
 *                                                                           *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS       *
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED *
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A           *
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER *
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,  *
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,       *
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR        *
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF    *
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING      *
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS        *
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.              *
*                                                                            *
\*===========================================================================*/

/*===========================================================================*\
*                                                                            *
*   $Revision$                                                       *
*   $LastChangedBy$                                                *
*   $Date$                    *
*                                                                            *
\*===========================================================================*/

//=============================================================================
//
//  CLASS DiffGeoT - IMPLEMENTATION
//
//=============================================================================

#define DIFFGEOT_C

//== INCLUDES =================================================================

//#include <mb_base/mb_base.hh>
#include "DiffGeoT.hh"

#ifdef WIN32
#undef min
#undef max
#endif


//== NAMESPACES ===============================================================

namespace Remeshing {

//== IMPLEMENTATION ==========================================================


template <class Mesh>
DiffGeoT<Mesh>::
DiffGeoT(Mesh& _mesh) 
  : mesh_(_mesh),
    weights_computed_(false),
    area_computed_(false)
{
}


//-----------------------------------------------------------------------------


template <class Mesh>
DiffGeoT<Mesh>::
~DiffGeoT()
{
  mesh_.remove_property(area_);
  mesh_.remove_property(gauss_curvature_);
  mesh_.remove_property(mean_curvature_);
  mesh_.remove_property(edge_weight_);
}


//-----------------------------------------------------------------------------


template <class Mesh>
void 
DiffGeoT<Mesh>::
compute(unsigned int _post_smoothing_iters)
{
  compute_edge_weights();
  compute_area();
  compute_gauss_curvature();
  compute_mean_curvature();
  post_smoothing(_post_smoothing_iters);
}


//-----------------------------------------------------------------------------


template <class Mesh>
void 
DiffGeoT<Mesh>::
compute_edge_weights()
{
  if (!edge_weight_.is_valid())
    mesh_.add_property(edge_weight_);


  typename Mesh::EdgeIter          e_it, e_end(mesh_.edges_end());
  typename Mesh::HalfedgeHandle    heh2;
  typename Mesh::Scalar            weight;


  for (e_it=mesh_.edges_begin(); e_it!=e_end; ++e_it)
  {
    const typename Mesh::HalfedgeHandle heh0   = mesh_.halfedge_handle(*e_it, 0);
    const typename Mesh::Point&         p0     = mesh_.point(mesh_.to_vertex_handle(heh0));

    const typename Mesh::HalfedgeHandle heh1   = mesh_.halfedge_handle(*e_it, 1);
    const typename Mesh::Point&         p1     = mesh_.point(mesh_.to_vertex_handle(heh1));

    heh2   = mesh_.next_halfedge_handle(heh0);
    const typename Mesh::Point&         p2     = mesh_.point(mesh_.to_vertex_handle(heh2));
    const typename Mesh::Point          d0     = (p0 - p2).normalize();
    const typename Mesh::Point          d1     = (p1 - p2).normalize();
    weight = 1.0 / tan(acos(d0|d1));

    heh2   = mesh_.next_halfedge_handle(heh1);
    const typename Mesh::Point&         p3     = mesh_.point(mesh_.to_vertex_handle(heh2));
    const typename Mesh::Point          d2     = (p0 - p3).normalize();
    const typename Mesh::Point          d3     = (p1 - p3).normalize();
    weight += 1.0 / tan(acos(d2|d3));

    mesh_.property(edge_weight_, *e_it) = weight;
  }


  weights_computed_ = true;
}


//-----------------------------------------------------------------------------


template <class Mesh>
void 
DiffGeoT<Mesh>::
compute_area()
{
  if (!area_.is_valid())
    mesh_.add_property(area_);


  typename Mesh::VertexIter  v_it, v_end(mesh_.vertices_end());


  for (v_it=mesh_.vertices_begin(); v_it!=v_end; ++v_it)
    mesh_.property(area_, *v_it) = compute_area(*v_it);


  area_computed_ = true;
}


//-----------------------------------------------------------------------------


template <class Mesh>
typename Mesh::Scalar
DiffGeoT<Mesh>::
compute_area(VertexHandle _vh) const
{
  typename Mesh::HalfedgeHandle    		heh0, heh1, heh2;
  typename Mesh::VertexOHalfedgeIter  voh_it;

  typename Mesh::Point   P, Q, R, PQ, QR, PR;
  double  normPQ, normQR, normPR;
  double  angleP, angleQ, angleR;
  double  area;
  double  ub(0.999), lb(-0.999);
  const double epsilon( 1e-12 );

  area = 0.0;

  for (voh_it=mesh_.voh_iter(_vh); voh_it.is_valid(); ++voh_it)
  {
    heh0 = *voh_it;
    heh1 = mesh_.next_halfedge_handle(heh0);
    heh2 = mesh_.next_halfedge_handle(heh1);

    if (mesh_.is_boundary(heh0)) continue;

    P = mesh_.point(mesh_.to_vertex_handle(heh2));
    Q = mesh_.point(mesh_.to_vertex_handle(heh0));
    R = mesh_.point(mesh_.to_vertex_handle(heh1));

    (PQ = Q) -= P;
    (QR = R) -= Q;
    (PR = R) -= P;
      
    normPQ = PQ.norm();
    normQR = QR.norm();
    normPR = PR.norm();

    if( normPQ <= epsilon || normQR <= epsilon || normPR <= epsilon )
	continue;

    angleP =  (PQ | PR) / normPQ / normPR;
    angleQ = -(PQ | QR) / normPQ / normQR;
    angleR =  (QR | PR) / normQR / normPR;

    if (angleP > ub) angleP = ub;  else if (angleP < lb) angleP = lb;
    if (angleQ > ub) angleQ = ub;  else if (angleQ < lb) angleQ = lb;
    if (angleR > ub) angleR = ub;  else if (angleR < lb) angleR = lb;

    angleP = acos(angleP);
    angleQ = acos(angleQ);
    angleR = acos(angleR);

    if (angleP >= M_PI_2 || angleQ >= M_PI_2 || angleR >= M_PI_2)
    {
      if (angleP >= M_PI_2)
	area += 0.25 * (PQ % PR).norm();
      else
	area += 0.125 * (PQ % PR).norm();
    }
    else
    {
      area += 0.125 * (normPR * normPR / tan(angleQ) + 
		       normPQ * normPQ / tan(angleR));
    }
  }
    
  return area;
}


//-----------------------------------------------------------------------------


template <class Mesh>
void 
DiffGeoT<Mesh>::
compute_gauss_curvature()
{
  if (!gauss_curvature_.is_valid())
    mesh_.add_property(gauss_curvature_);

  if (!area_computed_)
    compute_area();


  typename Mesh::VertexIter        v_it, v_end(mesh_.vertices_end());
  typename Mesh::VertexVertexIter  vv_it, vv_it2;
  typename Mesh::Point             d0, d1;
  typename Mesh::Scalar            curv, count;
  typename Mesh::Scalar            angles, cos_angle;
  typename Mesh::Scalar            lb(-1.0), ub(1.0);


  // compute for all non-boundary vertices
  for (v_it=mesh_.vertices_begin(); v_it!=v_end; ++v_it)
  {
    if (!mesh_.is_boundary(*v_it))
    {
      angles = 0.0;

      for (vv_it=mesh_.vv_iter(*v_it); vv_it.is_valid(); ++vv_it)
      {
        vv_it2 = vv_it; ++vv_it2;

        d0 = (mesh_.point(*vv_it)  - mesh_.point(*v_it)).normalize();
        d1 = (mesh_.point(*vv_it2) - mesh_.point(*v_it)).normalize();

        cos_angle = std::max(lb, std::min(ub, (d0 | d1)));
        angles += acos(cos_angle);
      }

      mesh_.property(gauss_curvature_, *v_it) = (2*M_PI-angles) / mesh_.property(area_, *v_it);
    }
  }


  // boundary vertices: get from neighbors
  for (v_it=mesh_.vertices_begin(); v_it!=v_end; ++v_it)
  {
    if (mesh_.is_boundary(*v_it))
    {
      curv = count = 0.0;

      for (vv_it=mesh_.vv_iter(*v_it); vv_it.is_valid(); ++vv_it)
        if (!mesh_.is_boundary(*vv_it))
        {
          curv += mesh_.property(gauss_curvature_, *vv_it);
          ++count;
        }

      if (count)
	mesh_.property(gauss_curvature_, *v_it) = curv / count;
    }
  }
}


//-----------------------------------------------------------------------------


template <class Mesh>
void 
DiffGeoT<Mesh>::
compute_mean_curvature()
{
  if (!mean_curvature_.is_valid())
    mesh_.add_property(mean_curvature_);

  if (!weights_computed_)
    compute_edge_weights();

  if (!area_computed_)
    compute_area();


  typename Mesh::VertexIter        		v_it, v_end(mesh_.vertices_end());
  typename Mesh::VertexVertexIter  		vv_it;
  typename Mesh::VertexOHalfedgeIter  voh_it;
  typename Mesh::Scalar            		weight;
  typename Mesh::Point             		umbrella;
  typename Mesh::EdgeHandle        		eh;
  typename Mesh::Scalar            		curv, count;


  // compute for all non-boundary vertices
  for (v_it=mesh_.vertices_begin(); v_it!=v_end; ++v_it)
  {
    if (!mesh_.is_boundary(*v_it))
    {
      umbrella[0] = umbrella[1] = umbrella[2] = 0.0;

      for (voh_it=mesh_.voh_iter(*v_it); voh_it.is_valid(); ++voh_it)
      {
        eh        = mesh_.edge_handle(*voh_it);
        weight    = mesh_.property(edge_weight_, eh);
        umbrella += (mesh_.point(*v_it) - mesh_.point(mesh_.to_vertex_handle(*voh_it))) * weight;
      }

      mesh_.property(mean_curvature_, *v_it) =
          umbrella.norm() / (4.0 * mesh_.property(area_, *v_it));
    }
  }


  // boundary vertices: get from neighbors
  for (v_it=mesh_.vertices_begin(); v_it!=v_end; ++v_it)
  {
    if (mesh_.is_boundary(*v_it))
    {
      curv = count = 0.0;

      for (vv_it=mesh_.vv_iter(*v_it); vv_it.is_valid(); ++vv_it)
        if (!mesh_.is_boundary(*vv_it))
        {
          curv += mesh_.property(mean_curvature_, *vv_it);
          ++count;
        }

      if (count)
        mesh_.property(mean_curvature_, *v_it) = curv / count;
    }
  }
}


//-----------------------------------------------------------------------------


template <class Mesh>
void 
DiffGeoT<Mesh>::
post_smoothing(unsigned int _iters)
{
  // early out
  if (!_iters) return;


  // something should already be there...
  if (! (weights_computed_ && 
	 area_computed_ &&
	 mean_curvature_.is_valid() &&
	 gauss_curvature_.is_valid()) )
  {
    std::cerr << "DiffGeoT::post_smoothing: something is wrong!\n";
    return;
  }



  typename Mesh::VertexIter        		v_it, v_end(mesh_.vertices_end());
  typename Mesh::VertexOHalfedgeIter  voh_it;
  typename Mesh::Scalar            		w, ww;
  typename Mesh::Scalar            		gc, mc;
  typename Mesh::EdgeHandle        		eh;



  // add helper properties to store new values
  OpenMesh::VPropHandleT<Scalar>   new_gauss, new_mean;
  mesh_.add_property(new_gauss);
  mesh_.add_property(new_mean);



  for (unsigned int i=0; i<_iters; ++i)
  {

    // compute new value
    for (v_it=mesh_.vertices_begin(); v_it!=v_end; ++v_it)
    {
      if (!mesh_.is_boundary(*v_it))
      {
        gc = mc = ww = 0.0;

        for (voh_it=mesh_.voh_iter(*v_it); voh_it.is_valid(); ++voh_it)
        {
          eh   = mesh_.edge_handle(*voh_it);
          ww  += (w = mesh_.property(edge_weight_, eh));
          mc  += w * mesh_.property(mean_curvature_,  mesh_.to_vertex_handle(*voh_it));
          gc  += w * mesh_.property(gauss_curvature_, mesh_.to_vertex_handle(*voh_it));
        }

        if (ww)
        {
          mesh_.property(new_mean,  *v_it) = mc / ww;
          mesh_.property(new_gauss, *v_it) = gc / ww;
        }
      }
    }


    // replace old by new value
    for (v_it=mesh_.vertices_begin(); v_it!=v_end; ++v_it)
    {
      if (!mesh_.is_boundary(*v_it))
      {
        mesh_.property(mean_curvature_,  *v_it) =  mesh_.property(new_mean, *v_it);

        mesh_.property(gauss_curvature_,  *v_it) = mesh_.property(new_gauss, *v_it);
      }
    }
  }


  // remove helper properties
  mesh_.remove_property(new_gauss);
  mesh_.remove_property(new_mean);
}


//=============================================================================
} // namespace Remeshing
//=============================================================================
