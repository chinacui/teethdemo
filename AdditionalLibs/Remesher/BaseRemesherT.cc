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
//  CLASS BaseRemesherT - IMPLEMENTATION
//
//=============================================================================

#define BASE_REMESHERT_C


//== INCLUDES =================================================================


#include "BaseRemesherT.hh"

// Remshing
#include "DiffGeoT.hh"

// OpenMesh
#include <OpenMesh/Core/Utils/Property.hh>
#include <OpenMesh/Core/System/omstream.hh>

// Selection Stuff
#include "MeshSelectionT.hh"


//== NAMESPACES ===============================================================

namespace Remeshing {


//== IMPLEMENTATION ==========================================================


template <class Mesh>
BaseRemesherT<Mesh>::
BaseRemesherT(Mesh& _mesh/*, ProgressEmitter* _progress*/)
  : mesh_(_mesh), refmesh_(0), bsp_(0)/*, progress_(_progress)*/ {

}


//-----------------------------------------------------------------------------


template <class Mesh>
BaseRemesherT<Mesh>::
~BaseRemesherT() {

    delete bsp_;
    delete refmesh_;
}


//-----------------------------------------------------------------------------


template <class Mesh>
void
BaseRemesherT<Mesh>::
init_reference()
{
  typename Mesh::VIter     v_it, v_end;
  typename Mesh::FIter     f_it, f_end;
  typename Mesh::CFVIter   fv_it;
  typename Mesh::VHandle   v0, v1, v2;


  // clear old stuff first
  if (bsp_) delete_reference();

  // tag non-locked vertices + one ring
  for (v_it=mesh_.vertices_begin(), v_end=mesh_.vertices_end();
      v_it!=v_end; ++v_it)
    mesh_.status(*v_it).set_tagged(!mesh_.status(*v_it).locked());

  MeshSelection::growVertexSelection(&mesh_);

  // create reference mesh
  refmesh_ = new Mesh();

  OpenMesh::VPropHandleT<typename Mesh::VHandle>  vhandles;
  mesh_.add_property(vhandles);

  for (v_it=mesh_.vertices_begin(), v_end=mesh_.vertices_end();
      v_it!=v_end; ++v_it)
    if (mesh_.status(*v_it).tagged())
      mesh_.property(vhandles, *v_it) = refmesh_->add_vertex(mesh_.point(*v_it));

  for (f_it=mesh_.faces_begin(), f_end=mesh_.faces_end();
      f_it!=f_end; ++f_it)
  {
    fv_it = mesh_.cfv_iter(*f_it);
    v0    = *fv_it;
    v1    = *(++fv_it);
    v2    = *(++fv_it);

    if (mesh_.status(v0).tagged() &&
        mesh_.status(v1).tagged() &&
        mesh_.status(v2).tagged())
    {
      refmesh_->add_face( mesh_.property(vhandles, v0),
          mesh_.property(vhandles, v1),
          mesh_.property(vhandles, v2) );
    }
  }

  MeshSelection::clearVertexSelection(&mesh_);
  mesh_.remove_property(vhandles);

  // need vertex normals
  refmesh_->request_face_normals();
  refmesh_->request_vertex_normals();
  refmesh_->update_normals(/*Mesh::ANGLE_WEIGHTED*/);


  // build BSP
  bsp_ = new BSP(*refmesh_);
  bsp_->reserve(refmesh_->n_faces());
  for (f_it=refmesh_->faces_begin(), f_end=refmesh_->faces_end();
       f_it!=f_end; ++f_it)
    bsp_->push_back(*f_it);
  bsp_->build(10, 100);
}


//-----------------------------------------------------------------------------


template <class Mesh>
void
BaseRemesherT<Mesh>::
delete_reference()
{
  delete bsp_;      bsp_     = 0;
  delete refmesh_;  refmesh_ = 0;
}


//-----------------------------------------------------------------------------


template <class Mesh>
void
BaseRemesherT<Mesh>::
project_to_reference(typename Mesh::VertexHandle _vh) const
{
  typename Mesh::Point          p = mesh_.point(_vh);
  typename Mesh::FaceHandle     fh = bsp_->nearest(p).handle;

  if ( ! fh.is_valid() )
    return;

  typename Mesh::CFVIter        fv_it = refmesh_->cfv_iter(fh);


  const typename Mesh::Point&   p0 = refmesh_->point(*fv_it);
  typename Mesh::Normal         n0 = refmesh_->normal(*fv_it);
  const typename Mesh::Point&   p1 = refmesh_->point(*(++fv_it));
  typename Mesh::Normal         n1 = refmesh_->normal(*fv_it);
  const typename Mesh::Point&   p2 = refmesh_->point(*(++fv_it));
  typename Mesh::Normal         n2 = refmesh_->normal(*fv_it);


  // project
  ACG::Geometry::distPointTriangle(p, p0, p1, p2, p);
  mesh_.set_point(_vh, p);


  // get barycentric coordinates
  if (!ACG::Geometry::baryCoord(p, p0, p1, p2, p))
    p[0] = p[1] = p[2] = 1.0/3.0;


  // interpolate normal
  typename Mesh::Normal n;
  n  = (n0 *= p[0]);
  n += (n1 *= p[1]);
  n += (n2 *= p[2]);
  n.normalize();
  mesh_.set_normal(_vh, n);
}


//-----------------------------------------------------------------------------


template <class Mesh>
void
BaseRemesherT<Mesh>::
remesh(unsigned int           _iters,
       unsigned int           _area_iters,
       bool                   _use_projection,
       Selection              _selection) {

    try
    {
        if (_selection == VERTEX_SELECTION)
          prepare_vertex_selection();
        else if (_selection == FACE_SELECTION)
          prepare_face_selection();
        remeshh(_iters, _area_iters, _use_projection);
    }
    catch (std::bad_alloc&)
    {
    mesh_.clear();
    omerr() << "Remeshig: Out of memory\n";
    }

    cleanup();
}

//-----------------------------------------------------------------------------

template <class Mesh>
void
BaseRemesherT<Mesh>::
prepare_vertex_selection()
{
  typename Mesh::EIter     e_it, e_end;
  typename Mesh::VIter     v_it, v_end;
  typename Mesh::FIter     f_it, f_end;
  typename Mesh::CFVIter   fv_it;
  typename Mesh::CVOHIter  vh_it;
  typename Mesh::VHandle   v0, v1;
  typename Mesh::FHandle   f0, f1;


  // need vertex and edge status
  mesh_.request_vertex_status();
  mesh_.request_edge_status();
  mesh_.request_face_status();



  // if nothing selected -> select all
  nothing_selected_ = true;
  for (v_it=mesh_.vertices_begin(), v_end=mesh_.vertices_end();
       v_it!=v_end; ++v_it)
    if (mesh_.status(*v_it).selected())
    { nothing_selected_ = false; break; }

  if (nothing_selected_)
    for (v_it=mesh_.vertices_begin(), v_end=mesh_.vertices_end();
   v_it!=v_end; ++v_it)
      mesh_.status(*v_it).set_selected(true);



  // lock un-selected vertices & edges
  for (v_it=mesh_.vertices_begin(), v_end=mesh_.vertices_end();
       v_it!=v_end; ++v_it)
    mesh_.status(*v_it).set_locked(!mesh_.status(*v_it).selected());

  for (e_it=mesh_.edges_begin(), e_end=mesh_.edges_end();
       e_it!=e_end; ++e_it)
  {
    v0 = mesh_.to_vertex_handle(mesh_.halfedge_handle(*e_it, 0));
    v1 = mesh_.to_vertex_handle(mesh_.halfedge_handle(*e_it, 1));
    mesh_.status(*e_it).set_locked(mesh_.status(v0).locked() ||
          mesh_.status(v1).locked());
  }



  // handle feature corners:
  // lock corner vertices (>2 feature edges) and endpoints
  // of feature lines (1 feature edge)
  for (v_it=mesh_.vertices_begin(), v_end=mesh_.vertices_end();
      v_it!=v_end; ++v_it)
  {
    if (mesh_.status(*v_it).feature())
    {
      int c=0;
      for (vh_it=mesh_.cvoh_iter(*v_it); vh_it.is_valid(); ++vh_it)
        if (mesh_.status(mesh_.edge_handle(*vh_it)).feature())
          ++c;
      if (c!=2) mesh_.status(*v_it).set_locked(true);
    }
  }



  // build reference mesh
  init_reference();
//   if (emit_progress_)  Progress().step(5);


  // add properties
  mesh_.add_property(valences_);
  mesh_.add_property(update_);
  mesh_.add_property(area_);
}

//-----------------------------------------------------------------------------


template <class Mesh>
void
BaseRemesherT<Mesh>::
prepare_face_selection()
{
  typename Mesh::EIter     e_it, e_end;
  typename Mesh::VIter     v_it, v_end;
  typename Mesh::FIter     f_it, f_end;
  typename Mesh::CFVIter   fv_it;
  typename Mesh::CVOHIter  vh_it;
  typename Mesh::VHandle   v0, v1;
  typename Mesh::FHandle   f0, f1;


  // need vertex and edge status
  mesh_.request_vertex_status();
  mesh_.request_edge_status();
  mesh_.request_face_status();

  // if nothing selected -> select all
  nothing_selected_ = true;
  for (f_it = mesh_.faces_begin(), f_end = mesh_.faces_end(); f_it != f_end;
      ++f_it)
  {
    if (mesh_.status(*f_it).selected())
    {
      nothing_selected_ = false;
      break;
    }
  }

  if (nothing_selected_)
    MeshSelection::selectAllFaces(&mesh_);



  // lock un-selected vertices & edges
  for (v_it = mesh_.vertices_begin(), v_end = mesh_.vertices_end();
      v_it != v_end; ++v_it)
  {
    bool all_faces_selected = true;

    for (typename Mesh::ConstVertexFaceIter vf_it = mesh_.cvf_iter(*v_it);
        vf_it.is_valid(); ++vf_it)
    {
      if (!mesh_.status(*vf_it).selected())
      {
        all_faces_selected = false;
        break;
      }
    }
    mesh_.status(*v_it).set_locked(!all_faces_selected);
  }

  for (e_it = mesh_.edges_begin(), e_end = mesh_.edges_end();
      e_it != e_end; ++e_it)
  {
    if (mesh_.is_boundary(*e_it))
    {
      mesh_.status(*e_it).set_locked(true);
    }
    else
    {
      f0 = mesh_.face_handle(mesh_.halfedge_handle(*e_it, 0));
      f1 = mesh_.face_handle(mesh_.halfedge_handle(*e_it, 1));

      mesh_.status(*e_it).set_locked(!(mesh_.status(f0).selected() && mesh_.status(f1).selected()));
    }
  }

  MeshSelection::clearFaceSelection(&mesh_);

  // handle feature corners:
  // lock corner vertices (>2 feature edges) and endpoints
  // of feature lines (1 feature edge)
  for (v_it=mesh_.vertices_begin(), v_end=mesh_.vertices_end();
      v_it!=v_end; ++v_it)
  {
    if (mesh_.status(*v_it).feature())
    {
      int c=0;
      for (vh_it=mesh_.cvoh_iter(*v_it); vh_it.is_valid(); ++vh_it)
        if (mesh_.status(mesh_.edge_handle(*vh_it)).feature())
          ++c;
      if (c!=2) mesh_.status(*v_it).set_locked(true);
    }
  }



  // build reference mesh
  init_reference();
//   if (emit_progress_)  Progress().step(5);


  // add properties
  mesh_.add_property(valences_);
  mesh_.add_property(update_);
  mesh_.add_property(area_);
}


//-----------------------------------------------------------------------------


template <class Mesh>
void
BaseRemesherT<Mesh>::
remeshh(unsigned int           _iters,
	unsigned int _area_iters, bool _use_projection) {

    double progress_step = 100.0 / (((double)_iters/4.0 + (double)_area_iters));
    double total_progress = 0.0;

    for (unsigned int i = 0; i < _iters; ++i) {

        split_long_edges();
        //if(progress_ != NULL) {
        //    total_progress += progress_step;
        //    progress_->sendProgressSignal(total_progress);
        //}

        collapse_short_edges();
        //if(progress_ != NULL) {
        //    total_progress += progress_step;
        //    progress_->sendProgressSignal(total_progress);
        //}

        flip_edges();
        //if(progress_ != NULL) {
        //    total_progress += progress_step;
        //    progress_->sendProgressSignal(total_progress);
        //}

        tangential_smoothing(_use_projection);
        //if(progress_ != NULL) {
        //    total_progress += progress_step;
        //    progress_->sendProgressSignal(total_progress);
        //}
    }

    if (_area_iters) {
        balanace_area(_area_iters, _use_projection);
        //if(progress_ != NULL) {
        //    total_progress += progress_step;
        //    progress_->sendProgressSignal(total_progress);
        //}
    }

    // feature edges block flips -> might lead to caps
    remove_caps();
}


//-----------------------------------------------------------------------------


template <class Mesh>
void
BaseRemesherT<Mesh>::
cleanup()
{
  typename Mesh::EIter     e_it, e_end;
  typename Mesh::VIter     v_it, v_end;
  typename Mesh::FIter     f_it, f_end;
  typename Mesh::CFVIter   fv_it;
  typename Mesh::CVOHIter  vh_it;
  typename Mesh::VHandle   v0, v1;
  typename Mesh::FHandle   f0, f1;


  // remove properties
  mesh_.remove_property(valences_);
  mesh_.remove_property(update_);
  mesh_.remove_property(area_);


  // remove reference again
  delete_reference();


  // unlock all vertices & edges
  for (v_it=mesh_.vertices_begin(), v_end=mesh_.vertices_end();
       v_it!=v_end; ++v_it)
  {
    mesh_.status(*v_it).set_locked(false);
  }
  for (e_it=mesh_.edges_begin(), e_end=mesh_.edges_end();
       e_it!=e_end; ++e_it)
  {
    mesh_.status(*e_it).set_locked(false);
  }


  // de-select if nothing was selected before
  if (nothing_selected_)
    for (v_it=mesh_.vertices_begin(), v_end=mesh_.vertices_end();
	 v_it!=v_end; ++v_it)
      mesh_.status(*v_it).set_selected(false);



  // free vertex and edge status
  mesh_.release_vertex_status();
  mesh_.release_edge_status();
  mesh_.release_face_status();
}


//-----------------------------------------------------------------------------


template <class Mesh>
void
BaseRemesherT<Mesh>::
split_long_edges()
{
  typename Mesh::VIter v_it, v_end;
  typename Mesh::EIter e_it, e_end;
  typename Mesh::VHandle v0, v1, vh;
  typename Mesh::EHandle eh, e0, e1;
  typename Mesh::FHandle f0, f1, f2, f3;

  bool ok, is_feature, is_boundary;
  int i;

  // un-tagg
  for (v_it = mesh_.vertices_begin(), v_end = mesh_.vertices_end(); v_it != v_end; ++v_it)
    mesh_.status(*v_it).set_tagged(false);

  // handle Nastran PIDs during refinement
  OpenMesh::FPropHandleT<unsigned int> pids;
  mesh_.get_property_handle(pids, "Nastran PIDs");

  // split long edges
  for (ok = false, i = 0; !ok && i < 100; ++i) {
    ok = true;

    for (e_it = mesh_.edges_begin(), e_end = mesh_.edges_end(); e_it != e_end; ++e_it) {
      v0 = mesh_.to_vertex_handle(mesh_.halfedge_handle(*e_it, 0));
      v1 = mesh_.to_vertex_handle(mesh_.halfedge_handle(*e_it, 1));

      if (!mesh_.status(*e_it).locked() && is_too_long(v0, v1)) {
        const typename Mesh::Point& p0 = mesh_.point(v0);
        const typename Mesh::Point& p1 = mesh_.point(v1);

        is_feature = mesh_.status(*e_it).feature();
        is_boundary = mesh_.is_boundary(*e_it);

        vh = mesh_.add_vertex((p0 + p1) * 0.5);
        mesh_.split(*e_it, vh);

        mesh_.status(vh).set_selected(true);

        if (is_feature) {
          eh = (is_boundary ? mesh_.edge_handle(mesh_.n_edges() - 2) : mesh_.edge_handle(mesh_.n_edges() - 3));

          mesh_.status(eh).set_feature(true);
          mesh_.status(vh).set_feature(true);
        } else {
          project_to_reference(vh);
        }

        if (pids.is_valid()) {
          e0 = *e_it;
          e1 = (is_boundary ? mesh_.edge_handle(mesh_.n_edges() - 2) : mesh_.edge_handle(mesh_.n_edges() - 3));

          f0 = mesh_.face_handle(mesh_.halfedge_handle(e0, 0));
          f1 = mesh_.face_handle(mesh_.halfedge_handle(e0, 1));
          f2 = mesh_.face_handle(mesh_.halfedge_handle(e1, 0));
          f3 = mesh_.face_handle(mesh_.halfedge_handle(e1, 1));

          if (f0.is_valid() && f3.is_valid()) mesh_.property(pids, f3) = mesh_.property(pids, f0);

          if (f1.is_valid() && f2.is_valid()) mesh_.property(pids, f2) = mesh_.property(pids, f1);
        }

        ok = false;
      }
    }
  }

  if (i == 100) omlog() << "split break\n";
}


//-----------------------------------------------------------------------------


template <class Mesh>
void
BaseRemesherT<Mesh>::
collapse_short_edges()
{
  typename Mesh::EIter e_it, e_end;
  typename Mesh::CVVIter vv_it;
  typename Mesh::VHandle v0, v1;
  typename Mesh::HHandle h0, h1, h01, h10;
  typename Mesh::Point p;
  bool ok, skip, b0, b1, l0, l1, f0, f1;
  int i;
  bool hcol01, hcol10;

  for (ok = false, i = 0; !ok && i < 100; ++i) {
    ok = true;

    for (e_it = mesh_.edges_begin(), e_end = mesh_.edges_end(); e_it != e_end; ++e_it) {
      if (!mesh_.status(*e_it).deleted() && !mesh_.status(*e_it).locked()) {
        h10 = mesh_.halfedge_handle(*e_it, 0);
        h01 = mesh_.halfedge_handle(*e_it, 1);
        v0 = mesh_.to_vertex_handle(h10);
        v1 = mesh_.to_vertex_handle(h01);

        if (is_too_short(v0, v1)) {
          // get status
          b0 = mesh_.is_boundary(v0);
          b1 = mesh_.is_boundary(v1);
          l0 = mesh_.status(v0).locked();
          l1 = mesh_.status(v1).locked();
          f0 = mesh_.status(v0).feature();
          f1 = mesh_.status(v1).feature();
          hcol01 = hcol10 = true;

          if (mesh_.status(*e_it).feature() && !f0 && !f1)
            std::cerr << "Bad luck" << std::endl;

          // boundary rules
          if (b0 && b1) {
            if (!mesh_.is_boundary(*e_it))
              continue;
          } else if (b0)
            hcol01 = false;
          else if (b1)
            hcol10 = false;

          // locked rules
          if (l0 && l1)
            continue;
          else if (l0)
            hcol01 = false;
          else if (l1)
            hcol10 = false;

          // feature rules

          // the other two edges removed by collapse must not be features
          h0 = mesh_.prev_halfedge_handle(h01);
          h1 = mesh_.next_halfedge_handle(h10);
          if (mesh_.status(mesh_.edge_handle(h0)).feature() || mesh_.status(mesh_.edge_handle(h1)).feature())
            hcol01 = false;
          h0 = mesh_.prev_halfedge_handle(h10);
          h1 = mesh_.next_halfedge_handle(h01);
          if (mesh_.status(mesh_.edge_handle(h0)).feature() || mesh_.status(mesh_.edge_handle(h1)).feature())
            hcol10 = false;

          if (f0 && f1) {
            // edge must be feature
            if (!mesh_.status(*e_it).feature())
              continue;

            // the other two edges removed by collapse must not be features
            h0 = mesh_.prev_halfedge_handle(h01);
            h1 = mesh_.next_halfedge_handle(h10);
            if (mesh_.status(mesh_.edge_handle(h0)).feature() || mesh_.status(mesh_.edge_handle(h1)).feature())
              hcol01 = false;
            h0 = mesh_.prev_halfedge_handle(h10);
            h1 = mesh_.next_halfedge_handle(h01);
            if (mesh_.status(mesh_.edge_handle(h0)).feature() || mesh_.status(mesh_.edge_handle(h1)).feature())
              hcol10 = false;
          } else if (f0)
            hcol01 = false;
          else if (f1)
            hcol10 = false;

          // topological rules
          if (hcol01)
            hcol01 = mesh_.is_collapse_ok(h01);
          if (hcol10)
            hcol10 = mesh_.is_collapse_ok(h10);

          // both collapses possible: collapse into vertex w/ higher valence
          if (hcol01 && hcol10) {
            if (mesh_.valence(v0) < mesh_.valence(v1))
              hcol10 = false;
            else
              hcol01 = false;
          }

          // try v1 -> v0
          if (hcol10) {
            // don't create too long edges
            skip = false;
            for (vv_it = mesh_.cvv_iter(v1); vv_it.is_valid() && !skip; ++vv_it)
              if (is_too_long(v0, *vv_it))
                skip = true;

            if (!skip) {
              mesh_.collapse(h10);
              ok = false;
            }
          }

          // try v0 -> v1
          else if (hcol01) {
            // don't create too long edges
            skip = false;
            for (vv_it = mesh_.cvv_iter(v0); vv_it.is_valid() && !skip; ++vv_it)
              if (is_too_long(v1, *vv_it))
                skip = true;

            if (!skip) {
              mesh_.collapse(h01);
              ok = false;
            }
          }
        }
      }
    }
  }

  mesh_.garbage_collection();

  if (i == 100)
    omlog() << "collapse break\n";
}


//-----------------------------------------------------------------------------


template <class Mesh>
void
BaseRemesherT<Mesh>::
flip_edges()
{
  typename Mesh::EIter e_it, e_end;
  typename Mesh::VIter v_it, v_end;
  typename Mesh::VHandle v0, v1, v2, v3, vh;
  typename Mesh::HHandle hh;
  typename Mesh::Point p;
  typename Mesh::FHandle fh;

  int val0, val1, val2, val3;
  int val_opt0, val_opt1, val_opt2, val_opt3;
  int ve0, ve1, ve2, ve3, ve_before, ve_after;
  bool ok;
  int i;

  // compute vertex valences
  for (v_it = mesh_.vertices_begin(), v_end = mesh_.vertices_end(); v_it != v_end; ++v_it)
    mesh_.property(valences_, *v_it) = mesh_.valence(*v_it);

  // flip all edges
  for (ok = false, i = 0; !ok && i < 100; ++i) {
    ok = true;

    for (e_it = mesh_.edges_begin(), e_end = mesh_.edges_end(); e_it != e_end; ++e_it) {
      if (!mesh_.status(*e_it).locked() && !mesh_.status(*e_it).feature()) {
        hh = mesh_.halfedge_handle(*e_it, 0);
        v0 = mesh_.to_vertex_handle(hh);
        v2 = mesh_.to_vertex_handle(mesh_.next_halfedge_handle(hh));
        if ( !mesh_.next_halfedge_handle(hh).is_valid() ) {
          std::cerr << "Error v2" << std::endl;
          continue;
        }
        hh = mesh_.halfedge_handle(*e_it, 1);
        v1 = mesh_.to_vertex_handle(hh);
        v3 = mesh_.to_vertex_handle(mesh_.next_halfedge_handle(hh));
        if ( !mesh_.next_halfedge_handle(hh).is_valid() ) {
          std::cerr << "Error v3" << std::endl;
          continue;
        }

        if ( !v2.is_valid())
          continue;

        if (!mesh_.status(v0).locked() && !mesh_.status(v1).locked() && !mesh_.status(v2).locked()
            && !mesh_.status(v3).locked()) {
          val0 = mesh_.property(valences_, v0);
          val1 = mesh_.property(valences_, v1);
          val2 = mesh_.property(valences_, v2);
          val3 = mesh_.property(valences_, v3);

          val_opt0 = (mesh_.is_boundary(v0) ? 4 : 6);
          val_opt1 = (mesh_.is_boundary(v1) ? 4 : 6);
          val_opt2 = (mesh_.is_boundary(v2) ? 4 : 6);
          val_opt3 = (mesh_.is_boundary(v3) ? 4 : 6);

          ve0 = (val0 - val_opt0);
          ve0 *= ve0;
          ve1 = (val1 - val_opt1);
          ve1 *= ve1;
          ve2 = (val2 - val_opt2);
          ve2 *= ve2;
          ve3 = (val3 - val_opt3);
          ve3 *= ve3;

          ve_before = ve0 + ve1 + ve2 + ve3;

          --val0;
          --val1;
          ++val2;
          ++val3;

          ve0 = (val0 - val_opt0);
          ve0 *= ve0;
          ve1 = (val1 - val_opt1);
          ve1 *= ve1;
          ve2 = (val2 - val_opt2);
          ve2 *= ve2;
          ve3 = (val3 - val_opt3);
          ve3 *= ve3;

          ve_after = ve0 + ve1 + ve2 + ve3;

          if (ve_before > ve_after && mesh_.is_flip_ok(*e_it)) {
            mesh_.flip(*e_it);
            --mesh_.property(valences_, v0);
            --mesh_.property(valences_, v1);
            ++mesh_.property(valences_, v2);
            ++mesh_.property(valences_, v3);
            ok = false;
          }
        }
      }
    }
  }

  if (i == 100)
    omlog() << "flip break\n";
}


//-----------------------------------------------------------------------------


template <class Mesh>
void
BaseRemesherT<Mesh>::
tangential_smoothing(bool _use_projection)
{
  typename Mesh::VIter     v_it, v_end(mesh_.vertices_end());
  typename Mesh::CVVIter   vv_it;
  typename Mesh::Scalar    valence;
  typename Mesh::Point     u, n;


  // tag active vertices
  for (v_it=mesh_.vertices_begin(); v_it!=v_end; ++v_it)
    mesh_.status(*v_it).set_tagged( !mesh_.status(*v_it).locked() &&
				   !mesh_.status(*v_it).feature() &&
				   !mesh_.is_boundary(*v_it) );


  // smooth
  for (int iters=0; iters<10; ++iters)
  {
    for (v_it=mesh_.vertices_begin(); v_it!=v_end; ++v_it)
    {
      if (mesh_.status(*v_it).tagged())
      {
	u.vectorize(0.0);
	valence = 0;

	for (vv_it=mesh_.cvv_iter(*v_it); vv_it.is_valid(); ++vv_it)
	{
	  u += mesh_.point(*vv_it);
	  ++valence;
	}

	if (valence)
	{
	  u *= (1.0/valence);
	  u -= mesh_.point(*v_it);
	  n  = mesh_.normal(*v_it);
	  n *= (u | n);
	  u -= n;
	}

	mesh_.property(update_, *v_it) = u;
      }
    }

    for (v_it=mesh_.vertices_begin(); v_it!=v_end; ++v_it)
      if (mesh_.status(*v_it).tagged())
	mesh_.point(*v_it) += mesh_.property(update_, *v_it);
  }


  // reset tagged status
  for (v_it=mesh_.vertices_begin(); v_it!=v_end; ++v_it)
    mesh_.status(*v_it).set_tagged(false);


  // project
  if (_use_projection)
    for (v_it=mesh_.vertices_begin(); v_it!=v_end; ++v_it)
      if (!mesh_.status(*v_it).locked() && !mesh_.status(*v_it).feature())
	project_to_reference(*v_it);
}


//-----------------------------------------------------------------------------


template <class Mesh>
void
BaseRemesherT<Mesh>::
balanace_area(unsigned int _iters, bool _use_projection)
{
  typename Mesh::VIter     v_it, v_end(mesh_.vertices_end());
  typename Mesh::CVVIter   vv_it;
  typename Mesh::Scalar    w, ww;
  typename Mesh::Point     u, n;


  DiffGeoT<Mesh>  diffgeo(mesh_);


  // tag active vertices
  for (v_it=mesh_.vertices_begin(); v_it!=v_end; ++v_it)
  {
    bool active = ( !mesh_.status(*v_it).locked() &&
        !mesh_.status(*v_it).feature() &&
        !mesh_.is_boundary(*v_it) );

    // don't move neighbors of boundary vertices
    for (vv_it=mesh_.cvv_iter(*v_it); active && vv_it.is_valid(); ++vv_it)
      if (mesh_.is_boundary(*vv_it))
        active = false;

    mesh_.status(*v_it).set_tagged( active );
  }


  // tag2 vertices for which area has to be computed
  for (v_it=mesh_.vertices_begin(); v_it!=v_end; ++v_it)
  {
    mesh_.status(*v_it).set_tagged2(false);
    for (vv_it=mesh_.cvv_iter(*v_it); vv_it.is_valid(); ++vv_it)
    {
      if (mesh_.status(*vv_it).tagged())
      {
        mesh_.status(*v_it).set_tagged2(true);
        break;
      }
    }
  }



  for (unsigned int bla=0; bla<_iters; ++bla)
  {
    // compute vertex areas
    for (v_it=mesh_.vertices_begin(); v_it!=v_end; ++v_it)
      if (mesh_.status(*v_it).tagged2())
        mesh_.property(area_, *v_it) = pow(diffgeo.compute_area(*v_it), 2);



    // smooth
    for (int iters=0; iters<10; ++iters)
    {
      for (v_it=mesh_.vertices_begin(); v_it!=v_end; ++v_it)
      {
        if (mesh_.status(*v_it).tagged())
        {
          u.vectorize(0.0);
          ww   = 0;

          for (vv_it=mesh_.cvv_iter(*v_it); vv_it.is_valid(); ++vv_it)
          {
            w   = mesh_.property(area_, *vv_it);
            u  += mesh_.point(*vv_it) * w;
            ww += w;
          }

          if (ww)
          {
            u *= (1.0/ww);
            u -= mesh_.point(*v_it);
            n  = mesh_.normal(*v_it);
            n *= (u | n);
            u -= n;

            u *= 0.1;
          }

          mesh_.property(update_, *v_it) = u;
        }
      }

      for (v_it=mesh_.vertices_begin(); v_it!=v_end; ++v_it)
        if (!mesh_.status(*v_it).locked() &&
            !mesh_.status(*v_it).feature() &&
            !mesh_.is_boundary(*v_it))
          mesh_.point(*v_it) += mesh_.property(update_, *v_it);
    }


    // project
    if (_use_projection)
      for (v_it=mesh_.vertices_begin(); v_it!=v_end; ++v_it)
        if (!mesh_.status(*v_it).locked() && !mesh_.status(*v_it).feature())
          project_to_reference(*v_it);


//     if (emit_progress_)  
//       if (!Progress().step(3))
// 	break;
  }
}


//-----------------------------------------------------------------------------


template <class Mesh>
void
BaseRemesherT<Mesh>::
remove_caps()
{
  typename Mesh::EdgeIter        e_it, e_end(mesh_.edges_end());
  typename Mesh::HalfedgeHandle  h;
  typename Mesh::Scalar          a0, a1, amin, aa(cos(170.0 * M_PI / 180.0));
  typename Mesh::VHandle         v, vb, vd;
  typename Mesh::FHandle         fb, fd;
  typename Mesh::Point           a, b, c, d;


  // handle Nastran PIDs for edge flips
  OpenMesh::FPropHandleT<unsigned int> pids;
  mesh_.get_property_handle(pids, "Nastran PIDs");


  for (e_it=mesh_.edges_begin(); e_it!=e_end; ++e_it)
  {
    if (!mesh_.status(*e_it).locked() &&
        mesh_.is_flip_ok(*e_it))
    {
      h  = mesh_.halfedge_handle(*e_it, 0);
      a  = mesh_.point(mesh_.to_vertex_handle(h));

      h  = mesh_.next_halfedge_handle(h);
      b  = mesh_.point(vb=mesh_.to_vertex_handle(h));

      h  = mesh_.halfedge_handle(*e_it, 1);
      c  = mesh_.point(mesh_.to_vertex_handle(h));

      h  = mesh_.next_halfedge_handle(h);
      d  = mesh_.point(vd=mesh_.to_vertex_handle(h));

      a0 = ( (a-b).normalize() | (c-b).normalize() );
      a1 = ( (a-d).normalize() | (c-d).normalize() );

      if (a0 < a1)  { amin = a0; v = vb; }
      else          { amin = a1; v = vd; }

      // is it a cap?
      if (amin < aa)
      {
        // feature edge and feature vertex -> seems to be intended
        if (mesh_.status(*e_it).feature() && mesh_.status(v).feature())
          continue;

        // handle PIDs: flip = split + collapse
        if (pids.is_valid())
        {
          fb = mesh_.face_handle(mesh_.halfedge_handle(*e_it, 0));
          fd = mesh_.face_handle(mesh_.halfedge_handle(*e_it, 1));

          if (v == vb)
            mesh_.property(pids, fb) = mesh_.property(pids, fd);
          else
            mesh_.property(pids, fd) = mesh_.property(pids, fb);
        }

        // project v onto feature edge
        if (mesh_.status(*e_it).feature())
          mesh_.set_point(v, (a+c)*0.5);

        // flip
        mesh_.flip(*e_it);
      }
    }
  }
}


//=============================================================================
} // namespace Remeshing
//=============================================================================
