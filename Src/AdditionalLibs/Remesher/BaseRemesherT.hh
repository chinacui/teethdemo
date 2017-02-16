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
//  CLASS BaseRemesherT
//
//=============================================================================


#ifndef BASE_REMESHERT_HH
#define BASE_REMESHERT_HH

/**
 *  BaseRemesher implements a modified version of a remesher approach presented by Mario Botsch and Leif Kobbelt in
 *  "A Remeshing Approach to Multiresolution Modeling" published at Symposium on Geometry Processing 2004 p. 189-196
 *  Paper can found at: http://www.graphics.rwth-aachen.de/publications/0000/37/
 */


//== INCLUDES =================================================================

#include <OpenMesh/Core/Utils/Property.hh>
#include "Geometry/bsp/TriangleBSPT.hh"

//#include "../ProgressEmitter.hh"

//== FORWARDDECLARATIONS ======================================================


//== NAMESPACES ===============================================================

namespace Remeshing {


//== CLASS DEFINITION =========================================================


template <class Mesh>
class BaseRemesherT
{
public:

  enum Selection
  {
    VERTEX_SELECTION,
    FACE_SELECTION
  };

  typedef typename Mesh::Scalar        Scalar;
  typedef typename Mesh::Point         Point;
  typedef typename Mesh::EdgeHandle    EdgeHandle;
  typedef typename Mesh::VertexHandle  VertexHandle;


  BaseRemesherT(Mesh& _mesh/*, ProgressEmitter* _progress = NULL*/);
  virtual ~BaseRemesherT();

  void remesh(unsigned int  _iters,
              unsigned int  _area_iters,
              bool          _use_projection = true,
              Selection     _selection=VERTEX_SELECTION);



protected:

  /// prepare for remeshing only selected vertices (if no vertex was selected, remesh whole mesh)
  void prepare_vertex_selection();
  /// prepare for remeshing only vertices which are fully surrounded by selected faces (if no face was selected, remesh whole mesh)
  void prepare_face_selection();
  void remeshh(unsigned int _iters, unsigned int _aiters, bool _proj);
  void cleanup();

  virtual void init_reference();
  virtual void delete_reference();
  virtual void project_to_reference(VertexHandle _vh) const;

  void split_long_edges();
  void collapse_short_edges();
  void flip_edges();
  void tangential_smoothing(bool _use_projection);
  void balanace_area(unsigned int _iters, bool _use_projection);
  void remove_caps();

  virtual bool is_too_long  (VertexHandle _v0, VertexHandle _v1) const = 0;
  virtual bool is_too_short (VertexHandle _v0, VertexHandle _v1) const = 0;


  
protected:

  typedef OpenMeshTriangleBSPT<Mesh>  BSP;

  Mesh&  mesh_;
  Mesh*  refmesh_;
  BSP*   bsp_;
  bool   nothing_selected_;

  OpenMesh::VPropHandleT<int>     valences_;
  OpenMesh::VPropHandleT<Point>   update_;
  OpenMesh::VPropHandleT<Scalar>  area_;

  //ProgressEmitter* progress_;
};


//=============================================================================
} // namespace Remeshing
//=============================================================================
#if defined(INCLUDE_TEMPLATES) && !defined(BASE_REMESHERT_C)
#define BASE_REMESHERT_TEMPLATES
#include "BaseRemesherT.cc"
#endif
//=============================================================================
#endif // BASE_REMESHERT_HH defined
//=============================================================================

