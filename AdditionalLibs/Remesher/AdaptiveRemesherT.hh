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
//  CLASS AdaptiveRemesherT
//
//=============================================================================


#ifndef ADAPTIVE_REMESHERT_HH
#define ADAPTIVE_REMESHERT_HH


//== INCLUDES =================================================================

#include "BaseRemesherT.hh"
#include "DiffGeoT.hh"
#include <OpenMesh/Core/Utils/Property.hh>

//#include "../ProgressEmitter.hh"

//== NAMESPACES ===============================================================

namespace Remeshing {


//== CLASS DEFINITION =========================================================


template <class Mesh>
class AdaptiveRemesherT : public BaseRemesherT<Mesh>
{
  typedef typename BaseRemesherT<Mesh>::Selection Selection;
public:

  typedef BaseRemesherT<Mesh>          Base;
  typedef typename Mesh::Scalar        Scalar;
  typedef typename Mesh::Point         Point;
  typedef typename Mesh::EdgeHandle    EdgeHandle;
  typedef typename Mesh::VertexHandle  VertexHandle;


  AdaptiveRemesherT(Mesh& _mesh/*, ProgressEmitter* _progress = NULL*/) : Base(_mesh/*, _progress*/) {}

  void remesh(Scalar        _error,
              Scalar        _min_edge_length,
              Scalar        _max_edge_length,
              unsigned int  _iters,
              bool          _use_projection = true,
              Selection     _selection=BaseRemesherT<Mesh>::VERTEX_SELECTION);



protected:

  virtual void init_reference();
  virtual void project_to_reference(VertexHandle _vh) const;

  void compute_curvature(Mesh& _mesh, OpenMesh::VPropHandleT<Scalar> _ph);

  virtual bool is_too_long  (VertexHandle _v0, VertexHandle _v1) const;
  virtual bool is_too_short (VertexHandle _v0, VertexHandle _v1) const;


protected:

  Scalar  error_, emax_, emin_;

  OpenMesh::VPropHandleT<Scalar>  refcurv_, curvature_;
};


//=============================================================================
} // namespace Remeshing
//=============================================================================
#if defined(INCLUDE_TEMPLATES) && !defined(ADAPTIVE_REMESHERT_C)
#define ADAPTIVE_REMESHERT_TEMPLATES
#include "AdaptiveRemesherT.cc"
#endif
//=============================================================================
#endif // ADAPTIVE_REMESHERT_HH defined
//=============================================================================

