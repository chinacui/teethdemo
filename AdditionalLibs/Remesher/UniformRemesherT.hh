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
//  CLASS UniformRemesherT
//
//=============================================================================


#ifndef UNIFORM_REMESHERT_HH
#define UNIFORM_REMESHERT_HH


//== INCLUDES =================================================================

#include "BaseRemesherT.hh"

//#include "../ProgressEmitter.hh"

//== NAMESPACES ===============================================================

namespace Remeshing {


//== CLASS DEFINITION =========================================================


template <class Mesh>
class UniformRemesherT : public BaseRemesherT<Mesh>
{
  typedef typename BaseRemesherT<Mesh>::Selection Selection;
public:

  typedef BaseRemesherT<Mesh>          Base;
  typedef typename Mesh::Scalar        Scalar;
  typedef typename Mesh::Point         Point;
  typedef typename Mesh::EdgeHandle    EdgeHandle;
  typedef typename Mesh::VertexHandle  VertexHandle;


  UniformRemesherT(Mesh& _mesh/*, ProgressEmitter* _progress = NULL*/) : Base(_mesh/*, _progress*/) {}


  void remesh(Scalar        _edge_length,
	      unsigned int  _iters,
	      unsigned int  _area_iters,
	      bool          _use_projection = true,
	      Selection     _selection = BaseRemesherT<Mesh>::VERTEX_SELECTION)
  {
//     emin_ = 4.0/5.0 * _edge_length;  sqr_emin_ = emin_ * emin_;
//     emax_ = 4.0/3.0 * _edge_length;  sqr_emax_ = emax_ * emax_;
    emin_ = 0.7 * _edge_length;  sqr_emin_ = emin_ * emin_;
    emax_ = 1.4 * _edge_length;  sqr_emax_ = emax_ * emax_;
    Base::remesh(_iters, _area_iters, _use_projection,_selection);
  }



protected:

  virtual bool is_too_long  (VertexHandle _v0, VertexHandle _v1) const
  {
    return (Base::mesh_.point(_v0) - 
	    Base::mesh_.point(_v1)).sqrnorm() > sqr_emax_;
  }

  virtual bool is_too_short (VertexHandle _v0, VertexHandle _v1) const
  {
    return (Base::mesh_.point(_v0) - 
	    Base::mesh_.point(_v1)).sqrnorm() < sqr_emin_;
  }


protected:

  Scalar  emax_, sqr_emax_, emin_, sqr_emin_;
};


//=============================================================================
} // namespace Remeshing
//=============================================================================
#endif // UNIFORM_REMESHERT_HH defined
//=============================================================================

