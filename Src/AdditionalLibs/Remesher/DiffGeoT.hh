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
//  CLASS DiffGeoT
//
//=============================================================================


#ifndef DIFFGEO_HH
#define DIFFGEO_HH


//== INCLUDES =================================================================

#include <OpenMesh/Core/Utils/Property.hh>


//== NAMESPACES ===============================================================

namespace Remeshing {

//== CLASS DEFINITION =========================================================


template <class Mesh>
class DiffGeoT
{
public:

  typedef typename Mesh::Scalar        Scalar;
  typedef typename Mesh::VertexHandle  VertexHandle;

   
  DiffGeoT(Mesh& _mesh);
  ~DiffGeoT();


  void compute(unsigned int _post_smoothing_iters=0);

  void compute_edge_weights();
  void compute_area();
  void compute_gauss_curvature();
  void compute_mean_curvature();
  void post_smoothing(unsigned int _iters);


  Scalar compute_area(VertexHandle _vh) const;


  Scalar area(VertexHandle _vh) const { 
    return mesh_.property(area_, _vh);
  }

  Scalar gauss_curvature(VertexHandle _vh) const {
    return mesh_.property(gauss_curvature_, _vh);
  }

  Scalar mean_curvature(VertexHandle _vh) const { 
    return mesh_.property(mean_curvature_, _vh);
  }

  Scalar min_curvature(VertexHandle _vh) const { 
    const Scalar zero(0.0);
    Scalar H = mean_curvature(_vh);
    Scalar K = gauss_curvature(_vh);
    return H - sqrt(std::max(zero, H*H-K));
  }
  
  Scalar max_curvature(VertexHandle _vh) const { 
    const Scalar zero(0.0);
    Scalar H = mean_curvature(_vh);
    Scalar K = gauss_curvature(_vh);
    return H + sqrt(std::max(zero, H*H-K));
  }


private:

  Mesh&  mesh_;

  OpenMesh::VPropHandleT<Scalar>   area_;
  OpenMesh::VPropHandleT<Scalar>   gauss_curvature_;
  OpenMesh::VPropHandleT<Scalar>   mean_curvature_;
  OpenMesh::EPropHandleT<Scalar>   edge_weight_;

  bool weights_computed_, area_computed_;
};


//=============================================================================
} // namespace Remeshing
//=============================================================================
#if defined(INCLUDE_TEMPLATES) && !defined(DIFFGEOT_C)
#define DIFFGEO_TEMPLATES
#include "DiffGeoT.cc"
#endif
//=============================================================================
#endif // DIFFGEO_HH defined
//=============================================================================

