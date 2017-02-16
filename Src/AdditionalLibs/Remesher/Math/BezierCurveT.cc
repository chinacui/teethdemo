/*===========================================================================*\
 *                                                                           *
 *                              OpenFlipper                                  *
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
 *                                                                           *
\*===========================================================================*/

/*===========================================================================*\
 *                                                                           *
 *   $Revision$                                                       *
 *   $Author$                                                      *
 *   $Date$                   *
 *                                                                           *
\*===========================================================================*/



//=============================================================================
//
//  CLASS BezierCurve - IMPLEMENTATION
//
//=============================================================================

#define ACG_BEZIERCURVE_C

//== INCLUDES =================================================================


#include "BezierCurveT.hh"


//== IMPLEMENTATION ========================================================== 


namespace ACG {


//-----------------------------------------------------------------------------


template <class Point>
Point
BezierCurveT<Point>::
operator()(Scalar _t) const
{
  // copy controll points
  std::vector<Point> b(*this);

  unsigned int n = b.size()-1, k;
  Scalar t0(1.0-_t), t1(_t);


  // de Casteljau
  unsigned int i, j;
  for (i=0; i<n; ++i)
    for (j=0, k=n-i; j<k; ++j)
      b[j] = t0*b[j] + t1*b[j+1];


  return b[0];
}


//-----------------------------------------------------------------------------


template <class Point>
void
BezierCurveT<Point>::
subdivide(Scalar _t, Self& _curve0, Self& _curve1) const
{
  // copy controll points
  std::vector<Point> b(*this);

  int n = degree();
  Scalar t0(1.0-_t), t1(_t);


  _curve0.clear();
  _curve0.reserve(n+1);
  _curve1.clear();
  _curve1.reserve(n+1);

  std::vector<Point>  tmp;
  tmp.reserve(n+1);


  // de Casteljau
  int i, j, k;
  for (i=0; i<n; ++i)
  {
    _curve0.push_back(b[0]);
    tmp.push_back(b[n-i]);

    for (j=0, k=n-i; j<k; ++j)
      b[j] = t0*b[j] + t1*b[j+1];
  }

  _curve0.push_back(b[0]);
  tmp.push_back(b[0]);


  for (i=n; i>=0; --i)
    _curve1.push_back(tmp[i]);


  assert(_curve0.degree() == n);
  assert(_curve1.degree() == n);
}


//=============================================================================
} // namespace ACG
//=============================================================================
