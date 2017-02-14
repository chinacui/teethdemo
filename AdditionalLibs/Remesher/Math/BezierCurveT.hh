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
//  CLASS BezierCurveT
//
//=============================================================================


#ifndef ACG_BEZIERCURVE_HH
#define ACG_BEZIERCURVE_HH


//== INCLUDES =================================================================


#include "VectorT.hh"
#include <vector>


//== NAMESPACES  ==============================================================


namespace ACG {


//== CLASS DEFINITION =========================================================


/** Bezier curve. Derived from std::vector<Vector<Scalar, Dimension>>.
 */

template <class Point>
class BezierCurveT : public std::vector<Point>
{
public:

  typedef typename Point::value_type  Scalar;
  typedef BezierCurveT<Point>         Self;
  typedef std::vector<Point>          Base;


  /// constructor
  BezierCurveT() {}

  /// destructor
  ~BezierCurveT() {}


  /// return degree (= size()-1)
  unsigned int degree() const { return Base::size()-1; }


  /// evaluate curve at parameter _t using deCasteljau
  Point operator()(Scalar _t) const;


  /** subdivide curve at parameter _t, store the two resulting
      curves in _curve0, _curve1
  */
  void subdivide(Scalar _t, Self& _curve0, Self& _curve1) const;

};


//=============================================================================
} // namespace ACG
//=============================================================================
#if defined(INCLUDE_TEMPLATES) && !defined(ACG_BEZIERCURVE_C)
#define ACG_BEZIERCURVE_TEMPLATES
#include "BezierCurveT.cc"
#endif
//=============================================================================
#endif // ACG_BEZIERCURVE_HH defined
//=============================================================================

