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
 *   $Revision$                                                      *
 *   $Author$                                                      *
 *   $Date$                    *
 *                                                                           *
\*===========================================================================*/


#ifndef ACG_BSPLINEBASIS_HH_
#define ACG_BSPLINEBASIS_HH_

#include <vector>
#include <ACG/Math/VectorT.hh>

namespace ACG {


/** \brief Find the span of a parameter value.
  *
  * @param _t parameter value
  * @param _degree spline degree
  * @param _knots knotvector
*/
template<typename Scalar>
Vec2i
bsplineSpan(Scalar _t,
  int _degree,
  const std::vector<Scalar>& _knots);


/** \brief Evaluate basis functions in a span
  *
  * Uses a fast algorithm to compute all basis functions in a span.
  * The output vector _N must have enough space to hold _span[1] - _span[0] values.
  *
  * @param _N output basis function values
  * @param _span span of the parameter _t
  * @param _t parameter value
  * @param _knots knotvector
*/
template<typename Scalar>
void
bsplineBasisFunctions( std::vector<Scalar>& _N,
  const Vec2i& _span,
  Scalar _t,
  const std::vector<Scalar>& _knots);



/** \brief Compute derivatives of basis functions in a span
  *
  * Uses a fast algorithm to compute all derivatives of basis functions in a span.
  * The output vector _ders must have enough space to hold _span[1] - _span[0] values.
  *
  * @param _ders output derivative values of all basis functions in the span
  * @param _span span of the parameter _t
  * @param _t parameter value
  * @param _der order of derivative (ie. first, second, third order...)
  * @param _knots knotvector
  * @param _functionVals output basis function values, null is accepted (optional)
  */
template<typename Scalar>
void
bsplineBasisDerivatives( std::vector<Scalar>& _ders,
  const Vec2i& _span,
  Scalar _t,
  int _der,
  const std::vector<Scalar>& _knots,
  std::vector<Scalar>* _functionVals);


/** \brief Evaluate a single basis function
  *
  * Uses the recursive definition of the basis function.
  * This is slower than using bsplineBasisFunctions() for all values in a span.
  *
  * @param _i index of basis function
  * @param _degree spline degree
  * @param _t parameter value
  * @param _knots knotvector
  * @return basis function value
*/
template<typename Scalar>
Scalar
bsplineBasisFunction(int _i,
  int _degree,
  Scalar _t,
  const std::vector<Scalar>& _knots);

/** \brief Compute derivative of a single basis function
  *
  * Uses the recursive definition of the basis function.
  * This is slower than using bsplineBasisDerivatives() for all values in a span.
  *
  * @param _i index of basis function
  * @param _degree spline degree
  * @param _t parameter value
  * @param _der order of derivative
  * @param _knots knotvector
  * @return derivative of basis function value
*/
template<typename Scalar>
Scalar
bsplineBasisDerivative(int _i,
  int _degree,
  Scalar _t,
  int _der,
  const std::vector<Scalar>& _knots);



} /* namespace ACG */

#if defined(INCLUDE_TEMPLATES) && !defined(ACG_BSPLINEBASIS_C)
#define ACG_BSPLINEBASIS_TEMPLATES
#include "BSplineBasis.cc"
#endif

#endif /* ACG_BSPLINEBASIS_HH_ */