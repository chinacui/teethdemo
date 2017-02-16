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
//  CLASS Geometry - IMPLEMENTATION
//
//=============================================================================

#define ACG_BSPLINEBASIS_C

//== INCLUDES =================================================================

#include "BSplineBasis.hh"

//----------------------------------------------------------------------------


namespace ACG {


//== IMPLEMENTATION ========================================================== 



template<typename Scalar>
Vec2i 
bsplineSpan(Scalar _t, 
  int _degree,
  const std::vector<Scalar>& _knots)
{
  // binary search

  int lo = _degree;
  int hi = _knots.size() - 1 - _degree;

  Scalar upperBound = _knots[hi];

  if (_t >= upperBound)
    return Vec2i(hi-1-_degree, hi - 1);

  int mid = (lo + hi) >> 1;

  while (_t < _knots[mid] || _t >= _knots[mid + 1])
  {
    if (_t < _knots[mid])
      hi = mid;
    else
      lo = mid;

    mid = (lo + hi) >> 1;
  }

  return Vec2i(mid - _degree, mid);


  
  // linear search:
//   
//   int i = 0;
// 
//   if (_t >= upperBound)
//     i = _knots.size() - _degree - 2;
//   else
//   {
//     while (_t >= _knots[i]) i++;
//     while (_t < _knots[i])  i--;
//   }
// 
//   return Vec2i(i - _degree, i);
}


template<typename Scalar>
void
bsplineBasisFunctions( std::vector<Scalar>& _N,
    const Vec2i& _span,
    Scalar _t,
    const std::vector<Scalar>& _knots)
{
  // inverted triangular scheme
  // "The NURBS Book" : Algorithm A2.2 p70


  // degree
  int p = _span[1] - _span[0];

  int i = _span[1];


  _N[0] = 1.0;

  // alloc temp buffer
  static std::vector<Scalar> left(p+1);
  static std::vector<Scalar> right(p+1);

  if (left.size() < size_t(p+1))
  {
    left.resize(p+1);
    right.resize(p+1);
  }

  // compute
  for (int j = 1; j <= p; ++j)
  {
    left[j] = _t - _knots[i + 1 - j];
    right[j] = _knots[i + j] - _t;

    Scalar saved = 0.0;

    for (int r = 0; r < j; ++r)
    {
      Scalar tmp = _N[r] / (right[r + 1] + left[j - r]);
      _N[r] = saved + right[r + 1] * tmp;
      saved = left[j - r] * tmp;
    }
    _N[j] = saved;
  }
}


template<typename Scalar>
void 
bsplineBasisDerivatives( std::vector<Scalar>& _ders,
    const Vec2i& _span,
    Scalar _t, 
    int _der, 
    const std::vector<Scalar>& _knots,
    std::vector<Scalar>* _functionVals )
{
  // The NURBS Book  p72  algorithm A2.3

  int p = _span[1] - _span[0];
  int p1 = p+1;
  int i = _span[1];


  // alloc temp arrays
  static std::vector< std::vector<Scalar> > ndu(p1);
  static std::vector<Scalar> left(p1);
  static std::vector<Scalar> right(p1);
  static std::vector<Scalar> a(2*p1);

  if (ndu[0].size() < size_t(p1))
  {
    ndu.resize(p1);
    for (int i = 0; i < p1; ++i)
      ndu[i].resize(p1);

    left.resize(p1);
    right.resize(p1);

    a.resize(2*p1);
  }




  ndu[0][0] = 1.0;

  for (int j = 1; j <= p; ++j)
  {
    left[j]= _t - _knots[i + 1 - j];
    right[j]= _knots[i + j] - _t;
    Scalar saved = 0.0;

    for (int r = 0; r < j; ++r)
    {
      // lower triangle
      ndu[j][r] = right[r+1] + left[j-r];
      Scalar tmp = ndu[r][j-1] / ndu[j][r];

      ndu[r][j] = saved + right[r+1] * tmp;
      saved = left[j-r] * tmp;
    }
    ndu[j][j] = saved;
  }

  // load the basis functions
  if (_functionVals)
  {
    for (int j = 0; j <= p; ++j)
      (*_functionVals)[j] = ndu[j][p];
  }

  // compute derivatives



  for (int r = 0; r <= p; ++r)
  {
    int s1 = 0, s2 = p1; // s1, s2: row offsets of linearized 2d array a[2][p+1]
    a[0] = 1.0;

    for (int k = 1; k <= _der; ++k)
    {
      Scalar d = 0.0;
      int rk = r - k, pk = p - k;

      if (r >= k)
      {
        a[s2 + 0] = a[s1 + 0] / ndu[pk+1][rk];
        d = a[s2] * ndu[rk][pk];
      }

      int j1 = (rk >= -1) ? 1 : -rk;
      int j2 = (r - 1 <= pk) ? k-1 : p-r;

      for (int j = j1; j <= j2; ++j)
      {
        a[s2 + j] = (a[s1 + j] - a[s1 + j-1]) / ndu[pk+1][rk+j];
        d += a[s2 + j] * ndu[rk+j][pk];
      }

      if (r <= pk)
      {
        a[s2 + k] = -a[s1 + k-1] / ndu[pk+1][r];
        d += a[s2 + k] * ndu[r][pk];
      }

      if (k == _der)
        _ders[r] = d;

      std::swap(s1, s2); // switch rows
    }
  }

  // correct factors

  int r = p;

  for (int k = 1; k <= _der; ++k)
  {
    Scalar rf = Scalar(r);
    for (int j = 0; j <= p; ++j)
      _ders[j] *= rf;

    r *= (p - k);
  }
}



template<typename Scalar>
Scalar 
bsplineBasisFunction(int _i,
  int _degree, 
  Scalar _t, 
  const std::vector<Scalar>& _knots)
{
  int m = _knots.size() - 1;

  // Mansfield Cox deBoor recursion
  if ((_i == 0 && _t == _knots[0]) || (_i == m - _degree - 1 && _t == _knots[m]))
    return 1.0;

  if (_degree == 0) {
    if (_t >= _knots[_i] && _t < _knots[_i + 1])
      return 1.0;
    else
      return 0.0;
  }

  double Nin1 = basisFunction(_i, _degree - 1, _t, _knots);
  double Nin2 = basisFunction(_i + 1, _degree - 1, _t, _knots);

  double fac1 = 0;
  //   if ((_knotvector(_i+_n) - _knotvector(_i)) > 0.000001 )
  if ((_knots[_i + _degree] - _knots[_i]) != 0)
    fac1 = (_t - _knots[_i]) / (_knots[_i + _degree] - _knots[_i]);

  double fac2 = 0;
  //   if ( (_knotvector(_i+1+_n) - _knotvector(_i+1)) > 0.000001 )
  if ((_knots[_i + 1 + _degree] - _knots[_i + 1]) != 0)
    fac2 = (_knots[_i + 1 + _degree] - _t) / (_knots[_i + 1 + _degree] - _knots[_i + 1]);

  //   std::cout << "Nin1 = " << Nin1 << ", Nin2 = " << Nin2 << ", fac1 = " << fac1 << ", fac2 = " << fac2 << std::endl;

  return (fac1*Nin1 + fac2*Nin2);
}


template<typename Scalar>
Scalar
bsplineBasisDerivative(int _i,
  int _degree,
  Scalar _t, 
  int _der, 
  const std::vector<Scalar>& _knots)
{
  assert(_degree >= 0);
  assert(_i >= 0);


  if (_der == 0)
    return basisFunction(_i, _degree, _t, _knots);

  Scalar Nin1 = derivativeBasisFunction(_i, _degree - 1, _t, _der - 1, _knots);
  Scalar Nin2 = derivativeBasisFunction(_i + 1, _degree - 1, _t, _der - 1, _knots);

  Scalar fac1 = 0;
  if (std::abs(_knots[_i + _degree] - _knots[_i]) > 1e-6)
    fac1 = Scalar(_degree) / (_knots[_i + _degree] - _knots[_i]);

  Scalar fac2 = 0;
  if (std::abs(_knots[_i + _degree + 1] - _knots[_i + 1]) > 1e-6)
    fac2 = Scalar(_degree) / (_knots[_i + _degree + 1] - _knots[_i + 1]);

  return (fac1*Nin1 - fac2*Nin2);
}

//=============================================================================
} // namespace ACG
//=============================================================================
