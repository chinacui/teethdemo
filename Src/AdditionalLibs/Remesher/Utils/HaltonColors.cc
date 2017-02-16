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
 *   $Date$                   *
 *                                                                           *
\*===========================================================================*/

//=============================================================================
//
//  CLASS HaltonColors (by Marcel Campen)
//
//=============================================================================


//== INCLUDES =================================================================


#include <ACG/Utils/HaltonColors.hh>


//== NAMESPACES ===============================================================

namespace ACG {

//== CLASS IMPLEMENTATION =========================================================


HaltonColors::HaltonColors(int skip)
{
  // skip first 250 sequence elements to lower discrepancy even further.
  current[0] = skip;
  current[1] = skip;
  current[2] = skip;

  // initialize prime bases for H,S,L. Especially the first should be small such that already
  // small numbers of generated colors are distributed over the whole color circle.
  bases[0] = 5;
  bases[1] = 13;
  bases[2] = 17;

  inverse_bases[0] = 1.0f / bases[0];
  inverse_bases[1] = 1.0f / bases[1];
  inverse_bases[2] = 1.0f / bases[2];
}

float HaltonColors::halton(int index)
{
  int base = bases[index];
  float inverse_base = inverse_bases[index];
  float H = 0;
  float half = inverse_base;
  int I = current[index];
  current[index] += 1;
  while (I > 0) {
    int digit = I % base;
    H = H + half * digit;
    I = (int)(inverse_base * (I - digit));
    half *= inverse_base;
  }
  return H;
}

float HaltonColors::random_interval(int index, float min, float max)
{
  return halton(index) * (max - min) + min;
}

ACG::Vec4f HaltonColors::HSL2RGB(double h, double sl, double l)
{
  double v;
  double r, g, b;

  r = l;
  g = l;
  b = l;

  v = (l <= 0.5) ? (l * (1.0 + sl)) : (l + sl - l * sl);

  if (v > 0) {
    double m;
    double sv;
    int sextant;
    double fract, vsf, mid1, mid2;

    m = l + l - v;
    sv = (v - m) / v;
    h *= 6.0;
    sextant = (int) h;
    fract = h - sextant;
    vsf = v * sv * fract;
    mid1 = m + vsf;
    mid2 = v - vsf;

    switch (sextant) {
      case 0:
        r = v;
        g = mid1;
        b = m;
        break;
      case 1:
        r = mid2;
        g = v;
        b = m;
        break;
      case 2:
        r = m;
        g = v;
        b = mid1;
        break;
      case 3:
        r = m;
        g = mid2;
        b = v;
        break;
      case 4:
        r = mid1;
        g = m;
        b = v;
        break;
      case 5:
        r = v;
        g = m;
        b = mid2;
        break;
    }
  }

  return Vec4f((float)r, (float)g, (float)b, 1.0f);
}

ACG::Vec4f HaltonColors::generateNextColor() {
    float h = random_interval(0, 0.0f , 0.9f ); // 0.9 instead of 1.0 to suppress natural bias towards red
    float s = random_interval(1, 0.40f, 0.80f); // saturation between 40% and 80%
    float l = random_interval(2, 0.30f, 0.60f); // lightness between 30% and 60%
    return HSL2RGB(h, s, l);
}

ACG::Vec4f HaltonColors::get_next_color()
{
    return generateNextColor();
}


//=============================================================================
} // ACG namespace end
//=============================================================================
//=============================================================================

