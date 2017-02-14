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
//  CLASS ColorCoder
//
//=============================================================================

//== INCLUDES =================================================================

#include <ACG/Utils/ColorCoder.hh>

//== NAMESPACES ===============================================================

namespace ACG {

//== CLASS DEFINITION =========================================================


  /// Default constructor.
ColorCoder::ColorCoder(float _min, float _max, bool _signed)
{
  set_range(_min, _max, _signed);
}

/// set the color coding range for unsigned coding
void ColorCoder::set_range(float _min, float _max, bool _signed)
{
  if (_min == _max) {
    val0_ = val1_ = val2_ = val3_ = val4_ = _min;
  } else {
    if (_min > _max)
      std::swap(_min, _max);
    val0_ = _min;
    val4_ = _max;
    val2_ = 0.5 * (val0_ + val4_);
    val1_ = 0.5 * (val0_ + val2_);
    val3_ = 0.5 * (val2_ + val4_);
  }
  signed_mode_ = _signed;
}

/// color coding
ACG::Vec4uc ColorCoder::color4(float _v) const
{
  return signed_mode_ ? color_signed(_v) : color_unsigned(_v);
}

/// color coding
ACG::Vec3uc ColorCoder::color(float _v) const
{
  ACG::Vec4uc c;
  if (signed_mode_)
    c = color_signed(_v);
  else
    c = color_unsigned(_v);
  return (ACG::Vec3uc(c[0], c[1], c[2]) / 255.f);
}

/// color coding
ACG::Vec3f ColorCoder::color_float(float _v) const
{
  ACG::Vec4uc c;
  if (signed_mode_)
    c = color_signed(_v);
  else
    c = color_unsigned(_v);
  return (ACG::Vec3f(c[0], c[1], c[2]) / 255.f);
}

/// color coding
ACG::Vec4f ColorCoder::color_float4(float _v) const
{
  ACG::Vec4uc c;
  if (signed_mode_)
    c = color_signed(_v);
  else
    c = color_unsigned(_v);
  return (ACG::Vec4f(c[0], c[1], c[2], c[3]) / 255.f);
}

/// color coding
QColor ColorCoder::color_qcolor(float _v) const
{
  ACG::Vec4uc c;
  if (signed_mode_)
    c = color_signed(_v);
  else
    c = color_unsigned(_v);

  return(QColor(c[0], c[1], c[2], c[3]));
}

/// min scalar value
float ColorCoder::min() const
{
  return val0_;
}
/// max scalar value
float ColorCoder::max() const
{
  return val4_;
}

ACG::Vec4uc ColorCoder::color_unsigned(float _v) const
{
  if (val4_ <= val0_)
    return ACG::Vec4uc(0, 0, 255, 255);

  unsigned char u;

  if (_v < val0_)
    return ACG::Vec4uc(0, 0, 255, 255);
  if (_v > val4_)
    return ACG::Vec4uc(255, 0, 0, 255);

  if (_v <= val2_) {
    // [v0, v1]
    if (_v <= val1_) {
      u = (unsigned char) (255.0 * (_v - val0_) / (val1_ - val0_));
      return ACG::Vec4uc(0, u, 255, 255);
    }
    // ]v1, v2]
    else {
      u = (unsigned char) (255.0 * (_v - val1_) / (val2_ - val1_));
      return ACG::Vec4uc(0, 255, 255 - u, 255);
    }
  } else {
    // ]v2, v3]
    if (_v <= val3_) {
      u = (unsigned char) (255.0 * (_v - val2_) / (val3_ - val2_));
      return ACG::Vec4uc(u, 255, 0, 255);
    }
    // ]v3, v4]
    else {
      u = (unsigned char) (255.0 * (_v - val3_) / (val4_ - val3_));
      return ACG::Vec4uc(255, 255 - u, 0, 255);
    }
  }
}

ACG::Vec4uc ColorCoder::color_signed(float _v) const
{
  if (val4_ <= val0_)
    return ACG::Vec4uc(0, 255, 0, 255);

  unsigned char r, g, b;

  if (_v < val0_)
    _v = val0_;
  else if (_v > val4_)
    _v = val4_;

  if (_v < 0.0) {
    r = val0_ ? (unsigned char) (255.0 * _v / val0_) : 0;
    b = 0;
  } else {
    r = 0;
    b = val4_ ? (unsigned char) (255.0 * _v / val4_) : 0;
  }
  g = 255 - r - b;

  return ACG::Vec4uc(r, g, b, 255);
}


//=============================================================================
}// namespace ACG
//=============================================================================
//=============================================================================

