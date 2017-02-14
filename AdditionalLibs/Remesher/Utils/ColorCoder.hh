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

#ifndef ACG_COLORCODER_HH
#define ACG_COLORCODER_HH

//== INCLUDES =================================================================

#include <ACG/Math/VectorT.hh>
#include <ACG/Config/ACGDefines.hh>
#include <QColor>

//== NAMESPACES ===============================================================

namespace ACG {

//== CLASS DEFINITION =========================================================

/** \brief Class for generating nice colors for doubles
 *
 *
 */
class ACGDLLEXPORT ColorCoder {
public:

  /// Default constructor.
  ColorCoder(float _min = 0.0, float _max = 1.0, bool _signed = false);

  /// set the color coding range for unsigned coding
  void set_range(float _min, float _max, bool _signed);

  /// color coding
  ACG::Vec4uc color4(float _v) const;

  /// color coding
  ACG::Vec3uc color(float _v) const;

  /// color coding
  ACG::Vec3f color_float(float _v) const;

  /// color coding
  ACG::Vec4f color_float4(float _v) const;

  /// color coding
  QColor color_qcolor(float _v) const;

  /// min scalar value
  float min() const;

  /// max scalar value
  float max() const;

  // Make the color coder usable as a function operator.
  inline ACG::Vec4f operator() (float _v) const {
      return color_float4(_v);
  }

private:

  ACG::Vec4uc color_unsigned(float _v) const;
  ACG::Vec4uc color_signed(float _v) const;

  float val0_, val1_, val2_, val3_, val4_;
  bool signed_mode_;
};

//=============================================================================
}// namespace ACG
//=============================================================================
#endif // ACG_COLORCODER_HH defined
//=============================================================================

