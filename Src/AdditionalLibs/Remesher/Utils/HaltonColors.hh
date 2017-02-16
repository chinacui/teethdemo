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


#ifndef ACG_HALTONCOLORS_HH
#define ACG_HALTONCOLORS_HH

//== INCLUDES =================================================================


#include <ACG/Math/VectorT.hh>
#include "../Config/ACGDefines.hh"

#include "ColorGenerator.hh"

//== NAMESPACES ===============================================================

namespace ACG {

  //== CLASS DEFINITION =========================================================

/** \brief Implementation of halton colors
 Provides deterministic pseudo-random low-discrepancy colors with a
 uniform distribution over a visually pleasing subset of HSL space,
 independent of the number of colors required.
 Simply instantiate and use get_next_color().
*/
class ACGDLLEXPORT HaltonColors : public ColorGenerator {

public:

  /// Default constructor
  HaltonColors(int skip = 250);

  /// Generate the next color (legacy method)
  ACG::Vec4f get_next_color();

  /// Generate the next color
  virtual ACG::Vec4f generateNextColor();

private:

  float halton(int index);
  float random_interval(int index, float min, float max);
  ACG::Vec4f HSL2RGB(double h, double sl, double l);

  int current[3]; // current Halton index
  int bases[3];   // Halton prime bases
  float inverse_bases[3];
};


  //=============================================================================
} 
//=============================================================================
#endif // ACG_HALTONCOLORS_HH defined
//=============================================================================

