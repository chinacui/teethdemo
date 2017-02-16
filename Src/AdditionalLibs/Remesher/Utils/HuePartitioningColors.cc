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

#include "HuePartitioningColors.hh"

#include <QColor>

namespace ACG {

#define DEFAULT_SATURATION 1.0f
#define DEFAULT_VALUE      1.0f

HuePartitioningColors::HuePartitioningColors(float _alpha, float _baseHue) :
        currentSubdiv_(2), currentIt_(0), currentTriadIt_(0),
        alpha_(_alpha), baseHue_(_baseHue),
        defaultSaturation_(DEFAULT_SATURATION), defaultValue_(DEFAULT_VALUE)  {
}

static inline float wrap01(float v) {
    float dummy;
    if (v > 1.0f) v = modff(v, &dummy);
    if (v < 0.0f) v = 1.0f - modff(v, &dummy);
    return v;
}

Vec4f HuePartitioningColors::generateNextColor() {
    const float resultHue =
            baseHue_
            + .33333333f / currentSubdiv_ * currentIt_
            + .33333333f * currentTriadIt_;

    // Convert color to RGB and store result.
    double r, g, b, a;
    QColor::fromHsvF(wrap01(resultHue), defaultSaturation_,
            defaultValue_, alpha_).getRgbF(
            &r, &g, &b, &a);
    const Vec4f result(r, g, b, alpha_);

    /*
     * Increment the three level iterator state.
     */

    // Level 1: Increment triad it.
    if (++currentTriadIt_ <= 2)
        return result;

    // Level 2: Increment subdiv it.
    currentTriadIt_ = 0;

    // Starting from the second subdivison, only visit odd numbers.
    currentIt_ += (currentSubdiv_ == 2 ? 1 : 2);
    if (currentIt_ < currentSubdiv_)
        return result;

    // Level 3: Increment subdivision level.
    currentIt_ = 1;
    currentSubdiv_ *= 2;

    return result;
}

} /* namespace ACG */
