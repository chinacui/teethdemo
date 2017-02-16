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

#ifndef HUEPARTITIONINGCOLORS_HH_
#define HUEPARTITIONINGCOLORS_HH_

#include <ACG/Math/VectorT.hh>
#include "../Config/ACGDefines.hh"
#include "ColorGenerator.hh"

namespace ACG {

/**
 * The HuePartitioningColors generator tries to generate a set of well
 * distinguishable and esthetically somewhat pleasing colors.
 *
 * Note that it intentionally behaves totally deterministic. (Because
 * reproducibility rocks.)
 */
class ACGDLLEXPORT HuePartitioningColors : public ColorGenerator {
    public:
        /**
         * Constructor
         *
         * @param _alpha The alpha value for all the colors.
         *
         * @param _baseHue The HSV-hue from which to start. This is the
         * hue of the first requested color. Default is an utterly
         * delighting shade of blue awesomeness.
         */
        HuePartitioningColors(
                float _alpha = 1.0f,
                float _baseHue = 0.5694f);

        /**
         * @return A new color.
         */
        virtual Vec4f generateNextColor();

        /**
         * Convenience method if you just need a bunch of
         * colors and don't need to instantiate a ColorGenerator.
         *
         * See description of generateNextNColors() for details.
         */
        template <class OUTPUT_ITERATOR>
        static void generateNColors(int n, OUTPUT_ITERATOR oit) {
            HuePartitioningColors cg;
            cg.generateNextNColors(n, oit);
        }

    private:
        int currentSubdiv_;
        int currentIt_;
        int currentTriadIt_;
        float alpha_;
        float baseHue_;

        const float defaultSaturation_, defaultValue_;
};

} /* namespace ACG */
#endif /* HUEPARTITIONINGCOLORS_HH_ */
