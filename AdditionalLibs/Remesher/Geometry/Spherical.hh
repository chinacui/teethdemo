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


/*
 * Spherical.hh
 *
 *  Created on: Sep 18, 2012
 *      Author: ebke
 */

#ifndef ACG_GEOMETRY_SPHERICAL_HH_
#define ACG_GEOMETRY_SPHERICAL_HH_

#include <ACG/Math/GLMatrixT.hh>
#include <cmath>
#include <stdexcept>

namespace ACG {
namespace Geometry {

static const double epsilon = 1e-6;

namespace Spherical_Private {

template<class Vec>
static inline typename Vec::value_type angleBetween(const Vec &v1, const Vec &v2, const Vec &normal) {
    typedef typename Vec::value_type Scalar;

    /*
     * Handle unstable 0 degree special case.
     * (0 degrees can easily become 360 degrees.)
     */
    const Scalar dotProd = v1 | v2;
    if (std::abs(dotProd - 1.0) < epsilon) return 0;

    /*
     * General case: Use determinant to check >/< 180 degree.
     */
    const Scalar det = GLMatrixT<Scalar>(normal, v1, v2).determinant();

    /*
     * Interpret arc cos accrodingly.
     */
    const Scalar arcos_angle = std::acos(std::max(-1.0, std::min(1.0, dotProd)));

    if (det >= -1e-6)
        return arcos_angle;
    else
        return 2 * M_PI - arcos_angle;
}

template<class Vec>
static inline typename Vec::value_type angleBetween(const Vec &v1, const Vec &v2) {
    const Vec normal = (v1 % v2).normalized();
    return angleBetween(v1, v2, normal);
}

} /* namespace Spherical_Private */

/**
 * Compute inner angle sum of the spherical triangle specified by
 * the three supplied unit vectors.
 *
 * @param n0 Must be unit length.
 * @param n1 Must be unit length.
 * @param n2 Must be unit length.
 *
 * @return Inner angle sum of the specified spherical triangle.
 */
template<class Vec>
static inline typename Vec::value_type sphericalInnerAngleSum(const Vec &n0, const Vec &n1, const Vec &n2) {
    typedef typename Vec::value_type Scalar;

    const Scalar a = Spherical_Private::angleBetween(n0, n1);
    const Scalar b = Spherical_Private::angleBetween(n1, n2);
    const Scalar c = Spherical_Private::angleBetween(n2, n0);
    if (a < epsilon || b < epsilon || c < epsilon) return M_PI;

    const Scalar s = .5 * (a + b + c);
    const Scalar sin_s = std::sin(s);
    const Scalar sin_a = std::sin(a);
    const Scalar sin_b = std::sin(b);
    const Scalar sin_c = std::sin(c);

#ifndef NDEBUG
    if (std::sin(s - a) < -1e-4) throw std::logic_error("ACG::Geometry::sphericalInnerAngleSum(): 0 > s-a");
    if (std::sin(s - b) < -1e-4) throw std::logic_error("ACG::Geometry::sphericalInnerAngleSum(): 0 > s-b");
    if (std::sin(s - c) < -1e-4) throw std::logic_error("ACG::Geometry::sphericalInnerAngleSum(): 0 > s-c");
#endif

    const Scalar alpha_2 = std::acos(std::min(1.0, std::sqrt(sin_s * std::max(.0, std::sin(s - a)) / (sin_b * sin_c))));
    const Scalar beta_2  = std::acos(std::min(1.0, std::sqrt(sin_s * std::max(.0, std::sin(s - b)) / (sin_c * sin_a))));
    const Scalar gamma_2 = std::acos(std::min(1.0, std::sqrt(sin_s * std::max(.0, std::sin(s - c)) / (sin_a * sin_b))));

    return 2 * (alpha_2 + beta_2 + gamma_2);
}

/**
 * Compute gauss curvature of spherical polyhedral spanned by the
 * supplied unit length normals. Normals have to be supplied in
 * CCW order.
 *
 * @return Gauss curvature of the supplied spherical polyhedral.
 */
template<class Vec, class INPUT_ITERATOR>
static inline typename Vec::value_type sphericalPolyhedralGaussCurv(INPUT_ITERATOR normals_begin, INPUT_ITERATOR normals_end) {
    typedef typename Vec::value_type Scalar;

    if (normals_begin == normals_end) return 0;
    Vec n0 = *(normals_begin++);

    if (normals_begin == normals_end) return 0;
    Vec n2 = *(normals_begin++);

    Scalar result = 0;
    while (normals_begin != normals_end) {
        /*
         * Next triangle.
         */
        const Vec n1 = n2;
        n2 = *(normals_begin++);

        const Scalar sign = GLMatrixT<Scalar>(n0, n1, n2).determinant() >= 0 ? 1 : -1;
        result += sign * (sphericalInnerAngleSum(n0, n1, n2) - M_PI);
    }

    return result;
}

} /* namespace Geometry */
} /* namespace ACG */

#endif /* ACG_GEOMETRY_SPHERICAL_HH_ */
