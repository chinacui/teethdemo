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
 *   $Revision: 9595 $                                                       *
 *   $Author: wilden $                                                       *
 *   $Date: 2010-06-17 12:48:23 +0200 (Thu, 17 Jun 2010) $                   *
 *                                                                           *
\*===========================================================================*/



//=============================================================================
//
//  CLASS DualQuaternion
//
//=============================================================================

#ifndef ACG_DUALQUATERNION_HH
#define ACG_DUALQUATERNION_HH


//== INCLUDES =================================================================

#include "QuaternionT.hh"


//== NAMESPACES  ==============================================================

namespace ACG {


//== CLASS DEFINITION =========================================================


/**
    \brief DualQuaternion class for representing rigid motions in 3d

    This is an implementation of:

    techreport{kavan-06-dual,
    author = "Ladislav Kavan and Steven Collins and Carol O'Sullivan and Jiri Zara",
    series = "Technical report TCD-CS-2006-46, Trinity College Dublin",
    title = "{D}ual {Q}uaternions for {R}igid {T}ransformation {B}lending",
    url = "http://www.cgg.cvut.cz/~kavanl1/",
    year = "2006"
    }
*/

template <class Scalar>
class DualQuaternionT
{
public:

  typedef QuaternionT<Scalar>      Quaternion;
  typedef DualQuaternionT<Scalar>  DualQuaternion;
  typedef VectorT<Scalar,3>        Vec3;
  typedef VectorT<Scalar,4>        Vec4;
  typedef Matrix4x4T<Scalar>       Matrix;

  /// real and dual quaternion parts
  Quaternion real_;
  Quaternion dual_;
  
  
  // Constructors
  //
  
  /// Default constructor ( constructs an identity dual quaternion )
  DualQuaternionT();
  
  /// Copy constructor
  DualQuaternionT(const DualQuaternion& _other);
  
  /// Construct from given real,dual parts
  DualQuaternionT(const Quaternion& _real, const Quaternion& _dual);
  
  /// Construct from 8 scalars
  DualQuaternionT(Scalar _Rw, Scalar _Rx, Scalar _Ry, Scalar _Rz,
                  Scalar _Dw, Scalar _Dx, Scalar _Dy, Scalar _Dz);

  /// Construct from a rotation given as quaternion
  DualQuaternionT(Quaternion _rotation);
  
  /// Construct from a translatation given as a vector
  DualQuaternionT(const Vec3& _translation);
  
  /// Construct from a translation+rotation
  DualQuaternionT(const Vec3& _translation, const Quaternion& _rotation);
  
  /// Construct from a rigid transformation given as matrix
  DualQuaternionT(const Matrix& _transformation);
  
  // default quaternions
  
  /// identity dual quaternion [ R(1, 0, 0, 0), D(0,0,0,0) ]
  static DualQuaternion identity();

  /// zero dual quaternion [ R(0, 0, 0, 0), D(0,0,0,0) ]
  static DualQuaternion zero();
  
  // Operators
  //
  
  /// conjugate dual quaternion
  DualQuaternion conjugate() const;

  /// invert dual quaternion
  DualQuaternion invert() const;
  
  /// normalize dual quaternion
  void normalize();
  
  /// dual quaternion comparison
  bool operator==(const DualQuaternion& _other) const;
  bool operator!=(const DualQuaternion& _other) const;
  
  /// addition
  DualQuaternion operator+(const DualQuaternion& _other) const;
  DualQuaternion& operator+=(const DualQuaternion& _other);
  
  /// substraction 
  DualQuaternion operator-(const DualQuaternion& _other) const;
  DualQuaternion& operator-=(const DualQuaternion& _other);

  /// dualQuaternion * dualQuaternion
  DualQuaternion operator*(const DualQuaternion& _q) const;

  /// dualQuaternion * scalar
  DualQuaternion operator*(const Scalar& _scalar) const;
  
  /// dualQuaternion *= dualQuaternion
  DualQuaternion& operator*=(const DualQuaternion& _q);

  /// Access as one big vector
  Scalar& operator [](const unsigned int& b);

  /// linear interpolation of dual quaternions. Result is normalized afterwards
  template <typename VectorType>
  static DualQuaternion interpolate(VectorType& _weights, const std::vector<DualQuaternion>& _dualQuaternions);
  
  /// Transform a point with the dual quaternion
  Vec3 transform_point(const Vec3& _point) const;
  
  /// Transform a vector with the dual quaternion
  Vec3 transform_vector(const Vec3& _point) const;
  
  /// print some info about the DQ
  void printInfo();

};


typedef DualQuaternionT<float>  DualQuaternionf;
typedef DualQuaternionT<double> DualQuaterniond;


//=============================================================================
} // namespace ACG
//=============================================================================
#if defined(INCLUDE_TEMPLATES) && !defined(ACG_DUALQUATERNIONT_C)
#define ACG_QUATERNIONT_TEMPLATES
#include "DualQuaternionT.cc"
#endif
//=============================================================================
#endif // ACG_DUALQUATERNION_HH defined
//=============================================================================

