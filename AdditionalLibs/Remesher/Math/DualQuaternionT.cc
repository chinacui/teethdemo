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
//  CLASS DualQuaternionT - IMPLEMENTATION
//
//=============================================================================


#define ACG_DUALQUATERNIONT_C


//== INCLUDES =================================================================


#include "DualQuaternionT.hh"
#include <iostream>

//== IMPLEMENTATION ========================================================== 


namespace ACG {

  //-----------------------------------------------------------------------------
  
  /// Default constructor ( constructs an identity dual quaternion )
  template <typename Scalar>
  DualQuaternionT<Scalar>::DualQuaternionT(){
    *this = DualQuaternion::identity();
  }
  
  //-----------------------------------------------------------------------------
  
  /// Copy constructor
  template <typename Scalar>
  DualQuaternionT<Scalar>::DualQuaternionT(const DualQuaternion& _other){
    real_ = _other.real_;
    dual_ = _other.dual_;
  }
  
  //-----------------------------------------------------------------------------
  
  /// Construct from given real,dual parts
  template <typename Scalar>
  DualQuaternionT<Scalar>::DualQuaternionT(const Quaternion& _real, const Quaternion& _dual){
    real_ = _real;
    dual_ = _dual;
  }
  
  //-----------------------------------------------------------------------------
  
  /// Construct from 8 scalars
  template <typename Scalar>
  DualQuaternionT<Scalar>::DualQuaternionT(Scalar _Rw, Scalar _Rx, Scalar _Ry, Scalar _Rz,
                                           Scalar _Dw, Scalar _Dx, Scalar _Dy, Scalar _Dz){
    real_ = Quaternion(_Rw, _Rx, _Ry, _Rz);
    dual_ = Quaternion(_Dw, _Dx, _Dy, _Dz);
  }
  
  //-----------------------------------------------------------------------------

  /// Construct from a rotation given as quaternion
  template <typename Scalar>
  DualQuaternionT<Scalar>::DualQuaternionT(Quaternion _rotation){
    real_ = _rotation;
    dual_ = Quaternion(0.0,0.0,0.0,0.0);
  }
  
  //-----------------------------------------------------------------------------
  
  /// Construct from a translatation given as a vector
  template <typename Scalar>
  DualQuaternionT<Scalar>::DualQuaternionT(const Vec3& _translation){
    real_.identity();
    dual_ = Quaternion( 0.0, 0.5 * _translation[0], 0.5 * _translation[1], 0.5 * _translation[2] );
  }
  
  //-----------------------------------------------------------------------------
  
  /// Construct from a translation+rotation
  template <typename Scalar>
  DualQuaternionT<Scalar>::DualQuaternionT(const Vec3& _translation, const Quaternion& _rotation){

    real_ = _rotation;
    dual_ = Quaternion( 0.0, 0.5 * _translation[0], 0.5 * _translation[1], 0.5 * _translation[2] );
    
    dual_ *= real_;
  }
  
  //-----------------------------------------------------------------------------
  
  /// Construct from a rigid transformation given as matrix
  template <typename Scalar>
  DualQuaternionT<Scalar>::DualQuaternionT(const Matrix& _transformation){
    real_ = Quaternion(_transformation); //the quaternion constructor ignores the translation
    dual_ = Quaternion( 0.0, 0.5 * _transformation(0,3), 0.5 * _transformation(1,3), 0.5 * _transformation(2,3) );
    
    dual_ *= real_;
  }
  
  //-----------------------------------------------------------------------------
  
  /// identity dual quaternion [ R(1, 0, 0, 0), D(0,0,0,0) ]
  template <typename Scalar>
  DualQuaternionT<Scalar> DualQuaternionT<Scalar>::identity(){
  
    Quaternion real;
    real.identity();
    
    Quaternion dual = Quaternion(0.0, 0.0, 0.0, 0.0);
    
    return DualQuaternion( real, dual );
  }

  //-----------------------------------------------------------------------------

  /// zero dual quaternion [ R(0, 0, 0, 0), D(0,0,0,0) ]
  template <typename Scalar>
  DualQuaternionT<Scalar> DualQuaternionT<Scalar>::zero(){
  
    Quaternion real = Quaternion(0.0, 0.0, 0.0, 0.0);
    Quaternion dual = Quaternion(0.0, 0.0, 0.0, 0.0);
    
    return DualQuaternion( real, dual );
  }
  
  //-----------------------------------------------------------------------------

  /// conjugate dual quaternion
  template <typename Scalar>
  DualQuaternionT<Scalar> DualQuaternionT<Scalar>::conjugate() const {
    return DualQuaternion( real_.conjugate(), dual_.conjugate() );
  }

  //-----------------------------------------------------------------------------
  
  /// invert dual quaternion
  template <typename Scalar>
  DualQuaternionT<Scalar> DualQuaternionT<Scalar>::invert() const {

    double sqrLen0 = real_.sqrnorm();
    double sqrLenE = 2.0 * (real_ | dual_);

    if ( sqrLen0 > 0.0 ){
      
      double invSqrLen0 = 1.0/sqrLen0;
      double invSqrLenE = -sqrLenE/(sqrLen0*sqrLen0);

      DualQuaternion conj = conjugate();
      
      conj.real_ = invSqrLen0 * conj.real_;
      conj.dual_ = invSqrLen0 * conj.dual_ + invSqrLenE * conj.real_;
      
      return conj;
      
    } else
      return DualQuaternion::zero();
  }
  
  //-----------------------------------------------------------------------------
  
  /// normalize dual quaternion
  template <typename Scalar>
  void DualQuaternionT<Scalar>::normalize() {
    
    const double magn = 1.0/real_.norm();
    const double magnSqr = 1.0/real_.sqrnorm();
    
    // normalize rotation
    real_ *= magn;
    dual_ *= magn;

    // normalize the rest
    dual_ -= ((real_| dual_)* magnSqr) * real_;

  }
  
  //-----------------------------------------------------------------------------
  
  /// dual quaternion comparison
  template <typename Scalar>
  bool DualQuaternionT<Scalar>::operator==(const DualQuaternion& _other) const {
    return (_other.real_ == real_) && (_other.dual_ == dual_);
  }
  
  //-----------------------------------------------------------------------------
  
  /// dual quaternion comparison
  template <typename Scalar>
  bool DualQuaternionT<Scalar>::operator!=(const DualQuaternion& _other) const {
    return (_other.real_ != real_) || (_other.dual_ != dual_);
  }
  
  //-----------------------------------------------------------------------------
  
  /// addition
  template <typename Scalar>
  DualQuaternionT<Scalar> DualQuaternionT<Scalar>::operator+(const DualQuaternion& _other) const {
    return DualQuaternion( real_ + _other.real_, dual_ + _other.dual_ );
  }
  
  //-----------------------------------------------------------------------------
  
  /// addition
  template <typename Scalar>
  DualQuaternionT<Scalar>& DualQuaternionT<Scalar>::operator+=(const DualQuaternion& _other) {
    real_ = real_ + _other.real_;
    dual_ = dual_ + _other.dual_;
    
    return (*this);
  }
  
  //-----------------------------------------------------------------------------
  
  /// substraction 
  template <typename Scalar>
  DualQuaternionT<Scalar> DualQuaternionT<Scalar>::operator-(const DualQuaternion& _other) const {
    return DualQuaternion( real_ - _other.real_, dual_ - _other.dual_ );
  }

  //-----------------------------------------------------------------------------

  /// substraction
  template <typename Scalar>
  DualQuaternionT<Scalar>& DualQuaternionT<Scalar>::operator-=(const DualQuaternion& _other) {
    real_ -= _other.real_;
    dual_ -= _other.dual_;
    
    return (*this);
  }
  
  //-----------------------------------------------------------------------------
  
  /// dualQuaternion * dualQuaternion
  template <typename Scalar>
  DualQuaternionT<Scalar> DualQuaternionT<Scalar>::operator*(const DualQuaternion& _q) const {
    return DualQuaternion( real_ * _q.real_, real_ * _q.dual_ + dual_ * _q.real_ );
  }

  //-----------------------------------------------------------------------------
  
  /// dualQuaternion *= dualQuaternion
  template <typename Scalar>
  DualQuaternionT<Scalar>& DualQuaternionT<Scalar>::operator*=(const DualQuaternion& _q){
    dual_  = real_ * _q.dual_ + dual_ * _q.real_;
    real_ *= _q.real_;
    
    return (*this);
  }
  
  //-----------------------------------------------------------------------------
  
  /// dualQuaternion * scalar
  template <typename Scalar>
  DualQuaternionT<Scalar> DualQuaternionT<Scalar>::operator*(const Scalar& _scalar) const {
    DualQuaternion q;
    
    q.real_ = real_ * _scalar;
    q.dual_ = dual_ * _scalar;
    
    return q;
  }

  //-----------------------------------------------------------------------------

  template <typename Scalar>
  Scalar& DualQuaternionT<Scalar>::operator [](const unsigned int& b) {
    if ( b < 4 ) {
      return real_[b];
    } else if ( b < 8 ) {
      return dual_[b - 4];
    } else {
      // Invalid operation, write error and return anything.
      std::cerr << "Error in Dualquaternion operator[], index b out of range [0...7]" << std::endl;
      return real_[0];
    }
  }


  //-----------------------------------------------------------------------------

  /// linear interpolation of dual quaternions. Result is normalized afterwards.
  template <typename Scalar> template<typename VectorType>
  DualQuaternionT<Scalar> DualQuaternionT<Scalar>::interpolate(VectorType& _weights, const std::vector<DualQuaternion>& _dualQuaternions)
  {
    if ( (_weights.size() != _dualQuaternions.size()) || (_weights.size() == 0) ){
      std::cerr << "Cannot interpolate dual quaternions ( weights: " << _weights.size() << ", DQs: " << _dualQuaternions.size() << std::endl;
      return DualQuaternionT::zero();
    }

    // Find max weight for pivoting to that quaternion,
    // so shortest arc is taken (see: 'coping antipodality' in the paper )
    uint pivotID = 0;
    
    for (uint i=1; i<_dualQuaternions.size(); i++)
      if (_weights[pivotID] < _weights[i])
        pivotID = i;

    DualQuaternion pivotDQ = _dualQuaternions[ pivotID ];
    Quaternion pivotQ = pivotDQ.real_;

    DualQuaternion res = DualQuaternion::zero();
    
    for (uint i=0; i<_dualQuaternions.size(); i++){

      DualQuaternion currentDQ = _dualQuaternions[i];
      Quaternion currentQ = currentDQ.real_;
      
      // Make sure dot product is >= 0
      if ( ( currentQ | pivotQ ) < 0.0 )
        _weights[i] = -_weights[i];

      res += _dualQuaternions[i] * _weights[i];
    }

    res.normalize();

    return res;
  }

  //-----------------------------------------------------------------------------

  /// transform a given point with this dual quaternion
  template <typename Scalar>
  VectorT<Scalar,3> DualQuaternionT<Scalar>::transform_point(const Vec3& _point) const
  {
    ///TODO check if this is a unit dual quaternion

    // for details about the calculation see the paper (algorithm 1)

    Vec3 p(_point);
    
    double r  = real_[0];
    Vec3   rv = Vec3(real_[1], real_[2], real_[3]);

    double d  = dual_[0];
    Vec3   dv = Vec3(dual_[1], dual_[2], dual_[3]);
    
    Vec3 tempVec = (rv % p) + r * p;
    p+=2.0*(rv % tempVec);

    Vec3 t = dv % rv;
    t+= d * rv;
    t+= -r*dv;
    p+=-2.0*t;

    return p;
  }
  
  //-----------------------------------------------------------------------------

  /// transform a given point with this dual quaternion
  template <typename Scalar>
  VectorT<Scalar,3> DualQuaternionT<Scalar>::transform_vector(const Vec3& _point) const
  {
    ///TODO check if this is a unit dual quaternion

    // for details about the calculation see the paper (algorithm 1)

    Vec3 p(_point);
    
    double r  = real_[0];
    Vec3   rv = Vec3(real_[1], real_[2], real_[3]);

    Vec3 tempVec = (rv % p) + r * p;
    p+=2.0*(rv % tempVec);

    return p;
}

  //-----------------------------------------------------------------------------
  
  /// print some info about the DQ
  template <typename Scalar>
  void DualQuaternionT<Scalar>::printInfo(){
    
    std::cerr << "Real Part:" << std::endl;
    real_.print_info();
    std::cerr << "Dual Part:" << std::endl;
    dual_.print_info();
  }

  //-----------------------------------------------------------------------------
  
//=============================================================================
} // namespace ACG
//=============================================================================
