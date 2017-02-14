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
 *   $Revision$                                                         *
 *   $Author$                                                      *
 *   $Date$                   *
 *                                                                           *
\*===========================================================================*/




//=============================================================================
//
//  IMPLEMENTATION
//
//=============================================================================

#define ALGORITHMSANGLE_C

//== INCLUDES =================================================================

#include "AlgorithmsAngleT.hh"

#include <OpenMesh/Core/Geometry/MathDefs.hh>

#include <cmath>
#include <iostream>


//== NAMESPACES ===============================================================

namespace ACG {
namespace Geometry {

//== IMPLEMENTATION ==========================================================

/// Return false if x is not a number
inline bool isNan(double x) {
   return (x != x);
}


template < typename VectorT , typename ValueT >
ValueT
getFullangle( VectorT _vector1 , VectorT _vector2 , const VectorT& _normal , bool& _skip )
{
   //Project vectors into tangent plane defined by _normal
   _vector1 = _vector1 - _normal * ( _vector1 | _normal );
   _vector2 = _vector2 - _normal * ( _vector2 | _normal );
   _vector1.normalize();
   _vector2.normalize();

   //calculate projection onto right Vector (used to decide if vector2 is left or right of vector1
   const double right = ( ( _normal % _vector1 ) | _vector2 ) ;

   double sp    = ( _vector1 | _vector2 );
   
   //Catch some errors with scalar product and the following acos
   if (sp < -1.0) {
     sp = -1.0;
   }

   if (sp > 1.0) {
      sp = 1.0;
   }

   double angle = acos(sp);

   // catch some possible nans
   _skip = ( isNan(right) || isNan(angle) ) ;

   if ( right < 0 ) {
      angle = 2 * M_PI - angle;
   }

   return angle;
}

template < typename ValueT >
inline
ValueT
angleDist( const ValueT& angle0 , const ValueT& angle1 ) {
  ValueT dist = fabs( angle1 - angle0 );
  return ( std::min( dist , 2 * M_PI - dist) );
}

template < typename ValueT >
inline
ValueT
getAngle( const ValueT& _cos ,
                  const ValueT& _sin )
{
   const double angle_asin   = asin( OpenMesh::sane_aarg(_sin) );
   const double angle_acos  = acos( OpenMesh::sane_aarg(_cos) );

   if ( angle_asin >= 0 ) { //Quadrant 1,2
      if ( angle_acos >= 0 ) { // Quadrant 1
          return angle_asin;
      } else { //Quadrant 2
         return (M_PI - angle_asin);
      }
   } else {  //Quadrant 3,4
      if ( angle_acos >= 0 ) { // Quadrant 4
         return (2 * M_PI + angle_asin);
      } else { //Quadrant 3
         return (M_PI - angle_asin);
      }
   }
}

template < typename ValueT >
inline
ValueT
radToDeg( const ValueT& _angle ) {
  return ( _angle / M_PI * 180);
}

template < typename ValueT >
inline
ValueT
degToRad( const ValueT& _angle ) {
   return ( _angle / 180 * M_PI );
}



//=============================================================================
} // Geometry Namespace
} // ACG Namespace
//=============================================================================
