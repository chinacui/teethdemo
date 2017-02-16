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
//
//=============================================================================


#ifndef ALGORITHMSANGLE_HH
#define ALGORITHMSANGLE_HH


/*! \file AlgorithmsAngleT.hh
    \brief Functions for geometric operations related to angles
    
    General file with template functions handling angles
*/

//== INCLUDES =================================================================

//== FORWARDDECLARATIONS ======================================================

//== NAMESPACES ===============================================================

/// Namespace providing different geometric functions concerning angles
namespace ACG {
namespace Geometry {
   
//== CLASS DEFINITION =========================================================   
                       
/** Return a fully parametrized angle
   @param _vector1 vector pointing away from center, angle = 0
   @param _vector2 vector for which the angle should be calculated                    
   @param _normal  the normal used to compute if vector2 is left or right of vecor1
   @param _skip    Flag to catch nan. If true nan occurred and value should not be used
   @return the angle (0 .. 2 * PI)
*/
template < typename VectorT , typename ValueT >
ValueT
getFullangle( VectorT _vector1, 
              VectorT _vector2, 
              const   VectorT& _normal, 
              bool&   _skip );


/** Calculate the difference between two angles ( minimum distance )
   @param angle0 angle1
   @param angle1 angle2
   @return The difference between the angles (0..PI)
*/
template < typename ValueT >
inline
ValueT
angleDist( const ValueT& angle0 , 
           const ValueT& angle1 );
                  
/** Calculate the correct 2D angle if cos and sin of the angle are given
    This function calculates based on the signs of the acos and asin of the 
    given angles, in which quadrant the angle is and returns the full angle
    in radians
   @param _cos cos of angle
   @param _sin sin of angle
   @return angle
*/
template < typename ValueT >
inline
ValueT
getAngle( const ValueT& _cos , 
          const ValueT& _sin );      

/** Convert angle from radians to degree
   @param _angle in radians
   @return angle in degree
*/
template < typename ValueT >
inline
ValueT
radToDeg( const ValueT& _angle );      

/** Convert angle from degree to radians
   @param _angle in degree 
   @return angle in radians
*/
template < typename ValueT >
inline
ValueT
degToRad( const ValueT& _angle ); 


//=============================================================================
} // Geometry Namespace 
} // ACG Namespace 
//=============================================================================
#if defined(INCLUDE_TEMPLATES) && !defined(ALGORITHMSANGLE_C)
#define ALGORITHMSANGLE_TEMPLATES
#include "AlgorithmsAngleT.cc"
#endif
//=============================================================================
#endif // ALGORITHMSANGLE_HH defined
//=============================================================================

