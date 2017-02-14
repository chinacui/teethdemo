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
*   $Revision: 21022 $                                                       *
*   $Author: moebius $                                                      *
*   $Date: 2015-07-17 08:23:03 +0200 (Fri, 17 Jul 2015) $                   *
*                                                                           *
\*===========================================================================*/


#ifndef ACG_TRIANGULATOR_HH
#define ACG_TRIANGULATOR_HH



#include <ACG/Math/VectorT.hh>
#include <ACG/Config/ACGDefines.hh>

#include <vector>
#include <list>


namespace ACG{



class ACGDLLEXPORT Triangulator
{
public:
  /** \brief Execute the triangulation algorithm on a polygon
   *
   * @param _pos polygon vertex positions (ccw)
  */
  Triangulator(const std::vector<Vec3f>& _pos);

  /** \brief Destructor
  */
  virtual ~Triangulator();

  /** \brief Get number of triangles
   * @return number of triangles after triangulation
  */
  int numTriangles() const { return numTris_; }

  /** \brief Get local vertex index 
   *
   * @param _i index to index buffer in range [0, .., numTris*3 - 1]
   * @return vertex index in range [0, .., _pos.size()-1], where _pos is the input vector of positions as passed to the constructor
  */
  int index(int _i) const { return tris_[_i]; }

  /** \brief Get local index buffer
  *
  * @return index buffer
  */
  const std::vector<int>& indices() const { return tris_; }

  /** \brief Is the polygon convex?
  *
  * @return convex true or false
  */
  bool convex() const { return convex_; }

  /** \brief Get number of reflex vertices
  *
  * A reflex vertex is a vertex with an inner angle larger than 180 deg.
  * @return number of reflex vertices
  */
  int numReflexVertices() const { return numReflexVertices_; }

  /** \brief Check if a vertex is reflex
  *
  * A reflex vertex is a vertex with an inner angle larger than 180 deg.
  * @param _i local vertex index in range [0, .., _pos.size()-1]
  * @return reflex vertex true or false
  */
  bool isReflexVertex(int _i) const;

  /** \brief Check if the triangulation was successful
  *
  * @return success true or false
  */
  bool success() const { return !status_; }


private:

  void initVertexList();

  // ear clipping algorithm in O(n^2)
  int earClippingN2();

  // ear clipping algorithm in O(n^3)
  int earClippingN3();


  // disable assignment operator
  Triangulator& operator=(const Triangulator&){ return *this;  }


  float triangleAreaSign(const Vec2f& v0, const Vec2f& v1, const Vec2f& v2) const;
  float distancePointToSegmentSq(const Vec2f& v0, const Vec2f& v1, const Vec2f& pt) const;

  // test if vertex v1 is reflex in triangle v0-v1-v2  (inner angle at v1 is greater than 180 deg)
  bool isReflexVertex(const Vec2f& v0, const Vec2f& v1, const Vec2f& v2) const;



  bool pointInTriangle(const Vec2f& v0, const Vec2f& v1, const Vec2f& v2, const Vec2f& pt) const;


  struct RingVertex
  {
    RingVertex() {}

    RingVertex(int i, bool r, const Vec2f& x, RingVertex* p, RingVertex* n)
      : id(i), reflex(r), pos(x), prev(p), next(n) {  }

    int id;
    bool reflex;
    Vec2f pos;

    RingVertex* prev;
    RingVertex* next;
  };


  bool updateReflexVertex(RingVertex* v);
  void addEar(RingVertex* _earTip);


  const int polySize_;
  int numRemaningVertices_;
  int numTris_;
  int numReflexVertices_;

  int status_;
  bool convex_;


  std::vector<Vec2f> pos_;

  std::vector<RingVertex> vertices_;
  std::list<RingVertex*> reflexVertices_;

  std::vector<int> tris_;

};




//=============================================================================
}


//=============================================================================
#endif // ACG_TRIANGULATOR_HH defined
//=============================================================================