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


#include "Triangulator.hh"


#include <iostream>



namespace ACG {


Triangulator::Triangulator(const std::vector<Vec3f>& _pos)
  : polySize_(_pos.size()), numRemaningVertices_(_pos.size()), numTris_(0), 
  numReflexVertices_(0),
  status_(-1), convex_(false)
{
  if (polySize_ < 3)
    return;


  if (polySize_ == 3)
  {
    numTris_ = 1;
    tris_.resize(3);
    tris_[0] = 0;
    tris_[1] = 1;
    tris_[2] = 2;

    numRemaningVertices_ = 0;
    convex_ = true;
    status_ = 0;
  }
  else
  {
    // project vertices onto the 2d plane of the polygon.
    // the projection plane is chosen orthogonal to the polygon surface normal

    // use Newell's Method to compute the surface normal
    // http://www.opengl.org/wiki/Calculating_a_Surface_Normal

    Vec3f n(0.0f, 0.0f, 0.0f);

    for (int i = 0; i < polySize_; ++i)
    {
      int next = (i + 1) % polySize_;

      Vec3f a = _pos[i] - _pos[next];
      Vec3f b = _pos[i] + _pos[next];

      n[0] += a[1] * b[2];
      n[1] += a[2] * b[0];
      n[2] += a[0] * b[1];
    }

    // project to 2d
    pos_.resize(polySize_);

    Vec3f axis[3] = { Vec3f(1.0f, 0.0f, 0.0f), Vec3f(0.0f, 1.0f, 0.0f), n };

    // orthonormalize projection axes
    axis[2].normalize();

    // make sure first axis is linearly independent from the normal
    while (std::abs(axis[0] | axis[2]) > 0.95f || (axis[0].sqrnorm() < 0.001f))
    {
      for (int i = 0; i < 3; ++i)
        axis[0][i] = float(rand()) / float(RAND_MAX) * 2.0f - 1.0f;

      axis[0].normalize();
    }

    // make axis[0] orthogonal to normal
    axis[0] = axis[0] - axis[2] * (axis[0] | axis[2]);
    axis[0].normalize();
    axis[1] = axis[2] % axis[0];


    for (int i = 0; i < polySize_; ++i)
    {
      // project onto polygon plane
      pos_[i][0] = axis[0] | _pos[i];
      pos_[i][1] = axis[1] | _pos[i];
    }


    // create triangle fans if there is at most one concave vertex
    int reflexVertexID = 0;

    for (int i = 0; i < polySize_; ++i)
    {
      // test vertex (i+1)
      if (isReflexVertex(pos_[i], pos_[(i + 1) % polySize_], pos_[(i + 2) % polySize_]))
      {
        ++numReflexVertices_;
        reflexVertexID = (i + 1) % polySize_;
      }
    }

    convex_ = !numReflexVertices_;


    if (numReflexVertices_ <= 1)
    {
      // create triangle fans
      numTris_ = polySize_ - 2;
      tris_.resize(numTris_ * 3);
      numRemaningVertices_ = 0;
      status_ = 0;

      for (int i = 0; i < numTris_; ++i)
      {
        tris_[i * 3] = reflexVertexID;
        tris_[i * 3 + 1] = (reflexVertexID + i + 1) % polySize_;
        tris_[i * 3 + 2] = (reflexVertexID + i + 2) % polySize_;
      }

    }
    else
    {
      // use the ear clipping algorithm

      earClippingN2();
//      earClippingN3();

//      triangulateExternal();
    }
  }
}


Triangulator::~Triangulator()
{
}


bool Triangulator::isReflexVertex(const Vec2f& v0, const Vec2f& v1, const Vec2f& v2) const
{
  // compute the sign of the cross product of the edges sharing vertex v1
  // <0 : inner angle greater than 180 deg  (reflex)
  // >0 : inner angle less than 180 deg (convex)

  Vec2f u = v2 - v1;
  Vec2f v = v0 - v1;

  return u[0] * v[1] - u[1] * v[0] < 0.0f;
}

bool Triangulator::isReflexVertex(int i) const
{
  int p = (i + polySize_ - 1) % polySize_;
  int n = (i + 1) % polySize_;

  return isReflexVertex(pos_[p], pos_[i], pos_[n]);
}

float Triangulator::triangleAreaSign(const Vec2f& v0, const Vec2f& v1, const Vec2f& v2) const
{
  // cross product
  return (v0[0] - v2[0]) * (v1[1] - v2[1]) - (v1[0] - v2[0]) * (v0[1] - v2[1]);
}

float Triangulator::distancePointToSegmentSq(const Vec2f& v0, const Vec2f& v1, const Vec2f& pt) const
{
  Vec2f segment = v1 - v0;
  Vec2f vec = pt - v0;
  float dp = vec | segment;

  if (dp < 0.0f)
    return vec.sqrnorm();
  
  float segSq = segment.sqrnorm();
  dp /= segSq;

  if (dp < 1.0f)
    return vec.sqrnorm() - dp * dp * segSq;

  vec = pt - v1;
  return vec.sqrnorm();
}

bool Triangulator::pointInTriangle(const Vec2f& v0, const Vec2f& v1, const Vec2f& v2, const Vec2f& pt) const
{
  // ACG implementation based on barycentric coordinates (slow)
//  return Geometry::isInTriangle(pt, v0, v1, v2);


  // fast implementation based on triangle areas
  // http://www.gamedev.net/topic/295943-is-this-a-better-point-in-triangle-test-2d/

  return triangleAreaSign(pt, v0, v1) >= 0.0f &&
    triangleAreaSign(pt, v1, v2) >= 0.0f &&
    triangleAreaSign(pt, v2, v0) >= 0.0f;


//   // more accurate algorithm:
//   // http://totologic.blogspot.de/2014/01/accurate-point-in-triangle-test.html
//     // note: didn't improve accuracy at all for problematic polygons
// 
//   static const float eps = 1e-4f;
//   static const float eps2 = eps*eps;
// 
// 
//   // point in aabb of triangle
//   Vec2f aabbMin = v0;
// 
//   aabbMin.minimize(v1);
//   aabbMin.minimize(v2);
// 
//   if (pt[0] + eps < aabbMin[0] || pt[1] + eps < aabbMin[1])
//     return false;
// 
//   Vec2f aabbMax = v0;
//   aabbMax.maximize(v1);
//   aabbMax.maximize(v2);
// 
//   if (pt[0] > aabbMax[0] + eps || pt[1] > aabbMax[1] + eps)
//     return false;
//   
// 
//   if (triangleAreaSign(pt, v0, v1) >= 0.0f &&
//     triangleAreaSign(pt, v1, v2) >= 0.0f &&
//     triangleAreaSign(pt, v2, v0) >= 0.0f)
//     return true;
// 
//   
//   return (distancePointToSegmentSq(v0, v1, pt) <= eps2 ||
//     distancePointToSegmentSq(v1, v2, pt) <= eps2 ||
//     distancePointToSegmentSq(v2, v0, pt) <= eps2);
}


void Triangulator::initVertexList()
{
  vertices_.resize(polySize_);
  reflexVertices_.clear();
  for (int i = 0; i < polySize_; ++i)
  {
    int p = (i + polySize_ - 1) % polySize_;
    int n = (i + 1) % polySize_;

    vertices_[i] = RingVertex(i, isReflexVertex(pos_[p], pos_[i], pos_[n]), pos_[i], &vertices_[p], &vertices_[n]);

    if (vertices_[i].reflex)
      reflexVertices_.push_back(&vertices_[i]);
  }
}



int Triangulator::earClippingN3()
{
  // http://www.geometrictools.com/Documentation/TriangulationByEarClipping.pdf
  // O(n^3)

  numTris_ = 0;
  status_ = 0;

  initVertexList();

  tris_.resize((polySize_ - 2) * 3);

  int numIterations = 0;
  RingVertex* firstVertex = &vertices_[0];

  while (numTris_ < polySize_ - 2)
  {
    // find an ear in the remaining polygon
    bool hasEars = false;

    RingVertex* curVertex = firstVertex;
    do 
    {
      curVertex->reflex = isReflexVertex(curVertex->prev->pos, curVertex->pos, curVertex->next->pos);

      if (!curVertex->reflex)
      {
        // triangle containment test
        bool isEar = true;

        // check all remaining vertices r for containment
        for (RingVertex* r = curVertex->next->next; r != curVertex->prev; r = r->next)
        {
          if (pointInTriangle(curVertex->prev->pos, curVertex->pos, curVertex->next->pos, r->pos))
          {
            isEar = false;
            break;
          }
        }

        // found an ear
        if (isEar)
        {
          // triangulate ear
          hasEars = true;
          tris_[numTris_ * 3] = curVertex->prev->id;
          tris_[numTris_ * 3 + 1] = curVertex->id;
          tris_[numTris_ * 3 + 2] = curVertex->next->id;
          ++numTris_;

          // remove vertex from linked list
          curVertex->prev->next = curVertex->next;
          curVertex->next->prev = curVertex->prev;

          break;
        }
      }

      curVertex = curVertex->next;
      ++numIterations;
    } while (curVertex != firstVertex);

    firstVertex = firstVertex->next;

    // create triangle fans and hope for good result if there are no more ears
    if (!hasEars && (numTris_ + 2 < polySize_))
    {
      for (RingVertex* iteratorVertex = firstVertex->next; iteratorVertex != firstVertex->prev; iteratorVertex = iteratorVertex->next)
      {
        tris_[numTris_ * 3] = firstVertex->id;
        tris_[numTris_ * 3 + 1] = iteratorVertex->id;
        tris_[numTris_ * 3 + 2] = iteratorVertex->next->id;
        ++numTris_;
      }

      assert(numTris_ == polySize_ - 2);
    }
  }

  return numTris_;
}



int Triangulator::earClippingN2()
{
  // http://www.geometrictools.com/Documentation/TriangulationByEarClipping.pdf
  // O(n^2)

  numTris_ = 0;
  status_ = 0;

  initVertexList();

  // triangulate

  int numTries = 0; // # checked vertices per iteration that aren't ears
  numRemaningVertices_ = polySize_; // size of currently remaining polygon
  int numIterations = 0;

  tris_.resize((polySize_ - 2) * 3);

  RingVertex* curVertex = &vertices_[0];
  while (numRemaningVertices_ > 3)
  {
    // check if the current vertex is an ear tip
    bool isEar = false;

    if (!curVertex->reflex)
    {
      // test current vertex for ear property
      isEar = true;
      for (std::list<RingVertex*>::iterator it = reflexVertices_.begin(); isEar && it != reflexVertices_.end(); ++it)
      {
        // skip direct neighbors
        if (*it == curVertex->prev || *it == curVertex->next)
          continue;

        // if any remaining vertex is inside the triangle, the current vertex is not an ear
        if (pointInTriangle(curVertex->prev->pos, curVertex->pos, curVertex->next->pos, (*it)->pos))
          isEar = false;
      }


      // found an ear
      if (isEar)
      {
        addEar(curVertex);
        numTries = 0;
      }
    }

    if (!isEar)
      ++numTries;

    if (numTries > numRemaningVertices_)
    {
      // something went wrong
      // create a triangle anyway and hope the result is ok

      addEar(curVertex);
      numTries = 0;

      status_ = 1;
    }

    curVertex = curVertex->next;
    ++numIterations;
  }


  // add the last remaining triangle
  if (numRemaningVertices_ == 3)
  {
    tris_[numTris_ * 3 + 0] = curVertex->prev->id;
    tris_[numTris_ * 3 + 1] = curVertex->id;
    tris_[numTris_ * 3 + 2] = curVertex->next->id;
    ++numTris_;
  }


  return numTris_;
}

bool Triangulator::updateReflexVertex(RingVertex* v)
{
  if (v->reflex)
  {
    // check reflex property
    v->reflex = isReflexVertex(v->prev->pos, v->pos, v->next->pos);

    // update list of reflex vertices
    if (!v->reflex)
      reflexVertices_.remove(v);
  }

  return v->reflex;
}

void Triangulator::addEar(RingVertex* _earTip)
{
  // add ear triangle
  tris_[numTris_ * 3 + 0] = _earTip->prev->id;
  tris_[numTris_ * 3 + 1] = _earTip->id;
  tris_[numTris_ * 3 + 2] = _earTip->next->id;

  // remove ear vertex from the linked list
  _earTip->prev->next = _earTip->next;
  _earTip->next->prev = _earTip->prev;

  // update reflex vertices list by checking the neighboring vertices
  updateReflexVertex(_earTip->prev);
  updateReflexVertex(_earTip->next);

  --numRemaningVertices_;
  ++numTris_;
}



}