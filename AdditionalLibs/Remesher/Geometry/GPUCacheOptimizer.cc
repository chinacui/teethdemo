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
 *   $Revision$															 *
 *   $Author$														 *
 *   $Date$													 *
 *                                                                           *
\*===========================================================================*/

//=============================================================================

#include "GPUCacheOptimizer.hh"
#include <cassert>
#include <cmath>
#include <vector>
#include <cstring>

//=============================================================================

namespace ACG
{

//=============================================================================

GPUCacheOptimizer::GPUCacheOptimizer( unsigned int NumTris, unsigned int NumVerts, unsigned int IndexSize, const void* pIndices) :
        m_NumVerts(NumVerts),
        m_NumTris(NumTris),
        m_IndexSize(IndexSize),
        m_pIndices(pIndices),
        m_NumTransformations(0)
{
  m_pTriMap = new unsigned int[m_NumTris];
}

GPUCacheOptimizer::~GPUCacheOptimizer(void)
{
	delete [] m_pTriMap;
}

//=============================================================================

unsigned int GPUCacheOptimizer::GetIndex(unsigned int i) const
{
	assert(i < m_NumTris * 3);

	return GetIndex(i, m_IndexSize, m_pIndices);
}

unsigned int GPUCacheOptimizer::GetIndex(unsigned int i, unsigned int IndexSize, const void* pIB)
{
	switch (IndexSize)
	{
	case 4: return ((const unsigned int*)pIB)[i]; break;
	case 2: return ((const unsigned short*)pIB)[i]; break;
	case 1: return ((const unsigned char*)pIB)[i]; break;
	default:
		assert(i == 1 || i == 2 || i == 4); // throw error
	}
	return 0xFFFFFFFF;
}

void GPUCacheOptimizer::SetIndex(unsigned int i, unsigned int val, unsigned int IndexSize, void* pIB)
{
	switch (IndexSize)
	{
	case 4: ((unsigned int*)pIB)[i] = val; break;
	case 2: ((unsigned short*)pIB)[i] = val; break;
	case 1: ((unsigned char*)pIB)[i] = val; break;
	default:
		assert(i == 1 || i == 2 || i == 4); // throw error
	}
}

//=============================================================================

void GPUCacheOptimizer::WriteIndexBuffer(unsigned int DstIndexSize, void* pDst)
{
	assert(DstIndexSize == 1 ||DstIndexSize == 2 || DstIndexSize == 4);
	// TODO: warning log, if DstIndexSize < m_IndexSize

	// support for 'in-place' operation via tmpbuf copy
	char* pSrc = (char*)m_pIndices;

	int bTmpCopy = 0;
	if (pDst == pSrc)
	{
		pSrc = new char[m_IndexSize * m_NumTris * 3];
		memcpy(pSrc, m_pIndices, m_IndexSize * m_NumTris * 3);

		bTmpCopy = 1;
	}

	for (unsigned int i = 0; i < m_NumTris; ++i)
	{
		for (int k = 0; k < 3; ++k)	
		{
			unsigned int TriVertex = GetIndex(m_pTriMap[i] * 3 + k, m_IndexSize, pSrc);

			// copy remapped tri indices
			SetIndex(i * 3 + k, TriVertex, DstIndexSize, pDst);
		}
	}

	if (bTmpCopy) delete [] pSrc;
}

//=============================================================================

void GPUCacheOptimizer::RemapVertices(unsigned int NumTris, unsigned int NumVerts,
									  const unsigned int* pVertMap,
									  unsigned int IndexSize, void* pInOutIndices,
									  unsigned int VertexStride, void* pInOutVertices)
{
	if (pVertMap && pInOutIndices && pInOutVertices && VertexStride)
	{
		// make tmp vertex buffer copy
		char* pTmpBuf = new char[VertexStride * NumVerts];
		memcpy(pTmpBuf, pInOutVertices, VertexStride * NumVerts);

		char* pVertexOut = (char*)pInOutVertices;

		// apply on vertex buffer

		for (unsigned int i = 0; i < NumVerts; ++i)
		{
			// some mapping destinations might be invalid
			//  this vertex is unused,  ignore then
			if (pVertMap[i] < NumVerts)
				memcpy(pVertexOut + pVertMap[i] * VertexStride,
					pTmpBuf + i * VertexStride, VertexStride);
		}

		// apply on index buffer

		for (unsigned int i = 0; i < NumTris * 3; ++i)
		{
			// IndexBuffer[i] = VertMap[IndexBuffer[i]]

			unsigned int v = GetIndex(i, IndexSize, pInOutIndices);
			SetIndex(i, pVertMap[v], IndexSize, pInOutIndices);
		}

		delete [] pTmpBuf;
	}
}

//=============================================================================

void GPUCacheOptimizer::OptimizeVertices(unsigned int NumTris, unsigned int NumVerts,
										 unsigned int IndexSize, const void* pIndices,
										 unsigned int* pVertMap)
{
	// straight forward algorithm
	// simply iterate over indices and increment vertex location if unvisited vertex found
	unsigned int uCounter = 0; // vertex counter

	memset(pVertMap, 0xFF, NumVerts * sizeof(unsigned int));

	for (unsigned int i = 0; i < NumTris * 3; ++i)
	{
		unsigned int vertex;

		if (IndexSize == 2) vertex = ((const unsigned short*)pIndices)[i];
		else vertex = ((const unsigned int*)pIndices)[i];

		if (pVertMap[vertex] == 0xFFFFFFFF)
			pVertMap[vertex] = uCounter++;
	}
}

//=============================================================================

unsigned int GPUCacheOptimizer::ComputeNumberOfVertexTransformations(unsigned int VertexCacheSize)
{
	if (m_NumTransformations) return m_NumTransformations;

	unsigned int NumIndices = 3 * m_NumTris;
	if (!NumIndices) return 0;

	unsigned int* Cache = new unsigned int[VertexCacheSize];
	unsigned int NumSlotsInUse = 0;
	unsigned int LastSlot = 0;
	m_NumTransformations = 0;

	for (unsigned int i = 0; i < m_NumTris; ++i)
	{
		unsigned int t = m_pTriMap[i];

		// for each vertex of triangle t:
		for (int k = 0; k < 3; ++k)
		{
			unsigned int Idx = GetIndex(t * 3 + k); // vertex index

			int bInCache = 0;
			for (unsigned int k = 0; k < NumSlotsInUse && !bInCache; ++k)
			{
				if (Cache[k] == Idx) bInCache = 1;
			}

			if (!bInCache)
			{
				++m_NumTransformations;
				if (NumSlotsInUse < VertexCacheSize)
				{
					Cache[NumSlotsInUse++] = Idx; 
					++LastSlot;
				}
				else
				{
					if (LastSlot == VertexCacheSize) LastSlot = 0;
					Cache[LastSlot++] = Idx;
				}
			}
		}
		
	}

	delete [] Cache;

	return m_NumTransformations;
}

//=============================================================================

float GPUCacheOptimizer::ComputeACMR(unsigned int VertexCacheSize)
{
	unsigned int NumT = ComputeNumberOfVertexTransformations(VertexCacheSize);
	return float(NumT) / float(m_NumTris);
}

//=============================================================================

float GPUCacheOptimizer::ComputeATVR(unsigned int VertexCacheSize)
{
	unsigned int NumT = ComputeNumberOfVertexTransformations(VertexCacheSize);
	return float(NumT) / float(m_NumVerts);
}

//=============================================================================

// forsyth's score function
void GPUCacheOptimizer::Opt_Vertex::FindScore(unsigned int MaxSizeVertexCache)
{



	float fNewScore = -1.0f; // -1 : vertex unused
	if (iNumTrisLeft > 0)
	{
		if (iCachePos < 0) fNewScore = 0.0f; // not in FIFO
		else
		{

			if (iCachePos < 3){ // last tri => fixed score

			  const float LastTriScore = 0.75f;
				fNewScore = LastTriScore;

			} else
			{
			  const float CacheDecayPower = 1.5f;

				// check for cache_pos < MaxSizeCachePos here..
				// Points for being high in the cache.
				const float Scaler = 1.0f / (MaxSizeVertexCache - 3);
				fNewScore = 1.0f - ( iCachePos - 3 ) * Scaler;
				fNewScore = powf( fNewScore, CacheDecayPower);
			}
		}

		// Bonus points for having a low number of tris still to
		// use the vert, so we get rid of lone verts quickly.

	  const float ValenceBoostScale = 2.0f;
	  const float ValenceBoostPower = 0.5f;

		float ValenceBoost = powf( float(iNumTrisLeft), -float(ValenceBoostPower));
		fNewScore += ValenceBoostScale * ValenceBoost;
	}

	fScore = fNewScore;
}

void GPUCacheOptimizer::Opt_Vertex::RemoveTriFromList(unsigned int tri)
{
	for (int k = 0; k < iNumTrisLeft; ++k)
	{
		// replace tri with last tri in this list
		if (pTris[k] == tri)
		{
			pTris[k] = pTris[iNumTrisLeft-1];
			break;
		}
	}
	--iNumTrisLeft;
}

//=============================================================================
// tipsify

GPUCacheOptimizerTipsify::GPUCacheOptimizerTipsify(unsigned int CacheSize, unsigned int NumTris, unsigned int NumVerts, unsigned int IndexSize, const void *pIndices)
: GPUCacheOptimizer(NumTris, NumVerts, IndexSize, pIndices)
{
  // optimization notes:
  //  - cache unfriendly layout of Opt_Vertex: high read/write access cost to Opt_vertex data

	if (NumVerts < 3 || !NumTris) return;

	Opt_Vertex* pVerts = new Opt_Vertex[NumVerts];
	Opt_Tris*   pTris  = new Opt_Tris[NumTris];

	// OPTIMIZATION: memset vs constructor initialization: memset - 400 ms,  constructor - 770 ms (at 10 mil vertices)
	memset(pVerts, 0, sizeof(Opt_Vertex) * NumVerts);

	// build adjacency, same start as in forsyth class
	for (unsigned int i = 0; i < NumTris; ++i)
	{
    // copy vertex indices of this tri
    Opt_Tris* pThisTri = pTris + i;

		// copy vertex indices of this tri
    for (unsigned int k = 0; k < 3; ++k)
    {

      const unsigned int idx = GetIndex(i * 3 + k);
      pThisTri->v[k] = idx;

      // count # tris per vertex
      ++pVerts[idx].iNumTrisTotal;
    }
	}

	// OPTIMIZATION: allocate one buffer for the complete vertex adjacency list (instead of allocating a new buffer for each vertex)
	const unsigned int vertexTriAdjBufSize   = NumTris * 3;
	unsigned int vertexTriAdjBufOffset       = 0;
	unsigned int* vertexTriAdjBuf            = new unsigned int[vertexTriAdjBufSize];

	// create list of tris per vertex
	for (unsigned int i = 0; i < NumTris; ++i)
	{
		// add this tri to per vertex tri list
		for (int k = 0; k < 3; ++k)
		{
			Opt_Vertex* pV = pVerts + pTris[i].v[k];
			if (!pV->pTris) {
			  pV->pTris = vertexTriAdjBuf + vertexTriAdjBufOffset;
			  vertexTriAdjBufOffset += pV->iNumTrisTotal;
			}

			// abuse <numTrisLeft> as temporal up counter 
			// (automatically sums to numTris, exactly what we want)
			pV->pTris[pV->iNumTrisLeft++] = i;

			pV->iCachePos = 0;
		}
	}

	// use the cache_pos of the OptFaces_Vertex as the time stamp

	//=============================================================================
	// OPTIMIZATION:
	//  push and pop on DeadEndVertexStack greatly increases processing time
	// -> replace with fixed size ring stack
	//	std::vector<unsigned int> DeadEndVertexStack;
	//	DeadEndVertexStack.reserve(2048);
	RingStack DeadEndVertexStack(128);


	int f          = 0; // arbitrary starting index (vertex)
	int iTimeStamp = CacheSize + 1;
	unsigned int i = 1; // cursor

	unsigned int numTrisAdded = 0;

	std::vector<unsigned int> N; // 1-ring of next candidates
	N.reserve(2048);

	while (f >= 0)
	{
		N.clear();

		// this vertex
		Opt_Vertex* pV = pVerts + f;

		// for each adjacent tri of this vertex
		for (int m = 0; m < pV->iNumTrisTotal; ++m)
		{
			Opt_Tris* pT = pTris + pV->pTris[m];

			if (!pT->bAdded)
			{
				// append
				m_pTriMap[numTrisAdded++] = pV->pTris[m];

				for (int k = 0; k < 3; ++k)
				{
					// push to cache
				  const unsigned int v = pT->v[k];
				  //         DeadEndVertexStack.push_back(pT->v[k]);
				  DeadEndVertexStack.push(v);

					// insert
				  N.push_back(v);

				  Opt_Vertex* adjV = pVerts + v;

				  --adjV->iNumTrisLeft;

				  if (iTimeStamp - adjV->iCachePos > (int)CacheSize)
				    adjV->iCachePos = iTimeStamp++;
				}
				pT->bAdded = 1;
			}
		}


		// select next fanning vertex
		// Get-Next-Vertex
		{
			int n = -1, p = -1; // best candidate and priority
			for (unsigned int k = 0; k < N.size(); ++k)
			{
				// for each vertex in N
				Opt_Vertex* pV = pVerts + N[k];

				if (pV->iNumTrisLeft > 0)
				{
					// error here in pseudo code:
					//  literal p should be named m here
					//  to find the best vertex
					int m = 0;
					if (iTimeStamp - pV->iCachePos + 2 * pV->iNumTrisLeft <= (int)CacheSize)
						m = iTimeStamp - pV->iCachePos;

					if (m > p)
					{
						p = m;
						n = N[k];
					}
				}
			}

			if (n == -1)
			{
				// Skip-Dead-End
				while (DeadEndVertexStack.length() && (n == -1))
				{
					//					unsigned int d = DeadEndVertexStack.back();
					//					DeadEndVertexStack.pop_back();
					const unsigned int d = DeadEndVertexStack.pop();

					if (pVerts[d].iNumTrisLeft > 0)
						n = d;
				}

				while (i+1 < NumVerts && (n == -1))
				{
					++i;
					if (pVerts[i].iNumTrisLeft > 0)
						n = i;
				}
			}

			f = n;
		}
	}

	// debugging purpose
	// 	int capac = N.capacity();
	// 	capac = DeadEndVertexStack.capacity();

	delete [] vertexTriAdjBuf;

	delete [] pVerts;
	delete [] pTris;
}

//=============================================================================

GPUCacheEfficiencyTester::GPUCacheEfficiencyTester(unsigned int NumTris, unsigned int NumVerts,
												   unsigned int IndexSize, const void* pIndices)
: GPUCacheOptimizer(NumTris, NumVerts, IndexSize, pIndices)
{
	for (unsigned int i = 0; i < NumTris; ++i) m_pTriMap[i] = i;
}

//=============================================================================


}
