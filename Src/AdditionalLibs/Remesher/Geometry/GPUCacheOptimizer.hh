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

#ifndef ACG_GPU_CACHE_OPT_HH
#define ACG_GPU_CACHE_OPT_HH


//== INCLUDES =================================================================

#include <ACG/Config/ACGDefines.hh>

//== FORWARDDECLARATIONS ======================================================

//== NAMESPACES ===============================================================

namespace ACG {

//== CLASS DEFINITION =========================================================


/** \class GPUCacheOptimizer GPUCacheOptimizer.hh <ACG/Scenegraph/MeshNodeT.hh>

    This class optimizes a triangle list for efficient vertex cache usage.
*/

class ACGDLLEXPORT GPUCacheOptimizer
{
private:
  // copy ops are private to prevent copying
  GPUCacheOptimizer(const GPUCacheOptimizer&);            // no implementation
  GPUCacheOptimizer& operator=(const GPUCacheOptimizer&); // no implementation
public:

	/** \brief constructor
	*     
	* The constructor needs a mesh on which this node will work.
	*
	* @param NumTris number of triangles
	*	@param NumVerts number of vertices
	*	@param IndexSize size in bytes of one index: 1, 2, 4 supported
	*	@param pIndices [in] index buffer
	*/
	GPUCacheOptimizer(unsigned int NumTris, unsigned int NumVerts, unsigned int IndexSize, const void* pIndices);
	virtual ~GPUCacheOptimizer(void);

	
	/** \brief Retrieves the map from dst triangle to src triangle.
	* how to remap:
	*  for each triangle t in DstTriBuffer
	*    DstTriBuffer[t] = SrcTriBuffer[TriangleMap[t]]
	* 
	* you can also use WriteIndexBuffer() to get the result
	*/
	const unsigned int* GetTriangleMap() const {return m_pTriMap;}

	/** \brief Applies the remapping on the initial pIndices constructor's param
	* and stores the result in the given destination buffer.
	*
	* @param DstIndexSize size in bytes of one index in the dst buffer
	*	@param pDst [out] index buffer to store the result may also be the same memory as pIndices of constructor's call
	* NOTE:
	* make sure pIndices is not modified/deleted between constructor's and this call
	*/
	void WriteIndexBuffer(unsigned int DstIndexSize, void* pDst);


	/** \brief Reorders vertex buffer to minimize memory address jumps.
	 *
	 * for best performace, use AFTER triangle remapping
	 * see description on RemapVertices() on how to apply this map
	 *
	 *
	 * @param pIndices [in] index buffer
	 * @param pVertMap [out] vertex remap, allocate NumVerts unsigned ints before calling
	 *	      dst vertex index = pVertMap[src vertex index]
	 *				NOTE: if a vertex is not referenced by any triangle, the value 0xFFFFFFFF
	 *				will be stored as the destination index!
	 * @param NumTris   Number of triangles
	 * @param NumVerts  Number of vertices
	 * @param IndexSize Size of the index
	*/
	static void OptimizeVertices(unsigned int NumTris, unsigned int NumVerts, unsigned int IndexSize,
								 const void* pIndices, unsigned int* pVertMap);
	// this function is declared static to be able to operate on the whole model
	// instead of just a subset of it
	// example use:
	// - optimize triangle list per material group
  // - optimize vertex buffer on whole mesh independently of material subsets


	/** \brief Applies the remap table of OptimizeVertices to a vertex and index buffer
	 *
	 * Pseudo code for manual remapping
	 * \code
	 *   for each index i in IndexBuffer:
	 *     IndexBuffer[i] = VertMap[IndexBuffer[i]]
	 *   TmpBuf = VertexBuffer
	 *   for each vertex v in TmpBuf
	 *     if (VertMap[v] != 0xFFFFFFFF)
	 *       VertexBuffer[VertMap[v]] = TmpBuf[v]
	 * \endcode
	 *
   * @param NumVerts  Number of vertices
	 * @param	pVertMap  vertex remap, result from OptimizeVertices() (input)
	 * @param	IndexSize size in bytes of one index: 1, 2, 4 supported
	 * @param pInOutIndices (triangle list) index buffer, remapped after call (input/output)
	 * @param VertexStride  size in bytes of one vertex
	 * @param pInOutVertices vertex buffer, remapped after call (input/output)
	 * @param NumTris   Number of triangles
	 *
	 */
	static void RemapVertices(unsigned int NumTris, unsigned int NumVerts, const unsigned int* pVertMap,
		unsigned int IndexSize, void* pInOutIndices, unsigned int VertexStride, void* pInOutVertices);


	/** \brief 
	* Returns the total number of vertex transforms performed with a certain VertexCache.
	*/
	unsigned int ComputeNumberOfVertexTransformations(unsigned int VertexCacheSize = 16);

	/** \brief Measures the efficiency use of the vertex cache.
	*  ACMR: Average Cache Miss Ratio
	* @return ratio: # vertex transformations / # tris
	*/
	float ComputeACMR(unsigned int VertexCacheSize = 16);

	/** \brief Measures the efficiency use of the vertex cache.
	* ATVR: Average Transform to Vertex Ratio
	* similar to ACMR, but easier to interpret
	* the optimal value is 1.0 given a mesh w/o duplicate vertices
	* @return ratio: # vertex transformations / # verts
	*/
	float ComputeATVR(unsigned int VertexCacheSize = 16);

protected:
	// look up  m_pIndices w.r.t. index size at location 'i'
	unsigned int GetIndex(unsigned int i) const;

	static unsigned int GetIndex(unsigned int i, unsigned int IndexSize, const void* pIB);
	static void SetIndex(unsigned int i, unsigned int val, unsigned int IndexSize, void* pIB);

	virtual void MakeAbstract() = 0;

protected:


	unsigned int m_NumVerts;
	unsigned int m_NumTris;
	
	/**
	TriMap[new tri index] = old tri index
	allocated in base class, computed in child class
	*/
	unsigned int* m_pTriMap;

private:
	unsigned int m_IndexSize;
	const void* m_pIndices;


	unsigned int m_NumTransformations;


protected:
	/** internal data structures used
	 forsyth and tipsify implementation
	*/

	struct Opt_Vertex
	{
		// Opt_Vertex(): iCachePos(-1), fScore(0.0f), iNumTrisTotal(0), iNumTrisLeft(0), pTris(0) {}
		// ~Opt_Vertex() {delete [] pTris;}

		int iCachePos;
		float fScore;
		/// # tris using this vertex
		int iNumTrisTotal;
		/// # tris left to add to final result
		int iNumTrisLeft;
		unsigned int* pTris;

		/// forsyth's score function
		void FindScore(unsigned int MaxSizeVertexCache);

		void RemoveTriFromList(unsigned int tri);
	};

	struct Opt_Tris
	{
		Opt_Tris() : bAdded(0), fScore(0.0f) 
		{ v[0] = v[1] = v[2] = 0xDEADBEEF;}

		int bAdded;
		/// sum of scores of vertices
		float fScore;
		/// vertices of this triangle
		unsigned int v[3];

		inline void FindScore(const Opt_Vertex* pVertices)
		{
			fScore = 0.0f;
			for (int i = 0; i < 3; ++i)
				fScore += pVertices[v[i]].fScore;
		}
	};
};

//////////////////////////////////////////////////////////////////////////

/** \class GPUCacheOptimizerTipsify GPUCacheOptimizer.hh

    Implementation of "Fast Triangle Reordering for Vertex Locality and Reduced Overdraw" by Sander et. al.
	 http://www.cs.princeton.edu/gfx/pubs/Sander_2007_%3ETR/index.php
*/

class ACGDLLEXPORT GPUCacheOptimizerTipsify : public GPUCacheOptimizer
{
public:

	/** \brief The actual computation happens here in this constructor
	 *
   * @param NumVerts  Number of vertices
   * @param IndexSize size in bytes of one index: 1, 2, 4 supported
   * @param pIndices  index buffer
   * @param CacheSize number of entries in the vertex cache
   * @param NumTris   Number of triangles
	*/
	GPUCacheOptimizerTipsify(unsigned int CacheSize,
	                         unsigned int NumTris,
	                         unsigned int NumVerts,
	                         unsigned int IndexSize,
	                         const void* pIndices);

private:

	void MakeAbstract(){}

	/// Simple and fast fixed size stack used in tipsify implementation
	struct RingStack
	{
	private:
		unsigned int* pStack;
		unsigned int uiStart, uiLen;
		unsigned int uiSize;

		inline unsigned int pos(unsigned int i) const
		{
			unsigned int t = uiStart + i;
			if (t >= uiLen) t -= uiLen;
			return t;
		}

	public:
		RingStack(unsigned int _uiSize) :
		  uiStart(0),
		  uiLen(0),
		  uiSize(_uiSize)
    {
      pStack = new unsigned int[uiSize];
    }

    RingStack(const RingStack& _other)
    {
      // Copy meta data
      uiSize  = _other.uiSize;
      uiLen   = _other.uiLen;
      uiStart = _other.uiStart;

      // Create empty storage
      pStack = new unsigned int[uiSize];

      // Copy storage from original to current storage
      for (unsigned int  i = 0 ; i < uiSize; ++i )
        pStack[i] = _other.pStack[i];

    }

		~RingStack() {delete [] pStack;}

		unsigned int length() const {return uiLen;} ///< current stack length
		unsigned int size() const {return uiSize;}  ///< reserved stack size i.e. maximum length

		inline void push(unsigned int v)
		{
			if (uiLen == uiSize)
			{
				pStack[uiStart++] = v;
				if (uiStart == uiSize) uiStart = 0;
			}
			else
				pStack[pos(uiLen++)] = v; // pos(uiLen) gives the index of the last element + 1
		}

		inline unsigned int pop()
		{
			if (uiSize && uiLen) return pStack[pos(--uiLen)];
			return 0xFFFFFFFF;
		}
	};
};

/** \class GPUCacheEfficiencyTester GPUCacheOptimizer.hh

    simple class providing ATVR and ACMR computations w/o any optimizing
*/
class ACGDLLEXPORT GPUCacheEfficiencyTester : public GPUCacheOptimizer
{
public:
	GPUCacheEfficiencyTester(unsigned int NumTris, unsigned int NumVerts, unsigned int IndexSize, const void* pIndices);

private:
	void MakeAbstract(){}
};



//=============================================================================
} // namespace ACG
//=============================================================================
#endif // ACG_GPU_CACHE_OPT_HH defined
//=============================================================================

