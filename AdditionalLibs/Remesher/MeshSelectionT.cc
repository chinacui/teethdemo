/*===========================================================================*\
*                                                                            *
*                              OpenFlipper                                   *
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
*                                                                            *
\*===========================================================================*/

/*===========================================================================*\
*                                                                            *
*   $Revision$                                                       *
*   $LastChangedBy$                                                *
*   $Date$                     *
*                                                                            *
\*===========================================================================*/




//=============================================================================
//
//  IMPLEMENTATION
//
//=============================================================================

#define MESHSELECTION_C

//== INCLUDES =================================================================

#include "MeshSelectionT.hh"
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>

#include <stack>
#include <set>
//== NAMESPACES ===============================================================

namespace MeshSelection {

//== IMPLEMENTATION ==========================================================

//=========================================================
//== Vertex Selection =====================================
//=========================================================

template< typename MeshT >
inline
void selectVertices(MeshT* _mesh, const std::vector< int >& _vertices) {
  const int n_vertices = (int)_mesh->n_vertices();

  for ( uint i = 0 ; i < _vertices.size() ; ++i )
    if ( (_vertices[i] >= 0) && ( _vertices[i] < n_vertices ) )  {
      typename MeshT::VertexHandle vh(_vertices[i]);
      _mesh->status(vh).set_selected(true);
    }
}

//=========================================================

template< typename MeshT >
inline
void unselectVertices(MeshT* _mesh, const std::vector< int >& _vertices) {
  const int n_vertices = (int)_mesh->n_vertices();

  for ( uint i = 0 ; i < _vertices.size() ; ++i )
    if ( (_vertices[i] >= 0) && ( _vertices[i] < n_vertices ) )  {
      typename MeshT::VertexHandle vh(_vertices[i]);
      _mesh->status(vh).set_selected(false);
    }
}

//=========================================================

template< typename MeshT >
inline
void selectAllVertices(MeshT* _mesh) {
   typename MeshT::VertexIter v_it, v_end=_mesh->vertices_end();

   for (v_it = _mesh->vertices_begin(); v_it != v_end ; ++v_it)
      _mesh->status(*v_it).set_selected(true);

}

//=========================================================

template< typename MeshT >
inline
void clearVertexSelection(MeshT* _mesh) {
   typename MeshT::VertexIter v_it, v_end=_mesh->vertices_end();

   for (v_it = _mesh->vertices_begin(); v_it != v_end ; ++v_it)
      _mesh->status(*v_it).set_selected(false);
}

//=========================================================

template< typename MeshT >
inline
void invertVertexSelection(MeshT* _mesh) {
   typename MeshT::VertexIter v_it, v_end=_mesh->vertices_end();

   for (v_it = _mesh->vertices_begin(); v_it != v_end ; ++v_it)
      _mesh->status(*v_it).set_selected( ! _mesh->status(*v_it).selected());
}

//=========================================================


template< typename MeshT >
inline
void selectBoundaryVertices(MeshT* _mesh) {
   typename MeshT::HalfedgeIter he_it, he_end=_mesh->halfedges_end();

   for (he_it = _mesh->halfedges_begin(); he_it != he_end ; ++he_it)
      if (_mesh->is_boundary(*he_it) ) {
         _mesh->status(_mesh->to_vertex_handle(*he_it)).set_selected(true);
         _mesh->status(_mesh->from_vertex_handle(*he_it)).set_selected(true);
      }
}

//-----------------------------------------------------------------------------

template< typename MeshT >
inline
void shrinkVertexSelection(MeshT* _mesh) {
   OpenMesh::VPropHandleT< bool > temp_shrink;

   _mesh->add_property( temp_shrink, "Temp property for Vertex selection shrinking" );

   typename MeshT::VertexIter v_it, v_end=_mesh->vertices_end();

   // initialize property ( copy status to new property )
   for (v_it = _mesh->vertices_begin(); v_it != v_end ; ++v_it)
      _mesh->property(temp_shrink,*v_it) = _mesh->status(*v_it).selected();

   // update selection
   for (v_it = _mesh->vertices_begin(); v_it != v_end ; ++v_it)
      if ( _mesh->property(temp_shrink,*v_it) ) {
         _mesh->status(*v_it).set_selected( true );

         for ( typename MeshT::VertexVertexIter vv_it(*_mesh,*v_it); vv_it.is_valid(); ++vv_it)
            if ( ! _mesh->property(temp_shrink,*vv_it) ){
                _mesh->status(*v_it).set_selected( false );
                break;
            }
      }

   _mesh->remove_property(temp_shrink);
}

//=========================================================

template< typename MeshT >
inline
void growVertexSelection(MeshT* _mesh) {
   OpenMesh::VPropHandleT< bool > temp_grow;

   _mesh->add_property( temp_grow, "Temp property for Vertex selection growing" );

   // initialize property ( copy status to new property )
   typename MeshT::VertexIter v_it, v_end=_mesh->vertices_end();
   for (v_it = _mesh->vertices_begin(); v_it != v_end ; ++v_it)
      _mesh->property(temp_grow,*v_it) = _mesh->status(*v_it).selected();

   // update selection
   for (v_it = _mesh->vertices_begin(); v_it != v_end ; ++v_it)
      if ( _mesh->property(temp_grow,*v_it) )
         for ( typename MeshT::VertexVertexIter vv_it(*_mesh,*v_it); vv_it.is_valid(); ++vv_it)
            _mesh->status(*vv_it).set_selected( true );

   _mesh->remove_property(temp_grow);
}

//=========================================================

template< typename MeshT >
inline
std::vector< int > getVertexSelection(MeshT* _mesh) {
  std::vector< int > selection;

  for ( typename MeshT::VertexIter v_it= _mesh->vertices_begin() ; v_it != _mesh->vertices_end() ; ++v_it )
    if ( _mesh->status(*v_it).selected() )
      selection.push_back( v_it->idx() );

  return selection;
}

//=========================================================

template< typename MeshT >
inline
std::vector< int > getVertexSelection(MeshT* _mesh, bool& _invert) {
  std::vector< int > selection;

  int count = 0;

  for ( typename MeshT::VertexIter v_it= _mesh->vertices_begin() ; v_it != _mesh->vertices_end() ; ++v_it )
    if ( _mesh->status(*v_it).selected() )
      ++count;

  if ( count > (int)( _mesh->n_vertices() / 2) )
    _invert = true;
  else
    _invert = false;

  for ( typename MeshT::VertexIter v_it= _mesh->vertices_begin() ; v_it != _mesh->vertices_end() ; ++v_it )
    if ( _mesh->status(*v_it).selected() ^ _invert )
      selection.push_back( v_it->idx() );

  return selection;
}


template< typename MeshT >
inline
void selectBoundaryVertices(MeshT* _mesh, const typename MeshT::VertexHandle& _vh){

  OpenMesh::VPropHandleT< bool > visited;
  _mesh->add_property(visited, "Visited Vertices");

  typename MeshT::VertexIter v_it, v_end = _mesh->vertices_end();
  for (v_it = _mesh->vertices_begin(); v_it != v_end; ++v_it)
    _mesh->property(visited, *v_it) = false;

  std::stack< typename MeshT::VertexHandle > stack;
  stack.push( _vh );

  while (!stack.empty()){

    typename MeshT::VertexHandle vh = stack.top();
    stack.pop();

    if (_mesh->property(visited,vh))
      continue;

    //find outgoing boundary-edges
    for (typename MeshT::VertexOHalfedgeIter voh_it(*_mesh,vh); voh_it.is_valid(); ++voh_it)
      if ( _mesh->is_boundary( _mesh->edge_handle( *voh_it ) ) )
        stack.push( _mesh->to_vertex_handle(*voh_it) );

    //select vertex
    _mesh->property(visited,vh) = true;
    _mesh->status( vh ).set_selected(true);
  }
  _mesh->remove_property(visited);
}

template< typename MeshT >
inline
void convertVertexToEdgeSelection(MeshT* _mesh, const std::vector< int >& _vertices) {

  for ( std::vector<int>::const_iterator v = _vertices.begin(); v != _vertices.end(); ++v) {

    typename MeshT::VertexHandle vh(*v);
    typename MeshT::VertexOHalfedgeIter ohe_iter = _mesh->voh_iter(vh);

    for (; ohe_iter.is_valid(); ++ohe_iter) {
      // test if both incident vertices are in _vertices
      typename MeshT::VertexHandle ovh = _mesh->to_vertex_handle(*ohe_iter);
      // search for ovh in _vertices
      for(std::vector<int>::const_iterator it = _vertices.begin(); it != _vertices.end(); ++it) {
        if((*it) == ovh.idx()) {
          _mesh->status(_mesh->edge_handle(*ohe_iter)).set_selected(true);
          break;
        }
      }
    }
  }
}

template< typename MeshT >
inline
void convertVertexToEdgeSelection(MeshT* _mesh) {

  typename MeshT::VertexIter v_it, v_end = _mesh->vertices_end();
  for (v_it = _mesh->vertices_begin(); v_it != v_end; ++v_it) {

    if ( _mesh->status( *v_it ).selected() ) {
      typename MeshT::VertexOHalfedgeIter ohe_iter = _mesh->voh_iter(*v_it);

      for (; ohe_iter.is_valid(); ++ohe_iter) {
        // test if both incident vertices are in _vertices
        typename MeshT::VertexHandle ovh = _mesh->to_vertex_handle(*ohe_iter);
        if (_mesh->status(ovh).selected())
          _mesh->status(_mesh->edge_handle(*ohe_iter)).set_selected(true);
      }
    }
  }
}

template< typename MeshT >
inline
void convertVertexToHalfedgeSelection(MeshT* _mesh, const std::vector< int >& _vertices) {

  for (std::vector<int>::const_iterator v = _vertices.begin(); v != _vertices.end(); ++v) {

    typename MeshT::VertexHandle vh(*v);
    typename MeshT::VertexOHalfedgeIter ohe_iter = _mesh->voh_iter(vh);

    for (; ohe_iter.is_valid(); ++ohe_iter) {
      // test if both incident vertices are in _vertices
      typename MeshT::VertexHandle ovh = _mesh->to_vertex_handle(*ohe_iter);
      // search for ovh in _vertices
      for(std::vector<int>::const_iterator it = _vertices.begin(); it != _vertices.end(); ++it) {
        if((*it) == ovh.idx()) {
          _mesh->status(*ohe_iter).set_selected(true);
          _mesh->status(_mesh->opposite_halfedge_handle(*ohe_iter)).set_selected(true);
          break;
        }
      }
    }
  }
}

template< typename MeshT >
inline
void convertVertexToHalfedgeSelection(MeshT* _mesh) {

  typename MeshT::VertexIter v_it, v_end = _mesh->vertices_end();
  
  for (v_it = _mesh->vertices_begin(); v_it != v_end; ++v_it) {

    if ( _mesh->status( *v_it ).selected() ) {

      typename MeshT::VertexOHalfedgeIter ohe_iter = _mesh->voh_iter(*v_it);

      for (; ohe_iter.is_valid(); ++ohe_iter) {
        // test if both incident vertices are in _vertices
        typename MeshT::VertexHandle ovh = _mesh->to_vertex_handle(*ohe_iter);
        if (_mesh->status(ovh).selected()) {
          _mesh->status(*ohe_iter).set_selected(true);
          _mesh->status(_mesh->opposite_halfedge_handle(*ohe_iter)).set_selected(true);
        }
      }
    }
  }
}

template< typename MeshT >
inline
void convertVertexToFaceSelection(MeshT* _mesh, const std::vector< int >& _vertices) {

  for(typename MeshT::FaceIter f_it = _mesh->faces_begin(); f_it != _mesh->faces_end(); ++f_it) {
    typename MeshT::FaceVertexIter fv_it = _mesh->fv_iter(*f_it);
    // go over each vertex of each face and test if it's selected
    bool allfound = true;
    for(; fv_it.is_valid(); ++fv_it) {
      // search fv_it in _vertices
      bool onefound = false;
      for(std::vector<int>::const_iterator it = _vertices.begin(); it != _vertices.end(); ++it) {
        if((*it) == fv_it->idx()) { onefound = true; break; }
      }
      if(!onefound) {
        allfound = false;
        break;
      }
    }
    if(allfound) {
      // all incident vertices are selected -> select face
      _mesh->status(*f_it).set_selected(true);
    }
  }
}

template< typename MeshT >
inline
void convertVertexToFaceSelection(MeshT* _mesh) {

  typename MeshT::FaceIter f_it, f_end = _mesh->faces_end();
  
  for (f_it = _mesh->faces_begin(); f_it != f_end; ++f_it) {

    typename MeshT::FaceVertexIter fv_it = _mesh->fv_iter(*f_it);
    // test if all incident vertices are selected
    bool allfound = true;
    for(; fv_it.is_valid(); ++fv_it) {
      if(!_mesh->status(*fv_it).selected()) {
        allfound = false;
        break;
      }
    }
    if(allfound)
      _mesh->status(*f_it).set_selected(true);
  }
}

template< typename MeshT >
inline
void convertVertexSelectionToFeatureVertices(MeshT* _mesh) {

    for (typename MeshT::VertexIter v_it = _mesh->vertices_begin(); v_it != _mesh->vertices_end(); ++v_it) {

        if (_mesh->status(*v_it).selected()) {

            _mesh->status(*v_it).set_feature(true);
        } else {
            _mesh->status(*v_it).set_feature(false);
        }
    }
}

template< typename MeshT >
inline
void convertFeatureVerticesToVertexSelection(MeshT* _mesh) {

    for (typename MeshT::VertexIter v_it = _mesh->vertices_begin(); v_it != _mesh->vertices_end(); ++v_it) {

        if (_mesh->status(*v_it).feature()) {

            _mesh->status(*v_it).set_selected(true);
        } else {
            _mesh->status(*v_it).set_selected(false);
        }
    }
}

template< typename MeshT >
inline
void clearFeatureVertices(MeshT* _mesh) {

    for (typename MeshT::VertexIter v_it = _mesh->vertices_begin(); v_it != _mesh->vertices_end(); ++v_it) {

        _mesh->status(*v_it).set_feature(false);
    }
}

//=========================================================
//== Modeling Regions =====================================
//=========================================================

template< typename MeshT >
inline
void setArea(MeshT* _mesh, const std::vector< int >& _vertices , unsigned int _type, bool _state) {
  for ( uint i = 0 ; i < _vertices.size() ; ++i ) {
    if ( _vertices[i] > (int)_mesh->n_vertices() )
      continue;

    typename MeshT::VertexHandle vh(_vertices[i]);
    _mesh->status(vh).change_bit(_type, _state);
  }
}

template< typename MeshT >
inline
void setArea(MeshT* _mesh , unsigned int _type, bool _state) {
  for ( typename MeshT::VertexIter v_it= _mesh->vertices_begin() ; v_it != _mesh->vertices_end() ; ++v_it )
    _mesh->status(*v_it).change_bit(_type,  _state);
}

template< typename MeshT >
inline
std::vector< int > getArea(MeshT* _mesh, unsigned int _type) {
  std::vector< int > selection;

  for ( typename MeshT::VertexIter v_it= _mesh->vertices_begin() ; v_it != _mesh->vertices_end() ; ++v_it )
    if ( _mesh->status(*v_it).is_bit_set( _type ) )
      selection.push_back( v_it->idx() );

  return selection;
}

template< typename MeshT >
inline
std::vector< int > getArea(MeshT* _mesh, unsigned int _type , bool& _invert) {
  std::vector< int > selection;

  int count = 0;

  for ( typename MeshT::VertexIter v_it= _mesh->vertices_begin() ; v_it != _mesh->vertices_end() ; ++v_it )
    if ( _mesh->status(*v_it).is_bit_set( _type ) )
      ++count;

  if ( count > (int)( _mesh->n_vertices() / 2) )
    _invert = true;
  else
    _invert = false;

  for ( typename MeshT::VertexIter v_it= _mesh->vertices_begin() ; v_it != _mesh->vertices_end() ; ++v_it )
    if ( _mesh->status(*v_it).is_bit_set( _type ) ^ _invert )
      selection.push_back( v_it->idx() );

  return selection;
}


//=========================================================
//== Edge Selection =====================================
//=========================================================

template< typename MeshT >
inline
void selectEdges(MeshT* _mesh, const std::vector< int >& _edges) {
  const int n_edges = (int)_mesh->n_edges();

  for ( uint i = 0 ; i < _edges.size() ; ++i )
    if ( (_edges[i] >= 0) && ( _edges[i] < n_edges ) )  {
      typename MeshT::EdgeHandle eh(_edges[i]);
      _mesh->status(eh).set_selected(true);
    }
}

//=========================================================

template< typename MeshT >
inline
void unselectEdges(MeshT* _mesh, const std::vector< int >& _edges) {
  const int n_edges = (int)_mesh->n_edges();

  for ( uint i = 0 ; i < _edges.size() ; ++i )
    if ( (_edges[i] >= 0) && ( _edges[i] < n_edges ) )  {
      typename MeshT::EdgeHandle eh(_edges[i]);
      _mesh->status(eh).set_selected(false);
    }
}

//=========================================================

template< typename MeshT >
inline
void selectAllEdges(MeshT* _mesh) {
  typename MeshT::EdgeIter e_it, e_end=_mesh->edges_end();

  for (e_it = _mesh->edges_begin(); e_it != e_end ; ++e_it)
  _mesh->status(*e_it).set_selected(true);
}

//=========================================================

template< typename MeshT >
inline
void clearEdgeSelection(MeshT* _mesh) {
  typename MeshT::EdgeIter e_it, e_end=_mesh->edges_end();

  for (e_it = _mesh->edges_begin(); e_it != e_end ; ++e_it)
    _mesh->status(*e_it).set_selected(false);
}

//=========================================================

template< typename MeshT >
inline
void invertEdgeSelection(MeshT* _mesh) {
  typename MeshT::EdgeIter e_it, e_end=_mesh->edges_end();

  for (e_it = _mesh->edges_begin(); e_it != e_end ; ++e_it)
    _mesh->status(*e_it).set_selected( ! _mesh->status(*e_it).selected());
}

//=========================================================

template<typename MeshT>
inline
void growEdgeSelection(MeshT* _mesh) {
    std::set<typename MeshT::EdgeHandle> selectedEhs;
    for (typename MeshT::EdgeIter e_it = _mesh->edges_begin(), e_end = _mesh->edges_end();
            e_it != e_end; ++e_it) {

        if (!_mesh->status(*e_it).selected()) continue;

        const typename MeshT::HalfedgeHandle he = _mesh->halfedge_handle(*e_it, 0);
        const typename MeshT::VertexHandle vhs[] = { _mesh->from_vertex_handle(he),
                                                     _mesh->to_vertex_handle(he) };

        for (int i = 0; i < 2; ++i) {
            for (typename MeshT::VertexEdgeIter ve_it = _mesh->ve_begin(vhs[i]), ve_end = _mesh->ve_end(vhs[i]);
                    ve_it != ve_end; ++ve_it) {

                selectedEhs.insert(*ve_it);
            }
        }

    }

    for (typename std::set<typename MeshT::EdgeHandle>::const_iterator it = selectedEhs.begin(); it != selectedEhs.end(); ++it)
        _mesh->status(*it).set_selected(true);
}

//=========================================================


template< typename MeshT >
inline
void selectBoundaryEdges(MeshT* _mesh) {
  typename MeshT::EdgeIter e_it, e_end=_mesh->edges_end();

  for (e_it = _mesh->edges_begin(); e_it != e_end ; ++e_it)
    if ( _mesh->is_boundary( _mesh->halfedge_handle(*e_it,0) ) ||
         _mesh->is_boundary( _mesh->halfedge_handle(*e_it,1) ) )
      _mesh->status(*e_it).set_selected( true );
}

//=========================================================

template< typename MeshT >
inline
std::vector< int > getEdgeSelection(MeshT* _mesh) {
  std::vector< int > selection;

  for ( typename MeshT::EdgeIter e_it= _mesh->edges_begin() ; e_it != _mesh->edges_end() ; ++e_it )
    if ( _mesh->status(*e_it).selected() )
      selection.push_back( e_it->idx() );

  return selection;
}

//=========================================================

template< typename MeshT >
inline
std::vector< int > getEdgeSelection(MeshT* _mesh, bool& _invert) {
  std::vector< int > selection;

  int count = 0;

  for ( typename MeshT::VertexIter e_it= _mesh->edges_begin() ; e_it != _mesh->edges_end() ; ++e_it )
    if ( _mesh->status(*e_it).selected() )
      ++count;

  if ( count > (int)( _mesh->n_vertices() / 2) )
    _invert = true;
  else
    _invert = false;

  for ( typename MeshT::VertexIter e_it= _mesh->edges_begin() ; e_it != _mesh->edges_end() ; ++e_it )
    if ( _mesh->status(*e_it).selected() ^ _invert )
      selection.push_back( e_it->idx() );

  return selection;
}

template< typename MeshT >
inline
void convertEdgeToVertexSelection(MeshT* _mesh, const std::vector< int >& _edges) {

	for (std::vector<int>::const_iterator e = _edges.begin(); e != _edges.end(); ++e) {

		typename MeshT::EdgeHandle eh(*e);
		typename MeshT::HalfedgeHandle heh0 = _mesh->halfedge_handle(eh, 0);

		typename MeshT::VertexHandle vh0 = _mesh->to_vertex_handle(heh0);
		typename MeshT::VertexHandle vh1 = _mesh->from_vertex_handle(heh0);

		_mesh->status(vh0).set_selected(true);
		_mesh->status(vh1).set_selected(true);
	}
}

template< typename MeshT >
inline
void convertEdgeToVertexSelection(MeshT* _mesh) {

  for ( typename MeshT::EdgeIter e_it= _mesh->edges_begin() ; e_it != _mesh->edges_end() ; ++e_it )
    
    if ( _mesh->status(*e_it).selected() ){

      typename MeshT::HalfedgeHandle heh0 = _mesh->halfedge_handle(*e_it, 0);

      typename MeshT::VertexHandle vh0 = _mesh->to_vertex_handle(heh0);
      typename MeshT::VertexHandle vh1 = _mesh->from_vertex_handle(heh0);

      _mesh->status(vh0).set_selected(true);
      _mesh->status(vh1).set_selected(true);
    }
}

template< typename MeshT >
inline
void convertEdgeToFaceSelection(MeshT* _mesh, const std::vector< int >& _edges) {

  for(typename MeshT::FaceIter f_it = _mesh->faces_begin(); f_it != _mesh->faces_end(); ++f_it) {
    typename MeshT::FaceEdgeIter fe_it = _mesh->fe_iter(*f_it);
    // go over each edge of each face and test if it's selected
    bool allfound = true;
    for(; fe_it.is_valid(); ++fe_it) {
      // search fe_it in _edges
      bool onefound = false;
      for(std::vector<int>::const_iterator it = _edges.begin(); it != _edges.end(); ++it) {
        if((*it) == fe_it->idx()) { onefound = true; break; }
      }
      if(!onefound) {
        allfound = false;
        break;
      }
    }
    if(allfound) {
      // all incident vertices are selected -> select face
      _mesh->status(*f_it).set_selected(true);
    }
  }
}

template< typename MeshT >
inline
void convertEdgeToFaceSelection(MeshT* _mesh) {

  typename MeshT::FaceIter f_it, f_end = _mesh->faces_end();
  
  for (f_it = _mesh->faces_begin(); f_it != f_end; ++f_it) {

    typename MeshT::FaceEdgeIter fe_it = _mesh->fe_iter(*f_it);
    // test if all incident edges are selected
    bool allfound = true;
    for(; fe_it.is_valid(); ++fe_it) {
      if(!_mesh->status(*fe_it).selected()) {
        allfound = false;
        break;
      }
    }
    if(allfound)
      _mesh->status(*f_it).set_selected(true);
  }
}

template< typename MeshT >
inline
void convertEdgeToHalfedgeSelection(MeshT* _mesh) {

  for ( typename MeshT::EdgeIter e_it= _mesh->edges_begin() ; e_it != _mesh->edges_end() ; ++e_it )
    
    if ( _mesh->status(*e_it).selected() ){

      _mesh->status(_mesh->halfedge_handle(*e_it, 0)).set_selected(true);
      _mesh->status(_mesh->halfedge_handle(*e_it, 1)).set_selected(true);
    }
}

template< typename MeshT >
inline
void convertEdgeSelectionToFeatureEdges(MeshT* _mesh) {

    for (typename MeshT::EdgeIter e_it = _mesh->edges_begin(); e_it != _mesh->edges_end(); ++e_it) {

        if (_mesh->status(*e_it).selected()) {

            _mesh->status(*e_it).set_feature(true);
        } else {
            _mesh->status(*e_it).set_feature(false);
        }
    }
}

template< typename MeshT >
inline
void convertFeatureEdgesToEdgeSelection(MeshT* _mesh) {

    for (typename MeshT::EdgeIter e_it = _mesh->edges_begin(); e_it != _mesh->edges_end(); ++e_it) {

        if (_mesh->status(*e_it).feature()) {

            _mesh->status(*e_it).set_selected(true);
        } else {
            _mesh->status(*e_it).set_selected(false);
        }
    }
}

template< typename MeshT >
inline
void clearFeatureEdges(MeshT* _mesh) {

    for (typename MeshT::EdgeIter e_it = _mesh->edges_begin(); e_it != _mesh->edges_end(); ++e_it) {

        _mesh->status(*e_it).set_feature(false);
    }
}

//=========================================================
//== Halfedge Selection =====================================
//=========================================================

template< typename MeshT >
inline
void selectHalfedges(MeshT* _mesh, const std::vector< int >& _halfedges) {
  const int n_halfedges = (int)_mesh->n_halfedges();

  for ( uint i = 0 ; i < _halfedges.size() ; ++i )
    if ( (_halfedges[i] >= 0) && ( _halfedges[i] < n_halfedges ) )  {
      typename MeshT::HalfedgeHandle heh(_halfedges[i]);
      _mesh->status(heh).set_selected(true);
    }
}

//=========================================================

template< typename MeshT >
inline
void unselectHalfedges(MeshT* _mesh, const std::vector< int >& _halfedges) {
  const int n_halfedges = (int)_mesh->n_halfedges();

  for ( uint i = 0 ; i < _halfedges.size() ; ++i )
    if ( (_halfedges[i] >= 0) && ( _halfedges[i] < n_halfedges ) )  {
      typename MeshT::HalfedgeHandle heh(_halfedges[i]);
      _mesh->status(heh).set_selected(false);
    }
}

//=========================================================

template< typename MeshT >
inline
void selectAllHalfedges(MeshT* _mesh) {
  typename MeshT::HalfedgeIter he_it, he_end=_mesh->halfedges_end();

  for (he_it = _mesh->halfedges_begin(); he_it != he_end ; ++he_it)
  _mesh->status(*he_it).set_selected(true);
}

//=========================================================

template< typename MeshT >
inline
void clearHalfedgeSelection(MeshT* _mesh) {
  typename MeshT::HalfedgeIter he_it, he_end=_mesh->halfedges_end();

  for (he_it = _mesh->halfedges_begin(); he_it != he_end ; ++he_it)
    _mesh->status(*he_it).set_selected(false);
}

//=========================================================

template< typename MeshT >
inline
void invertHalfedgeSelection(MeshT* _mesh) {
  typename MeshT::HalfedgeIter he_it, he_end=_mesh->halfedges_end();

  for (he_it = _mesh->halfedges_begin(); he_it != he_end ; ++he_it)
    _mesh->status(*he_it).set_selected( ! _mesh->status(*he_it).selected());
}

//=========================================================


template< typename MeshT >
inline
void selectBoundaryHalfedges(MeshT* _mesh) {
  typename MeshT::HalfedgeIter he_it, he_end=_mesh->halfedges_end();

  for (he_it = _mesh->halfedges_begin(); he_it != he_end ; ++he_it)
    if ( _mesh->is_boundary( *he_it))
      _mesh->status(*he_it).set_selected( true );
}

//=========================================================

template< typename MeshT >
inline
std::vector< int > getHalfedgeSelection(MeshT* _mesh) {
  std::vector< int > selection;

  for ( typename MeshT::HalfedgeIter he_it= _mesh->halfedges_begin() ; he_it != _mesh->halfedges_end() ; ++he_it )
    if ( _mesh->status(*he_it).selected() )
      selection.push_back( he_it->idx() );

  return selection;
}

template< typename MeshT >
inline
void convertHalfedgeToVertexSelection(MeshT* _mesh) {
    
    for ( typename MeshT::HalfedgeIter he_it= _mesh->halfedges_begin() ; he_it != _mesh->halfedges_end() ; ++he_it ) {
        
        if(_mesh->status(*he_it).selected()) {
            _mesh->status(_mesh->to_vertex_handle(*he_it)).set_selected(true);
            _mesh->status(_mesh->from_vertex_handle(*he_it)).set_selected(true);
        }
    }
}

template< typename MeshT >
inline
void convertHalfedgeToEdgeSelection(MeshT* _mesh) {
    
    for ( typename MeshT::HalfedgeIter he_it= _mesh->halfedges_begin() ; he_it != _mesh->halfedges_end() ; ++he_it ) {
        
        if(_mesh->status(*he_it).selected()) {
            _mesh->status(_mesh->edge_handle(*he_it)).set_selected(true);
        }
    }
}

template< typename MeshT >
inline
void convertHalfedgeToFaceSelection(MeshT* _mesh) {
    // Note: A face is not only selected
    // iff all incident halfedges are selected but
    // at least one of them. This is, however,
    // desired in some cases.
    for ( typename MeshT::HalfedgeIter he_it= _mesh->halfedges_begin() ; he_it != _mesh->halfedges_end() ; ++he_it ) {
        
        if(_mesh->status(*he_it).selected()) {
            _mesh->status(_mesh->face_handle(*he_it)).set_selected(true);
        }
    }
}

//=========================================================
//== Face Selection =======================================
//=========================================================

template< typename MeshT >
inline
void selectFaces(MeshT* _mesh, const std::vector< int >& _faces) {
  const int n_faces = (int)_mesh->n_faces();

  for ( uint i = 0 ; i < _faces.size() ; ++i )
    if ( (_faces[i] >= 0) && ( _faces[i] < n_faces ) )  {
      typename MeshT::FaceHandle fh(_faces[i]);
      _mesh->status(fh).set_selected(true);
    }
}

//=========================================================

template< typename MeshT >
inline
void unselectFaces(MeshT* _mesh, const std::vector< int >& _faces) {
  const int n_faces = (int)_mesh->n_faces();

  for ( uint i = 0 ; i < _faces.size() ; ++i )
    if ( (_faces[i] >= 0) && ( _faces[i] < n_faces ) )  {
      typename MeshT::FaceHandle fh(_faces[i]);
      _mesh->status(fh).set_selected(false);
    }
}

//=========================================================

template< typename MeshT >
inline
void selectAllFaces(MeshT* _mesh) {
   typename MeshT::FaceIter f_it, f_end=_mesh->faces_end();

   for (f_it = _mesh->faces_begin(); f_it != f_end ; ++f_it)
      _mesh->status(*f_it).set_selected(true);
}

//=========================================================


template< typename MeshT >
inline
void clearFaceSelection(MeshT* _mesh) {
   typename MeshT::FaceIter f_it, f_end=_mesh->faces_end();

   for (f_it = _mesh->faces_begin(); f_it != f_end ; ++f_it)
      _mesh->status(*f_it).set_selected(false);
}

//-----------------------------------------------------------------------------


template< typename MeshT >
inline
void invertFaceSelection(MeshT* _mesh) {
   typename MeshT::FaceIter f_it, f_end=_mesh->faces_end();

   for (f_it = _mesh->faces_begin(); f_it != f_end ; ++f_it)
      _mesh->status(*f_it).set_selected( ! _mesh->status(*f_it).selected());
}

//=========================================================

template< typename MeshT >
inline
void selectBoundaryFaces(MeshT* _mesh) {
  typename MeshT::HalfedgeIter he_it, he_end=_mesh->halfedges_end();

  for (he_it = _mesh->halfedges_begin(); he_it != he_end ; ++he_it)
    if (_mesh->is_boundary(*he_it) ) {
        for (typename MeshT::VertexFaceIter vf_it(*_mesh ,_mesh->to_vertex_handle(*he_it) ) ; vf_it.is_valid() ; ++vf_it)
          _mesh->status(*vf_it).set_selected(true);
        for (typename MeshT::VertexFaceIter vf_it(*_mesh ,_mesh->from_vertex_handle(*he_it) ) ; vf_it.is_valid() ; ++vf_it)
          _mesh->status(*vf_it).set_selected(true);
    }
}

//=========================================================


template< typename MeshT >
inline
void shrinkFaceSelection(MeshT* _mesh) {
   OpenMesh::FPropHandleT< bool > temp_shrink;

   _mesh->add_property( temp_shrink, "Temp property for Face selection shrinking" );

   typename MeshT::FaceIter f_it, f_end=_mesh->faces_end();

   // initialize property ( copy status to new property )
   for (f_it = _mesh->faces_begin(); f_it != f_end ; ++f_it)
      _mesh->property(temp_shrink,*f_it) = _mesh->status(*f_it).selected();

   // Shrink selection ( deselects all faces which are adjacent to a boundary vertex of the original selection)
   for (f_it = _mesh->faces_begin(); f_it != f_end ; ++f_it)
      if ( _mesh->property(temp_shrink,*f_it) ) {
         bool boundary = false;
         for ( typename MeshT::FaceVertexIter fv_it(*_mesh,*f_it); fv_it.is_valid() ; ++fv_it) {
            for ( typename MeshT::VertexFaceIter vf_it(*_mesh,*fv_it); vf_it.is_valid() ; ++vf_it) {
               if ( ! _mesh->property(temp_shrink,*vf_it) ) {
                 boundary = true;
               }
            }
            if ( boundary )
               break;
         }

         _mesh->status(*f_it).set_selected( !boundary );
      }

   _mesh->remove_property(temp_shrink);
}

//=========================================================

template< typename MeshT >
inline
void growFaceSelection(MeshT* _mesh) {
  OpenMesh::FPropHandleT< bool > temp_grow;

  _mesh->add_property( temp_grow, "Temp property for Face selection growing" );

  typename MeshT::FaceIter f_it, f_end=_mesh->faces_end();

  // initialize property ( copy status to new property )
  for (f_it = _mesh->faces_begin(); f_it != f_end ; ++f_it)
    _mesh->property(temp_grow,*f_it) = _mesh->status(*f_it).selected();

  // Grow selection ( selects all faces which are adjacent to a vertex of a already selected face)
  for (f_it = _mesh->faces_begin(); f_it != f_end ; ++f_it)
    if ( _mesh->property(temp_grow,*f_it) )
        for ( typename MeshT::FaceVertexIter fv_it(*_mesh,*f_it); fv_it.is_valid() ; ++fv_it)
          for ( typename MeshT::VertexFaceIter vf_it(*_mesh,*fv_it); vf_it.is_valid() ; ++vf_it)
              _mesh->status(*vf_it).set_selected( true );

  _mesh->remove_property(temp_grow);
}

//=========================================================

template< typename MeshT >
inline
std::vector< int > getFaceSelection(MeshT* _mesh) {
  std::vector< int > selection;

  for ( typename MeshT::FaceIter f_it= _mesh->faces_begin() ; f_it != _mesh->faces_end() ; ++f_it )
    if ( _mesh->status(*f_it).selected() )
      selection.push_back( f_it->idx() );

  return selection;
}

//=========================================================

template< typename MeshT >
inline
std::vector< int > getFaceSelection(MeshT* _mesh, bool& _invert) {
  std::vector< int > selection;

  int count = 0;

  for ( typename MeshT::FaceIter f_it= _mesh->faces_begin() ; f_it != _mesh->faces_end() ; ++f_it )
    if ( _mesh->status(*f_it).selected() )
      ++count;

  if ( count > (int)( _mesh->n_vertices() / 2) )
    _invert = true;
  else
    _invert = false;

  for ( typename MeshT::FaceIter f_it= _mesh->faces_begin() ; f_it != _mesh->faces_end() ; ++f_it )
    if ( _mesh->status(*f_it).selected() ^ _invert )
      selection.push_back( f_it->idx() );

  return selection;
}

template< typename MeshT >
inline
void convertFaceToVertexSelection(MeshT* _mesh, const std::vector< int >& _faces) {

	for (std::vector<int>::const_iterator f = _faces.begin(); f != _faces.end(); ++f) {

		typename MeshT::FaceHandle fh(*f);

		typename MeshT::FaceVertexIter v_iter = _mesh->fv_iter(fh);

		for (; v_iter.is_valid(); ++v_iter) {
			_mesh->status(*v_iter).set_selected(true);
		}
	}
}

template< typename MeshT >
inline
void convertFaceToVertexSelection(MeshT* _mesh) {

  for ( typename MeshT::FaceIter f_it= _mesh->faces_begin() ; f_it != _mesh->faces_end() ; ++f_it )
    
    if ( _mesh->status(*f_it).selected() ){

      typename MeshT::FaceVertexIter v_iter = _mesh->fv_iter(*f_it);

      for (; v_iter.is_valid(); ++v_iter)
        _mesh->status(*v_iter).set_selected(true);
    }
}

template< typename MeshT >
inline
void convertFaceToEdgeSelection(MeshT* _mesh, const std::vector< int >& _faces) {

	for (std::vector<int>::const_iterator f = _faces.begin(); f != _faces.end(); ++f) {

		typename MeshT::FaceHandle fh(*f);

		typename MeshT::FaceEdgeIter e_iter = _mesh->fe_iter(fh);

		for (; e_iter.is_valid(); ++e_iter) {
			_mesh->status(*e_iter).set_selected(true);
		}
	}
}

template< typename MeshT >
inline
void convertFaceToEdgeSelection(MeshT* _mesh) {

  for ( typename MeshT::FaceIter f_it= _mesh->faces_begin() ; f_it != _mesh->faces_end() ; ++f_it )
    
    if ( _mesh->status(*f_it).selected() ){

      typename MeshT::FaceEdgeIter e_iter = _mesh->fe_iter(*f_it);

      for (; e_iter.is_valid(); ++e_iter)
        _mesh->status(*e_iter).set_selected(true);
    }
}

template< typename MeshT >
inline
void convertFaceToHalfedgeSelection(MeshT* _mesh) {

  for ( typename MeshT::FaceIter f_it= _mesh->faces_begin() ; f_it != _mesh->faces_end() ; ++f_it )
    
    if ( _mesh->status(*f_it).selected() ){

      typename MeshT::FaceHalfedgeIter fh_iter = _mesh->fh_iter(*f_it);

      for (; fh_iter.is_valid(); ++fh_iter)
        _mesh->status(*fh_iter).set_selected(true);
    }
}

template< typename MeshT >
inline
void convertFaceSelectionToFeatureFaces(MeshT* _mesh) {

    for (typename MeshT::FaceIter f_it = _mesh->faces_begin(); f_it != _mesh->faces_end(); ++f_it) {

        if (_mesh->status(*f_it).selected()) {

            _mesh->status(*f_it).set_feature(true);
        } else {
            _mesh->status(*f_it).set_feature(false);
        }
    }
}

template< typename MeshT >
inline
void convertFeatureFacesToFaceSelection(MeshT* _mesh) {

    for (typename MeshT::FaceIter f_it = _mesh->faces_begin(); f_it != _mesh->faces_end(); ++f_it) {

        if (_mesh->status(*f_it).feature()) {

            _mesh->status(*f_it).set_selected(true);
        } else {
            _mesh->status(*f_it).set_selected(false);
        }
    }
}

template< typename MeshT >
inline
void clearFeatureFaces(MeshT* _mesh) {

    for (typename MeshT::FaceIter f_it = _mesh->faces_begin(); f_it != _mesh->faces_end(); ++f_it) {

        _mesh->status(*f_it).set_feature(false);
    }
}

//=============================================================================
} // MeshSelection Namespace
//=============================================================================
