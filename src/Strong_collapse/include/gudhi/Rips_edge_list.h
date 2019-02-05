/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Clément Maria, Pawel Dlotko, Vincent Rouvreau Siddharth Pritam
 *
 *    Copyright (C) 2016  INRIA
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef RIPS_edge_list_H_
#define RIPS_edge_list_H_

#include <gudhi/Debug_utils.h>
#include <gudhi/graph_simplicial_complex.h>
#include <boost/graph/adjacency_list.hpp>

#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <limits>  // for numeric_limits
#include <utility>  // for pair<>


namespace Gudhi {

namespace rips_edge_list {

/**
 * \class Rips_complex
 * \brief Rips complex data structure.
 * 
 * \ingroup rips_complex
 * 
 * \details
 * The data structure is a one skeleton graph, or Rips graph, containing edges when the edge length is less or equal
 * to a given threshold. Edge length is computed from a user given point cloud with a given distance function, or a
 * distance matrix.
 * 
 * \tparam Filtration_value is the type used to store the filtration values of the simplicial complex.
 */
template<typename Filtration_value>
class Rips_edge_list {
 public:
  /**
   * \brief Type of the one skeleton graph stored inside the Rips complex structure.
   */
  // typedef typename boost::adjacency_list < boost::vecS, boost::vecS, boost::undirectedS
  // , boost::property < vertex_filtration_t, Filtration_value >
  // , boost::property < edge_filtration_t, Filtration_value >> OneSkeletonGraph;

 private:
  typedef int Vertex_handle;

 public:
  /** \brief Rips_complex constructor from a list of points.
   *
   * @param[in] points Range of points.
   * @param[in] threshold Rips value.
   * @param[in] distance distance function that returns a `Filtration_value` from 2 given points.
   * 
   * \tparam ForwardPointRange must be a range for which `std::begin` and `std::end` return input iterators on a
   * point.
   *
   * \tparam Distance furnishes `operator()(const Point& p1, const Point& p2)`, where
   * `Point` is a point from the `ForwardPointRange`, and that returns a `Filtration_value`.
   */
  template<typename ForwardPointRange, typename Distance >
  Rips_edge_list(const ForwardPointRange& points, Filtration_value threshold, Distance distance) {
    compute_proximity_graph(points, threshold, distance);
  }

  /** \brief Rips_complex constructor from a distance matrix.
   *
   * @param[in] distance_matrix Range of distances.
   * @param[in] threshold Rips value.
   * 
   * \tparam DistanceMatrix must have a `size()` method and on which `distance_matrix[i][j]` returns
   * the distance between points \f$i\f$ and \f$j\f$ as long as \f$ 0 \leqslant i < j \leqslant
   * distance\_matrix.size().\f$
   */
  template<typename DistanceMatrix>
  Rips_edge_list(const DistanceMatrix& distance_matrix, Filtration_value threshold) {
    compute_proximity_graph(boost::irange((size_t)0, distance_matrix.size()), threshold,
                            [&](size_t i, size_t j){return distance_matrix[j][i];});
  }

  /** \brief Initializes the egde list (one skeleton) complex from the Rips graph 
   *
   * \tparam EdgeListForRips must meet `EdgeListForRips` concept.
   * 
   * @param[in] edges EdgeListForRips to be created.
   * @param[in] dim_max graph expansion for Rips until this given maximal dimension.
   * @exception std::invalid_argument In debug mode, if `edges.num_vertices()` does not return 0.
   *
   */
  template <typename EdgeListForRips>
  void create_edges(EdgeListForRips& edge_list) {
    GUDHI_CHECK(edges.num_vertices() == 0,
                std::invalid_argument("Rips_complex::create_complex - edge list is not empty"));

    // sort the tuple (filteration_valuem, (v1,v2){edge}) 
    //By default the sort is done on the first element, so here it's filteration value.
    std::sort(edges.begin(),edges.end());
    for(size_t i = 0; i < edges.size(); i++ )
      edge_list.emplace_back(edges.at(i));

  }

 private:
  /** \brief Computes the proximity graph of the points.
   *
   * If points contains n elements, the proximity graph is the graph with n vertices, and an edge [u,v] iff the
   * distance function between points u and v is smaller than threshold.
   *
   * \tparam ForwardPointRange furnishes `.begin()` and `.end()`
   * methods.
   *
   * \tparam Distance furnishes `operator()(const Point& p1, const Point& p2)`, where
   * `Point` is a point from the `ForwardPointRange`, and that returns a `Filtration_value`.
   */
  template< typename ForwardPointRange, typename Distance >
  void compute_proximity_graph(const ForwardPointRange& points, Filtration_value threshold,
               Distance distance) {
    edges.clear();
    // Compute the proximity graph of the points.
    // If points contains n elements, the proximity graph is the graph with n vertices, and an edge [u,v] iff the
    // distance function between points u and v is smaller than threshold.
    // --------------------------------------------------------------------------------------------
    // Creates the vector of edges and its filtration values (returned by distance function)
    Vertex_handle idx_u = 0;
    for (auto it_u = std::begin(points); it_u != std::end(points); ++it_u, ++idx_u) {
      Vertex_handle idx_v = idx_u + 1;
      for (auto it_v = it_u + 1; it_v != std::end(points); ++it_v, ++idx_v) {
        Filtration_value fil = distance(*it_u, *it_v);
        if (fil <= threshold) {
          edges.emplace_back(fil, idx_u, idx_v);
        }
      }
    }

    // --------------------------------------------------------------------------------------------
    // Creates the proximity graph from edges and sets the property with the filtration value.
    // Number of points is labeled from 0 to idx_u-1
    // --------------------------------------------------------------------------------------------
    // Do not use : rips_skeleton_graph_ = OneSkeletonGraph(...) -> deep copy of the graph (boost graph is not
    // move-enabled)
    // rips_skeleton_graph_.~OneSkeletonGraph();
    // new(&rips_skeleton_graph_)OneSkeletonGraph(edges.begin(), edges.end(), edges_fil.begin(), idx_u);

    // auto vertex_prop = boost::get(vertex_filtration_t(), rips_skeleton_graph_);

    // using vertex_iterator = typename boost::graph_traits<OneSkeletonGraph>::vertex_iterator;
    // vertex_iterator vi, vi_end;
    // for (std::tie(vi, vi_end) = boost::vertices(rips_skeleton_graph_);
    //      vi != vi_end; ++vi) {
    //   boost::put(vertex_prop, *vi, 0.);
    // }
  }

 private:
  // OneSkeletonGraph rips_skeleton_graph_;
  std::vector< std::tuple < Filtration_value, Vertex_handle, Vertex_handle > > edges;
  // std::vector< Filtration_value > edges_fil;
};

}  // namespace rips_complex

}  // namespace Gudhi

#endif  // RIPS_COMPLEX_H_