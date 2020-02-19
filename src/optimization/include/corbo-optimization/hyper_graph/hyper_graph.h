/*********************************************************************
 *
 *  Software License Agreement
 *
 *  Copyright (c) 2020,
 *  TU Dortmund - Institute of Control Theory and Systems Engineering.
 *  All rights reserved.
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 *  Authors: Christoph Rösmann
 *********************************************************************/

#ifndef SRC_OPTIMIZATION_INCLUDE_CORBO_OPTIMIZATION_HYPER_GRAPH_HYPER_GRAPH_H_
#define SRC_OPTIMIZATION_INCLUDE_CORBO_OPTIMIZATION_HYPER_GRAPH_HYPER_GRAPH_H_

#include <corbo-core/types.h>
#include <corbo-optimization/hyper_graph/edge_set.h>
#include <corbo-optimization/hyper_graph/vertex_set.h>

#include <memory>
#include <vector>

namespace corbo {

/**
 * @brief hyper-graph representation
 *
 * @ingroup optimization hyper-graph
 *
 * @see VertexInterface EdgeInterface BaseEdge DiscretizationGrid
 *
 * @author Christoph Rösmann (christoph.roesmann@tu-dortmund.de)
 */
class HyperGraph
{
 public:
    using Ptr = std::shared_ptr<HyperGraph>;

    // enum class EdgeType { Objective, Equality, Inequality };

    HyperGraph()
    {
        setEdgeSet(std::make_shared<OptimizationEdgeSet>());
        setVertexSet(std::make_shared<VertexSet>());
    }

    HyperGraph(OptimizationEdgeSet::Ptr edges, VertexSetInterface::Ptr vertices)
    {
        setEdgeSet(edges);
        setVertexSet(vertices);
    }

    virtual ~HyperGraph() {}

    void setEdgeSet(OptimizationEdgeSet::Ptr edges) { _edges = edges; }
    void setVertexSet(VertexSetInterface::Ptr vertices) { _vertices = vertices; }

    bool hasEdgeSet() const { return (bool)_edges; }
    bool hasVertexSet() const { return (bool)_vertices; }

    OptimizationEdgeSet::Ptr getEdgeSet() const { return _edges; }
    VertexSetInterface::Ptr getVertexSet() const { return _vertices; }

    OptimizationEdgeSet* getEdgeSetRaw() const { return _edges.get(); }
    VertexSetInterface* getVertexSetRaw() const { return _vertices.get(); }

    //! Return number of objective edges
    // int numObjectiveEdges() { return (int)_edges_objective.size(); }
    //! Return number of equality constraint edges
    // int numEqualityEdges() { return (int)_edges_equalities.size(); }
    //! Return number of inequality constraint edges
    // int numInequalityEdges() { return (int)_edges_inequalities.size(); }
    //! Return number of vertices
    // int numVertices() { return (int)_vertices.size(); }

    //! Add a vertex to the graph (without taking ownership)
    // void addVertex(VertexInterface* vtx);
    //! Add objective edge to the graph (and taking ownership)
    // void addObjectiveEdge(EdgeInterface::UPtr edge);
    //! Add equality constraint edge to the graph (and taking ownership)
    // void addEqualityConstraintEdge(EdgeInterface::UPtr edge);
    //! Add inequality constraint edge to the graph (and taking ownership)
    // void addInequalityConstraintEdge(EdgeInterface::UPtr edge);

    //! Precompute the structure of the graph to speed-up later computations

    /**
     * @brief Reset hyper-graph structure
     *
     * Calling this method does only reset flags that capture the state of structure precomputation.
     * Edges and vertices are not touched, otherwise refer to clear() or clearConnectedEdges().
     */
    // void resetStructure();

    //! Determine, if all vertices that are attached to an edge are part of the hyper graph
    bool checkGraphConsistency();

 protected:
    OptimizationEdgeSet::Ptr _edges;  // TODO(roesmann): generalize to EdgeSetInterface?
    VertexSetInterface::Ptr _vertices;

    // internal states
    // bool _edge_obj_idx_precomputed  = false;
    // bool _edge_eq_idx_precomputed   = false;
    // bool _edge_ineq_idx_precomputed = false;
    // bool _vertex_idx_precomputed    = false;
    // int _finite_bounds_precomputed  = false;
};

}  // namespace corbo

#endif  // SRC_OPTIMIZATION_INCLUDE_CORBO_OPTIMIZATION_HYPER_GRAPH_HYPER_GRAPH_H_
