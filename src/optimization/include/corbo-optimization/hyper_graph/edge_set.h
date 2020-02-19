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

#ifndef SRC_OPTIMIZATION_INCLUDE_CORBO_OPTIMIZATION_HYPER_GRAPH_EDGE_SET_H_
#define SRC_OPTIMIZATION_INCLUDE_CORBO_OPTIMIZATION_HYPER_GRAPH_EDGE_SET_H_

#include <corbo-optimization/hyper_graph/edge_interface.h>

#include <corbo-core/console.h>

#include <memory>
#include <vector>

namespace corbo {

/**
 * @brief Abstract class representing a set of edges
 *
 * @ingroup hyper-graph
 *
 * @see HyperGraph EdgeInterface VertexSet
 *
 * @author Christoph Rösmann (christoph.roesmann@tu-dortmund.de)
 */
class EdgeSetInterface
{
 public:
    using Ptr = std::shared_ptr<EdgeSetInterface>;

    //! Virtual destructor
    virtual ~EdgeSetInterface() {}

    //! Precompute edge indices in the hyper-graph (e.g. for the Jacobian structure)
    virtual void computeEdgeIndices() = 0;
    virtual void reserveEdgeCacheMemory(int est_value_cache_size, int est_jacobians_cache_size) = 0;

    virtual void clear() = 0;

    void setModified(bool modified) { _modified = modified; }
    bool isModified() const { return _modified; }

 protected:
    void setEdgeIdx(BaseEdge& edge, int idx) { edge._edge_idx = idx; }
    void setEdgeIdx(BaseMixedEdge& edge, int obj_idx, int eq_idx, int ineq_idx)
    {
        edge._edge_idx_obj  = obj_idx;
        edge._edge_idx_eq   = eq_idx;
        edge._edge_idx_ineq = ineq_idx;
    }
    bool _modified = true;
};

class OptimizationEdgeSet : public EdgeSetInterface
{
 public:
    using Ptr = std::shared_ptr<OptimizationEdgeSet>;

    void clear() override;

    void computeEdgeIndices() override;
    void registerEdgesAtVertices(VertexSetInterface& vertices);
    void registerEdgesAtVertices();

    void getDimensions(int& non_lsq_obj_dim, int& lsq_obj_dim, int& eq_dim, int& ineq_dim);

    void addEdges(std::initializer_list<BaseEdge::Ptr> objective_edges, std::initializer_list<BaseEdge::Ptr> lsq_objective_edges,
                  std::initializer_list<BaseEdge::Ptr> equality_edges, std::initializer_list<BaseEdge::Ptr> inequality_edges,
                  std::initializer_list<BaseMixedEdge::Ptr> mixed_edges)
    {
        setModified(true);
        if (objective_edges.size() > 0) _objectives.insert(_objectives.end(), objective_edges.begin(), objective_edges.end());
        if (lsq_objective_edges.size() > 0) _lsq_objectives.insert(_lsq_objectives.end(), lsq_objective_edges.begin(), lsq_objective_edges.end());
        if (equality_edges.size() > 0) _equalities.insert(_equalities.end(), equality_edges.begin(), equality_edges.end());
        if (inequality_edges.size() > 0) _inequalities.insert(_inequalities.end(), inequality_edges.begin(), inequality_edges.end());
        if (mixed_edges.size() > 0) _mixed.insert(_mixed.end(), mixed_edges.begin(), mixed_edges.end());
    }

    void addEdges(std::initializer_list<BaseEdge::Ptr> objective_edges, std::initializer_list<BaseEdge::Ptr> equality_edges,
                  std::initializer_list<BaseEdge::Ptr> inequality_edges, std::initializer_list<BaseMixedEdge::Ptr> mixed_edges)
    {
        setModified(true);
        if (objective_edges.size() > 0)
        {
            for (auto& edge : objective_edges)
            {
                if (edge->isLeastSquaresForm())
                    _lsq_objectives.push_back(edge);
                else
                    _objectives.push_back(edge);
            }
        }
        if (equality_edges.size() > 0) _equalities.insert(_equalities.end(), equality_edges.begin(), equality_edges.end());
        if (inequality_edges.size() > 0) _inequalities.insert(_inequalities.end(), inequality_edges.begin(), inequality_edges.end());
        if (mixed_edges.size() > 0) _mixed.insert(_mixed.end(), mixed_edges.begin(), mixed_edges.end());
    }

    void addObjectiveEdge(BaseEdge::Ptr edge)
    {
        if (edge->isLeastSquaresForm())
            addLsqObjectiveEdge(edge);
        else
        {
            getObjectiveEdgesRef().push_back(edge);
        }

        // PRINT_ERROR_COND(edge->isLeastSquaresForm(),
        //                 "OptimizationEdgeSet::addLsqObjectiveEdge(): The added edge does return isLeastSquaresForm() == true. You should add it
        //                 "
        //                 "using addLsqObjectiveEdge.");
    }
    void addLsqObjectiveEdge(BaseEdge::Ptr edge)
    {
        PRINT_ERROR_COND(!edge->isLeastSquaresForm(),
                         "OptimizationEdgeSet::addLsqObjectiveEdge(): The added edge does not return isLeastSquaresForm() == true.");
        getLsqObjectiveEdgesRef().push_back(edge);
    }
    void addEqualityEdge(BaseEdge::Ptr edge) { getEqualityEdgesRef().push_back(edge); }
    void addInequalityEdge(BaseEdge::Ptr edge) { getInequalityEdgesRef().push_back(edge); }
    void addMixedEdge(BaseMixedEdge::Ptr edge) { getMixedEdgesRef().push_back(edge); }

    std::vector<BaseEdge::Ptr>& getObjectiveEdgesRef()
    {
        setModified(true);
        return _objectives;
    }
    const std::vector<BaseEdge::Ptr>& getObjectiveEdges() const { return _objectives; }
    std::vector<BaseEdge::Ptr>& getLsqObjectiveEdgesRef()
    {
        setModified(true);
        return _lsq_objectives;
    }
    const std::vector<BaseEdge::Ptr>& getLsqObjectiveEdges() const { return _lsq_objectives; }
    std::vector<BaseEdge::Ptr>& getEqualityEdgesRef()
    {
        setModified(true);
        return _equalities;
    }
    const std::vector<BaseEdge::Ptr>& getEqualityEdges() const { return _equalities; }
    std::vector<BaseEdge::Ptr>& getInequalityEdgesRef()
    {
        setModified(true);
        return _inequalities;
    }
    const std::vector<BaseEdge::Ptr>& getInequalityEdges() const { return _inequalities; }
    std::vector<BaseMixedEdge::Ptr>& getMixedEdgesRef()
    {
        setModified(true);
        return _mixed;
    }
    const std::vector<BaseMixedEdge::Ptr>& getMixedEdges() const { return _mixed; }

    bool hasOnlyLeastSquaresObjectives() const { return _objectives.empty(); }

    void reserveEdgeCacheMemory(int est_value_cache_size, int est_jacobians_cache_size) override;

    bool isEdgeCacheEmpty();
    void clearEdgeCache();

 protected:
    //! Precompute overall edge indices in the hyper-graph (e.g. for the Jacobian structure)
    void computeObjectiveEdgeIndices(std::vector<BaseEdge::Ptr>& edges, int& idx, bool lsq_edges);
    void computeEdgeIndices(std::vector<BaseEdge::Ptr>& edges, int& idx);
    void computeEdgeIndices(std::vector<BaseMixedEdge::Ptr>& edges, int& idx_obj, int& idx_lsq_obj, int& idx_eq, int& idx_ineq);

 private:
    std::vector<BaseEdge::Ptr> _objectives;
    std::vector<BaseEdge::Ptr> _lsq_objectives;
    std::vector<BaseEdge::Ptr> _equalities;
    std::vector<BaseEdge::Ptr> _inequalities;
    std::vector<BaseMixedEdge::Ptr> _mixed;
};

}  // namespace corbo

#endif  // SRC_OPTIMIZATION_INCLUDE_CORBO_OPTIMIZATION_HYPER_GRAPH_EDGE_SET_H_
