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

#include <corbo-optimization/hyper_graph/edge_set.h>

#include <corbo-optimization/hyper_graph/vertex_set.h>

namespace corbo {

void OptimizationEdgeSet::computeEdgeIndices()
{
    int idx_non_lsq_obj = 0;
    int idx_lsq_obj     = 0;
    int idx_eq          = 0;
    int idx_ineq        = 0;
    computeObjectiveEdgeIndices(_objectives, idx_non_lsq_obj, false);
    computeObjectiveEdgeIndices(_lsq_objectives, idx_lsq_obj, true);
    computeEdgeIndices(_equalities, idx_eq);
    computeEdgeIndices(_inequalities, idx_ineq);
    computeEdgeIndices(_mixed, idx_non_lsq_obj, idx_lsq_obj, idx_eq, idx_ineq);
}

void OptimizationEdgeSet::reserveEdgeCacheMemory(int est_value_cache_size, int est_jacobians_cache_size)
{
    for (BaseEdge::Ptr& edge : _objectives) edge->reserveCacheMemory(est_value_cache_size, est_jacobians_cache_size);
    for (BaseEdge::Ptr& edge : _lsq_objectives) edge->reserveCacheMemory(est_value_cache_size, est_jacobians_cache_size);
    for (BaseEdge::Ptr& edge : _equalities) edge->reserveCacheMemory(est_value_cache_size, est_jacobians_cache_size);
    for (BaseEdge::Ptr& edge : _inequalities) edge->reserveCacheMemory(est_value_cache_size, est_jacobians_cache_size);
    for (BaseMixedEdge::Ptr& edge : _mixed) edge->reserveCacheMemory(est_value_cache_size, est_jacobians_cache_size);
}

bool OptimizationEdgeSet::isEdgeCacheEmpty()
{
    for (BaseEdge::Ptr& edge : _objectives)
    {
        if (edge->getCache().sizeValues() > 0) return false;
        if (edge->getCache().sizeJacobians() > 0) return false;
    }
    for (BaseEdge::Ptr& edge : _lsq_objectives)
    {
        if (edge->getCache().sizeValues() > 0) return false;
        if (edge->getCache().sizeJacobians() > 0) return false;
    }
    for (BaseEdge::Ptr& edge : _equalities)
    {
        if (edge->getCache().sizeValues() > 0) return false;
        if (edge->getCache().sizeJacobians() > 0) return false;
    }
    for (BaseEdge::Ptr& edge : _inequalities)
    {
        if (edge->getCache().sizeValues() > 0) return false;
        if (edge->getCache().sizeJacobians() > 0) return false;
    }
    for (BaseMixedEdge::Ptr& edge : _mixed)
    {
        if (edge->getObjectiveCache().sizeValues() > 0) return false;
        if (edge->getObjectiveCache().sizeJacobians() > 0) return false;
        if (edge->getEqualityCache().sizeValues() > 0) return false;
        if (edge->getEqualityCache().sizeJacobians() > 0) return false;
        if (edge->getInequalityCache().sizeValues() > 0) return false;
        if (edge->getInequalityCache().sizeJacobians() > 0) return false;
    }
    return true;
}

void OptimizationEdgeSet::clearEdgeCache()
{
    for (BaseEdge::Ptr& edge : _objectives) edge->getCache().clear();
    for (BaseEdge::Ptr& edge : _lsq_objectives) edge->getCache().clear();
    for (BaseEdge::Ptr& edge : _equalities) edge->getCache().clear();
    for (BaseEdge::Ptr& edge : _inequalities) edge->getCache().clear();
    for (BaseMixedEdge::Ptr& edge : _mixed)
    {
        edge->getObjectiveCache().clear();
        edge->getEqualityCache().clear();
        edge->getInequalityCache().clear();
    }
}

void OptimizationEdgeSet::computeObjectiveEdgeIndices(std::vector<BaseEdge::Ptr>& edges, int& idx, bool lsq_edges)
{
    if (edges.empty()) return;

    setEdgeIdx(*edges[0], idx);
    int n = edges.size();
    for (int i = 0; i < n; ++i)  // we include the last edge to return the correct idx
    {
        if (lsq_edges)  // (edges[i]->isLeastSquaresForm())
        {
            idx = edges[i]->getEdgeIdx() + edges[i]->getDimension();
        }
        else  // non-lsq objectives can only have a single dimension
        {
            idx = edges[i]->getEdgeIdx() + 1;
        }
        if (i < n - 1) setEdgeIdx(*edges[i + 1], idx);
    }
}

void OptimizationEdgeSet::computeEdgeIndices(std::vector<BaseEdge::Ptr>& edges, int& idx)
{
    if (edges.empty()) return;

    setEdgeIdx(*edges[0], idx);
    int n = edges.size();
    for (int i = 0; i < n; ++i)  // we include the last edge to return the correct idx
    {
        idx = edges[i]->getEdgeIdx() + edges[i]->getDimension();
        if (i < n - 1) setEdgeIdx(*edges[i + 1], idx);
    }
}

void OptimizationEdgeSet::computeEdgeIndices(std::vector<BaseMixedEdge::Ptr>& edges, int& idx_obj, int& idx_lsq_obj, int& idx_eq, int& idx_ineq)
{
    if (edges.empty()) return;

    if (edges[0]->isObjectiveLeastSquaresForm())
    {
        setEdgeIdx(*edges[0], idx_lsq_obj, idx_eq, idx_ineq);
    }
    else
    {
        setEdgeIdx(*edges[0], idx_obj, idx_eq, idx_ineq);
    }

    int n = edges.size();
    for (int i = 0; i < n; ++i)  // we include the last edge to return the correct idx
    {
        if (edges[i]->isObjectiveLeastSquaresForm())
        {
            idx_obj = edges[i]->getEdgeObjectiveIdx() + edges[i]->getObjectiveDimension();  // resulting dimension is always 1 -> f(x)^T * f(x)
        }
        else
        {
            idx_obj = edges[i]->getEdgeObjectiveIdx() + 1;
        }
        idx_eq   = edges[i]->getEdgeEqualityIdx() + edges[i]->getEqualityDimension();
        idx_ineq = edges[i]->getEdgeInequalityIdx() + edges[i]->getInequalityDimension();

        if (i < n - 1)
        {
            setEdgeIdx(*edges[i + 1], edges[i]->isObjectiveLeastSquaresForm() ? idx_lsq_obj : idx_obj, idx_eq, idx_ineq);
        }
    }
}

void OptimizationEdgeSet::registerEdgesAtVertices(VertexSetInterface& vertices)
{
    vertices.clearConnectedEdges();
    registerEdgesAtVertices();
}

void OptimizationEdgeSet::registerEdgesAtVertices()
{
    for (BaseEdge::Ptr& edge : _objectives)
    {
        for (int i = 0; i < edge->getNumVertices(); ++i) edge->getVertexRaw(i)->registerObjectiveEdge(edge.get());
    }
    for (BaseEdge::Ptr& edge : _lsq_objectives)
    {
        for (int i = 0; i < edge->getNumVertices(); ++i) edge->getVertexRaw(i)->registerLsqObjectiveEdge(edge.get());
    }
    for (BaseEdge::Ptr& edge : _equalities)
    {
        for (int i = 0; i < edge->getNumVertices(); ++i) edge->getVertexRaw(i)->registerEqualityEdge(edge.get());
    }
    for (BaseEdge::Ptr& edge : _inequalities)
    {
        for (int i = 0; i < edge->getNumVertices(); ++i) edge->getVertexRaw(i)->registerInequalityEdge(edge.get());
    }
    for (BaseMixedEdge::Ptr& edge : _mixed)
    {
        for (int i = 0; i < edge->getNumVertices(); ++i) edge->getVertexRaw(i)->registerMixedEdge(edge.get());
    }
}

void OptimizationEdgeSet::getDimensions(int& non_lsq_obj_dim, int& lsq_obj_dim, int& eq_dim, int& ineq_dim)
{
    non_lsq_obj_dim = 0;
    lsq_obj_dim     = 0;
    eq_dim          = 0;
    ineq_dim        = 0;

    // plain objectives
    for (BaseEdge::Ptr& edge : _objectives)
    {
        if (edge->getDimension() > 0)
        {
            non_lsq_obj_dim = 1;
            break;
        }
    }
    for (BaseEdge::Ptr& edge : _lsq_objectives)
    {
        lsq_obj_dim += edge->getDimension();
    }
    for (BaseEdge::Ptr& edge : _equalities)
    {
        eq_dim += edge->getDimension();
    }
    for (BaseEdge::Ptr& edge : _inequalities)
    {
        ineq_dim += edge->getDimension();
    }
    for (BaseMixedEdge::Ptr& edge : _mixed)
    {
        if (edge->isObjectiveLeastSquaresForm())
            lsq_obj_dim += edge->getObjectiveDimension();
        else
            non_lsq_obj_dim = 1;

        eq_dim += edge->getEqualityDimension();
        ineq_dim += edge->getInequalityDimension();
    }
}

void OptimizationEdgeSet::clear()
{
    setModified(true);
    _objectives.clear();
    _lsq_objectives.clear();
    _equalities.clear();
    _inequalities.clear();
    _mixed.clear();
}

}  // namespace corbo
