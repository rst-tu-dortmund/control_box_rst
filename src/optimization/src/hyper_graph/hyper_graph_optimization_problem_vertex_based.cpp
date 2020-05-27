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

#include <corbo-optimization/hyper_graph/hyper_graph_optimization_problem_vertex_based.h>

namespace corbo {

static constexpr const double cd_delta     = 1e-9;
static constexpr const double cd_neg2delta = -2 * cd_delta;
static constexpr const double cd_scalar    = 1.0 / (2 * cd_delta);

void HyperGraphOptimizationProblemVertexBased::precomputeEdgeQuantities()
{
    assert(_graph.hasEdgeSet());

    OptimizationEdgeSet::Ptr edges = _graph.getEdgeSet();
    if (edges->isModified())
    {
        edges->getDimensions(_dim_non_lsq_obj, _dim_lsq_obj, _dim_eq, _dim_ineq);
        edges->computeEdgeIndices();
        edges->reserveEdgeCacheMemory(3, 2);  // max 3 values (central differences plus current values; max 2 jacobians (for hessian)
        // we make use of the internal edge cache here; TODO(roesmann) we might not iterate over previously allocated
        // edges
        edges->registerEdgesAtVertices(*_graph.getVertexSet());  // we also need to register edges to later iterate vertices
        edges->setModified(false);
    }
}

void HyperGraphOptimizationProblemVertexBased::precomputeConnectedMixedEdges(const VertexInterface* vertex, bool objective, bool equality,
                                                                             bool inequality)
{
    for (BaseMixedEdge* edge : vertex->getConnectedMixedEdgesRef())
    {
        if (!equality && !inequality && edge->getObjectiveDimension() == 0) continue;
        if (!objective && !inequality && edge->getEqualityDimension() == 0) continue;
        if (!equality && !objective && edge->getInequalityDimension() == 0) continue;
        // TODO(roesmann): these are not all cases

        edge->precompute();
    }
}

void HyperGraphOptimizationProblemVertexBased::computeObjectiveValuesCached(const VertexInterface* vertex, bool include_lsq_edges,
                                                                            bool include_nonmixed, bool include_mixed)
{
    if (include_nonmixed)
    {
        for (BaseEdge* edge : vertex->getConnectedObjectiveEdgesRef())
        {
            edge->computeValuesCached();
        }
        if (include_lsq_edges)
        {
            for (BaseEdge* edge : vertex->getConnectedLsqObjectiveEdgesRef())
            {
                edge->computeSquaredNormOfValuesCached();
            }
        }
    }
    if (include_mixed)
    {
        for (BaseMixedEdge* edge : vertex->getConnectedMixedEdgesRef())
        {
            if (edge->getObjectiveDimension() == 0) continue;
            if (!include_lsq_edges && edge->isObjectiveLeastSquaresForm()) continue;

            edge->precompute();

            if (edge->isObjectiveLeastSquaresForm())
                edge->computeSquaredNormOfObjectiveValuesCached();
            else
                edge->computeObjectiveValuesCached();
        }
    }
}

void HyperGraphOptimizationProblemVertexBased::finalizeObjectiveGradient(const VertexInterface* vertex, double& gradient_coeff,
                                                                         bool include_lsq_edges, bool include_nonmixed, bool include_mixed,
                                                                         bool precompute_mixed)
{
    if (include_nonmixed)
    {
        for (BaseEdge* edge : vertex->getConnectedObjectiveEdgesRef())
        {
            Eigen::VectorXd values_minus(edge->getDimension());
            edge->computeValues(values_minus);
            gradient_coeff += cd_scalar * (edge->getCache().topValues() - values_minus).sum();
            edge->getCache().popValues();
        }
        if (include_lsq_edges)
        {
            for (BaseEdge* edge : vertex->getConnectedLsqObjectiveEdgesRef())
            {
                gradient_coeff += cd_scalar * (edge->getCache().topValues()[0] - edge->computeSquaredNormOfValues());
                edge->getCache().popValues();
            }
        }
    }
    if (include_mixed)
    {
        for (BaseMixedEdge* edge : vertex->getConnectedMixedEdgesRef())
        {
            if (edge->getObjectiveDimension() == 0) continue;
            if (!include_lsq_edges && edge->isObjectiveLeastSquaresForm()) continue;

            if (precompute_mixed) edge->precompute();

            if (edge->isObjectiveLeastSquaresForm())
            {
                gradient_coeff += cd_scalar * (edge->getObjectiveCache().topValues()[0] - edge->computeSquaredNormOfObjectiveValues());
            }
            else
            {
                Eigen::VectorXd values_minus(edge->getObjectiveDimension());
                edge->computeObjectiveValues(values_minus);
                gradient_coeff += cd_scalar * (edge->getObjectiveCache().topValues() - values_minus).sum();
            }
            edge->getObjectiveCache().popValues();
        }
    }
}

void HyperGraphOptimizationProblemVertexBased::computeLsqObjectiveValuesCached(const VertexInterface* vertex, bool include_nonmixed,
                                                                               bool include_mixed, bool precompute_mixed)
{
    if (include_nonmixed)
    {
        for (BaseEdge* edge : vertex->getConnectedLsqObjectiveEdgesRef())
        {
            edge->computeValuesCached();
        }
    }
    if (include_mixed)
    {
        for (BaseMixedEdge* edge : vertex->getConnectedMixedEdgesRef())
        {
            if (edge->getObjectiveDimension() == 0) continue;
            if (!edge->isObjectiveLeastSquaresForm()) continue;

            if (precompute_mixed) edge->precompute();

            edge->computeObjectiveValuesCached();
        }
    }
}

void HyperGraphOptimizationProblemVertexBased::finalizeLsqObjectiveJacobian(const VertexInterface* vertex, int vtx_idx,
                                                                            Eigen::Ref<Eigen::MatrixXd>& jacobian, const double* multipliers,
                                                                            bool include_nonmixed, bool include_mixed)
{
    if (include_nonmixed)
    {
        for (BaseEdge* edge : vertex->getConnectedLsqObjectiveEdgesRef())
        {
            Eigen::VectorXd values_minus(edge->getDimension());
            edge->computeValues(values_minus);

            if (multipliers)
            {
                Eigen::Map<const Eigen::VectorXd> multipliers_map(multipliers + edge->getEdgeIdx(), edge->getDimension());
                jacobian.block(edge->getEdgeIdx(), vtx_idx, edge->getDimension(), 1) =
                    cd_scalar * multipliers_map.cwiseProduct(edge->getCache().topValues() - values_minus);
            }
            else
            {
                jacobian.block(edge->getEdgeIdx(), vtx_idx, edge->getDimension(), 1) = cd_scalar * (edge->getCache().topValues() - values_minus);
            }
            edge->getCache().popValues();
        }
    }

    if (include_mixed)
    {
        for (BaseMixedEdge* edge : vertex->getConnectedMixedEdgesRef())
        {
            if (edge->getObjectiveDimension() == 0) continue;
            if (!edge->isObjectiveLeastSquaresForm()) continue;

            edge->precompute();

            Eigen::VectorXd values_minus(edge->getObjectiveDimension());
            edge->computeObjectiveValues(values_minus);
            if (multipliers)
            {
                Eigen::Map<const Eigen::VectorXd> multipliers_map(multipliers + edge->getEdgeObjectiveIdx(), edge->getObjectiveDimension());
                jacobian.block(edge->getEdgeObjectiveIdx(), vtx_idx, edge->getObjectiveDimension(), 1).noalias() =
                    cd_scalar * multipliers_map.cwiseProduct(edge->getObjectiveCache().topValues() - values_minus);
            }
            else
            {
                jacobian.block(edge->getEdgeObjectiveIdx(), vtx_idx, edge->getObjectiveDimension(), 1).noalias() =
                    cd_scalar * (edge->getObjectiveCache().topValues() - values_minus);
            }
            edge->getObjectiveCache().popValues();
        }
    }
}

void HyperGraphOptimizationProblemVertexBased::computeLsqObjectiveJacobianStructureForVertex(const VertexInterface* vertex, int vtx_idx,
                                                                                             Eigen::Ref<Eigen::VectorXi> i_row,
                                                                                             Eigen::Ref<Eigen::VectorXi> j_col, int& nnz_idx,
                                                                                             int row_offset)
{
    for (BaseEdge* edge : vertex->getConnectedLsqObjectiveEdgesRef())
    {
        const int edge_idx = edge->getEdgeIdx() + row_offset;

        for (int i = 0; i < edge->getDimension(); ++i)
        {
            i_row[nnz_idx] = edge_idx + i;
            j_col[nnz_idx] = vtx_idx;
            ++nnz_idx;
        }
    }

    for (BaseMixedEdge* edge : vertex->getConnectedMixedEdgesRef())
    {
        if (edge->getObjectiveDimension() == 0) continue;
        if (!edge->isObjectiveLeastSquaresForm()) continue;

        const int edge_idx = edge->getEdgeObjectiveIdx() + row_offset;

        for (int i = 0; i < edge->getObjectiveDimension(); ++i)
        {
            i_row[nnz_idx] = edge_idx + i;
            j_col[nnz_idx] = vtx_idx;
            ++nnz_idx;
        }
    }
}

void HyperGraphOptimizationProblemVertexBased::finalizeLsqObjectiveJacobianSparseValues(const VertexInterface* vertex, int& nnz_idx,
                                                                                        Eigen::Ref<Eigen::VectorXd>& values,
                                                                                        const double* multipliers, bool precompute_mixed)
{
    for (BaseEdge* edge : vertex->getConnectedLsqObjectiveEdgesRef())
    {
        Eigen::VectorXd values_minus(edge->getDimension());
        edge->computeValues(values_minus);

        if (multipliers)
        {
            Eigen::Map<const Eigen::VectorXd> multipliers_map(multipliers + edge->getEdgeIdx(), edge->getDimension());
            values.segment(nnz_idx, edge->getDimension()) = cd_scalar * multipliers_map.cwiseProduct(edge->getCache().topValues() - values_minus);
            nnz_idx += edge->getDimension();
        }
        else
        {
            values.segment(nnz_idx, edge->getDimension()) = cd_scalar * (edge->getCache().topValues() - values_minus);
            nnz_idx += edge->getDimension();
        }
        edge->getCache().popValues();
    }

    for (BaseMixedEdge* edge : vertex->getConnectedMixedEdgesRef())
    {
        if (edge->getObjectiveDimension() == 0) continue;
        if (!edge->isObjectiveLeastSquaresForm()) continue;

        if (precompute_mixed) edge->precompute();

        Eigen::VectorXd values_minus(edge->getObjectiveDimension());
        edge->computeObjectiveValues(values_minus);
        if (multipliers)
        {
            Eigen::Map<const Eigen::VectorXd> multipliers_map(multipliers + edge->getEdgeObjectiveIdx(), edge->getObjectiveDimension());
            values.segment(nnz_idx, edge->getObjectiveDimension()).noalias() =
                cd_scalar * multipliers_map.cwiseProduct(edge->getObjectiveCache().topValues() - values_minus);
            nnz_idx += edge->getObjectiveDimension();
        }
        else
        {
            values.segment(nnz_idx, edge->getObjectiveDimension()).noalias() = cd_scalar * (edge->getObjectiveCache().topValues() - values_minus);
            nnz_idx += edge->getObjectiveDimension();
        }
        edge->getObjectiveCache().popValues();
    }
}

void HyperGraphOptimizationProblemVertexBased::computeEqualitiesValuesCached(const VertexInterface* vertex, bool include_nonmixed, bool include_mixed,
                                                                             bool precompute_mixed)
{
    if (include_nonmixed)
    {
        for (BaseEdge* edge : vertex->getConnectedEqualityEdgesRef())
        {
            edge->computeValuesCached();
        }
    }
    if (include_mixed)
    {
        for (BaseMixedEdge* edge : vertex->getConnectedMixedEdgesRef())
        {
            if (edge->getEqualityDimension() == 0) continue;

            if (precompute_mixed) edge->precompute();

            edge->computeEqualityValuesCached();
        }
    }
}

void HyperGraphOptimizationProblemVertexBased::finalizeEqualitiesJacobian(const VertexInterface* vertex, int vtx_idx,
                                                                          Eigen::Ref<Eigen::MatrixXd>& jacobian, const double* multipliers,
                                                                          bool include_nonmixed, bool include_mixed, bool precompute_mixed)
{
    if (include_nonmixed)
    {
        for (BaseEdge* edge : vertex->getConnectedEqualityEdgesRef())
        {
            Eigen::VectorXd values_minus(edge->getDimension());
            edge->computeValues(values_minus);

            if (multipliers)
            {
                Eigen::Map<const Eigen::VectorXd> multipliers_map(multipliers + edge->getEdgeIdx(), edge->getDimension());
                jacobian.block(edge->getEdgeIdx(), vtx_idx, edge->getDimension(), 1) =
                    cd_scalar * multipliers_map.cwiseProduct(edge->getCache().topValues() - values_minus);
            }
            else
            {
                jacobian.block(edge->getEdgeIdx(), vtx_idx, edge->getDimension(), 1) = cd_scalar * (edge->getCache().topValues() - values_minus);
            }
            edge->getCache().popValues();
        }
    }

    if (include_mixed)
    {
        for (BaseMixedEdge* edge : vertex->getConnectedMixedEdgesRef())
        {
            if (edge->getEqualityDimension() == 0) continue;

            if (precompute_mixed) edge->precompute();

            Eigen::VectorXd values_minus(edge->getEqualityDimension());
            edge->computeEqualityValues(values_minus);
            if (multipliers)
            {
                Eigen::Map<const Eigen::VectorXd> multipliers_map(multipliers + edge->getEdgeEqualityIdx(), edge->getEqualityDimension());
                jacobian.block(edge->getEdgeEqualityIdx(), vtx_idx, edge->getEqualityDimension(), 1).noalias() =
                    cd_scalar * multipliers_map.cwiseProduct(edge->getEqualityCache().topValues() - values_minus);
            }
            else
            {
                jacobian.block(edge->getEdgeEqualityIdx(), vtx_idx, edge->getEqualityDimension(), 1).noalias() =
                    cd_scalar * (edge->getEqualityCache().topValues() - values_minus);
            }
            edge->getEqualityCache().popValues();
        }
    }
}

void HyperGraphOptimizationProblemVertexBased::computeEqualitiesJacobianStructureForVertex(const VertexInterface* vertex, int vtx_idx,
                                                                                           Eigen::Ref<Eigen::VectorXi> i_row,
                                                                                           Eigen::Ref<Eigen::VectorXi> j_col, int& nnz_idx,
                                                                                           int row_offset)
{
    for (BaseEdge* edge : vertex->getConnectedEqualityEdgesRef())
    {
        const int edge_idx = edge->getEdgeIdx() + row_offset;

        for (int i = 0; i < edge->getDimension(); ++i)
        {
            i_row[nnz_idx] = edge_idx + i;
            j_col[nnz_idx] = vtx_idx;
            ++nnz_idx;
        }
    }

    for (BaseMixedEdge* edge : vertex->getConnectedMixedEdgesRef())
    {
        if (edge->getEqualityDimension() == 0) continue;

        const int edge_idx = edge->getEdgeEqualityIdx() + row_offset;

        for (int i = 0; i < edge->getEqualityDimension(); ++i)
        {
            i_row[nnz_idx] = edge_idx + i;
            j_col[nnz_idx] = vtx_idx;
            ++nnz_idx;
        }
    }
}

void HyperGraphOptimizationProblemVertexBased::finalizeEqualitiesJacobianSparseValues(const VertexInterface* vertex, int& nnz_idx,
                                                                                      Eigen::Ref<Eigen::VectorXd>& values, const double* multipliers,
                                                                                      bool precompute_mixed)
{
    for (BaseEdge* edge : vertex->getConnectedEqualityEdgesRef())
    {
        Eigen::VectorXd values_minus(edge->getDimension());
        edge->computeValues(values_minus);

        if (multipliers)
        {
            Eigen::Map<const Eigen::VectorXd> multipliers_map(multipliers + edge->getEdgeIdx(), edge->getDimension());
            values.segment(nnz_idx, edge->getDimension()) = cd_scalar * multipliers_map.cwiseProduct(edge->getCache().topValues() - values_minus);
            nnz_idx += edge->getDimension();
        }
        else
        {
            values.segment(nnz_idx, edge->getDimension()) = cd_scalar * (edge->getCache().topValues() - values_minus);
            nnz_idx += edge->getDimension();
        }
        edge->getCache().popValues();
    }

    for (BaseMixedEdge* edge : vertex->getConnectedMixedEdgesRef())
    {
        if (edge->getEqualityDimension() == 0) continue;

        if (precompute_mixed) edge->precompute();

        Eigen::VectorXd values_minus(edge->getEqualityDimension());
        edge->computeEqualityValues(values_minus);
        if (multipliers)
        {
            Eigen::Map<const Eigen::VectorXd> multipliers_map(multipliers + edge->getEdgeEqualityIdx(), edge->getEqualityDimension());
            values.segment(nnz_idx, edge->getEqualityDimension()).noalias() =
                cd_scalar * multipliers_map.cwiseProduct(edge->getEqualityCache().topValues() - values_minus);
            nnz_idx += edge->getEqualityDimension();
        }
        else
        {
            values.segment(nnz_idx, edge->getEqualityDimension()).noalias() = cd_scalar * (edge->getEqualityCache().topValues() - values_minus);
            nnz_idx += edge->getEqualityDimension();
        }
        edge->getEqualityCache().popValues();
    }
}

void HyperGraphOptimizationProblemVertexBased::computeInequalitiesValuesCached(const VertexInterface* vertex, bool include_nonmixed,
                                                                               bool include_mixed, bool precompute_mixed)
{
    if (include_nonmixed)
    {
        for (BaseEdge* edge : vertex->getConnectedInequalityEdgesRef())
        {
            edge->computeValuesCached();
        }
    }
    if (include_mixed)
    {
        for (BaseMixedEdge* edge : vertex->getConnectedMixedEdgesRef())
        {
            if (edge->getInequalityDimension() == 0) continue;

            if (precompute_mixed) edge->precompute();

            edge->computeInequalityValuesCached();
        }
    }
}

void HyperGraphOptimizationProblemVertexBased::finalizeInequalitiesJacobian(const VertexInterface* vertex, int vtx_idx,
                                                                            Eigen::Ref<Eigen::MatrixXd>& jacobian, const double* multipliers,
                                                                            bool include_nonmixed, bool include_mixed, bool precompute_mixed)
{
    if (include_nonmixed)
    {
        for (BaseEdge* edge : vertex->getConnectedInequalityEdgesRef())
        {
            Eigen::VectorXd values_minus(edge->getDimension());
            edge->computeValues(values_minus);

            if (multipliers)
            {
                Eigen::Map<const Eigen::VectorXd> multipliers_map(multipliers + edge->getEdgeIdx(), edge->getDimension());
                jacobian.block(edge->getEdgeIdx(), vtx_idx, edge->getDimension(), 1) =
                    cd_scalar * multipliers_map.cwiseProduct(edge->getCache().topValues() - values_minus);
            }
            else
            {
                jacobian.block(edge->getEdgeIdx(), vtx_idx, edge->getDimension(), 1) = cd_scalar * (edge->getCache().topValues() - values_minus);
            }
            edge->getCache().popValues();
        }
    }

    if (include_mixed)
    {
        for (BaseMixedEdge* edge : vertex->getConnectedMixedEdgesRef())
        {
            if (edge->getInequalityDimension() == 0) continue;

            if (precompute_mixed) edge->precompute();

            Eigen::VectorXd values_minus(edge->getInequalityDimension());
            edge->computeInequalityValues(values_minus);
            if (multipliers)
            {
                Eigen::Map<const Eigen::VectorXd> multipliers_map(multipliers + edge->getEdgeInequalityIdx(), edge->getInequalityDimension());
                jacobian.block(edge->getEdgeInequalityIdx(), vtx_idx, edge->getInequalityDimension(), 1).noalias() =
                    cd_scalar * multipliers_map.cwiseProduct(edge->getInequalityCache().topValues() - values_minus);
            }
            else
            {
                jacobian.block(edge->getEdgeInequalityIdx(), vtx_idx, edge->getInequalityDimension(), 1).noalias() =
                    cd_scalar * (edge->getInequalityCache().topValues() - values_minus);
            }
            edge->getInequalityCache().popValues();
        }
    }
}

void HyperGraphOptimizationProblemVertexBased::finalizeActiveInequalitiesJacobian(const VertexInterface* vertex, int vtx_idx,
                                                                                  Eigen::Ref<Eigen::MatrixXd>& jacobian, double weight,
                                                                                  bool include_nonmixed, bool include_mixed)
{
    const double weighted_scalar = weight * cd_scalar;

    if (include_nonmixed)
    {
        for (BaseEdge* edge : vertex->getConnectedInequalityEdgesRef())
        {
            assert(edge->getCache().sizeValues() > 1);  // we need both the actual values and the intermediate result from finite differences

            Eigen::Array<bool, -1, 1> active = edge->getCache().recentValues(EdgeCache::Previous).array() > 0.0;

            if (active.any())
            {
                Eigen::VectorXd values_minus(edge->getDimension());
                edge->computeValues(values_minus);

                // now set values to zero if inactive or multiply with weight otherwise
                for (int j = 0; j < edge->getDimension(); ++j)
                {
                    if (active[j]) jacobian(edge->getEdgeIdx() + j, vtx_idx) = weighted_scalar * (edge->getCache().topValues()(j) - values_minus(j));
                    // else if (!active[j])
                    // jacobian(edge->getEdgeIdx() + j, vtx_idx) = 0.0;  // should be zero already
                }
            }
            edge->getCache().popValues();  // intermediate result
            edge->getCache().popValues();  // values before finite differences
        }
    }

    if (include_mixed)
    {
        for (BaseMixedEdge* edge : vertex->getConnectedMixedEdgesRef())
        {
            if (edge->getInequalityDimension() == 0) continue;

            assert(edge->getInequalityCache().sizeValues() >
                   1);  // we need both the actual values and the intermediate result from finite differences

            Eigen::Array<bool, -1, 1> active = edge->getInequalityCache().recentValues(EdgeCache::Previous).array() > 0.0;

            if (active.any())
            {
                edge->precompute();

                Eigen::VectorXd values_minus(edge->getInequalityDimension());
                edge->computeInequalityValues(values_minus);

                // now set values to zero if inactive or multiply with weight otherwise
                for (int j = 0; j < edge->getInequalityDimension(); ++j)
                {
                    if (active[j])
                        jacobian(edge->getEdgeInequalityIdx() + j, vtx_idx) =
                            weighted_scalar * (edge->getInequalityCache().topValues()(j) - values_minus(j));
                    // else if (!active[j])
                    // jacobian(edge->getEdgeIdx() + j, vtx_idx) = 0.0;  // should be zero already
                }
            }
            edge->getInequalityCache().popValues();  // intermediate result
            edge->getInequalityCache().popValues();  // values before finite differences
        }
    }
}

void HyperGraphOptimizationProblemVertexBased::computeInequalitiesJacobianStructureForVertex(const VertexInterface* vertex, int vtx_idx,
                                                                                             Eigen::Ref<Eigen::VectorXi> i_row,
                                                                                             Eigen::Ref<Eigen::VectorXi> j_col, int& nnz_idx,
                                                                                             int row_offset)
{
    for (BaseEdge* edge : vertex->getConnectedInequalityEdgesRef())
    {
        const int edge_idx = edge->getEdgeIdx() + row_offset;

        for (int i = 0; i < edge->getDimension(); ++i)
        {
            i_row[nnz_idx] = edge_idx + i;
            j_col[nnz_idx] = vtx_idx;
            ++nnz_idx;
        }
    }

    for (BaseMixedEdge* edge : vertex->getConnectedMixedEdgesRef())
    {
        if (edge->getInequalityDimension() == 0) continue;

        const int edge_idx = edge->getEdgeInequalityIdx() + row_offset;

        for (int i = 0; i < edge->getInequalityDimension(); ++i)
        {
            i_row[nnz_idx] = edge_idx + i;
            j_col[nnz_idx] = vtx_idx;
            ++nnz_idx;
        }
    }
}

void HyperGraphOptimizationProblemVertexBased::finalizeInequalitiesJacobianSparseValues(const VertexInterface* vertex, int& nnz_idx,
                                                                                        Eigen::Ref<Eigen::VectorXd>& values,
                                                                                        const double* multipliers, bool precompute_mixed)
{
    for (BaseEdge* edge : vertex->getConnectedInequalityEdgesRef())
    {
        Eigen::VectorXd values_minus(edge->getDimension());
        edge->computeValues(values_minus);

        if (multipliers)
        {
            Eigen::Map<const Eigen::VectorXd> multipliers_map(multipliers + edge->getEdgeIdx(), edge->getDimension());
            values.segment(nnz_idx, edge->getDimension()) = cd_scalar * multipliers_map.cwiseProduct(edge->getCache().topValues() - values_minus);
            edge->getCache().popValues();
            nnz_idx += edge->getDimension();
        }
        else
        {
            values.segment(nnz_idx, edge->getDimension()) = cd_scalar * (edge->getCache().topValues() - values_minus);
            edge->getCache().popValues();
            nnz_idx += edge->getDimension();
        }
    }

    for (BaseMixedEdge* edge : vertex->getConnectedMixedEdgesRef())
    {
        if (edge->getInequalityDimension() == 0) continue;

        if (precompute_mixed) edge->precompute();

        Eigen::VectorXd values_minus(edge->getInequalityDimension());
        edge->computeInequalityValues(values_minus);
        if (multipliers)
        {
            Eigen::Map<const Eigen::VectorXd> multipliers_map(multipliers + edge->getEdgeInequalityIdx(), edge->getInequalityDimension());
            values.segment(nnz_idx, edge->getInequalityDimension()).noalias() =
                cd_scalar * multipliers_map.cwiseProduct(edge->getInequalityCache().topValues() - values_minus);
            edge->getInequalityCache().popValues();
            nnz_idx += edge->getInequalityDimension();
        }
        else
        {
            values.segment(nnz_idx, edge->getInequalityDimension()).noalias() = cd_scalar * (edge->getInequalityCache().topValues() - values_minus);
            edge->getInequalityCache().popValues();
            nnz_idx += edge->getInequalityDimension();
        }
    }
}

void HyperGraphOptimizationProblemVertexBased::finalizeActiveInequalitiesJacobianSparseValues(const VertexInterface* vertex, int& nnz_idx,
                                                                                              Eigen::Ref<Eigen::VectorXd>& values, double weight)
{
    const double weighted_scalar = weight * cd_scalar;

    for (BaseEdge* edge : vertex->getConnectedInequalityEdgesRef())
    {
        assert(edge->getCache().sizeValues() > 1);  // we need both the actual values and the intermediate result from finite differences

        Eigen::Array<bool, -1, 1> active = edge->getCache().recentValues(EdgeCache::Previous).array() > 0.0;

        if (active.any())
        {
            Eigen::VectorXd values_minus(edge->getDimension());
            edge->computeValues(values_minus);

            // now set values to zero if inactive or multiply with weight otherwise
            for (int j = 0; j < edge->getDimension(); ++j)
            {
                if (active[j]) values(nnz_idx) = weighted_scalar * (edge->getCache().topValues()(j) - values_minus(j));
                ++nnz_idx;
            }
        }
        else
        {
            nnz_idx += edge->getDimension();
        }
        edge->getCache().popValues();  // intermediate result
        edge->getCache().popValues();  // values before finite differences
    }

    for (BaseMixedEdge* edge : vertex->getConnectedMixedEdgesRef())
    {
        if (edge->getInequalityDimension() == 0) continue;

        assert(edge->getInequalityCache().sizeValues() > 1);  // we need both the actual values and the intermediate result from finite differences

        Eigen::Array<bool, -1, 1> active = edge->getInequalityCache().recentValues(EdgeCache::Previous).array() > 0.0;

        if (active.any())
        {
            edge->precompute();

            Eigen::VectorXd values_minus(edge->getInequalityDimension());
            edge->computeInequalityValues(values_minus);

            // now set values to zero if inactive or multiply with weight otherwise
            for (int j = 0; j < edge->getInequalityDimension(); ++j)
            {
                if (active[j]) values(nnz_idx) = weighted_scalar * (edge->getInequalityCache().topValues()(j) - values_minus(j));
                ++nnz_idx;
            }
        }
        else
        {
            nnz_idx += edge->getInequalityDimension();
        }
        edge->getInequalityCache().popValues();  // intermediate result
        edge->getInequalityCache().popValues();  // values before finite differences
    }
}

void HyperGraphOptimizationProblemVertexBased::computeGradientObjective(Eigen::Ref<Eigen::VectorXd> gradient)
{
    assert(_graph.hasVertexSet() && _graph.hasEdgeSet());
    assert(getGraph().getEdgeSet()->isEdgeCacheEmpty());
    if (!_graph_precomputed) precomputeGraphQuantities();

    assert(gradient.size() == getParameterDimension());

    gradient.setZero();

    for (VertexInterface* vertex : _graph.getVertexSet()->getActiveVertices())
    {
        if (vertex->isFixed()) continue;

        int free_idx = 0;
        for (int i = 0; i < vertex->getDimension(); ++i)
        {
            if (vertex->isFixedComponent(i)) continue;

            vertex->plus(i, cd_delta);

            computeObjectiveValuesCached(vertex, true);

            vertex->plus(i, cd_neg2delta);

            finalizeObjectiveGradient(vertex, gradient(vertex->getVertexIdx() + free_idx), true);

            vertex->plus(i, cd_delta);  // revert offset
            ++free_idx;
        }
    }
    assert(getGraph().getEdgeSet()->isEdgeCacheEmpty());
}

void HyperGraphOptimizationProblemVertexBased::computeGradientNonLsqObjective(Eigen::Ref<Eigen::VectorXd> gradient)
{
    assert(_graph.hasVertexSet() && _graph.hasEdgeSet());
    assert(getGraph().getEdgeSet()->isEdgeCacheEmpty());
    if (!_graph_precomputed) precomputeGraphQuantities();

    assert(gradient.size() == getParameterDimension());

    for (VertexInterface* vertex : _graph.getVertexSet()->getActiveVertices())
    {
        if (vertex->isFixed()) continue;

        int free_idx = 0;
        for (int i = 0; i < vertex->getDimension(); ++i)
        {
            if (vertex->isFixedComponent(i)) continue;

            vertex->plus(i, cd_delta);

            computeObjectiveValuesCached(vertex, false);

            vertex->plus(i, cd_neg2delta);

            finalizeObjectiveGradient(vertex, gradient(vertex->getVertexIdx() + free_idx), false);

            vertex->plus(i, cd_delta);  // revert offset
            ++free_idx;
        }
    }
    assert(getGraph().getEdgeSet()->isEdgeCacheEmpty());
}

void HyperGraphOptimizationProblemVertexBased::computeDenseJacobianLsqObjective(Eigen::Ref<Eigen::MatrixXd> jacobian, const double* multipliers)
{
    assert(_graph.hasVertexSet() && _graph.hasEdgeSet());
    assert(getGraph().getEdgeSet()->isEdgeCacheEmpty());
    if (!_graph_precomputed) precomputeGraphQuantities();

    assert(jacobian.rows() == getLsqObjectiveDimension());
    assert(jacobian.cols() == getParameterDimension());

    // we need to set to zero everything first
    jacobian.setZero();

    for (VertexInterface* vertex : _graph.getVertexSet()->getActiveVertices())
    {
        if (vertex->isFixed()) continue;

        int free_idx = 0;
        for (int i = 0; i < vertex->getDimension(); ++i)
        {
            if (vertex->isFixedComponent(i)) continue;

            vertex->plus(i, cd_delta);

            computeLsqObjectiveValuesCached(vertex);

            vertex->plus(i, cd_neg2delta);

            finalizeLsqObjectiveJacobian(vertex, vertex->getVertexIdx() + free_idx, jacobian, multipliers);

            vertex->plus(i, cd_delta);  // revert offset
            ++free_idx;
        }
    }
    assert(getGraph().getEdgeSet()->isEdgeCacheEmpty());
}

int HyperGraphOptimizationProblemVertexBased::computeSparseJacobianLsqObjectiveNNZ()
{
    int nnz = 0;

    assert(_graph.hasEdgeSet());
    if (!_graph_precomputed) precomputeGraphQuantities();

    OptimizationEdgeSet::Ptr edges = _graph.getEdgeSet();

    // Iterate lsq objective edges
    for (BaseEdge::Ptr& edge : edges->getLsqObjectiveEdgesRef())
    {
        for (int i = 0; i < edge->getNumVertices(); ++i)
        {
            nnz += edge->getDimension() * edge->getVertexRaw(i)->getDimensionUnfixed();  // block size
        }
    }
    // Iterate mixed edges
    for (BaseMixedEdge::Ptr& edge : edges->getMixedEdgesRef())
    {
        if (edge->getObjectiveDimension() == 0 || !edge->isObjectiveLeastSquaresForm()) continue;

        for (int i = 0; i < edge->getNumVertices(); ++i)
        {
            nnz += edge->getObjectiveDimension() * edge->getVertexRaw(i)->getDimensionUnfixed();  // block size
        }
    }
    return nnz;
}

void HyperGraphOptimizationProblemVertexBased::computeSparseJacobianLsqObjectiveStructure(Eigen::Ref<Eigen::VectorXi> i_row,
                                                                                          Eigen::Ref<Eigen::VectorXi> j_col)
{
    assert(i_row.size() == computeSparseJacobianLsqObjectiveNNZ());
    assert(j_col.size() == i_row.size());

    int nnz_idx = 0;

    for (VertexInterface* vertex : _graph.getVertexSet()->getActiveVertices())
    {
        if (vertex->isFixed()) continue;

        int free_idx = 0;
        for (int i = 0; i < vertex->getDimension(); ++i)
        {
            if (vertex->isFixedComponent(i)) continue;

            computeLsqObjectiveJacobianStructureForVertex(vertex, vertex->getVertexIdx() + free_idx, i_row, j_col, nnz_idx);

            ++free_idx;
        }
    }
}

void HyperGraphOptimizationProblemVertexBased::computeSparseJacobianLsqObjectiveValues(Eigen::Ref<Eigen::VectorXd> values, const double* multipliers)
{
    assert(values.size() == computeSparseJacobianLsqObjectiveNNZ());
    assert(getGraph().getEdgeSet()->isEdgeCacheEmpty());

    // we need to set to zero everything first
    values.setZero();

    int nnz_idx = 0;

    for (VertexInterface* vertex : _graph.getVertexSet()->getActiveVertices())
    {
        if (vertex->isFixed()) continue;

        int free_idx = 0;
        for (int i = 0; i < vertex->getDimension(); ++i)
        {
            if (vertex->isFixedComponent(i)) continue;

            vertex->plus(i, cd_delta);

            computeLsqObjectiveValuesCached(vertex);

            vertex->plus(i, cd_neg2delta);

            finalizeLsqObjectiveJacobianSparseValues(vertex, nnz_idx, values, multipliers);

            vertex->plus(i, cd_delta);  // revert offset
            ++free_idx;
        }
    }
    assert(getGraph().getEdgeSet()->isEdgeCacheEmpty());
}

void HyperGraphOptimizationProblemVertexBased::computeDenseJacobianEqualities(Eigen::Ref<Eigen::MatrixXd> jacobian, const double* multipliers)
{
    assert(_graph.hasVertexSet() && _graph.hasEdgeSet());
    assert(getGraph().getEdgeSet()->isEdgeCacheEmpty());
    if (!_graph_precomputed) precomputeGraphQuantities();

    assert(jacobian.rows() == getEqualityDimension());
    assert(jacobian.cols() == getParameterDimension());

    // we need to set to zero everything first
    jacobian.setZero();

    for (VertexInterface* vertex : _graph.getVertexSet()->getActiveVertices())
    {
        if (vertex->isFixed()) continue;

        int free_idx = 0;
        for (int i = 0; i < vertex->getDimension(); ++i)
        {
            if (vertex->isFixedComponent(i)) continue;

            vertex->plus(i, cd_delta);

            computeEqualitiesValuesCached(vertex);

            vertex->plus(i, cd_neg2delta);

            finalizeEqualitiesJacobian(vertex, vertex->getVertexIdx() + free_idx, jacobian, multipliers);

            vertex->plus(i, cd_delta);  // revert offset
            ++free_idx;
        }
    }
    assert(getGraph().getEdgeSet()->isEdgeCacheEmpty());
}

int HyperGraphOptimizationProblemVertexBased::computeSparseJacobianEqualitiesNNZ()
{
    int nnz = 0;

    assert(_graph.hasEdgeSet());
    if (!_graph_precomputed) precomputeGraphQuantities();

    OptimizationEdgeSet::Ptr edges = _graph.getEdgeSet();

    // Iterate equality edges
    for (BaseEdge::Ptr& edge : edges->getEqualityEdgesRef())
    {
        for (int i = 0; i < edge->getNumVertices(); ++i)
        {
            nnz += edge->getDimension() * edge->getVertexRaw(i)->getDimensionUnfixed();  // block size
        }
    }
    // Iterate mixed edges
    for (BaseMixedEdge::Ptr& edge : edges->getMixedEdgesRef())
    {
        if (edge->getEqualityDimension() == 0) continue;

        for (int i = 0; i < edge->getNumVertices(); ++i)
        {
            nnz += edge->getEqualityDimension() * edge->getVertexRaw(i)->getDimensionUnfixed();  // block size
        }
    }
    return nnz;
}

void HyperGraphOptimizationProblemVertexBased::computeSparseJacobianEqualitiesStructure(Eigen::Ref<Eigen::VectorXi> i_row,
                                                                                        Eigen::Ref<Eigen::VectorXi> j_col)
{
    assert(i_row.size() == computeSparseJacobianEqualitiesNNZ());
    assert(j_col.size() == i_row.size());

    int nnz_idx = 0;

    for (VertexInterface* vertex : _graph.getVertexSet()->getActiveVertices())
    {
        if (vertex->isFixed()) continue;

        int free_idx = 0;
        for (int i = 0; i < vertex->getDimension(); ++i)
        {
            if (vertex->isFixedComponent(i)) continue;

            computeEqualitiesJacobianStructureForVertex(vertex, vertex->getVertexIdx() + free_idx, i_row, j_col, nnz_idx);

            ++free_idx;
        }
    }
}

void HyperGraphOptimizationProblemVertexBased::computeSparseJacobianEqualitiesValues(Eigen::Ref<Eigen::VectorXd> values, const double* multipliers)
{
    assert(values.size() == computeSparseJacobianEqualitiesNNZ());
    assert(getGraph().getEdgeSet()->isEdgeCacheEmpty());

    // we need to set to zero everything first
    values.setZero();

    int nnz_idx = 0;

    for (VertexInterface* vertex : _graph.getVertexSet()->getActiveVertices())
    {
        if (vertex->isFixed()) continue;

        int free_idx = 0;
        for (int i = 0; i < vertex->getDimension(); ++i)
        {
            if (vertex->isFixedComponent(i)) continue;

            vertex->plus(i, cd_delta);

            computeEqualitiesValuesCached(vertex);

            vertex->plus(i, cd_neg2delta);

            finalizeEqualitiesJacobianSparseValues(vertex, nnz_idx, values, multipliers);

            vertex->plus(i, cd_delta);  // revert offset
            ++free_idx;
        }
    }
    assert(getGraph().getEdgeSet()->isEdgeCacheEmpty());
}

void HyperGraphOptimizationProblemVertexBased::computeDenseJacobianInequalities(Eigen::Ref<Eigen::MatrixXd> jacobian, const double* multipliers)
{
    assert(_graph.hasVertexSet() && _graph.hasEdgeSet());
    assert(getGraph().getEdgeSet()->isEdgeCacheEmpty());
    if (!_graph_precomputed) precomputeGraphQuantities();

    assert(jacobian.rows() == getInequalityDimension());
    assert(jacobian.cols() == getParameterDimension());

    // we need to set to zero everything first
    jacobian.setZero();

    for (VertexInterface* vertex : _graph.getVertexSet()->getActiveVertices())
    {
        if (vertex->isFixed()) continue;

        int free_idx = 0;
        for (int i = 0; i < vertex->getDimension(); ++i)
        {
            if (vertex->isFixedComponent(i)) continue;

            vertex->plus(i, cd_delta);

            computeInequalitiesValuesCached(vertex);

            vertex->plus(i, cd_neg2delta);

            finalizeInequalitiesJacobian(vertex, vertex->getVertexIdx() + free_idx, jacobian, multipliers);

            vertex->plus(i, cd_delta);  // revert offset
            ++free_idx;
        }
    }
    assert(getGraph().getEdgeSet()->isEdgeCacheEmpty());
}

int HyperGraphOptimizationProblemVertexBased::computeSparseJacobianInequalitiesNNZ()
{
    int nnz = 0;

    assert(_graph.hasEdgeSet());
    if (!_graph_precomputed) precomputeGraphQuantities();

    OptimizationEdgeSet::Ptr edges = _graph.getEdgeSet();

    // Iterate inequality edges
    for (BaseEdge::Ptr& edge : edges->getInequalityEdgesRef())
    {
        for (int i = 0; i < edge->getNumVertices(); ++i)
        {
            nnz += edge->getDimension() * edge->getVertexRaw(i)->getDimensionUnfixed();  // block size
        }
    }
    // Iterate mixed edges
    for (BaseMixedEdge::Ptr& edge : edges->getMixedEdgesRef())
    {
        if (edge->getInequalityDimension() == 0) continue;

        for (int i = 0; i < edge->getNumVertices(); ++i)
        {
            nnz += edge->getInequalityDimension() * edge->getVertexRaw(i)->getDimensionUnfixed();  // block size
        }
    }
    return nnz;
}

void HyperGraphOptimizationProblemVertexBased::computeSparseJacobianInequalitiesStructure(Eigen::Ref<Eigen::VectorXi> i_row,
                                                                                          Eigen::Ref<Eigen::VectorXi> j_col)
{
    assert(i_row.size() == computeSparseJacobianInequalitiesNNZ());
    assert(j_col.size() == i_row.size());

    int nnz_idx = 0;

    for (VertexInterface* vertex : _graph.getVertexSet()->getActiveVertices())
    {
        if (vertex->isFixed()) continue;

        int free_idx = 0;
        for (int i = 0; i < vertex->getDimension(); ++i)
        {
            if (vertex->isFixedComponent(i)) continue;

            computeInequalitiesJacobianStructureForVertex(vertex, vertex->getVertexIdx() + free_idx, i_row, j_col, nnz_idx);

            ++free_idx;
        }
    }
    assert(getGraph().getEdgeSet()->isEdgeCacheEmpty());
}

void HyperGraphOptimizationProblemVertexBased::computeSparseJacobianInequalitiesValues(Eigen::Ref<Eigen::VectorXd> values, const double* multipliers)
{
    assert(values.size() == computeSparseJacobianInequalitiesNNZ());
    assert(getGraph().getEdgeSet()->isEdgeCacheEmpty());

    // we need to set to zero everything first
    values.setZero();

    int nnz_idx = 0;

    for (VertexInterface* vertex : _graph.getVertexSet()->getActiveVertices())
    {
        if (vertex->isFixed()) continue;

        int free_idx = 0;
        for (int i = 0; i < vertex->getDimension(); ++i)
        {
            if (vertex->isFixedComponent(i)) continue;

            vertex->plus(i, cd_delta);

            computeInequalitiesValuesCached(vertex);

            vertex->plus(i, cd_neg2delta);

            finalizeInequalitiesJacobianSparseValues(vertex, nnz_idx, values, multipliers);

            vertex->plus(i, cd_delta);  // revert offset
            ++free_idx;
        }
    }
    assert(getGraph().getEdgeSet()->isEdgeCacheEmpty());
}

void HyperGraphOptimizationProblemVertexBased::computeDenseJacobianActiveInequalities(Eigen::Ref<Eigen::MatrixXd> jacobian, double weight)
{
    assert(_graph.hasVertexSet() && _graph.hasEdgeSet());
    assert(getGraph().getEdgeSet()->isEdgeCacheEmpty());
    if (!_graph_precomputed) precomputeGraphQuantities();

    assert(jacobian.rows() == getInequalityDimension());
    assert(jacobian.cols() == getParameterDimension());

    // we need to set to zero everything first
    jacobian.setZero();

    for (VertexInterface* vertex : _graph.getVertexSet()->getActiveVertices())
    {
        if (vertex->isFixed()) continue;

        int free_idx = 0;
        for (int i = 0; i < vertex->getDimension(); ++i)
        {
            if (vertex->isFixedComponent(i)) continue;

            computeInequalitiesValuesCached(vertex);  // we also need the current values

            vertex->plus(i, cd_delta);

            computeInequalitiesValuesCached(vertex);

            vertex->plus(i, cd_neg2delta);

            finalizeActiveInequalitiesJacobian(vertex, vertex->getVertexIdx() + free_idx, jacobian, weight);

            vertex->plus(i, cd_delta);  // revert offset
            ++free_idx;
        }
    }
    assert(getGraph().getEdgeSet()->isEdgeCacheEmpty());
}
void HyperGraphOptimizationProblemVertexBased::computeSparseJacobianActiveInequalitiesValues(Eigen::Ref<Eigen::VectorXd> values, double weight)
{
    assert(values.size() == computeSparseJacobianInequalitiesNNZ());
    assert(getGraph().getEdgeSet()->isEdgeCacheEmpty());

    // we need to set to zero everything first
    values.setZero();

    int nnz_idx = 0;

    for (VertexInterface* vertex : _graph.getVertexSet()->getActiveVertices())
    {
        if (vertex->isFixed()) continue;

        int free_idx = 0;
        for (int i = 0; i < vertex->getDimension(); ++i)
        {
            if (vertex->isFixedComponent(i)) continue;

            computeInequalitiesValuesCached(vertex);

            vertex->plus(i, cd_delta);

            computeInequalitiesValuesCached(vertex);

            vertex->plus(i, cd_neg2delta);

            finalizeActiveInequalitiesJacobianSparseValues(vertex, nnz_idx, values, weight);

            vertex->plus(i, cd_delta);  // revert offset
            ++free_idx;
        }
    }

    assert(getGraph().getEdgeSet()->isEdgeCacheEmpty());
}

int HyperGraphOptimizationProblemVertexBased::computeCombinedSparseJacobiansNNZ(bool objective_lsq, bool equality, bool inequality)
{
    int nnz = 0;
    if (objective_lsq) nnz += computeSparseJacobianLsqObjectiveNNZ();
    if (equality) nnz += computeSparseJacobianEqualitiesNNZ();
    if (inequality) nnz += computeSparseJacobianInequalitiesNNZ();
    return nnz;
}
void HyperGraphOptimizationProblemVertexBased::computeCombinedSparseJacobiansStructure(Eigen::Ref<Eigen::VectorXi> i_row,
                                                                                       Eigen::Ref<Eigen::VectorXi> j_col, bool objective_lsq,
                                                                                       bool equality, bool inequality)
{
    assert(i_row.size() == computeCombinedSparseJacobiansNNZ(objective_lsq, equality, inequality));
    assert(j_col.size() == i_row.size());

    int nnz_idx = 0;

    int eq_row_offset   = objective_lsq ? getLsqObjectiveDimension() : 0;
    int ineq_row_offset = equality ? eq_row_offset + getEqualityDimension() : eq_row_offset;

    for (VertexInterface* vertex : _graph.getVertexSet()->getActiveVertices())
    {
        if (vertex->isFixed()) continue;

        int free_idx = 0;
        for (int i = 0; i < vertex->getDimension(); ++i)
        {
            if (vertex->isFixedComponent(i)) continue;

            if (objective_lsq) computeLsqObjectiveJacobianStructureForVertex(vertex, vertex->getVertexIdx() + free_idx, i_row, j_col, nnz_idx, 0);
            if (equality)
                computeEqualitiesJacobianStructureForVertex(vertex, vertex->getVertexIdx() + free_idx, i_row, j_col, nnz_idx, eq_row_offset);
            if (inequality)
                computeInequalitiesJacobianStructureForVertex(vertex, vertex->getVertexIdx() + free_idx, i_row, j_col, nnz_idx, ineq_row_offset);

            ++free_idx;
        }
    }
}
void HyperGraphOptimizationProblemVertexBased::computeCombinedSparseJacobiansValues(Eigen::Ref<Eigen::VectorXd> values, bool objective_lsq,
                                                                                    bool equality, bool inequality, const double* multipliers_obj,
                                                                                    const double* multipliers_eq, const double* multipliers_ineq)
{
    assert(values.size() == computeCombinedSparseJacobiansNNZ(objective_lsq, equality, inequality));
    assert(getGraph().getEdgeSet()->isEdgeCacheEmpty());

    // we need to set to zero everything first
    values.setZero();

    int nnz_idx = 0;

    for (VertexInterface* vertex : _graph.getVertexSet()->getActiveVertices())
    {
        if (vertex->isFixed()) continue;

        int free_idx = 0;
        for (int i = 0; i < vertex->getDimension(); ++i)
        {
            if (vertex->isFixedComponent(i)) continue;

            vertex->plus(i, cd_delta);

            precomputeConnectedMixedEdges(vertex, objective_lsq, equality, inequality);

            if (objective_lsq) computeLsqObjectiveValuesCached(vertex, true, true, false);
            if (equality) computeEqualitiesValuesCached(vertex, true, true, false);
            if (inequality) computeInequalitiesValuesCached(vertex, true, true, false);

            vertex->plus(i, cd_neg2delta);

            precomputeConnectedMixedEdges(vertex, objective_lsq, equality, inequality);

            if (objective_lsq) finalizeLsqObjectiveJacobianSparseValues(vertex, nnz_idx, values, multipliers_obj, false);
            if (equality) finalizeEqualitiesJacobianSparseValues(vertex, nnz_idx, values, multipliers_eq, false);
            if (inequality) finalizeInequalitiesJacobianSparseValues(vertex, nnz_idx, values, multipliers_ineq, false);

            vertex->plus(i, cd_delta);  // revert offset
            ++free_idx;
        }
    }

    assert(getGraph().getEdgeSet()->isEdgeCacheEmpty());
}

void HyperGraphOptimizationProblemVertexBased::computeCombinedSparseJacobian(Eigen::SparseMatrix<double>& jacobian, bool objective_lsq, bool equality,
                                                                             bool inequality, bool finite_combined_bounds, bool active_ineq,
                                                                             double weight_eq, double weight_ineq, double weight_bounds,
                                                                             const Eigen::VectorXd* values, const Eigen::VectorXi* col_nnz)
{
    PRINT_ERROR_NAMED("Not yet implemented for the vertex based strategy");
    OptimizationProblemInterface::computeCombinedSparseJacobian(jacobian, objective_lsq, equality, inequality, finite_combined_bounds, active_ineq,
                                                                weight_eq, weight_ineq, weight_bounds, values, col_nnz);
}

void HyperGraphOptimizationProblemVertexBased::computeGradientObjectiveAndCombinedSparseJacobiansValues(Eigen::Ref<Eigen::VectorXd> gradient,
                                                                                                        Eigen::Ref<Eigen::VectorXd> jac_values,
                                                                                                        bool equality, bool inequality,
                                                                                                        const double* multipliers_eq,
                                                                                                        const double* multipliers_ineq)
{
    assert(jac_values.size() == computeCombinedSparseJacobiansNNZ(false, equality, inequality));
    assert(getGraph().getEdgeSet()->isEdgeCacheEmpty());

    // we need to set to zero everything first
    gradient.setZero();
    jac_values.setZero();

    int nnz_idx = 0;

    for (VertexInterface* vertex : _graph.getVertexSet()->getActiveVertices())
    {
        if (vertex->isFixed()) continue;

        int free_idx = 0;
        for (int i = 0; i < vertex->getDimension(); ++i)
        {
            if (vertex->isFixedComponent(i)) continue;

            vertex->plus(i, cd_delta);

            precomputeConnectedMixedEdges(vertex, true, equality, inequality);

            computeObjectiveValuesCached(vertex, true);
            if (equality) computeEqualitiesValuesCached(vertex, true, true, false);
            if (inequality) computeInequalitiesValuesCached(vertex, true, true, false);

            vertex->plus(i, cd_neg2delta);

            precomputeConnectedMixedEdges(vertex, true, equality, inequality);

            finalizeObjectiveGradient(vertex, gradient(vertex->getVertexIdx() + free_idx), true, true, true, false);
            if (equality) finalizeEqualitiesJacobianSparseValues(vertex, nnz_idx, jac_values, multipliers_eq, false);
            if (inequality) finalizeInequalitiesJacobianSparseValues(vertex, nnz_idx, jac_values, multipliers_ineq, false);

            vertex->plus(i, cd_delta);  // revert offset
            ++free_idx;
        }
    }

    assert(getGraph().getEdgeSet()->isEdgeCacheEmpty());
}

}  // namespace corbo
