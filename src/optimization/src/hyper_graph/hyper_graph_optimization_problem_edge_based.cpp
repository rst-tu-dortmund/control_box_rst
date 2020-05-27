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

#include <corbo-optimization/hyper_graph/hyper_graph_optimization_problem_edge_based.h>

namespace corbo {

// TODO(roesmann): reduce amount of code by introducing more general functions working on all edge containers!

void HyperGraphOptimizationProblemEdgeBased::computeGradientObjective(Eigen::Ref<Eigen::VectorXd> gradient)
{
    assert(_graph.hasEdgeSet());
    if (!_graph_precomputed) precomputeGraphQuantities();
    assert(gradient.size() == getParameterDimension());

    // we need to set to zero everything first
    gradient.setZero();

    OptimizationEdgeSet::Ptr edges = _graph.getEdgeSet();

    // Iterate plain objective edges
    for (BaseEdge::Ptr& edge : edges->getObjectiveEdgesRef())
    {
        // if (edge->providesJacobian())
        for (int i = 0; i < edge->getNumVertices(); ++i)
        {
            int vert_dim_unfixed = edge->getVertexRaw(i)->getDimensionUnfixed();
            if (vert_dim_unfixed == 0) continue;

            Eigen::MatrixXd block_jacobian(edge->getDimension(), vert_dim_unfixed);
            edge->computeJacobian(i, block_jacobian, nullptr);
            gradient.segment(edge->getVertexRaw(i)->getVertexIdx(), vert_dim_unfixed) +=
                block_jacobian.colwise().sum();  // sum up all values, we need a dimension  of one
        }
    }
    // Iterate lsq objective edges
    for (BaseEdge::Ptr& edge : edges->getLsqObjectiveEdgesRef())
    {
        for (int i = 0; i < edge->getNumVertices(); ++i)
        {
            int vert_dim_unfixed = edge->getVertexRaw(i)->getDimensionUnfixed();
            if (vert_dim_unfixed == 0) continue;

            Eigen::MatrixXd block_jacobian(edge->getDimension(), vert_dim_unfixed);
            edge->computeJacobian(i, block_jacobian, nullptr);

            // apply chain rule // TODO(roesmann): in case of forward differences we could use this for computing the jacobian
            Eigen::VectorXd values(edge->getDimension());
            edge->computeValues(values);

            gradient.segment(edge->getVertexRaw(i)->getVertexIdx(), vert_dim_unfixed) += 2.0 * values.transpose() * block_jacobian;
        }
    }
    // Iterate mixed edges
    for (BaseMixedEdge::Ptr& edge : edges->getMixedEdgesRef())
    {
        if (edge->getObjectiveDimension() == 0) continue;

        for (int i = 0; i < edge->getNumVertices(); ++i)
        {
            int vert_dim_unfixed = edge->getVertexRaw(i)->getDimensionUnfixed();
            if (vert_dim_unfixed == 0) continue;

            Eigen::MatrixXd block_jacobian(edge->getObjectiveDimension(), vert_dim_unfixed);
            edge->computeObjectiveJacobian(i, block_jacobian, nullptr);

            if (edge->isObjectiveLeastSquaresForm())
            {
                // apply chain rule // TODO(roesmann): in case of forward differences we could use this for computing the jacobian
                Eigen::VectorXd values(edge->getObjectiveDimension());
                edge->computeObjectiveValues(values);
                gradient.segment(edge->getVertexRaw(i)->getVertexIdx(), vert_dim_unfixed) += 2.0 * values.transpose() * block_jacobian;
            }
            else
            {
                gradient.segment(edge->getVertexRaw(i)->getVertexIdx(), vert_dim_unfixed) +=
                    block_jacobian.colwise().sum();  // sum up all values, we need a dimension of one
            }
        }
    }
}

void HyperGraphOptimizationProblemEdgeBased::computeGradientNonLsqObjective(Eigen::Ref<Eigen::VectorXd> gradient)
{
    assert(_graph.hasEdgeSet());
    if (!_graph_precomputed) precomputeGraphQuantities();
    assert(gradient.size() == getParameterDimension());

    // we need to set to zero everything first
    gradient.setZero();

    OptimizationEdgeSet::Ptr edges = _graph.getEdgeSet();

    // Iterate plain objective edges
    for (BaseEdge::Ptr& edge : edges->getObjectiveEdgesRef())
    {
        // if (edge->providesJacobian())
        for (int i = 0; i < edge->getNumVertices(); ++i)
        {
            int vert_dim_unfixed = edge->getVertexRaw(i)->getDimensionUnfixed();
            if (vert_dim_unfixed == 0) continue;

            Eigen::MatrixXd block_jacobian(edge->getDimension(), vert_dim_unfixed);
            edge->computeJacobian(i, block_jacobian, nullptr);
            gradient.segment(edge->getVertexRaw(i)->getVertexIdx(), vert_dim_unfixed) +=
                block_jacobian.colwise().sum();  // sum up all values, we need a dimension  of one
        }
    }

    // Iterate mixed edges
    for (BaseMixedEdge::Ptr& edge : edges->getMixedEdgesRef())
    {
        if (edge->getObjectiveDimension() == 0 || edge->isObjectiveLeastSquaresForm()) continue;

        for (int i = 0; i < edge->getNumVertices(); ++i)
        {
            int vert_dim_unfixed = edge->getVertexRaw(i)->getDimensionUnfixed();
            if (vert_dim_unfixed == 0) continue;

            Eigen::MatrixXd block_jacobian(edge->getObjectiveDimension(), vert_dim_unfixed);
            edge->computeObjectiveJacobian(i, block_jacobian, nullptr);

            gradient.segment(edge->getVertexRaw(i)->getVertexIdx(), vert_dim_unfixed) +=
                block_jacobian.colwise().sum();  // sum up all values, we need a dimension of one
        }
    }
}

void HyperGraphOptimizationProblemEdgeBased::computeDenseJacobianLsqObjective(Eigen::Ref<Eigen::MatrixXd> jacobian, const double* multipliers)
{
    assert(_graph.hasEdgeSet());
    if (!_graph_precomputed) precomputeGraphQuantities();
    assert(jacobian.rows() == getLsqObjectiveDimension());
    assert(jacobian.cols() == getParameterDimension());

    // we need to set to zero everything first
    jacobian.setZero();

    OptimizationEdgeSet::Ptr edges = _graph.getEdgeSet();

    // Iterate lsq objective edges
    for (BaseEdge::Ptr& edge : edges->getLsqObjectiveEdgesRef())
    {
        for (int i = 0; i < edge->getNumVertices(); ++i)
        {
            int vert_dim_unfixed = edge->getVertexRaw(i)->getDimensionUnfixed();
            if (vert_dim_unfixed == 0) continue;

            edge->computeJacobian(i,
                                  jacobian.block(edge->getEdgeIdx(), edge->getVertexRaw(i)->getVertexIdx(), edge->getDimension(), vert_dim_unfixed),
                                  multipliers ? multipliers + edge->getEdgeIdx() : nullptr);
        }
    }
    // Iterate mixed edges
    for (BaseMixedEdge::Ptr& edge : edges->getMixedEdgesRef())
    {
        if (edge->getObjectiveDimension() == 0 || !edge->isObjectiveLeastSquaresForm()) continue;

        for (int i = 0; i < edge->getNumVertices(); ++i)
        {
            int vert_dim_unfixed = edge->getVertexRaw(i)->getDimensionUnfixed();
            if (vert_dim_unfixed == 0) continue;

            edge->computeObjectiveJacobian(
                i,
                jacobian.block(edge->getEdgeObjectiveIdx(), edge->getVertexRaw(i)->getVertexIdx(), edge->getObjectiveDimension(), vert_dim_unfixed),
                multipliers ? multipliers + edge->getEdgeObjectiveIdx() : nullptr);
        }
    }
}

int HyperGraphOptimizationProblemEdgeBased::computeSparseJacobianLsqObjectiveNNZ()
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

void HyperGraphOptimizationProblemEdgeBased::computeSparseJacobianLsqObjectiveStructure(Eigen::Ref<Eigen::VectorXi> i_row,
                                                                                        Eigen::Ref<Eigen::VectorXi> j_col)
{
    assert(i_row.size() == computeSparseJacobianLsqObjectiveNNZ());
    assert(j_col.size() == i_row.size());

    OptimizationEdgeSet::Ptr edges = _graph.getEdgeSet();

    int nnz_idx = 0;

    // Iterate lsq objective edges
    for (BaseEdge::Ptr& edge : edges->getLsqObjectiveEdgesRef())
    {
        for (int vert_idx = 0; vert_idx < edge->getNumVertices(); ++vert_idx)
        {
            const VertexInterface* vertex = edge->getVertexRaw(vert_idx);
            // Iterate all free variables
            int idx_free = 0;
            for (int i = 0; i < vertex->getDimension(); ++i)
            {
                if (!vertex->isFixedComponent(i))
                {
                    // Iterate inner edge dimension
                    for (int j = 0; j < edge->getDimension(); ++j)
                    {
                        i_row[nnz_idx] = edge->getEdgeIdx() + j;
                        j_col[nnz_idx] = vertex->getVertexIdx() + idx_free;
                        ++nnz_idx;
                    }
                    ++idx_free;
                }
            }
        }
    }
    // Iterate mixed edges
    for (BaseMixedEdge::Ptr& edge : edges->getMixedEdgesRef())
    {
        if (edge->getObjectiveDimension() == 0 || !edge->isObjectiveLeastSquaresForm()) continue;

        for (int vert_idx = 0; vert_idx < edge->getNumVertices(); ++vert_idx)
        {
            const VertexInterface* vertex = edge->getVertexRaw(vert_idx);
            // Iterate all free variables
            int idx_free = 0;
            for (int i = 0; i < vertex->getDimension(); ++i)
            {
                if (!vertex->isFixedComponent(i))
                {
                    // Iterate inner edge dimension
                    for (int j = 0; j < edge->getObjectiveDimension(); ++j)
                    {
                        i_row[nnz_idx] = edge->getEdgeObjectiveIdx() + j;
                        j_col[nnz_idx] = vertex->getVertexIdx() + idx_free;
                        ++nnz_idx;
                    }
                    ++idx_free;
                }
            }
        }
    }
}

void HyperGraphOptimizationProblemEdgeBased::computeSparseJacobianLsqObjectiveValues(Eigen::Ref<Eigen::VectorXd> values, const double* multipliers)
{
    assert(values.size() == computeSparseJacobianLsqObjectiveNNZ());

    int nnz_idx = 0;

    OptimizationEdgeSet::Ptr edges = _graph.getEdgeSet();

    // Iterate lsq objective edges
    for (BaseEdge::Ptr& edge : edges->getLsqObjectiveEdgesRef())
    {
        for (int vert_idx = 0; vert_idx < edge->getNumVertices(); ++vert_idx)
        {
            int vert_dim_unfixed = edge->getVertexRaw(vert_idx)->getDimensionUnfixed();
            if (vert_dim_unfixed == 0) continue;

            Eigen::Map<Eigen::MatrixXd> block_jacobian(values.data() + nnz_idx, edge->getDimension(), vert_dim_unfixed);
            edge->computeJacobian(vert_idx, block_jacobian, multipliers ? multipliers + edge->getEdgeIdx() : nullptr);
            nnz_idx += block_jacobian.rows() * block_jacobian.cols();
        }
    }

    // Iterate mixed edges
    for (BaseMixedEdge::Ptr& edge : edges->getMixedEdgesRef())
    {
        if (edge->getObjectiveDimension() == 0 || !edge->isObjectiveLeastSquaresForm()) continue;

        for (int vert_idx = 0; vert_idx < edge->getNumVertices(); ++vert_idx)
        {
            int vert_dim_unfixed = edge->getVertexRaw(vert_idx)->getDimensionUnfixed();
            if (vert_dim_unfixed == 0) continue;

            Eigen::Map<Eigen::MatrixXd> block_jacobian(values.data() + nnz_idx, edge->getObjectiveDimension(), vert_dim_unfixed);
            edge->computeObjectiveJacobian(vert_idx, block_jacobian, multipliers ? multipliers + edge->getEdgeObjectiveIdx() : nullptr);
            nnz_idx += block_jacobian.rows() * block_jacobian.cols();
        }
    }
}

void HyperGraphOptimizationProblemEdgeBased::computeDenseJacobianEqualities(Eigen::Ref<Eigen::MatrixXd> jacobian, const double* multipliers)
{
    assert(_graph.hasEdgeSet());
    if (!_graph_precomputed) precomputeGraphQuantities();
    assert(jacobian.rows() == getEqualityDimension());
    assert(jacobian.cols() == getParameterDimension());

    // we need to set to zero everything first
    jacobian.setZero();

    OptimizationEdgeSet::Ptr edges = _graph.getEdgeSet();

    // Iterate equality edges
    for (BaseEdge::Ptr& edge : edges->getEqualityEdgesRef())
    {
        for (int i = 0; i < edge->getNumVertices(); ++i)
        {
            int vert_dim_unfixed = edge->getVertexRaw(i)->getDimensionUnfixed();
            if (vert_dim_unfixed == 0) continue;

            Eigen::MatrixXd block_jacobian(edge->getDimension(), vert_dim_unfixed);
            edge->computeJacobian(i,
                                  jacobian.block(edge->getEdgeIdx(), edge->getVertexRaw(i)->getVertexIdx(), edge->getDimension(), vert_dim_unfixed),
                                  multipliers ? multipliers + edge->getEdgeIdx() : nullptr);
        }
    }
    // Iterate mixed edges
    for (BaseMixedEdge::Ptr& edge : edges->getMixedEdgesRef())
    {
        if (edge->getEqualityDimension() == 0) continue;

        for (int i = 0; i < edge->getNumVertices(); ++i)
        {
            int vert_dim_unfixed = edge->getVertexRaw(i)->getDimensionUnfixed();
            if (vert_dim_unfixed == 0) continue;

            Eigen::MatrixXd block_jacobian(edge->getEqualityDimension(), vert_dim_unfixed);
            edge->computeEqualityJacobian(
                i, jacobian.block(edge->getEdgeEqualityIdx(), edge->getVertexRaw(i)->getVertexIdx(), edge->getEqualityDimension(), vert_dim_unfixed),
                multipliers ? multipliers + edge->getEdgeEqualityIdx() : nullptr);
        }
    }
}

int HyperGraphOptimizationProblemEdgeBased::computeSparseJacobianEqualitiesNNZ()
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

void HyperGraphOptimizationProblemEdgeBased::computeSparseJacobianEqualitiesStructure(Eigen::Ref<Eigen::VectorXi> i_row,
                                                                                      Eigen::Ref<Eigen::VectorXi> j_col)
{
    assert(i_row.size() == computeSparseJacobianEqualitiesNNZ());
    assert(j_col.size() == i_row.size());

    OptimizationEdgeSet::Ptr edges = _graph.getEdgeSet();

    int nnz_idx = 0;

    // Iterate equality edges
    for (BaseEdge::Ptr& edge : edges->getEqualityEdgesRef())
    {
        for (int vert_idx = 0; vert_idx < edge->getNumVertices(); ++vert_idx)
        {
            const VertexInterface* vertex = edge->getVertexRaw(vert_idx);
            // Iterate all free variables
            int idx_free = 0;
            for (int i = 0; i < vertex->getDimension(); ++i)
            {
                if (!vertex->isFixedComponent(i))
                {
                    // Iterate inner edge dimension
                    for (int j = 0; j < edge->getDimension(); ++j)
                    {
                        i_row[nnz_idx] = edge->getEdgeIdx() + j;
                        j_col[nnz_idx] = vertex->getVertexIdx() + idx_free;
                        ++nnz_idx;
                    }
                    ++idx_free;
                }
            }
        }
    }
    // Iterate mixed edges
    for (BaseMixedEdge::Ptr& edge : edges->getMixedEdgesRef())
    {
        if (edge->getEqualityDimension() == 0) continue;

        for (int vert_idx = 0; vert_idx < edge->getNumVertices(); ++vert_idx)
        {
            const VertexInterface* vertex = edge->getVertexRaw(vert_idx);
            // Iterate all free variables
            int idx_free = 0;
            for (int i = 0; i < vertex->getDimension(); ++i)
            {
                if (!vertex->isFixedComponent(i))
                {
                    // Iterate inner edge dimension
                    for (int j = 0; j < edge->getEqualityDimension(); ++j)
                    {
                        i_row[nnz_idx] = edge->getEdgeEqualityIdx() + j;
                        j_col[nnz_idx] = vertex->getVertexIdx() + idx_free;
                        ++nnz_idx;
                    }
                    ++idx_free;
                }
            }
        }
    }
}

void HyperGraphOptimizationProblemEdgeBased::computeSparseJacobianEqualitiesValues(Eigen::Ref<Eigen::VectorXd> values, const double* multipliers)
{
    assert(values.size() == computeSparseJacobianEqualitiesNNZ());

    int nnz_idx = 0;

    OptimizationEdgeSet::Ptr edges = _graph.getEdgeSet();

    // Iterate equality edges
    for (BaseEdge::Ptr& edge : edges->getEqualityEdgesRef())
    {
        for (int vert_idx = 0; vert_idx < edge->getNumVertices(); ++vert_idx)
        {
            int vert_dim_unfixed = edge->getVertexRaw(vert_idx)->getDimensionUnfixed();
            if (vert_dim_unfixed == 0) continue;

            Eigen::Map<Eigen::MatrixXd> block_jacobian(values.data() + nnz_idx, edge->getDimension(), vert_dim_unfixed);
            edge->computeJacobian(vert_idx, block_jacobian, multipliers ? multipliers + edge->getEdgeIdx() : nullptr);
            nnz_idx += block_jacobian.rows() * block_jacobian.cols();
        }
    }

    // Iterate mixed edges
    for (BaseMixedEdge::Ptr& edge : edges->getMixedEdgesRef())
    {
        if (edge->getEqualityDimension() == 0) continue;

        for (int vert_idx = 0; vert_idx < edge->getNumVertices(); ++vert_idx)
        {
            int vert_dim_unfixed = edge->getVertexRaw(vert_idx)->getDimensionUnfixed();
            if (vert_dim_unfixed == 0) continue;

            Eigen::Map<Eigen::MatrixXd> block_jacobian(values.data() + nnz_idx, edge->getEqualityDimension(), vert_dim_unfixed);
            edge->computeEqualityJacobian(vert_idx, block_jacobian, multipliers ? multipliers + edge->getEdgeEqualityIdx() : nullptr);
            nnz_idx += block_jacobian.rows() * block_jacobian.cols();
        }
    }
}

void HyperGraphOptimizationProblemEdgeBased::computeDenseJacobianInequalities(Eigen::Ref<Eigen::MatrixXd> jacobian, const double* multipliers)
{
    assert(_graph.hasEdgeSet());
    if (!_graph_precomputed) precomputeGraphQuantities();
    assert(jacobian.rows() == getInequalityDimension());
    assert(jacobian.cols() == getParameterDimension());

    // we need to set to zero everything first
    jacobian.setZero();

    OptimizationEdgeSet::Ptr edges = _graph.getEdgeSet();

    // Iterate inequality edges
    for (BaseEdge::Ptr& edge : edges->getInequalityEdgesRef())
    {
        for (int i = 0; i < edge->getNumVertices(); ++i)
        {
            int vert_dim_unfixed = edge->getVertexRaw(i)->getDimensionUnfixed();
            if (vert_dim_unfixed == 0) continue;

            Eigen::MatrixXd block_jacobian(edge->getDimension(), vert_dim_unfixed);
            edge->computeJacobian(i,
                                  jacobian.block(edge->getEdgeIdx(), edge->getVertexRaw(i)->getVertexIdx(), edge->getDimension(), vert_dim_unfixed),
                                  multipliers ? multipliers + edge->getEdgeIdx() : nullptr);
        }
    }
    // Iterate mixed edges
    for (BaseMixedEdge::Ptr& edge : edges->getMixedEdgesRef())
    {
        if (edge->getInequalityDimension() == 0) continue;

        for (int i = 0; i < edge->getNumVertices(); ++i)
        {
            int vert_dim_unfixed = edge->getVertexRaw(i)->getDimensionUnfixed();
            if (vert_dim_unfixed == 0) continue;

            Eigen::MatrixXd block_jacobian(edge->getInequalityDimension(), vert_dim_unfixed);
            edge->computeInequalityJacobian(
                i,
                jacobian.block(edge->getEdgeInequalityIdx(), edge->getVertexRaw(i)->getVertexIdx(), edge->getInequalityDimension(), vert_dim_unfixed),
                multipliers ? multipliers + edge->getEdgeInequalityIdx() : nullptr);
        }
    }
}

int HyperGraphOptimizationProblemEdgeBased::computeSparseJacobianInequalitiesNNZ()
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

void HyperGraphOptimizationProblemEdgeBased::computeSparseJacobianInequalitiesStructure(Eigen::Ref<Eigen::VectorXi> i_row,
                                                                                        Eigen::Ref<Eigen::VectorXi> j_col)
{
    assert(i_row.size() == computeSparseJacobianInequalitiesNNZ());
    assert(j_col.size() == i_row.size());

    OptimizationEdgeSet::Ptr edges = _graph.getEdgeSet();

    int nnz_idx = 0;

    // Iterate inequality edges
    for (BaseEdge::Ptr& edge : edges->getInequalityEdgesRef())
    {
        for (int vert_idx = 0; vert_idx < edge->getNumVertices(); ++vert_idx)
        {
            const VertexInterface* vertex = edge->getVertexRaw(vert_idx);
            // Iterate all free variables
            int idx_free = 0;
            for (int i = 0; i < vertex->getDimension(); ++i)
            {
                if (!vertex->isFixedComponent(i))
                {
                    // Iterate inner edge dimension
                    for (int j = 0; j < edge->getDimension(); ++j)
                    {
                        i_row[nnz_idx] = edge->getEdgeIdx() + j;
                        j_col[nnz_idx] = vertex->getVertexIdx() + idx_free;
                        ++nnz_idx;
                    }
                    ++idx_free;
                }
            }
        }
    }
    // Iterate mixed edges
    for (BaseMixedEdge::Ptr& edge : edges->getMixedEdgesRef())
    {
        if (edge->getInequalityDimension() == 0) continue;

        for (int vert_idx = 0; vert_idx < edge->getNumVertices(); ++vert_idx)
        {
            const VertexInterface* vertex = edge->getVertexRaw(vert_idx);
            // Iterate all free variables
            int idx_free = 0;
            for (int i = 0; i < vertex->getDimension(); ++i)
            {
                if (!vertex->isFixedComponent(i))
                {
                    // Iterate inner edge dimension
                    for (int j = 0; j < edge->getInequalityDimension(); ++j)
                    {
                        i_row[nnz_idx] = edge->getEdgeInequalityIdx() + j;
                        j_col[nnz_idx] = vertex->getVertexIdx() + idx_free;
                        ++nnz_idx;
                    }
                    ++idx_free;
                }
            }
        }
    }
}

void HyperGraphOptimizationProblemEdgeBased::computeSparseJacobianInequalitiesValues(Eigen::Ref<Eigen::VectorXd> values, const double* multipliers)
{
    assert(values.size() == computeSparseJacobianInequalitiesNNZ());

    int nnz_idx = 0;

    OptimizationEdgeSet::Ptr edges = _graph.getEdgeSet();

    // Iterate inequality edges
    for (BaseEdge::Ptr& edge : edges->getInequalityEdgesRef())
    {
        for (int vert_idx = 0; vert_idx < edge->getNumVertices(); ++vert_idx)
        {
            int vert_dim_unfixed = edge->getVertexRaw(vert_idx)->getDimensionUnfixed();
            if (vert_dim_unfixed == 0) continue;

            Eigen::Map<Eigen::MatrixXd> block_jacobian(values.data() + nnz_idx, edge->getDimension(), vert_dim_unfixed);
            edge->computeJacobian(vert_idx, block_jacobian, multipliers ? multipliers + edge->getEdgeIdx() : nullptr);
            nnz_idx += block_jacobian.rows() * block_jacobian.cols();
        }
    }

    // Iterate mixed edges
    for (BaseMixedEdge::Ptr& edge : edges->getMixedEdgesRef())
    {
        if (edge->getInequalityDimension() == 0) continue;

        for (int vert_idx = 0; vert_idx < edge->getNumVertices(); ++vert_idx)
        {
            int vert_dim_unfixed = edge->getVertexRaw(vert_idx)->getDimensionUnfixed();
            if (vert_dim_unfixed == 0) continue;

            Eigen::Map<Eigen::MatrixXd> block_jacobian(values.data() + nnz_idx, edge->getInequalityDimension(), vert_dim_unfixed);
            edge->computeInequalityJacobian(vert_idx, block_jacobian, multipliers ? multipliers + edge->getEdgeInequalityIdx() : nullptr);
            nnz_idx += block_jacobian.rows() * block_jacobian.cols();
        }
    }
}

void HyperGraphOptimizationProblemEdgeBased::computeDenseJacobianActiveInequalities(Eigen::Ref<Eigen::MatrixXd> jacobian, double weight)
{
    assert(_graph.hasEdgeSet());
    if (!_graph_precomputed) precomputeGraphQuantities();
    assert(jacobian.rows() == getInequalityDimension());
    assert(jacobian.cols() == getParameterDimension());

    // we need to set to zero everything first
    jacobian.setZero();

    OptimizationEdgeSet::Ptr edges = _graph.getEdgeSet();

    // Iterate inequality edges
    for (BaseEdge::Ptr& edge : edges->getInequalityEdgesRef())
    {
        for (int i = 0; i < edge->getNumVertices(); ++i)
        {
            const VertexInterface* vertex = edge->getVertexRaw(i);
            int vert_dim_unfixed          = vertex->getDimensionUnfixed();
            if (vert_dim_unfixed == 0) continue;

            // compute values and check which values are active
            // TODO(roesmann): we could obtain values computed before calling computeDenseJacobianActiveInequalities();
            Eigen::VectorXd values(edge->getDimension());
            edge->computeValues(values);
            Eigen::Array<bool, -1, 1> active = values.array() > 0.0;

            if (!active.any()) continue;

            // TODO(roesmann): we need to compute the complete block jacobian here (some rows are already zero according to array "active"
            Eigen::MatrixXd block_jacobian(edge->getDimension(), vert_dim_unfixed);
            edge->computeJacobian(i, jacobian.block(edge->getEdgeIdx(), vertex->getVertexIdx(), edge->getDimension(), vert_dim_unfixed), nullptr);

            // now set values to zero if inactive or multiply with weight otherwise
            for (int j = 0; j < edge->getDimension(); ++j)
            {
                if (active[j] && weight != 1.0)
                    jacobian.block(edge->getEdgeIdx() + j, edge->getVertexRaw(i)->getVertexIdx(), 1, vert_dim_unfixed) *= weight;
                else if (!active[j])
                    jacobian.block(edge->getEdgeIdx() + j, edge->getVertexRaw(i)->getVertexIdx(), 1, vert_dim_unfixed).setZero();  // should be zero
                // already
            }
        }
    }
    // Iterate mixed edges
    for (BaseMixedEdge::Ptr& edge : edges->getMixedEdgesRef())
    {
        if (edge->getInequalityDimension() == 0) continue;

        for (int i = 0; i < edge->getNumVertices(); ++i)
        {
            const VertexInterface* vertex = edge->getVertexRaw(i);
            int vert_dim_unfixed          = edge->getVertexRaw(i)->getDimensionUnfixed();
            if (vert_dim_unfixed == 0) continue;

            // compute values and check which values are active
            // TODO(roesmann): we could obtain values computed before calling computeDenseJacobianActiveInequalities();
            Eigen::VectorXd values(edge->getInequalityDimension());
            edge->precompute();
            edge->computeInequalityValues(values);
            Eigen::Array<bool, -1, 1> active = values.array() > 0.0;

            if (!active.any()) continue;

            Eigen::MatrixXd block_jacobian(edge->getInequalityDimension(), vert_dim_unfixed);
            edge->computeInequalityJacobian(
                i, jacobian.block(edge->getEdgeInequalityIdx(), vertex->getVertexIdx(), edge->getInequalityDimension(), vert_dim_unfixed), nullptr);

            // now set values to zero if inactive or multiply with weight otherwise
            for (int j = 0; j < edge->getInequalityDimension(); ++j)
            {
                if (active[j] && weight != 1.0)
                    jacobian.block(edge->getEdgeInequalityIdx() + j, vertex->getVertexIdx(), 1, vert_dim_unfixed) *= weight;
                else if (!active[j])
                    jacobian.block(edge->getEdgeInequalityIdx() + j, vertex->getVertexIdx(), 1, vert_dim_unfixed).setZero();  // should be zero
            }
        }
    }
}
void HyperGraphOptimizationProblemEdgeBased::computeSparseJacobianActiveInequalitiesValues(Eigen::Ref<Eigen::VectorXd> values, double weight)
{
    assert(values.size() == computeSparseJacobianInequalitiesNNZ());

    int jac_row_idx = 0;
    int nnz_idx     = 0;

    OptimizationEdgeSet::Ptr edges = _graph.getEdgeSet();

    // Iterate inequality edges
    for (BaseEdge::Ptr& edge : edges->getInequalityEdgesRef())
    {
        // compute values and check which values are active
        // TODO(roesmann): we could obtain values computed before calling computeDenseJacobianActiveInequalities();
        Eigen::VectorXd edge_values(edge->getDimension());
        edge->computeValues(edge_values);
        Eigen::Array<bool, -1, 1> active = edge_values.array() > 0.0;

        for (int vert_idx = 0; vert_idx < edge->getNumVertices(); ++vert_idx)
        {
            int vert_dim_unfixed = edge->getVertexRaw(vert_idx)->getDimensionUnfixed();
            if (vert_dim_unfixed == 0) continue;

            Eigen::Map<Eigen::MatrixXd> block_jacobian(values.data() + nnz_idx, edge->getDimension(), vert_dim_unfixed);
            nnz_idx += block_jacobian.rows() * block_jacobian.cols();

            if (!active.any())
            {
                block_jacobian.setZero();
                continue;
            }

            edge->computeJacobian(vert_idx, block_jacobian, nullptr);

            // now set values to zero if inactive or multiply with weight otherwise
            for (int j = 0; j < edge->getDimension(); ++j)
            {
                if (active[j] && weight != 1.0)
                    block_jacobian.block(j, 0, 1, vert_dim_unfixed) *= weight;
                else if (!active[j])
                    block_jacobian.block(j, 0, 1, vert_dim_unfixed).setZero();
            }
        }
        jac_row_idx += edge->getDimension();
    }

    // Iterate mixed edges
    for (BaseMixedEdge::Ptr& edge : edges->getMixedEdgesRef())
    {
        if (edge->getInequalityDimension() == 0) continue;

        // compute values and check which values are active
        // TODO(roesmann): we could obtain values computed before calling computeDenseJacobianActiveInequalities();
        Eigen::VectorXd edge_values(edge->getInequalityDimension());
        edge->precompute();
        edge->computeInequalityValues(edge_values);
        Eigen::Array<bool, -1, 1> active = edge_values.array() > 0.0;

        for (int vert_idx = 0; vert_idx < edge->getNumVertices(); ++vert_idx)
        {
            int vert_dim_unfixed = edge->getVertexRaw(vert_idx)->getDimensionUnfixed();
            if (vert_dim_unfixed == 0) continue;

            Eigen::Map<Eigen::MatrixXd> block_jacobian(values.data() + nnz_idx, edge->getInequalityDimension(), vert_dim_unfixed);
            nnz_idx += block_jacobian.rows() * block_jacobian.cols();

            if (!active.any())
            {
                block_jacobian.setZero();
                continue;
            }

            edge->computeInequalityJacobian(vert_idx, block_jacobian, nullptr);

            // now set values to zero if inactive or multiply with weight otherwise
            for (int j = 0; j < edge->getInequalityDimension(); ++j)
            {
                if (active[j] && weight != 1.0)
                    block_jacobian.block(j, 0, 1, vert_dim_unfixed) *= weight;
                else if (!active[j])
                    block_jacobian.block(j, 0, 1, vert_dim_unfixed).setZero();
            }
        }
        jac_row_idx += edge->getInequalityDimension();
    }
}

void HyperGraphOptimizationProblemEdgeBased::computeDenseJacobians(Eigen::Ref<Eigen::VectorXd> gradient_non_lsq_obj,
                                                                   Eigen::Ref<Eigen::MatrixXd> jacobian_lsq_obj,
                                                                   Eigen::Ref<Eigen::MatrixXd> jacobian_eq, Eigen::Ref<Eigen::MatrixXd> jacobian_ineq,
                                                                   const double* multipliers_lsq_obj, const double* multipliers_eq,
                                                                   const double* multipliers_ineq, bool active_ineq, double active_ineq_weight)
{
    assert(gradient_non_lsq_obj.size() == getParameterDimension());
    assert(jacobian_lsq_obj.rows() == getLsqObjectiveDimension());
    assert(jacobian_lsq_obj.cols() == getParameterDimension());
    assert(jacobian_eq.rows() == getEqualityDimension());
    assert(jacobian_eq.cols() == getParameterDimension());
    assert(jacobian_ineq.rows() == getInequalityDimension());
    assert(jacobian_ineq.cols() == getParameterDimension());

    PRINT_DEBUG_COND_ONCE(active_ineq && multipliers_ineq,
                          "HyperGraphOptimizationProblemEdgeBased::computeDenseJacobians(): multiplier_ineq is ignored if active_ineq is enabled");

    assert(_graph.hasEdgeSet());
    if (!_graph_precomputed) precomputeGraphQuantities();

    // int jac_row_idx = 0;
    // int nnz_idx     = 0;

    // set matrices to zero
    gradient_non_lsq_obj.setZero();
    jacobian_lsq_obj.setZero();
    jacobian_eq.setZero();
    jacobian_ineq.setZero();

    OptimizationEdgeSet::Ptr edges = _graph.getEdgeSet();

    // just copy and paste, the important thing is for MixedEdges -> precompute is called just once

    // Iterate plain objective edges (for gradient)
    for (BaseEdge::Ptr& edge : edges->getObjectiveEdgesRef())
    {
        for (int i = 0; i < edge->getNumVertices(); ++i)
        {
            int vert_dim_unfixed = edge->getVertexRaw(i)->getDimensionUnfixed();
            if (vert_dim_unfixed == 0) continue;

            Eigen::MatrixXd block_jacobian(edge->getDimension(), vert_dim_unfixed);
            edge->computeJacobian(i, block_jacobian, nullptr);
            gradient_non_lsq_obj.segment(edge->getVertexRaw(i)->getVertexIdx(), vert_dim_unfixed) +=
                block_jacobian.colwise().sum();  // sum up all values, we need a dimension  of one
        }
    }

    // Iterate lsq objective edges
    for (BaseEdge::Ptr& edge : edges->getLsqObjectiveEdgesRef())
    {
        for (int i = 0; i < edge->getNumVertices(); ++i)
        {
            int vert_dim_unfixed = edge->getVertexRaw(i)->getDimensionUnfixed();
            if (vert_dim_unfixed == 0) continue;

            Eigen::MatrixXd block_jacobian(edge->getDimension(), vert_dim_unfixed);
            edge->computeJacobian(
                i, jacobian_lsq_obj.block(edge->getEdgeIdx(), edge->getVertexRaw(i)->getVertexIdx(), edge->getDimension(), vert_dim_unfixed),
                multipliers_lsq_obj ? multipliers_lsq_obj + edge->getEdgeIdx() : nullptr);
        }
    }

    // Iterate equality edges
    for (BaseEdge::Ptr& edge : edges->getEqualityEdgesRef())
    {
        for (int i = 0; i < edge->getNumVertices(); ++i)
        {
            int vert_dim_unfixed = edge->getVertexRaw(i)->getDimensionUnfixed();
            if (vert_dim_unfixed == 0) continue;

            Eigen::MatrixXd block_jacobian(edge->getDimension(), vert_dim_unfixed);
            edge->computeJacobian(
                i, jacobian_eq.block(edge->getEdgeIdx(), edge->getVertexRaw(i)->getVertexIdx(), edge->getDimension(), vert_dim_unfixed),
                multipliers_eq ? multipliers_eq + edge->getEdgeIdx() : nullptr);
        }
    }

    // Iterate inequality edges
    for (BaseEdge::Ptr& edge : edges->getInequalityEdgesRef())
    {
        for (int i = 0; i < edge->getNumVertices(); ++i)
        {
            const VertexInterface* vertex = edge->getVertexRaw(i);
            int vert_dim_unfixed          = vertex->getDimensionUnfixed();
            if (vert_dim_unfixed == 0) continue;

            if (active_ineq)
            {
                // compute values and check which values are active
                // TODO(roesmann): we could obtain values computed before calling computeDenseJacobianActiveInequalities();
                Eigen::VectorXd edge_values(edge->getDimension());
                edge->computeValues(edge_values);
                Eigen::Array<bool, -1, 1> active = edge_values.array() > 0.0;

                if (!active.any()) continue;

                // TODO(roesmann): we need to compute the complete block jacobian here (some rows are already zero according to array "active"
                Eigen::MatrixXd block_jacobian(edge->getDimension(), vert_dim_unfixed);
                edge->computeJacobian(i, jacobian_ineq.block(edge->getEdgeIdx(), vertex->getVertexIdx(), edge->getDimension(), vert_dim_unfixed),
                                      nullptr);

                // now set values to zero if inactive or multiply with weight otherwise
                for (int j = 0; j < edge->getDimension(); ++j)
                {
                    if (active[j] && active_ineq_weight != 1.0)
                        jacobian_ineq.block(edge->getEdgeIdx() + j, vertex->getVertexIdx(), 1, vert_dim_unfixed) *= active_ineq_weight;
                    else if (!active[j])
                        jacobian_ineq.block(edge->getEdgeIdx() + j, vertex->getVertexIdx(), 1, vert_dim_unfixed).setZero();
                }
            }
            else
            {
                Eigen::MatrixXd block_jacobian(edge->getDimension(), vert_dim_unfixed);
                edge->computeJacobian(i, jacobian_ineq.block(edge->getEdgeIdx(), vertex->getVertexIdx(), edge->getDimension(), vert_dim_unfixed),
                                      multipliers_ineq ? multipliers_ineq + edge->getEdgeIdx() : nullptr);
            }
        }
    }

    // Iterate mixed edges
    for (BaseMixedEdge::Ptr& edge : edges->getMixedEdgesRef())
    {
        const double* mult_ineq_part    = (multipliers_ineq && !active_ineq) ? multipliers_ineq + edge->getEdgeInequalityIdx() : nullptr;
        const double* mult_eq_part      = multipliers_eq ? multipliers_eq + edge->getEdgeEqualityIdx() : nullptr;
        const double* mult_obj_lsq_part = multipliers_lsq_obj ? multipliers_lsq_obj + edge->getEdgeObjectiveIdx() : nullptr;

        for (int i = 0; i < edge->getNumVertices(); ++i)
        {
            const VertexInterface* vertex = edge->getVertexRaw(i);
            int vert_dim_unfixed          = vertex->getDimensionUnfixed();
            if (vert_dim_unfixed == 0) continue;

            Eigen::Ref<Eigen::MatrixXd> block_eq =
                jacobian_eq.block(edge->getEdgeEqualityIdx(), vertex->getVertexIdx(), edge->getEqualityDimension(), vert_dim_unfixed);

            Eigen::Ref<Eigen::MatrixXd> block_ineq =
                jacobian_ineq.block(edge->getEdgeInequalityIdx(), vertex->getVertexIdx(), edge->getInequalityDimension(), vert_dim_unfixed);

            if (edge->isObjectiveLeastSquaresForm())
            {
                Eigen::Ref<Eigen::MatrixXd> block_lsq_obj =
                    jacobian_lsq_obj.block(edge->getEdgeObjectiveIdx(), vertex->getVertexIdx(), edge->getObjectiveDimension(), vert_dim_unfixed);
                edge->computeJacobians(i, block_lsq_obj, block_eq, block_ineq, mult_obj_lsq_part, mult_eq_part, mult_ineq_part);
            }
            else
            {
                Eigen::MatrixXd block_lsq_obj(edge->getObjectiveDimension(), vert_dim_unfixed);
                edge->computeJacobians(i, block_lsq_obj, block_eq, block_ineq, mult_obj_lsq_part, mult_eq_part, mult_ineq_part);

                gradient_non_lsq_obj.segment(vertex->getVertexIdx(), vert_dim_unfixed) +=
                    block_lsq_obj.colwise().sum();  // sum up all values, we need a dimension of one
            }

            if (active_ineq)
            {
                // compute values and check which values are active
                // TODO(roesmann): we could obtain values computed before calling computeDenseJacobianActiveInequalities();
                Eigen::VectorXd ineq_values(edge->getInequalityDimension());
                edge->precompute();
                edge->computeInequalityValues(ineq_values);
                Eigen::Array<bool, -1, 1> active = ineq_values.array() > 0.0;

                // now set values to zero if inactive or multiply with weight otherwise
                for (int j = 0; j < edge->getInequalityDimension(); ++j)
                {
                    if (active[j] && active_ineq_weight != 1.0)
                        jacobian_ineq.block(edge->getEdgeInequalityIdx() + j, vertex->getVertexIdx(), 1, vert_dim_unfixed) *= active_ineq_weight;
                    else if (!active[j])
                        jacobian_ineq.block(edge->getEdgeInequalityIdx() + j, vertex->getVertexIdx(), 1, vert_dim_unfixed).setZero();
                }
            }
        }
    }
}

void HyperGraphOptimizationProblemEdgeBased::computeSparseJacobiansValues(Eigen::Ref<Eigen::VectorXd> values_obj,
                                                                          Eigen::Ref<Eigen::VectorXd> values_eq,
                                                                          Eigen::Ref<Eigen::VectorXd> values_ineq, const double* multipliers_obj,
                                                                          const double* multipliers_eq, const double* multipliers_ineq,
                                                                          bool active_ineq, double active_ineq_weight)
{
    assert(values_obj.size() == computeSparseJacobianLsqObjectiveNNZ());
    assert(values_eq.size() == computeSparseJacobianEqualitiesNNZ());
    assert(values_ineq.size() == computeSparseJacobianInequalitiesNNZ());

    PRINT_DEBUG_COND_ONCE(
        active_ineq && multipliers_ineq,
        "HyperGraphOptimizationProblemEdgeBased::computeSparseJacobiansValues(): multiplier_ineq is ignored if active_ineq is enabled");

    OptimizationEdgeSet::Ptr edges = _graph.getEdgeSet();

    int obj_nnz_idx  = 0;
    int eq_nnz_idx   = 0;
    int ineq_nnz_idx = 0;

    // Iterate lsq objective edges
    for (BaseEdge::Ptr& edge : edges->getLsqObjectiveEdgesRef())
    {
        const double* mult_obj_lsq_part = multipliers_obj ? multipliers_obj + edge->getEdgeIdx() : nullptr;

        for (int vert_idx = 0; vert_idx < edge->getNumVertices(); ++vert_idx)
        {
            int vert_dim_unfixed = edge->getVertexRaw(vert_idx)->getDimensionUnfixed();
            if (vert_dim_unfixed == 0) continue;

            Eigen::Map<Eigen::MatrixXd> block_jacobian(values_obj.data() + obj_nnz_idx, edge->getDimension(), vert_dim_unfixed);
            edge->computeJacobian(vert_idx, block_jacobian, mult_obj_lsq_part);
            obj_nnz_idx += block_jacobian.rows() * block_jacobian.cols();
        }
    }

    // Iterate equality edges
    for (BaseEdge::Ptr& edge : edges->getEqualityEdgesRef())
    {
        const double* mult_eq_part = multipliers_eq ? multipliers_eq + edge->getEdgeIdx() : nullptr;

        for (int vert_idx = 0; vert_idx < edge->getNumVertices(); ++vert_idx)
        {
            int vert_dim_unfixed = edge->getVertexRaw(vert_idx)->getDimensionUnfixed();
            if (vert_dim_unfixed == 0) continue;

            Eigen::Map<Eigen::MatrixXd> block_jacobian(values_eq.data() + eq_nnz_idx, edge->getDimension(), vert_dim_unfixed);
            edge->computeJacobian(vert_idx, block_jacobian, mult_eq_part);
            eq_nnz_idx += block_jacobian.rows() * block_jacobian.cols();
        }
    }

    // Iterate inequality edges
    for (BaseEdge::Ptr& edge : edges->getInequalityEdgesRef())
    {
        const double* mult_ineq_part = multipliers_ineq ? multipliers_ineq + edge->getEdgeIdx() : nullptr;

        for (int vert_idx = 0; vert_idx < edge->getNumVertices(); ++vert_idx)
        {
            int vert_dim_unfixed = edge->getVertexRaw(vert_idx)->getDimensionUnfixed();
            if (vert_dim_unfixed == 0) continue;

            if (active_ineq)
            {
                // compute values and check which values are active
                // TODO(roesmann): we could obtain values computed before calling computeDenseJacobianActiveInequalities();
                Eigen::VectorXd edge_values(edge->getDimension());
                edge->computeValues(edge_values);
                Eigen::Array<bool, -1, 1> active = edge_values.array() > 0.0;

                Eigen::Map<Eigen::MatrixXd> block_jacobian(values_ineq.data() + ineq_nnz_idx, edge->getDimension(), vert_dim_unfixed);
                ineq_nnz_idx += block_jacobian.rows() * block_jacobian.cols();

                if (!active.any())
                {
                    block_jacobian.setZero();
                    continue;
                }

                edge->computeJacobian(vert_idx, block_jacobian, nullptr);

                // now set values to zero if inactive or multiply with weight otherwise
                for (int j = 0; j < edge->getDimension(); ++j)
                {
                    if (active[j] && active_ineq_weight != 1.0)
                        block_jacobian.block(j, 0, 1, vert_dim_unfixed) *= active_ineq_weight;
                    else if (!active[j])
                        block_jacobian.block(j, 0, 1, vert_dim_unfixed).setZero();
                }
            }
            else
            {
                Eigen::Map<Eigen::MatrixXd> block_jacobian(values_ineq.data() + ineq_nnz_idx, edge->getDimension(), vert_dim_unfixed);
                edge->computeJacobian(vert_idx, block_jacobian, mult_ineq_part);
                ineq_nnz_idx += block_jacobian.rows() * block_jacobian.cols();
            }
        }
    }

    // Iterate mixed edges
    for (BaseMixedEdge::Ptr& edge : edges->getMixedEdgesRef())
    {
        // we currently need no gradient of the non-lsq terms
        if (!edge->isObjectiveLeastSquaresForm() && edge->getEqualityDimension() == 0 && edge->getInequalityDimension() == 0) continue;

        const double* mult_obj_lsq_part = multipliers_obj ? multipliers_obj + edge->getEdgeObjectiveIdx() : nullptr;
        const double* mult_eq_part      = multipliers_eq ? multipliers_eq + edge->getEdgeEqualityIdx() : nullptr;
        const double* mult_ineq_part    = (multipliers_ineq && !active_ineq) ? multipliers_ineq + edge->getEdgeInequalityIdx() : nullptr;

        for (int i = 0; i < edge->getNumVertices(); ++i)
        {
            const VertexInterface* vertex = edge->getVertexRaw(i);
            int vert_dim_unfixed          = vertex->getDimensionUnfixed();
            if (vert_dim_unfixed == 0) continue;

            Eigen::Map<Eigen::MatrixXd> block_eq(values_eq.data() + eq_nnz_idx, edge->getEqualityDimension(), vert_dim_unfixed);
            eq_nnz_idx += block_eq.rows() * block_eq.cols();

            Eigen::Map<Eigen::MatrixXd> block_ineq(values_ineq.data() + ineq_nnz_idx, edge->getInequalityDimension(), vert_dim_unfixed);
            ineq_nnz_idx += block_ineq.rows() * block_ineq.cols();

            if (edge->isObjectiveLeastSquaresForm())
            {
                Eigen::Map<Eigen::MatrixXd> block_lsq_obj(values_obj.data() + obj_nnz_idx, edge->getObjectiveDimension(), vert_dim_unfixed);
                obj_nnz_idx += block_lsq_obj.rows() * block_lsq_obj.cols();

                edge->computeJacobians(i, block_lsq_obj, block_eq, block_ineq, mult_obj_lsq_part, mult_eq_part, mult_ineq_part);
            }
            else
            {
                // Do nothing, the sparse version currently computes no gradient of the non-lsq parts
                // TODO(roesmann): the following matrix is useless and might be better removed (or the method should also compute the gradient)
                Eigen::MatrixXd block_lsq_obj(edge->getObjectiveDimension(), vert_dim_unfixed);
                edge->computeJacobians(i, block_lsq_obj, block_eq, block_ineq, mult_obj_lsq_part, mult_eq_part, mult_ineq_part);

                //                Eigen::MatrixXd block_lsq_obj(edge->getObjectiveDimension(), vert_dim_unfixed);
                //                edge->computeJacobians(i, block_lsq_obj, block_eq, block_ineq, mult_obj_lsq_part, mult_eq_part, mult_ineq_part);

                //                gradient_non_lsq_obj.segment(vertex->getVertexIdx(), vert_dim_unfixed) +=
                //                    block_lsq_obj.colwise().sum();  // sum up all values, we need a dimension of one
            }

            if (active_ineq)
            {
                // compute values and check which values are active
                // TODO(roesmann): we could obtain values computed before calling computeDenseJacobianActiveInequalities();
                Eigen::VectorXd ineq_values(edge->getInequalityDimension());
                edge->precompute();
                edge->computeInequalityValues(ineq_values);
                Eigen::Array<bool, -1, 1> active = ineq_values.array() > 0.0;

                // now set values to zero if inactive or multiply with weight otherwise
                for (int j = 0; j < edge->getInequalityDimension(); ++j)
                {
                    if (active[j] && active_ineq_weight != 1.0)
                        block_ineq.row(j) *= active_ineq_weight;
                    else if (!active[j])
                        block_ineq.row(j).setZero();
                }
            }
        }
    }
}

int HyperGraphOptimizationProblemEdgeBased::computeCombinedSparseJacobiansNNZ(bool objective_lsq, bool equality, bool inequality)
{
    int nnz = 0;
    if (objective_lsq) nnz += computeSparseJacobianLsqObjectiveNNZ();
    if (equality) nnz += computeSparseJacobianEqualitiesNNZ();
    if (inequality) nnz += computeSparseJacobianInequalitiesNNZ();
    return nnz;
}
void HyperGraphOptimizationProblemEdgeBased::computeCombinedSparseJacobiansStructure(Eigen::Ref<Eigen::VectorXi> i_row,
                                                                                     Eigen::Ref<Eigen::VectorXi> j_col, bool objective_lsq,
                                                                                     bool equality, bool inequality)
{
    assert(i_row.size() == computeCombinedSparseJacobiansNNZ(objective_lsq, equality, inequality));
    assert(j_col.size() == i_row.size());

    OptimizationEdgeSet::Ptr edges = _graph.getEdgeSet();

    int nnz_idx = 0;

    int equality_row_start   = objective_lsq ? getLsqObjectiveDimension() : 0;
    int inequality_row_start = equality ? equality_row_start + getEqualityDimension() : equality_row_start;

    if (objective_lsq)
    {
        // Iterate lsq objective edges
        for (BaseEdge::Ptr& edge : edges->getLsqObjectiveEdgesRef())
        {
            for (int vert_idx = 0; vert_idx < edge->getNumVertices(); ++vert_idx)
            {
                const VertexInterface* vertex = edge->getVertexRaw(vert_idx);
                // Iterate all free variables
                int idx_free = 0;
                for (int i = 0; i < vertex->getDimension(); ++i)
                {
                    if (!vertex->isFixedComponent(i))
                    {
                        // Iterate inner edge dimension
                        for (int j = 0; j < edge->getDimension(); ++j)
                        {
                            i_row[nnz_idx] = edge->getEdgeIdx() + j;
                            j_col[nnz_idx] = vertex->getVertexIdx() + idx_free;
                            ++nnz_idx;
                        }
                        ++idx_free;
                    }
                }
            }
        }
    }
    if (equality)
    {
        // Iterate equality edges
        for (BaseEdge::Ptr& edge : edges->getEqualityEdgesRef())
        {
            for (int vert_idx = 0; vert_idx < edge->getNumVertices(); ++vert_idx)
            {
                const VertexInterface* vertex = edge->getVertexRaw(vert_idx);
                // Iterate all free variables
                int idx_free = 0;
                for (int i = 0; i < vertex->getDimension(); ++i)
                {
                    if (!vertex->isFixedComponent(i))
                    {
                        // Iterate inner edge dimension
                        for (int j = 0; j < edge->getDimension(); ++j)
                        {
                            i_row[nnz_idx] = equality_row_start + edge->getEdgeIdx() + j;
                            j_col[nnz_idx] = vertex->getVertexIdx() + idx_free;
                            ++nnz_idx;
                        }
                        ++idx_free;
                    }
                }
            }
        }
    }
    if (inequality)
    {
        // Iterate equality edges
        for (BaseEdge::Ptr& edge : edges->getInequalityEdgesRef())
        {
            for (int vert_idx = 0; vert_idx < edge->getNumVertices(); ++vert_idx)
            {
                const VertexInterface* vertex = edge->getVertexRaw(vert_idx);
                // Iterate all free variables
                int idx_free = 0;
                for (int i = 0; i < vertex->getDimension(); ++i)
                {
                    if (!vertex->isFixedComponent(i))
                    {
                        // Iterate inner edge dimension
                        for (int j = 0; j < edge->getDimension(); ++j)
                        {
                            i_row[nnz_idx] = inequality_row_start + edge->getEdgeIdx() + j;
                            j_col[nnz_idx] = vertex->getVertexIdx() + idx_free;
                            ++nnz_idx;
                        }
                        ++idx_free;
                    }
                }
            }
        }
    }

    // Iterate mixed edges
    for (BaseMixedEdge::Ptr& edge : edges->getMixedEdgesRef())
    {
        for (int vert_idx = 0; vert_idx < edge->getNumVertices(); ++vert_idx)
        {
            const VertexInterface* vertex = edge->getVertexRaw(vert_idx);

            // we need to do the following in the same order as in computeCombinedSparseJacobiansValues()
            if (objective_lsq)
            {
                // Iterate all free variables
                int idx_free = 0;
                for (int i = 0; i < vertex->getDimension(); ++i)
                {
                    if (!vertex->isFixedComponent(i))
                    {
                        for (int j = 0; j < edge->getObjectiveDimension(); ++j)
                        {
                            i_row[nnz_idx] = edge->getEdgeObjectiveIdx() + j;
                            j_col[nnz_idx] = vertex->getVertexIdx() + idx_free;
                            ++nnz_idx;
                        }
                        ++idx_free;
                    }
                }
            }

            if (equality)
            {
                // Iterate all free variables
                int idx_free = 0;
                for (int i = 0; i < vertex->getDimension(); ++i)
                {
                    if (!vertex->isFixedComponent(i))
                    {
                        for (int j = 0; j < edge->getEqualityDimension(); ++j)
                        {
                            i_row[nnz_idx] = equality_row_start + edge->getEdgeEqualityIdx() + j;
                            j_col[nnz_idx] = vertex->getVertexIdx() + idx_free;
                            ++nnz_idx;
                        }
                        ++idx_free;
                    }
                }
            }

            if (inequality)
            {
                // Iterate all free variables
                int idx_free = 0;
                for (int i = 0; i < vertex->getDimension(); ++i)
                {
                    if (!vertex->isFixedComponent(i))
                    {
                        for (int j = 0; j < edge->getInequalityDimension(); ++j)
                        {
                            i_row[nnz_idx] = inequality_row_start + edge->getEdgeInequalityIdx() + j;
                            j_col[nnz_idx] = vertex->getVertexIdx() + idx_free;
                            ++nnz_idx;
                        }
                        ++idx_free;
                    }
                }
            }
        }
    }
}
void HyperGraphOptimizationProblemEdgeBased::computeCombinedSparseJacobiansValues(Eigen::Ref<Eigen::VectorXd> values, bool objective_lsq,
                                                                                  bool equality, bool inequality, const double* multipliers_obj,
                                                                                  const double* multipliers_eq, const double* multipliers_ineq)
{
    assert(values.size() == computeCombinedSparseJacobiansNNZ(objective_lsq, equality, inequality));

    int nnz_idx = 0;

    OptimizationEdgeSet::Ptr edges = _graph.getEdgeSet();

    if (objective_lsq)
    {
        // Iterate lsq objective edges
        for (BaseEdge::Ptr& edge : edges->getLsqObjectiveEdgesRef())
        {
            const double* mult_obj_lsq_part = multipliers_obj ? multipliers_obj + edge->getEdgeIdx() : nullptr;

            for (int vert_idx = 0; vert_idx < edge->getNumVertices(); ++vert_idx)
            {
                int vert_dim_unfixed = edge->getVertexRaw(vert_idx)->getDimensionUnfixed();
                if (vert_dim_unfixed == 0) continue;

                Eigen::Map<Eigen::MatrixXd> block_jacobian(values.data() + nnz_idx, edge->getDimension(), vert_dim_unfixed);
                edge->computeJacobian(vert_idx, block_jacobian, mult_obj_lsq_part);
                nnz_idx += block_jacobian.rows() * block_jacobian.cols();
            }
        }
    }

    if (equality)
    {
        // Iterate lsq objective edges
        for (BaseEdge::Ptr& edge : edges->getEqualityEdgesRef())
        {
            const double* mult_eq_part = multipliers_eq ? multipliers_eq + edge->getEdgeIdx() : nullptr;

            for (int vert_idx = 0; vert_idx < edge->getNumVertices(); ++vert_idx)
            {
                int vert_dim_unfixed = edge->getVertexRaw(vert_idx)->getDimensionUnfixed();
                if (vert_dim_unfixed == 0) continue;

                Eigen::Map<Eigen::MatrixXd> block_jacobian(values.data() + nnz_idx, edge->getDimension(), vert_dim_unfixed);
                edge->computeJacobian(vert_idx, block_jacobian, mult_eq_part);
                nnz_idx += block_jacobian.rows() * block_jacobian.cols();
            }
        }
    }

    if (inequality)
    {
        // Iterate lsq objective edges
        for (BaseEdge::Ptr& edge : edges->getInequalityEdgesRef())
        {
            const double* mult_ineq_part = multipliers_ineq ? multipliers_ineq + edge->getEdgeIdx() : nullptr;

            for (int vert_idx = 0; vert_idx < edge->getNumVertices(); ++vert_idx)
            {
                int vert_dim_unfixed = edge->getVertexRaw(vert_idx)->getDimensionUnfixed();
                if (vert_dim_unfixed == 0) continue;

                Eigen::Map<Eigen::MatrixXd> block_jacobian(values.data() + nnz_idx, edge->getDimension(), vert_dim_unfixed);
                edge->computeJacobian(vert_idx, block_jacobian, mult_ineq_part);
                nnz_idx += block_jacobian.rows() * block_jacobian.cols();
            }
        }
    }

    // Iterate mixed edges
    for (BaseMixedEdge::Ptr& edge : edges->getMixedEdgesRef())
    {
        if (!objective_lsq && edge->getEqualityDimension() == 0 && edge->getInequalityDimension() == 0) continue;
        // TODO(roesmann) implement the other cases as well

        const double* mult_obj_lsq_part = multipliers_obj ? multipliers_obj + edge->getEdgeObjectiveIdx() : nullptr;
        const double* mult_eq_part      = multipliers_eq ? multipliers_eq + edge->getEdgeEqualityIdx() : nullptr;
        const double* mult_ineq_part    = multipliers_ineq ? multipliers_ineq + edge->getEdgeInequalityIdx() : nullptr;

        for (int vert_idx = 0; vert_idx < edge->getNumVertices(); ++vert_idx)
        {
            int vert_dim_unfixed = edge->getVertexRaw(vert_idx)->getDimensionUnfixed();
            if (vert_dim_unfixed == 0) continue;

            if (objective_lsq && equality && inequality)
            {
                Eigen::Map<Eigen::MatrixXd> block_obj(values.data() + nnz_idx, edge->getObjectiveDimension(), vert_dim_unfixed);
                nnz_idx += block_obj.rows() * block_obj.cols();
                Eigen::Map<Eigen::MatrixXd> block_eq(values.data() + nnz_idx, edge->getEqualityDimension(), vert_dim_unfixed);
                nnz_idx += block_eq.rows() * block_eq.cols();
                Eigen::Map<Eigen::MatrixXd> block_ineq(values.data() + nnz_idx, edge->getInequalityDimension(), vert_dim_unfixed);
                nnz_idx += block_ineq.rows() * block_ineq.cols();

                edge->computeJacobians(vert_idx, block_obj, block_eq, block_ineq, mult_obj_lsq_part, mult_eq_part, mult_ineq_part);
            }
            else if (equality && inequality)
            {
                Eigen::Map<Eigen::MatrixXd> block_eq(values.data() + nnz_idx, edge->getEqualityDimension(), vert_dim_unfixed);
                nnz_idx += block_eq.rows() * block_eq.cols();
                Eigen::Map<Eigen::MatrixXd> block_ineq(values.data() + nnz_idx, edge->getInequalityDimension(), vert_dim_unfixed);
                nnz_idx += block_ineq.rows() * block_ineq.cols();

                edge->computeConstraintJacobians(vert_idx, block_eq, block_ineq, mult_eq_part, mult_ineq_part);
            }
            else
            {
                // TODO(roesmann) this is not really efficient yet
                if (objective_lsq)
                {
                    Eigen::Map<Eigen::MatrixXd> block_obj(values.data() + nnz_idx, edge->getObjectiveDimension(), vert_dim_unfixed);
                    nnz_idx += block_obj.rows() * block_obj.cols();
                    edge->computeObjectiveJacobian(vert_idx, block_obj, mult_obj_lsq_part);
                }
                if (equality)
                {
                    Eigen::Map<Eigen::MatrixXd> block_eq(values.data() + nnz_idx, edge->getEqualityDimension(), vert_dim_unfixed);
                    nnz_idx += block_eq.rows() * block_eq.cols();
                    edge->computeEqualityJacobian(vert_idx, block_eq, mult_eq_part);
                }
                if (inequality)
                {
                    Eigen::Map<Eigen::MatrixXd> block_ineq(values.data() + nnz_idx, edge->getInequalityDimension(), vert_dim_unfixed);
                    nnz_idx += block_ineq.rows() * block_ineq.cols();
                    edge->computeInequalityJacobian(vert_idx, block_ineq, mult_ineq_part);
                }
            }
        }
    }
}

void HyperGraphOptimizationProblemEdgeBased::computeCombinedSparseJacobian(Eigen::SparseMatrix<double>& jacobian, bool objective_lsq, bool equality,
                                                                           bool inequality, bool finite_combined_bounds, bool active_ineq,
                                                                           double weight_eq, double weight_ineq, double weight_bounds,
                                                                           const Eigen::VectorXd* values, const Eigen::VectorXi* col_nnz)
{
    jacobian.setZero();  // this might be obsolete as we explicitly set numeric zeros below.
    // TODO(roesmann): what is most efficient, if we just want to udpate the Jacobian?
    if (col_nnz) jacobian.reserve(*col_nnz);

    OptimizationEdgeSet::Ptr edges = _graph.getEdgeSet();

    int equality_row_start   = objective_lsq ? getLsqObjectiveDimension() : 0;
    int inequality_row_start = equality ? equality_row_start + getEqualityDimension() : equality_row_start;
    int bounds_start         = inequality ? inequality_row_start + getInequalityDimension() : inequality_row_start;

    if (objective_lsq && getLsqObjectiveDimension() > 0)
    {
        // Iterate lsq objective edges
        for (BaseEdge::Ptr& edge : edges->getLsqObjectiveEdgesRef())
        {
            for (int vert_idx = 0; vert_idx < edge->getNumVertices(); ++vert_idx)
            {
                const VertexInterface* vertex = edge->getVertexRaw(vert_idx);

                int vert_dim_unfixed = vertex->getDimensionUnfixed();
                if (vert_dim_unfixed == 0) continue;

                Eigen::MatrixXd block_jacobian(edge->getDimension(), vert_dim_unfixed);
                edge->computeJacobian(vert_idx, block_jacobian);

                // Iterate all free variables
                int idx_free = 0;
                for (int i = 0; i < vertex->getDimension(); ++i)
                {
                    if (!vertex->isFixedComponent(i))
                    {
                        // Iterate inner edge dimension
                        for (int j = 0; j < edge->getDimension(); ++j)
                        {
                            jacobian.coeffRef(edge->getEdgeIdx() + j, vertex->getVertexIdx() + idx_free) = block_jacobian(j, idx_free);
                        }
                        ++idx_free;
                    }
                }
            }
        }
    }

    if (equality && getEqualityDimension() > 0)
    {
        // Iterate equality edges
        for (BaseEdge::Ptr& edge : edges->getEqualityEdgesRef())
        {
            for (int vert_idx = 0; vert_idx < edge->getNumVertices(); ++vert_idx)
            {
                const VertexInterface* vertex = edge->getVertexRaw(vert_idx);

                int vert_dim_unfixed = vertex->getDimensionUnfixed();
                if (vert_dim_unfixed == 0) continue;

                Eigen::MatrixXd block_jacobian(edge->getDimension(), vert_dim_unfixed);
                edge->computeJacobian(vert_idx, block_jacobian);

                // Iterate all free variables
                int idx_free = 0;
                for (int i = 0; i < vertex->getDimension(); ++i)
                {
                    if (!vertex->isFixedComponent(i))
                    {
                        // Iterate inner edge dimension
                        for (int j = 0; j < edge->getDimension(); ++j)
                        {
                            jacobian.coeffRef(equality_row_start + edge->getEdgeIdx() + j, vertex->getVertexIdx() + idx_free) =
                                block_jacobian(j, idx_free) * weight_eq;
                        }
                        ++idx_free;
                    }
                }
            }
        }
    }

    if (inequality && getInequalityDimension() > 0)
    {
        // Iterate inequality edges
        for (BaseEdge::Ptr& edge : edges->getInequalityEdgesRef())
        {
            // check if active
            Eigen::Array<bool, -1, 1> active(edge->getDimension());
            if (active_ineq)
            {
                if (values)
                    active = values->segment(inequality_row_start + edge->getEdgeIdx(), edge->getDimension()).array() > 0.0;
                else
                {
                    Eigen::VectorXd values_tmp(edge->getDimension());
                    edge->computeValues(values_tmp);
                    active = values_tmp.array() > 0.0;
                }
            }
            else
                active.setConstant(true);

            for (int vert_idx = 0; vert_idx < edge->getNumVertices(); ++vert_idx)
            {
                const VertexInterface* vertex = edge->getVertexRaw(vert_idx);

                int vert_dim_unfixed = vertex->getDimensionUnfixed();
                if (vert_dim_unfixed == 0) continue;

                Eigen::MatrixXd block_jacobian(edge->getDimension(), vert_dim_unfixed);
                edge->computeJacobian(vert_idx, block_jacobian);

                // Iterate all free variables
                int idx_free = 0;
                for (int i = 0; i < vertex->getDimension(); ++i)
                {
                    if (!vertex->isFixedComponent(i))
                    {
                        // Iterate inner edge dimension
                        for (int j = 0; j < edge->getDimension(); ++j)
                        {
                            if (active[j])
                            {
                                jacobian.coeffRef(inequality_row_start + edge->getEdgeIdx() + j, vertex->getVertexIdx() + idx_free) =
                                    block_jacobian(j, idx_free) * weight_ineq;
                            }
                            else
                            {
                                jacobian.coeffRef(inequality_row_start + edge->getEdgeIdx() + j, vertex->getVertexIdx() + idx_free) = 0.0;
                            }
                        }
                        ++idx_free;
                    }
                }
            }
        }
    }

    // Iterate mixed edges
    for (BaseMixedEdge::Ptr& edge : edges->getMixedEdgesRef())
    {
        if (!objective_lsq && edge->getEqualityDimension() == 0 && edge->getInequalityDimension() == 0) continue;
        // TODO(roesmann) implement the other cases as well

        // check if active
        Eigen::Array<bool, -1, 1> active(edge->getInequalityDimension());
        if (active_ineq && inequality)
        {
            if (values)
                active = values->segment(inequality_row_start + edge->getEdgeInequalityIdx(), edge->getInequalityDimension()).array() > 0.0;
            else
            {
                Eigen::VectorXd values_tmp(edge->getInequalityDimension());
                edge->computeInequalityValues(values_tmp);
                active = values_tmp.array() > 0.0;
            }
        }
        else if (inequality)
            active.setConstant(true);

        for (int vert_idx = 0; vert_idx < edge->getNumVertices(); ++vert_idx)
        {
            const VertexInterface* vertex = edge->getVertexRaw(vert_idx);

            int vert_dim_unfixed = vertex->getDimensionUnfixed();
            if (vert_dim_unfixed == 0) continue;

            Eigen::MatrixXd block_obj(edge->getObjectiveDimension(), vert_dim_unfixed);
            Eigen::MatrixXd block_eq(edge->getEqualityDimension(), vert_dim_unfixed);
            Eigen::MatrixXd block_ineq(edge->getInequalityDimension(), vert_dim_unfixed);

            if (objective_lsq && equality && inequality)
            {
                edge->computeJacobians(vert_idx, block_obj, block_eq, block_ineq);
            }
            else if (equality && inequality)
            {
                edge->computeConstraintJacobians(vert_idx, block_eq, block_ineq);
            }
            else
            {
                // TODO(roesmann) this is not really efficient yet
                if (objective_lsq)
                {
                    edge->computeObjectiveJacobian(vert_idx, block_obj);
                }
                if (equality)
                {
                    edge->computeEqualityJacobian(vert_idx, block_eq);
                }
                if (inequality)
                {
                    edge->computeInequalityJacobian(vert_idx, block_ineq);
                }
            }

            // Iterate all free variables
            int idx_free = 0;

            for (int i = 0; i < vertex->getDimension(); ++i)
            {
                if (!vertex->isFixedComponent(i))
                {
                    // Iterate inner edge dimension
                    if (objective_lsq)
                    {
                        for (int j = 0; j < edge->getObjectiveDimension(); ++j)
                        {
                            jacobian.coeffRef(edge->getEdgeObjectiveIdx() + j, vertex->getVertexIdx() + idx_free) = block_obj(j, idx_free);
                        }
                    }
                    if (equality)
                    {
                        for (int j = 0; j < edge->getEqualityDimension(); ++j)
                        {
                            jacobian.coeffRef(equality_row_start + edge->getEdgeEqualityIdx() + j, vertex->getVertexIdx() + idx_free) =
                                block_eq(j, idx_free) * weight_eq;
                        }
                    }
                    if (inequality)
                    {
                        for (int j = 0; j < edge->getInequalityDimension(); ++j)
                        {
                            if (active[j])
                            {
                                jacobian.coeffRef(inequality_row_start + edge->getEdgeInequalityIdx() + j, vertex->getVertexIdx() + idx_free) =
                                    block_ineq(j, idx_free) * weight_ineq;
                            }
                            else
                            {
                                jacobian.coeffRef(inequality_row_start + edge->getEdgeInequalityIdx() + j, vertex->getVertexIdx() + idx_free) = 0.0;
                            }
                        }
                    }
                    ++idx_free;
                }
            }
        }
    }

    if (finite_combined_bounds && finiteCombinedBoundsDimension() > 0)
    {
        // we have a single value per row
        int jac_row_idx = bounds_start;
        // Iterate vertices
        for (const VertexInterface* vertex : _graph.getVertexSet()->getActiveVertices())
        {
            int vert_idx = vertex->getVertexIdx();
            int free_idx = 0;
            for (int i = 0; i < vertex->getDimension(); ++i)
            {
                if (vertex->isFixedComponent(i)) continue;

                if (vertex->hasFiniteLowerBound(i) || vertex->hasFiniteUpperBound(i))
                {
                    if (vertex->getData()[i] < vertex->getLowerBounds()[i])
                    {
                        jacobian.coeffRef(jac_row_idx, vert_idx + free_idx) = -weight_bounds;
                    }
                    else if (vertex->getData()[i] > vertex->getUpperBounds()[i])
                    {
                        jacobian.coeffRef(jac_row_idx, vert_idx + free_idx) = weight_bounds;
                    }
                    else
                        jacobian.coeffRef(jac_row_idx, vert_idx + free_idx) = 0.0;  // preserve structure

                    ++jac_row_idx;
                }
                ++free_idx;
            }
        }
    }
}

// useful for IPOPT (w/ hessian-approx)
void HyperGraphOptimizationProblemEdgeBased::computeGradientObjectiveAndCombinedSparseJacobiansValues(Eigen::Ref<Eigen::VectorXd> gradient,
                                                                                                      Eigen::Ref<Eigen::VectorXd> jac_values,
                                                                                                      bool equality, bool inequality,
                                                                                                      const double* multipliers_eq,
                                                                                                      const double* multipliers_ineq)
{
    assert(jac_values.size() == computeCombinedSparseJacobiansNNZ(false, equality, inequality));
    assert(gradient.size() == getParameterDimension());
    assert(equality == true);    // TODO(roesmann)
    assert(inequality == true);  // TODO(roesmann)

    // we need to set the gradient to zero first
    gradient.setZero();

    OptimizationEdgeSet::Ptr edges = _graph.getEdgeSet();

    // Iterate plain objective edges
    for (BaseEdge::Ptr& edge : edges->getObjectiveEdgesRef())
    {
        // if (edge->providesJacobian())
        for (int i = 0; i < edge->getNumVertices(); ++i)
        {
            int vert_dim_unfixed = edge->getVertexRaw(i)->getDimensionUnfixed();
            if (vert_dim_unfixed == 0) continue;

            Eigen::MatrixXd block_jacobian(edge->getDimension(), vert_dim_unfixed);
            edge->computeJacobian(i, block_jacobian, nullptr);
            gradient.segment(edge->getVertexRaw(i)->getVertexIdx(), vert_dim_unfixed) +=
                block_jacobian.colwise().sum();  // sum up all values, we need a dimension  of one
        }
    }
    // Iterate lsq objective edges
    for (BaseEdge::Ptr& edge : edges->getLsqObjectiveEdgesRef())
    {
        for (int i = 0; i < edge->getNumVertices(); ++i)
        {
            int vert_dim_unfixed = edge->getVertexRaw(i)->getDimensionUnfixed();
            if (vert_dim_unfixed == 0) continue;

            Eigen::MatrixXd block_jacobian(edge->getDimension(), vert_dim_unfixed);
            edge->computeJacobian(i, block_jacobian, nullptr);

            // apply chain rule // TODO(roesmann): in case of forward differences we could use this for computing the jacobian
            Eigen::VectorXd values(edge->getDimension());
            edge->computeValues(values);

            gradient.segment(edge->getVertexRaw(i)->getVertexIdx(), vert_dim_unfixed) += 2.0 * values.transpose() * block_jacobian;
        }
    }

    int nnz_idx = 0;  // jacobian

    if (equality)
    {
        // Iterate lsq objective edges
        for (BaseEdge::Ptr& edge : edges->getEqualityEdgesRef())
        {
            const double* mult_eq_part = multipliers_eq ? multipliers_eq + edge->getEdgeIdx() : nullptr;

            for (int vert_idx = 0; vert_idx < edge->getNumVertices(); ++vert_idx)
            {
                int vert_dim_unfixed = edge->getVertexRaw(vert_idx)->getDimensionUnfixed();
                if (vert_dim_unfixed == 0) continue;

                Eigen::Map<Eigen::MatrixXd> block_jacobian(jac_values.data() + nnz_idx, edge->getDimension(), vert_dim_unfixed);
                edge->computeJacobian(vert_idx, block_jacobian, mult_eq_part);
                nnz_idx += block_jacobian.rows() * block_jacobian.cols();
            }
        }
    }

    if (inequality)
    {
        // Iterate lsq objective edges
        for (BaseEdge::Ptr& edge : edges->getInequalityEdgesRef())
        {
            const double* mult_ineq_part = multipliers_ineq ? multipliers_ineq + edge->getEdgeIdx() : nullptr;

            for (int vert_idx = 0; vert_idx < edge->getNumVertices(); ++vert_idx)
            {
                int vert_dim_unfixed = edge->getVertexRaw(vert_idx)->getDimensionUnfixed();
                if (vert_dim_unfixed == 0) continue;

                Eigen::Map<Eigen::MatrixXd> block_jacobian(jac_values.data() + nnz_idx, edge->getDimension(), vert_dim_unfixed);
                edge->computeJacobian(vert_idx, block_jacobian, mult_ineq_part);
                nnz_idx += block_jacobian.rows() * block_jacobian.cols();
            }
        }
    }

    // Iterate mixed edges
    for (BaseMixedEdge::Ptr& edge : edges->getMixedEdgesRef())
    {
        // bool has_obj  = edge->getObjectiveDimension() > 0;
        // bool has_eq   = equality && edge->getEqualityDimension();
        // bool has_ineq = inequality && edge->getInequalityDimension();

        // ###############
        for (int i = 0; i < edge->getNumVertices(); ++i)
        {
            int vert_dim_unfixed = edge->getVertexRaw(i)->getDimensionUnfixed();
            if (vert_dim_unfixed == 0) continue;

            Eigen::MatrixXd block_jacobian(edge->getObjectiveDimension(), vert_dim_unfixed);
            edge->computeObjectiveJacobian(i, block_jacobian, nullptr);

            if (edge->isObjectiveLeastSquaresForm())
            {
                // apply chain rule // TODO(roesmann): in case of forward differences we could use this for computing the jacobian
                Eigen::VectorXd values(edge->getObjectiveDimension());
                edge->computeObjectiveValues(values);
                gradient.segment(edge->getVertexRaw(i)->getVertexIdx(), vert_dim_unfixed) += 2.0 * values.transpose() * block_jacobian;
            }
            else
            {
                gradient.segment(edge->getVertexRaw(i)->getVertexIdx(), vert_dim_unfixed) +=
                    block_jacobian.colwise().sum();  // sum up all values, we need a dimension of one
            }
        }
        // #########

        const double* mult_eq_part   = multipliers_eq ? multipliers_eq + edge->getEdgeEqualityIdx() : nullptr;
        const double* mult_ineq_part = multipliers_ineq ? multipliers_ineq + edge->getEdgeInequalityIdx() : nullptr;

        for (int vert_idx = 0; vert_idx < edge->getNumVertices(); ++vert_idx)
        {
            int vert_dim_unfixed = edge->getVertexRaw(vert_idx)->getDimensionUnfixed();
            if (vert_dim_unfixed == 0) continue;

            Eigen::MatrixXd block_obj(edge->getObjectiveDimension(), vert_dim_unfixed);
            Eigen::Map<Eigen::MatrixXd> block_eq(jac_values.data() + nnz_idx, edge->getEqualityDimension(), vert_dim_unfixed);
            nnz_idx += block_eq.rows() * block_eq.cols();
            Eigen::Map<Eigen::MatrixXd> block_ineq(jac_values.data() + nnz_idx, edge->getInequalityDimension(), vert_dim_unfixed);
            nnz_idx += block_ineq.rows() * block_ineq.cols();
            edge->computeJacobians(vert_idx, block_obj, block_eq, block_ineq, nullptr, mult_eq_part, mult_ineq_part);

            if (edge->isObjectiveLeastSquaresForm())
            {
                // apply chain rule // TODO(roesmann): in case of forward differences we could use this for computing the jacobian
                Eigen::VectorXd values(edge->getObjectiveDimension());
                edge->computeObjectiveValues(values);
                gradient.segment(edge->getVertexRaw(vert_idx)->getVertexIdx(), vert_dim_unfixed) += 2.0 * values.transpose() * block_obj;
            }
            else
            {
                gradient.segment(edge->getVertexRaw(vert_idx)->getVertexIdx(), vert_dim_unfixed) +=
                    block_obj.colwise().sum();  // sum up all values, we need a dimension of one
            }

            // if (has_obj && has_eq && has_ineq)
            // {
            // }

            //            if (equality && inequality)
            //            {
            //                Eigen::Map<Eigen::MatrixXd> block_eq(jac_values.data() + nnz_idx, edge->getEqualityDimension(), vert_dim_unfixed);
            //                nnz_idx += block_eq.rows() * block_eq.cols();
            //                Eigen::Map<Eigen::MatrixXd> block_ineq(jac_values.data() + nnz_idx, edge->getInequalityDimension(), vert_dim_unfixed);
            //                nnz_idx += block_ineq.rows() * block_ineq.cols();

            //                edge->computeConstraintJacobians(vert_idx, block_eq, block_ineq, mult_eq_part, mult_ineq_part);
            //            }
            //            else
            //            {
            //                if (equality)
            //                {
            //                    Eigen::Map<Eigen::MatrixXd> block_eq(jac_values.data() + nnz_idx, edge->getEqualityDimension(), vert_dim_unfixed);
            //                    nnz_idx += block_eq.rows() * block_eq.cols();
            //                    edge->computeEqualityJacobian(vert_idx, block_eq, mult_eq_part);
            //                }
            //                if (inequality)
            //                {
            //                    Eigen::Map<Eigen::MatrixXd> block_ineq(jac_values.data() + nnz_idx, edge->getInequalityDimension(),
            //                    vert_dim_unfixed);
            //                    nnz_idx += block_ineq.rows() * block_ineq.cols();
            //                    edge->computeEqualityJacobian(vert_idx, block_ineq, mult_ineq_part);
            //                }
            //            }
        }
    }
}

void HyperGraphOptimizationProblemEdgeBased::computeDenseHessianObjective(Eigen::Ref<Eigen::MatrixXd> hessian, double multiplier)
{
    assert(_graph.hasEdgeSet());
    if (!_graph_precomputed) precomputeGraphQuantities();
    assert(hessian.rows() == getParameterDimension());
    assert(hessian.rows() == hessian.cols());

    // we need to set to zero everything first
    // but we only compute the upper hessian via edges
    hessian.triangularView<Eigen::Upper>().setZero();

    if (getObjectiveDimension() == 0) return;

    OptimizationEdgeSet::Ptr edges = _graph.getEdgeSet();

    // Iterate plain objective edges
    for (BaseEdge::Ptr& edge : edges->getObjectiveEdgesRef())
    {
        if (edge->isLinear()) continue;

        for (int row_i = 0; row_i < edge->getNumVertices(); ++row_i)
        {
            const VertexInterface* vertex1 = edge->getVertexRaw(row_i);

            int vert_i_dim_unfixed = vertex1->getDimensionUnfixed();
            if (vert_i_dim_unfixed == 0) continue;

            // compute jacobian to speed up hessian computation
            Eigen::MatrixXd block_jacobian1(edge->getDimension(), vertex1->getDimensionUnfixed());
            edge->computeJacobian(row_i, block_jacobian1);

            for (int col_j = row_i; col_j < edge->getNumVertices(); ++col_j)  // hessian is symmetric
            {
                const VertexInterface* vertex2 = edge->getVertexRaw(col_j);

                int vert_j_dim_unfixed = vertex2->getDimensionUnfixed();
                if (vert_j_dim_unfixed == 0) continue;

                // warning: we require the hessian to be initialized with zeros
                // and that edge->computeHessian only adds (+=) values!!
                edge->computeHessianInc(row_i, col_j, block_jacobian1,
                                        hessian.block(vertex1->getVertexIdx(), vertex2->getVertexIdx(), vert_i_dim_unfixed, vert_j_dim_unfixed),
                                        nullptr, multiplier);
            }
        }
    }

    // Iterate lsq objective edges
    for (BaseEdge::Ptr& edge : edges->getLsqObjectiveEdgesRef())
    {
        // if (edge->isLinear()) continue; // we must not check this one, since Linear-squared is not linear ;-)

        for (int row_i = 0; row_i < edge->getNumVertices(); ++row_i)
        {
            const VertexInterface* vertex1 = edge->getVertexRaw(row_i);

            int vert_i_dim_unfixed = vertex1->getDimensionUnfixed();
            if (vert_i_dim_unfixed == 0) continue;

            // compute jacobian to speed up hessian computation
            Eigen::MatrixXd block_jacobian1(edge->getDimension(), vertex1->getDimensionUnfixed());
            edge->computeJacobian(row_i, block_jacobian1);

            for (int col_j = row_i; col_j < edge->getNumVertices(); ++col_j)  // hessian is symmetric
            {
                const VertexInterface* vertex2 = edge->getVertexRaw(col_j);

                int vert_j_dim_unfixed = vertex2->getDimensionUnfixed();
                if (vert_j_dim_unfixed == 0) continue;

                // we need the second jacobian for applying the chain rule
                Eigen::MatrixXd block_jacobian2(edge->getDimension(), vertex2->getDimensionUnfixed());
                edge->computeJacobian(col_j, block_jacobian2);

                // approximate H = 2*J^T*J
                // if (_hessian_lsq_obj_approx)
                //{
                hessian.block(vertex1->getVertexIdx(), vertex2->getVertexIdx(), vert_i_dim_unfixed, vert_j_dim_unfixed) +=
                    2.0 * multiplier * block_jacobian1.transpose() * block_jacobian2;
                //}
                //                else
                //                {
                //                    Eigen::MatrixXd block_hessian(vertex1->getDimensionUnfixed(), vertex2->getDimensionUnfixed());
                //                    edge->computeHessianInc(row_i, col_j, block_jacobian1, block_hessian, nullptr, 1.0);  // weight is already in
                //                    the next equation

                //                    // Eigen::VectorXd edge_values(edge->getDimension());
                //                    // edge->computeValues(edge_values);

                //                    hessian.block(vertex1->getVertexIdx(), vertex2->getVertexIdx(), vert_i_dim_unfixed, vert_j_dim_unfixed) +=
                //                        2.0 * (multiplier * block_jacobian1.transpose() * block_jacobian2 + block_hessian);  // the second part is
                //                        not
                //                        correct yet
                //                }
            }
        }
    }

    // Iterate mixed edges
    for (BaseMixedEdge::Ptr& edge : edges->getMixedEdgesRef())
    {
        if (edge->getObjectiveDimension() == 0 || edge->isObjectiveLinear()) continue;

        for (int row_i = 0; row_i < edge->getNumVertices(); ++row_i)
        {
            const VertexInterface* vertex1 = edge->getVertexRaw(row_i);

            int vert_i_dim_unfixed = vertex1->getDimensionUnfixed();
            if (vert_i_dim_unfixed == 0) continue;

            // compute jacobian to speed up hessian computation
            Eigen::MatrixXd block_jacobian1(edge->getObjectiveDimension(), vertex1->getDimensionUnfixed());
            edge->computeObjectiveJacobian(row_i, block_jacobian1);

            for (int col_j = row_i; col_j < edge->getNumVertices(); ++col_j)  // hessian is symmetric
            {
                const VertexInterface* vertex2 = edge->getVertexRaw(col_j);

                int vert_j_dim_unfixed = vertex2->getDimensionUnfixed();
                if (vert_j_dim_unfixed == 0) continue;

                if (edge->isObjectiveLeastSquaresForm())
                {
                    // we need the second jacobian for applying the chain rule
                    Eigen::MatrixXd block_jacobian2(edge->getObjectiveDimension(), vertex2->getDimensionUnfixed());
                    edge->computeObjectiveJacobian(col_j, block_jacobian2);

                    // approximate H = 2*J^T*J
                    // if (_hessian_lsq_obj_approx)
                    //{
                    hessian.block(vertex1->getVertexIdx(), vertex2->getVertexIdx(), vert_i_dim_unfixed, vert_j_dim_unfixed) +=
                        2.0 * multiplier * block_jacobian1.transpose() * block_jacobian2;
                }
                else
                {
                    // warning: we require the hessian to be initialized with zeros
                    // and that edge->computeHessian only adds (+=) values!!
                    edge->computeObjectiveHessianInc(
                        row_i, col_j, block_jacobian1,
                        hessian.block(vertex1->getVertexIdx(), vertex2->getVertexIdx(), vert_i_dim_unfixed, vert_j_dim_unfixed), nullptr, multiplier);
                }
            }
        }
    }

    // TODO(roesmann) at only_lower/only_upper flag to the function call in order to avoid the following copy
    hessian.triangularView<Eigen::Lower>() = hessian.adjoint().triangularView<Eigen::Lower>();
}

int HyperGraphOptimizationProblemEdgeBased::computeSparseHessianObjectiveNNZ(bool lower_part_only)
{
    int nnz = 0;

    assert(_graph.hasEdgeSet());
    if (!_graph_precomputed) precomputeGraphQuantities();

    OptimizationEdgeSet::Ptr edges = _graph.getEdgeSet();

    // Iterate plain objective edges
    for (BaseEdge::Ptr& edge : edges->getObjectiveEdgesRef())
    {
        if (edge->isLinear()) continue;

        for (int row_i = 0; row_i < edge->getNumVertices(); ++row_i)
        {
            const VertexInterface* vertex1 = edge->getVertexRaw(row_i);

            int vert_i_dim_unfixed = vertex1->getDimensionUnfixed();
            if (vert_i_dim_unfixed == 0) continue;

            int col_start = lower_part_only ? row_i : 0;
            for (int col_j = col_start; col_j < edge->getNumVertices(); ++col_j)  // hessian is symmetric
            {
                const VertexInterface* vertex2 = edge->getVertexRaw(col_j);

                int vert_j_dim_unfixed = vertex2->getDimensionUnfixed();
                if (vert_j_dim_unfixed == 0) continue;

                if (lower_part_only && vertex1 == vertex2)
                {
                    int n = vertex1->getDimensionUnfixed();
                    nnz += n * (n + 1) / 2;  // lower triangular including diagonal
                }
                else
                    nnz += vertex1->getDimensionUnfixed() * vertex2->getDimensionUnfixed();  // block size
            }
        }
    }

    // Iterate lsq objective edges
    for (BaseEdge::Ptr& edge : edges->getLsqObjectiveEdgesRef())
    {
        // if (edge->isLinear()) continue; // we must not check this one, since Linear-squared is not linear ;-)

        for (int row_i = 0; row_i < edge->getNumVertices(); ++row_i)
        {
            const VertexInterface* vertex1 = edge->getVertexRaw(row_i);

            int vert_i_dim_unfixed = vertex1->getDimensionUnfixed();
            if (vert_i_dim_unfixed == 0) continue;

            int col_start = lower_part_only ? row_i : 0;
            for (int col_j = col_start; col_j < edge->getNumVertices(); ++col_j)  // hessian is symmetric
            {
                const VertexInterface* vertex2 = edge->getVertexRaw(col_j);

                int vert_j_dim_unfixed = vertex2->getDimensionUnfixed();
                if (vert_j_dim_unfixed == 0) continue;

                if (lower_part_only && vertex1 == vertex2)
                {
                    int n = vertex1->getDimensionUnfixed();
                    nnz += n * (n + 1) / 2;  // lower triangular including diagonal
                }
                else
                    nnz += vertex1->getDimensionUnfixed() * vertex2->getDimensionUnfixed();  // block size
            }
        }
    }

    // Iterate mixed edges
    for (BaseMixedEdge::Ptr& edge : edges->getMixedEdgesRef())
    {
        if (edge->getObjectiveDimension() == 0 || edge->isObjectiveLinear()) continue;

        for (int row_i = 0; row_i < edge->getNumVertices(); ++row_i)
        {
            const VertexInterface* vertex1 = edge->getVertexRaw(row_i);

            int vert_i_dim_unfixed = vertex1->getDimensionUnfixed();
            if (vert_i_dim_unfixed == 0) continue;

            int col_start = lower_part_only ? row_i : 0;
            for (int col_j = col_start; col_j < edge->getNumVertices(); ++col_j)  // hessian is symmetric
            {
                const VertexInterface* vertex2 = edge->getVertexRaw(col_j);

                int vert_j_dim_unfixed = vertex2->getDimensionUnfixed();
                if (vert_j_dim_unfixed == 0) continue;

                if (lower_part_only && vertex1 == vertex2)
                {
                    int n = vertex1->getDimensionUnfixed();
                    nnz += n * (n + 1) / 2;  // lower triangular including diagonal
                }
                else
                    nnz += vertex1->getDimensionUnfixed() * vertex2->getDimensionUnfixed();  // block size
            }
        }
    }

    return nnz;
}

void HyperGraphOptimizationProblemEdgeBased::computeSparseHessianObjectiveStructure(Eigen::Ref<Eigen::VectorXi> i_row,
                                                                                    Eigen::Ref<Eigen::VectorXi> j_col, bool lower_part_only)
{
    assert(_graph.hasEdgeSet());
    if (!_graph_precomputed) precomputeGraphQuantities();

    assert(i_row.size() == computeSparseHessianObjectiveNNZ(lower_part_only));
    assert(i_row.size() == j_col.size());

    OptimizationEdgeSet::Ptr edges = _graph.getEdgeSet();

    int nnz_idx = 0;

    // Iterate plain objective edges
    for (BaseEdge::Ptr& edge : edges->getObjectiveEdgesRef())
    {
        if (edge->isLinear()) continue;

        for (int row_i = 0; row_i < edge->getNumVertices(); ++row_i)
        {
            const VertexInterface* vertex1 = edge->getVertexRaw(row_i);

            int vert_i_dim_unfixed = vertex1->getDimensionUnfixed();
            if (vert_i_dim_unfixed == 0) continue;

            int col_end = lower_part_only ? row_i + 1 : edge->getNumVertices();
            for (int col_j = 0; col_j < col_end; ++col_j)  // hessian is symmetric
            {
                const VertexInterface* vertex2 = edge->getVertexRaw(col_j);

                int vert_j_dim_unfixed = vertex2->getDimensionUnfixed();
                if (vert_j_dim_unfixed == 0) continue;

                if (lower_part_only && vertex1 == vertex2)
                {
                    for (int v1_idx = 0; v1_idx < vert_i_dim_unfixed; ++v1_idx)
                    {
                        for (int v2_idx = 0; v2_idx < v1_idx + 1; ++v2_idx)
                        {
                            i_row[nnz_idx] = vertex1->getVertexIdx() + v1_idx;
                            j_col[nnz_idx] = vertex2->getVertexIdx() + v2_idx;
                            ++nnz_idx;
                        }
                    }
                }
                else
                {
                    for (int v1_idx = 0; v1_idx < vert_i_dim_unfixed; ++v1_idx)
                    {
                        for (int v2_idx = 0; v2_idx < vert_j_dim_unfixed; ++v2_idx)
                        {
                            i_row[nnz_idx] = vertex1->getVertexIdx() + v1_idx;
                            j_col[nnz_idx] = vertex2->getVertexIdx() + v2_idx;
                            ++nnz_idx;
                        }
                    }
                }
            }
        }
    }

    // Iterate lsq objective edges
    for (BaseEdge::Ptr& edge : edges->getLsqObjectiveEdgesRef())
    {
        // if (edge->isLinear()) continue; // we must not check this one, since Linear-squared is not linear ;-)

        for (int row_i = 0; row_i < edge->getNumVertices(); ++row_i)
        {
            const VertexInterface* vertex1 = edge->getVertexRaw(row_i);

            int vert_i_dim_unfixed = vertex1->getDimensionUnfixed();
            if (vert_i_dim_unfixed == 0) continue;

            int col_end = lower_part_only ? row_i + 1 : edge->getNumVertices();
            for (int col_j = 0; col_j < col_end; ++col_j)  // hessian is symmetric
            {
                const VertexInterface* vertex2 = edge->getVertexRaw(col_j);

                int vert_j_dim_unfixed = vertex2->getDimensionUnfixed();
                if (vert_j_dim_unfixed == 0) continue;

                if (lower_part_only && vertex1 == vertex2)
                {
                    for (int v1_idx = 0; v1_idx < vert_i_dim_unfixed; ++v1_idx)
                    {
                        for (int v2_idx = 0; v2_idx < v1_idx + 1; ++v2_idx)
                        {
                            i_row[nnz_idx] = vertex1->getVertexIdx() + v1_idx;
                            j_col[nnz_idx] = vertex2->getVertexIdx() + v2_idx;
                            ++nnz_idx;
                        }
                    }
                }
                else
                {
                    for (int v1_idx = 0; v1_idx < vert_i_dim_unfixed; ++v1_idx)
                    {
                        for (int v2_idx = 0; v2_idx < vert_j_dim_unfixed; ++v2_idx)
                        {
                            i_row[nnz_idx] = vertex1->getVertexIdx() + v1_idx;
                            j_col[nnz_idx] = vertex2->getVertexIdx() + v2_idx;
                            ++nnz_idx;
                        }
                    }
                }
            }
        }
    }

    // Iterate mixed edges
    for (BaseMixedEdge::Ptr& edge : edges->getMixedEdgesRef())
    {
        if (edge->getObjectiveDimension() == 0 || edge->isObjectiveLinear()) continue;

        for (int row_i = 0; row_i < edge->getNumVertices(); ++row_i)
        {
            const VertexInterface* vertex1 = edge->getVertexRaw(row_i);

            int vert_i_dim_unfixed = vertex1->getDimensionUnfixed();
            if (vert_i_dim_unfixed == 0) continue;

            int col_end = lower_part_only ? row_i + 1 : edge->getNumVertices();
            for (int col_j = 0; col_j < col_end; ++col_j)  // hessian is symmetric
            {
                const VertexInterface* vertex2 = edge->getVertexRaw(col_j);

                int vert_j_dim_unfixed = vertex2->getDimensionUnfixed();
                if (vert_j_dim_unfixed == 0) continue;

                if (lower_part_only && vertex1 == vertex2)
                {
                    for (int v1_idx = 0; v1_idx < vert_i_dim_unfixed; ++v1_idx)
                    {
                        for (int v2_idx = 0; v2_idx < v1_idx + 1; ++v2_idx)
                        {
                            i_row[nnz_idx] = vertex1->getVertexIdx() + v1_idx;
                            j_col[nnz_idx] = vertex2->getVertexIdx() + v2_idx;
                            ++nnz_idx;
                        }
                    }
                }
                else
                {
                    for (int v1_idx = 0; v1_idx < vert_i_dim_unfixed; ++v1_idx)
                    {
                        for (int v2_idx = 0; v2_idx < vert_j_dim_unfixed; ++v2_idx)
                        {
                            i_row[nnz_idx] = vertex1->getVertexIdx() + v1_idx;
                            j_col[nnz_idx] = vertex2->getVertexIdx() + v2_idx;
                            ++nnz_idx;
                        }
                    }
                }
            }
        }
    }
}

void HyperGraphOptimizationProblemEdgeBased::computeSparseHessianObjectiveValues(Eigen::Ref<Eigen::VectorXd> values, double multiplier,
                                                                                 bool lower_part_only)
{
    assert(values.size() == computeSparseHessianObjectiveNNZ(lower_part_only));
    values.setZero();

    if (getObjectiveDimension() == 0) return;

    int nnz_idx = 0;

    OptimizationEdgeSet::Ptr edges = _graph.getEdgeSet();

    // Iterate plain objective edges
    for (BaseEdge::Ptr& edge : edges->getObjectiveEdgesRef())
    {
        if (edge->isLinear()) continue;

        for (int row_i = 0; row_i < edge->getNumVertices(); ++row_i)
        {
            const VertexInterface* vertex1 = edge->getVertexRaw(row_i);

            int vert_i_dim_unfixed = vertex1->getDimensionUnfixed();
            if (vert_i_dim_unfixed == 0) continue;

            // compute jacobian to speed up hessian computation
            Eigen::MatrixXd block_jacobian1(edge->getDimension(), vertex1->getDimensionUnfixed());
            edge->computeJacobian(row_i, block_jacobian1);

            int col_end = lower_part_only ? row_i + 1 : edge->getNumVertices();
            for (int col_j = 0; col_j < col_end; ++col_j)  // hessian is symmetric
            {
                const VertexInterface* vertex2 = edge->getVertexRaw(col_j);

                int vert_j_dim_unfixed = vertex2->getDimensionUnfixed();
                if (vert_j_dim_unfixed == 0) continue;

                // warning: we require the hessian to be initialized with zeros
                // and that edge->computeHessian only adds (+=) values!!
                if (lower_part_only && vertex1 == vertex2)
                {
                    // we only need the lower triangular part, so cache Hessian
                    Eigen::MatrixXd block_hessian_ii(vert_i_dim_unfixed, vert_j_dim_unfixed);
                    edge->computeHessian(row_i, col_j, block_jacobian1, block_hessian_ii, nullptr, multiplier);
                    for (int i = 0; i < vert_i_dim_unfixed; ++i)
                    {
                        for (int j = 0; j < i + 1; ++j)
                        {
                            values[nnz_idx] += block_hessian_ii.coeffRef(i, j);
                            ++nnz_idx;
                        }
                        // nnz_idx += vert_i_dim_unfixed * (vert_i_dim_unfixed + 1) / 2;
                    }
                }
                else
                {
                    Eigen::Map<Eigen::MatrixXd> block_hessian(values.data() + nnz_idx, vert_i_dim_unfixed, vert_j_dim_unfixed);
                    edge->computeHessianInc(row_i, col_j, block_jacobian1, block_hessian, nullptr, multiplier);
                    nnz_idx += vert_i_dim_unfixed * vert_j_dim_unfixed;
                }
            }
        }
    }

    // Iterate lsq objective edges
    for (BaseEdge::Ptr& edge : edges->getLsqObjectiveEdgesRef())
    {
        // if (edge->isLinear()) continue; // we must not check this one, since Linear-squared is not linear ;-)

        for (int row_i = 0; row_i < edge->getNumVertices(); ++row_i)
        {
            const VertexInterface* vertex1 = edge->getVertexRaw(row_i);

            int vert_i_dim_unfixed = vertex1->getDimensionUnfixed();
            if (vert_i_dim_unfixed == 0) continue;

            // compute jacobian to speed up hessian computation
            Eigen::MatrixXd block_jacobian1(edge->getDimension(), vertex1->getDimensionUnfixed());
            edge->computeJacobian(row_i, block_jacobian1);

            int col_end = lower_part_only ? row_i + 1 : edge->getNumVertices();
            for (int col_j = 0; col_j < col_end; ++col_j)  // hessian is symmetric
            {
                const VertexInterface* vertex2 = edge->getVertexRaw(col_j);

                int vert_j_dim_unfixed = vertex2->getDimensionUnfixed();
                if (vert_j_dim_unfixed == 0) continue;

                // we need the second jacobian for applying the chain rule
                Eigen::MatrixXd block_jacobian2(edge->getDimension(), vertex2->getDimensionUnfixed());
                edge->computeJacobian(col_j, block_jacobian2);

                // approximate H = 2*J^T*J
                if (lower_part_only && vertex1 == vertex2)
                {
                    // we only need the lower triangular part, so cache Hessian
                    Eigen::MatrixXd block_hessian_ii(vert_i_dim_unfixed, vert_j_dim_unfixed);
                    block_hessian_ii.noalias() = 2.0 * multiplier * block_jacobian1.transpose() * block_jacobian2;
                    for (int i = 0; i < vert_i_dim_unfixed; ++i)
                    {
                        for (int j = 0; j < i + 1; ++j)
                        {
                            values[nnz_idx] += block_hessian_ii.coeffRef(i, j);
                            ++nnz_idx;
                        }
                        // nnz_idx += vert_i_dim_unfixed * (vert_i_dim_unfixed + 1) / 2;
                    }
                }
                else
                {
                    Eigen::Map<Eigen::MatrixXd> block_hessian(values.data() + nnz_idx, vert_i_dim_unfixed, vert_j_dim_unfixed);
                    block_hessian += 2.0 * multiplier * block_jacobian1.transpose() * block_jacobian2;
                    nnz_idx += vert_i_dim_unfixed * vert_j_dim_unfixed;
                }
            }
        }
    }

    // Iterate mixed edges
    for (BaseMixedEdge::Ptr& edge : edges->getMixedEdgesRef())
    {
        if (edge->getObjectiveDimension() == 0 || edge->isObjectiveLinear()) continue;

        for (int row_i = 0; row_i < edge->getNumVertices(); ++row_i)
        {
            const VertexInterface* vertex1 = edge->getVertexRaw(row_i);

            int vert_i_dim_unfixed = vertex1->getDimensionUnfixed();
            if (vert_i_dim_unfixed == 0) continue;

            // compute jacobian to speed up hessian computation
            Eigen::MatrixXd block_jacobian1(edge->getObjectiveDimension(), vertex1->getDimensionUnfixed());
            edge->computeObjectiveJacobian(row_i, block_jacobian1);

            int col_end = lower_part_only ? row_i + 1 : edge->getNumVertices();
            for (int col_j = 0; col_j < col_end; ++col_j)  // hessian is symmetric
            {
                const VertexInterface* vertex2 = edge->getVertexRaw(col_j);

                int vert_j_dim_unfixed = vertex2->getDimensionUnfixed();
                if (vert_j_dim_unfixed == 0) continue;

                if (edge->isObjectiveLeastSquaresForm())
                {
                    // we need the second jacobian for applying the chain rule
                    Eigen::MatrixXd block_jacobian2(edge->getObjectiveDimension(), vertex2->getDimensionUnfixed());
                    edge->computeObjectiveJacobian(col_j, block_jacobian2);

                    // approximate H = 2*J^T*J
                    if (lower_part_only && vertex1 == vertex2)
                    {
                        // we only need the lower triangular part, so cache Hessian
                        Eigen::MatrixXd block_hessian_ii(vert_i_dim_unfixed, vert_j_dim_unfixed);
                        block_hessian_ii.noalias() = 2.0 * multiplier * block_jacobian1.transpose() * block_jacobian2;
                        for (int i = 0; i < vert_i_dim_unfixed; ++i)
                        {
                            for (int j = 0; j < i + 1; ++j)
                            {
                                values[nnz_idx] += block_hessian_ii.coeffRef(i, j);
                                ++nnz_idx;
                            }
                            // nnz_idx += vert_i_dim_unfixed * (vert_i_dim_unfixed + 1) / 2;
                        }
                    }
                    else
                    {
                        Eigen::Map<Eigen::MatrixXd> block_hessian(values.data() + nnz_idx, vert_i_dim_unfixed, vert_j_dim_unfixed);
                        block_hessian += 2.0 * multiplier * block_jacobian1.transpose() * block_jacobian2;
                        nnz_idx += vert_i_dim_unfixed * vert_j_dim_unfixed;
                    }
                }
                else
                {
                    // warning: we require the hessian to be initialized with zeros
                    // and that edge->computeHessian only adds (+=) values!!
                    if (lower_part_only && vertex1 == vertex2)
                    {
                        // we only need the lower triangular part, so cache Hessian
                        Eigen::MatrixXd block_hessian_ii(vert_i_dim_unfixed, vert_j_dim_unfixed);
                        edge->computeObjectiveHessian(row_i, col_j, block_jacobian1, block_hessian_ii, nullptr, multiplier);
                        for (int i = 0; i < vert_i_dim_unfixed; ++i)
                        {
                            for (int j = 0; j < i + 1; ++j)
                            {
                                values[nnz_idx] += block_hessian_ii.coeffRef(i, j);
                                ++nnz_idx;
                            }
                            // nnz_idx += vert_i_dim_unfixed * (vert_i_dim_unfixed + 1) / 2;
                        }
                    }
                    else
                    {
                        Eigen::Map<Eigen::MatrixXd> block_hessian(values.data() + nnz_idx, vert_i_dim_unfixed, vert_j_dim_unfixed);
                        edge->computeObjectiveHessianInc(row_i, col_j, block_jacobian1, block_hessian, nullptr, multiplier);
                        nnz_idx += vert_i_dim_unfixed * vert_j_dim_unfixed;
                    }
                }
            }
        }
    }
}

void HyperGraphOptimizationProblemEdgeBased::computeSparseHessianObjectiveLL(Eigen::SparseMatrix<double, Eigen::ColMajor, long long>& H,
                                                                             const Eigen::VectorXi* col_nnz, bool upper_part_only)
{
    H.setZero();
    // TODO(roesmann): what is most efficient, if we just want to udpate the Hessian?
    if (col_nnz) H.reserve(*col_nnz);

    OptimizationEdgeSet::Ptr edges = _graph.getEdgeSet();

    // Iterate plain objective edges
    for (BaseEdge::Ptr& edge : edges->getObjectiveEdgesRef())
    {
        if (edge->isLinear()) continue;

        for (int row_i = 0; row_i < edge->getNumVertices(); ++row_i)
        {
            const VertexInterface* vertex1 = edge->getVertexRaw(row_i);

            int vert_i_dim_unfixed = vertex1->getDimensionUnfixed();
            if (vert_i_dim_unfixed == 0) continue;

            // compute jacobian to speed up hessian computation
            Eigen::MatrixXd block_jacobian1(edge->getDimension(), vertex1->getDimensionUnfixed());
            edge->computeJacobian(row_i, block_jacobian1);

            int col_start = upper_part_only ? row_i : 0;
            for (int col_j = col_start; col_j < edge->getNumVertices(); ++col_j)  // hessian is symmetric
            {
                const VertexInterface* vertex2 = edge->getVertexRaw(col_j);

                int vert_j_dim_unfixed = vertex2->getDimensionUnfixed();
                if (vert_j_dim_unfixed == 0) continue;

                // warning: we require the hessian to be initialized with zeros
                // and that edge->computeHessian only adds (+=) values!!
                Eigen::MatrixXd block_hessian(vert_i_dim_unfixed, vert_j_dim_unfixed);
                edge->computeHessian(row_i, col_j, block_jacobian1, block_hessian);

                // insert values
                if (upper_part_only && vertex1 == vertex2)
                {
                    // block on the main diagonal, so use upper part only
                    for (int i = 0; i < block_hessian.rows(); ++i)
                    {
                        for (int j = i; j < block_hessian.cols(); ++j)
                        {
                            H.coeffRef(vertex1->getVertexIdx() + i, vertex2->getVertexIdx() + j) += block_hessian.coeffRef(i, j);
                        }
                    }
                }
                else
                {
                    for (int i = 0; i < block_hessian.rows(); ++i)
                    {
                        for (int j = 0; j < block_hessian.cols(); ++j)
                        {
                            H.coeffRef(vertex1->getVertexIdx() + i, vertex2->getVertexIdx() + j) += block_hessian.coeffRef(i, j);
                        }
                    }
                }
            }
        }
    }

    // Iterate lsq objective edges
    for (BaseEdge::Ptr& edge : edges->getLsqObjectiveEdgesRef())
    {
        // if (edge->isLinear()) continue; // we must not check this one, since Linear-squared is not linear ;-)

        for (int row_i = 0; row_i < edge->getNumVertices(); ++row_i)
        {
            const VertexInterface* vertex1 = edge->getVertexRaw(row_i);

            int vert_i_dim_unfixed = vertex1->getDimensionUnfixed();
            if (vert_i_dim_unfixed == 0) continue;

            // compute jacobian to speed up hessian computation
            Eigen::MatrixXd block_jacobian1(edge->getDimension(), vertex1->getDimensionUnfixed());
            edge->computeJacobian(row_i, block_jacobian1);

            int col_start = upper_part_only ? row_i : 0;
            for (int col_j = col_start; col_j < edge->getNumVertices(); ++col_j)  // hessian is symmetric
            {
                const VertexInterface* vertex2 = edge->getVertexRaw(col_j);

                int vert_j_dim_unfixed = vertex2->getDimensionUnfixed();
                if (vert_j_dim_unfixed == 0) continue;

                // we need the second jacobian for applying the chain rule
                Eigen::MatrixXd block_jacobian2(edge->getDimension(), vertex2->getDimensionUnfixed());
                edge->computeJacobian(col_j, block_jacobian2);

                // approximate H = 2*J^T*J
                Eigen::MatrixXd block_hessian(vert_i_dim_unfixed, vert_j_dim_unfixed);
                block_hessian = 2.0 * block_jacobian1.transpose() * block_jacobian2;

                // insert values
                if (upper_part_only && vertex1 == vertex2)
                {
                    // block on the main diagonal, so use upper part only
                    for (int i = 0; i < block_hessian.rows(); ++i)
                    {
                        for (int j = i; j < block_hessian.cols(); ++j)
                        {
                            H.coeffRef(vertex1->getVertexIdx() + i, vertex2->getVertexIdx() + j) += block_hessian.coeffRef(i, j);
                        }
                    }
                }
                else
                {
                    for (int i = 0; i < block_hessian.rows(); ++i)
                    {
                        for (int j = 0; j < block_hessian.cols(); ++j)
                        {
                            H.coeffRef(vertex1->getVertexIdx() + i, vertex2->getVertexIdx() + j) += block_hessian.coeffRef(i, j);
                        }
                    }
                }
            }
        }
    }

    // Iterate mixed edges
    for (BaseMixedEdge::Ptr& edge : edges->getMixedEdgesRef())
    {
        if (edge->getObjectiveDimension() == 0 || edge->isObjectiveLinear()) continue;

        for (int row_i = 0; row_i < edge->getNumVertices(); ++row_i)
        {
            const VertexInterface* vertex1 = edge->getVertexRaw(row_i);

            int vert_i_dim_unfixed = vertex1->getDimensionUnfixed();
            if (vert_i_dim_unfixed == 0) continue;

            // compute jacobian to speed up hessian computation
            Eigen::MatrixXd obj_block_jacobian1(edge->getObjectiveDimension(), vertex1->getDimensionUnfixed());
            Eigen::MatrixXd eq_block_jacobian1(edge->getEqualityDimension(), vertex1->getDimensionUnfixed());
            Eigen::MatrixXd ineq_block_jacobian1(edge->getInequalityDimension(), vertex1->getDimensionUnfixed());
            edge->computeJacobians(row_i, obj_block_jacobian1, eq_block_jacobian1, ineq_block_jacobian1);

            int col_start = upper_part_only ? row_i : 0;
            for (int col_j = col_start; col_j < edge->getNumVertices(); ++col_j)  // hessian is symmetric
            {
                const VertexInterface* vertex2 = edge->getVertexRaw(col_j);

                int vert_j_dim_unfixed = vertex2->getDimensionUnfixed();
                if (vert_j_dim_unfixed == 0) continue;

                Eigen::MatrixXd block_hessian(vert_i_dim_unfixed, vert_j_dim_unfixed);
                block_hessian.setZero();

                if (edge->isObjectiveLeastSquaresForm())
                {
                    if (!edge->isObjectiveLinear() && edge->getObjectiveDimension() != 0)
                    {
                        // we need the second jacobian for applying the chain rule
                        Eigen::MatrixXd obj_block_jacobian2(edge->getObjectiveDimension(), vertex2->getDimensionUnfixed());
                        edge->computeObjectiveJacobian(col_j, obj_block_jacobian2);

                        // approximate H = 2*J^T*J
                        block_hessian += 2.0 * obj_block_jacobian1.transpose() * obj_block_jacobian2;
                    }
                }
                else  // TODO(roesmann): check all reasonable cases:
                {
                    if (!edge->isObjectiveLinear() && edge->getObjectiveDimension() != 0)
                    {
                        edge->computeObjectiveHessianInc(row_i, col_j, obj_block_jacobian1, block_hessian, nullptr, 1.0);
                    }
                }
                // insert values
                if (upper_part_only && vertex1 == vertex2)
                {
                    // block on the main diagonal, so use upper part only
                    for (int i = 0; i < block_hessian.rows(); ++i)
                    {
                        for (int j = i; j < block_hessian.cols(); ++j)
                        {
                            H.coeffRef(vertex1->getVertexIdx() + i, vertex2->getVertexIdx() + j) += block_hessian.coeffRef(i, j);
                        }
                    }
                }
                else
                {
                    for (int i = 0; i < block_hessian.rows(); ++i)
                    {
                        for (int j = 0; j < block_hessian.cols(); ++j)
                        {
                            H.coeffRef(vertex1->getVertexIdx() + i, vertex2->getVertexIdx() + j) += block_hessian.coeffRef(i, j);
                        }
                    }
                }
            }
        }
    }
}

void HyperGraphOptimizationProblemEdgeBased::computeSparseHessianObjectiveNNZperCol(Eigen::Ref<Eigen::VectorXi> col_nnz, bool upper_part_only)
{
    assert(col_nnz.size() == getParameterDimension());

    col_nnz.setZero();
    // TODO(roesmann): what is most efficient, if we just want to udpate the Hessian?

    // TODO(roesmann): with the following strategy, we overestimate the number of nnz, since vertices are shared among multiple edges...
    // however, we do not want to create the complete sparsity pattern here by creating a matrix of zeros and ones (tagging).

    OptimizationEdgeSet::Ptr edges = _graph.getEdgeSet();

    // Iterate plain objective edges
    for (BaseEdge::Ptr& edge : edges->getObjectiveEdgesRef())
    {
        if (edge->isLinear()) continue;

        for (int row_i = 0; row_i < edge->getNumVertices(); ++row_i)
        {
            const VertexInterface* vertex1 = edge->getVertexRaw(row_i);

            int vert_i_dim_unfixed = vertex1->getDimensionUnfixed();
            if (vert_i_dim_unfixed == 0) continue;

            int col_start = upper_part_only ? row_i : 0;
            for (int col_j = col_start; col_j < edge->getNumVertices(); ++col_j)  // hessian is symmetric
            {
                const VertexInterface* vertex2 = edge->getVertexRaw(col_j);

                int vert_j_dim_unfixed = vertex2->getDimensionUnfixed();
                if (vert_j_dim_unfixed == 0) continue;

                if (upper_part_only && vertex1 == vertex2)
                {
                    // just the upper triangular part of the bock matrix including diagonal, so N*(N+1)/2 elements
                    for (int l = 0; l < vert_j_dim_unfixed; ++l)
                    {
                        col_nnz(vertex2->getVertexIdx() + l) += l + 1;  // first column 0, second 1, ...., vert_i_dim_unfixed
                    }
                }
                else
                {
                    col_nnz.segment(vertex2->getVertexIdx(), vert_j_dim_unfixed).array() += vert_i_dim_unfixed;
                }
            }
        }
    }

    // Iterate lsq objective edges
    for (BaseEdge::Ptr& edge : edges->getLsqObjectiveEdgesRef())
    {
        // if (edge->isLinear()) continue; // we must not check this one, since Linear-squared is not linear ;-)

        for (int row_i = 0; row_i < edge->getNumVertices(); ++row_i)
        {
            const VertexInterface* vertex1 = edge->getVertexRaw(row_i);

            int vert_i_dim_unfixed = vertex1->getDimensionUnfixed();
            if (vert_i_dim_unfixed == 0) continue;

            int col_start = upper_part_only ? row_i : 0;
            for (int col_j = col_start; col_j < edge->getNumVertices(); ++col_j)  // hessian is symmetric
            {
                const VertexInterface* vertex2 = edge->getVertexRaw(col_j);

                int vert_j_dim_unfixed = vertex2->getDimensionUnfixed();
                if (vert_j_dim_unfixed == 0) continue;

                if (upper_part_only && vertex1 == vertex2)
                {
                    // just the upper triangular part of the bock matrix including diagonal, so N*(N+1)/2 elements
                    for (int l = 0; l < vert_j_dim_unfixed; ++l)
                    {
                        col_nnz(vertex2->getVertexIdx() + l) += l + 1;  // first column 0, second 1, ...., vert_i_dim_unfixed
                    }
                }
                else
                {
                    col_nnz.segment(vertex2->getVertexIdx(), vert_j_dim_unfixed).array() += vert_i_dim_unfixed;
                }
            }
        }
    }

    // Iterate mixed edges
    for (BaseMixedEdge::Ptr& edge : edges->getMixedEdgesRef())
    {
        if (edge->getObjectiveDimension() == 0 || edge->isObjectiveLinear()) continue;

        for (int row_i = 0; row_i < edge->getNumVertices(); ++row_i)
        {
            const VertexInterface* vertex1 = edge->getVertexRaw(row_i);

            int vert_i_dim_unfixed = vertex1->getDimensionUnfixed();
            if (vert_i_dim_unfixed == 0) continue;

            int col_start = upper_part_only ? row_i : 0;
            for (int col_j = col_start; col_j < edge->getNumVertices(); ++col_j)  // hessian is symmetric
            {
                const VertexInterface* vertex2 = edge->getVertexRaw(col_j);

                int vert_j_dim_unfixed = vertex2->getDimensionUnfixed();
                if (vert_j_dim_unfixed == 0) continue;

                if (upper_part_only && vertex1 == vertex2)
                {
                    // just the upper triangular part of the bock matrix including diagonal, so N*(N+1)/2 elements
                    for (int l = 0; l < vert_j_dim_unfixed; ++l)
                    {
                        col_nnz(vertex2->getVertexIdx() + l) += l + 1;  // first column 0, second 1, ...., vert_i_dim_unfixed
                    }
                }
                else
                {
                    col_nnz.segment(vertex2->getVertexIdx(), vert_j_dim_unfixed).array() += vert_i_dim_unfixed;
                }
            }
        }
    }
}

int HyperGraphOptimizationProblemEdgeBased::computeSparseHessianEqualitiesNNZ(bool lower_part_only)
{
    int nnz = 0;

    assert(_graph.hasEdgeSet());
    if (!_graph_precomputed) precomputeGraphQuantities();

    if (getEqualityDimension() == 0) return nnz;

    OptimizationEdgeSet::Ptr edges = _graph.getEdgeSet();

    // Iterate equality edges
    for (BaseEdge::Ptr& edge : edges->getEqualityEdgesRef())
    {
        if (edge->isLinear()) continue;

        for (int row_i = 0; row_i < edge->getNumVertices(); ++row_i)
        {
            const VertexInterface* vertex1 = edge->getVertexRaw(row_i);

            int vert_i_dim_unfixed = vertex1->getDimensionUnfixed();
            if (vert_i_dim_unfixed == 0) continue;

            int col_end = lower_part_only ? row_i + 1 : edge->getNumVertices();
            for (int col_j = 0; col_j < col_end; ++col_j)  // hessian is symmetric
            {
                const VertexInterface* vertex2 = edge->getVertexRaw(col_j);

                int vert_j_dim_unfixed = vertex2->getDimensionUnfixed();
                if (vert_j_dim_unfixed == 0) continue;

                if (lower_part_only && vertex1 == vertex2)
                {
                    int n = vertex1->getDimensionUnfixed();
                    nnz += n * (n + 1) / 2;  // lower triangular including diagonal
                }
                else
                    nnz += vertex1->getDimensionUnfixed() * vertex2->getDimensionUnfixed();  // block size
            }
        }
    }
    // Iterate mixed edges
    for (BaseMixedEdge::Ptr& edge : edges->getMixedEdgesRef())
    {
        if (getEqualityDimension() == 0 || edge->isEqualityLinear()) continue;

        for (int row_i = 0; row_i < edge->getNumVertices(); ++row_i)
        {
            const VertexInterface* vertex1 = edge->getVertexRaw(row_i);

            int vert_i_dim_unfixed = vertex1->getDimensionUnfixed();
            if (vert_i_dim_unfixed == 0) continue;

            int col_end = lower_part_only ? row_i + 1 : edge->getNumVertices();
            for (int col_j = 0; col_j < col_end; ++col_j)  // hessian is symmetric
            {
                const VertexInterface* vertex2 = edge->getVertexRaw(col_j);

                int vert_j_dim_unfixed = vertex2->getDimensionUnfixed();
                if (vert_j_dim_unfixed == 0) continue;

                if (lower_part_only && vertex1 == vertex2)
                {
                    int n = vertex1->getDimensionUnfixed();
                    nnz += n * (n + 1) / 2;  // lower triangular including diagonal
                }
                else
                    nnz += vertex1->getDimensionUnfixed() * vertex2->getDimensionUnfixed();  // block size
            }
        }
    }

    return nnz;
}
void HyperGraphOptimizationProblemEdgeBased::computeSparseHessianEqualitiesStructure(Eigen::Ref<Eigen::VectorXi> i_row,
                                                                                     Eigen::Ref<Eigen::VectorXi> j_col, bool lower_part_only)
{
    assert(_graph.hasEdgeSet());
    if (!_graph_precomputed) precomputeGraphQuantities();

    assert(i_row.size() == computeSparseHessianEqualitiesNNZ(lower_part_only));
    assert(i_row.size() == j_col.size());

    OptimizationEdgeSet::Ptr edges = _graph.getEdgeSet();

    int nnz_idx = 0;

    // Iterate equality edges
    for (BaseEdge::Ptr& edge : edges->getEqualityEdgesRef())
    {
        if (edge->isLinear()) continue;

        for (int row_i = 0; row_i < edge->getNumVertices(); ++row_i)
        {
            const VertexInterface* vertex1 = edge->getVertexRaw(row_i);

            int vert_i_dim_unfixed = vertex1->getDimensionUnfixed();
            if (vert_i_dim_unfixed == 0) continue;

            int col_end = lower_part_only ? row_i + 1 : edge->getNumVertices();
            for (int col_j = 0; col_j < col_end; ++col_j)  // hessian is symmetric
            {
                const VertexInterface* vertex2 = edge->getVertexRaw(col_j);

                int vert_j_dim_unfixed = vertex2->getDimensionUnfixed();
                if (vert_j_dim_unfixed == 0) continue;

                if (lower_part_only && vertex1 == vertex2)
                {
                    for (int v1_idx = 0; v1_idx < vert_i_dim_unfixed; ++v1_idx)
                    {
                        for (int v2_idx = 0; v2_idx < v1_idx + 1; ++v2_idx)
                        {
                            i_row[nnz_idx] = vertex1->getVertexIdx() + v1_idx;
                            j_col[nnz_idx] = vertex2->getVertexIdx() + v2_idx;
                            ++nnz_idx;
                        }
                    }
                }
                else
                {
                    for (int v1_idx = 0; v1_idx < vert_i_dim_unfixed; ++v1_idx)
                    {
                        for (int v2_idx = 0; v2_idx < vert_j_dim_unfixed; ++v2_idx)
                        {
                            i_row[nnz_idx] = vertex1->getVertexIdx() + v1_idx;
                            j_col[nnz_idx] = vertex2->getVertexIdx() + v2_idx;
                            ++nnz_idx;
                        }
                    }
                }
            }
        }
    }

    // Iterate mixed edges
    for (BaseMixedEdge::Ptr& edge : edges->getMixedEdgesRef())
    {
        if (getEqualityDimension() == 0 || edge->isEqualityLinear()) continue;

        for (int row_i = 0; row_i < edge->getNumVertices(); ++row_i)
        {
            const VertexInterface* vertex1 = edge->getVertexRaw(row_i);

            int vert_i_dim_unfixed = vertex1->getDimensionUnfixed();
            if (vert_i_dim_unfixed == 0) continue;

            int col_end = lower_part_only ? row_i + 1 : edge->getNumVertices();
            for (int col_j = 0; col_j < col_end; ++col_j)  // hessian is symmetric
            {
                const VertexInterface* vertex2 = edge->getVertexRaw(col_j);

                int vert_j_dim_unfixed = vertex2->getDimensionUnfixed();
                if (vert_j_dim_unfixed == 0) continue;

                if (lower_part_only && vertex1 == vertex2)
                {
                    for (int v1_idx = 0; v1_idx < vert_i_dim_unfixed; ++v1_idx)
                    {
                        for (int v2_idx = 0; v2_idx < v1_idx + 1; ++v2_idx)
                        {
                            i_row[nnz_idx] = vertex1->getVertexIdx() + v1_idx;
                            j_col[nnz_idx] = vertex2->getVertexIdx() + v2_idx;
                            ++nnz_idx;
                        }
                    }
                }
                else
                {
                    for (int v1_idx = 0; v1_idx < vert_i_dim_unfixed; ++v1_idx)
                    {
                        for (int v2_idx = 0; v2_idx < vert_j_dim_unfixed; ++v2_idx)
                        {
                            i_row[nnz_idx] = vertex1->getVertexIdx() + v1_idx;
                            j_col[nnz_idx] = vertex2->getVertexIdx() + v2_idx;
                            ++nnz_idx;
                        }
                    }
                }
            }
        }
    }
}
void HyperGraphOptimizationProblemEdgeBased::computeSparseHessianEqualitiesValues(Eigen::Ref<Eigen::VectorXd> values, const double* multipliers,
                                                                                  bool lower_part_only)
{
    assert(values.size() == computeSparseHessianEqualitiesNNZ(lower_part_only));

    values.setZero();

    if (getEqualityDimension() == 0) return;

    int nnz_idx = 0;

    OptimizationEdgeSet::Ptr edges = _graph.getEdgeSet();

    // Iterate equality edges
    for (BaseEdge::Ptr& edge : edges->getEqualityEdgesRef())
    {
        if (edge->isLinear()) continue;

        const double* mult_eq_part = multipliers ? multipliers + edge->getEdgeIdx() : nullptr;

        for (int row_i = 0; row_i < edge->getNumVertices(); ++row_i)
        {
            const VertexInterface* vertex1 = edge->getVertexRaw(row_i);

            int vert_i_dim_unfixed = vertex1->getDimensionUnfixed();
            if (vert_i_dim_unfixed == 0) continue;

            // compute jacobian to speed up hessian computation
            Eigen::MatrixXd block_jacobian1(edge->getDimension(), vertex1->getDimensionUnfixed());
            edge->computeJacobian(row_i, block_jacobian1);

            int col_end = lower_part_only ? row_i + 1 : edge->getNumVertices();
            for (int col_j = 0; col_j < col_end; ++col_j)  // hessian is symmetric
            {
                const VertexInterface* vertex2 = edge->getVertexRaw(col_j);

                int vert_j_dim_unfixed = vertex2->getDimensionUnfixed();
                if (vert_j_dim_unfixed == 0) continue;

                // warning: we require the hessian to be initialized with zeros
                // and that edge->computeHessian only adds (+=) values!!
                if (lower_part_only && vertex1 == vertex2)
                {
                    // we only need the lower triangular part, so cache Hessian
                    Eigen::MatrixXd block_hessian_ii(vert_i_dim_unfixed, vert_j_dim_unfixed);
                    edge->computeHessian(row_i, col_j, block_jacobian1, block_hessian_ii, mult_eq_part);
                    for (int i = 0; i < vert_i_dim_unfixed; ++i)
                    {
                        for (int j = 0; j < i + 1; ++j)
                        {
                            values[nnz_idx] += block_hessian_ii.coeffRef(i, j);
                            ++nnz_idx;
                        }
                        // nnz_idx += vert_i_dim_unfixed * (vert_i_dim_unfixed + 1) / 2;
                    }
                }
                else
                {
                    Eigen::Map<Eigen::MatrixXd> block_hessian(values.data() + nnz_idx, vert_i_dim_unfixed, vert_j_dim_unfixed);
                    edge->computeHessianInc(row_i, col_j, block_jacobian1, block_hessian, mult_eq_part);
                    nnz_idx += vert_i_dim_unfixed * vert_j_dim_unfixed;
                }
            }
        }
    }

    // Iterate mixed edges
    for (BaseMixedEdge::Ptr& edge : edges->getMixedEdgesRef())
    {
        if (edge->getEqualityDimension() == 0 || edge->isEqualityLinear()) continue;

        const double* mult_eq_part = multipliers ? multipliers + edge->getEdgeEqualityIdx() : nullptr;

        for (int row_i = 0; row_i < edge->getNumVertices(); ++row_i)
        {
            const VertexInterface* vertex1 = edge->getVertexRaw(row_i);

            int vert_i_dim_unfixed = vertex1->getDimensionUnfixed();
            if (vert_i_dim_unfixed == 0) continue;

            // compute jacobian to speed up hessian computation
            Eigen::MatrixXd block_jacobian1(edge->getEqualityDimension(), vertex1->getDimensionUnfixed());
            edge->computeEqualityJacobian(row_i, block_jacobian1);

            int col_end = lower_part_only ? row_i + 1 : edge->getNumVertices();
            for (int col_j = 0; col_j < col_end; ++col_j)  // hessian is symmetric
            {
                const VertexInterface* vertex2 = edge->getVertexRaw(col_j);

                int vert_j_dim_unfixed = vertex2->getDimensionUnfixed();
                if (vert_j_dim_unfixed == 0) continue;

                // warning: we require the hessian to be initialized with zeros
                // and that edge->computeHessian only adds (+=) values!!
                if (lower_part_only && vertex1 == vertex2)
                {
                    // we only need the lower triangular part, so cache Hessian
                    Eigen::MatrixXd block_hessian_ii(vert_i_dim_unfixed, vert_j_dim_unfixed);
                    edge->computeEqualityHessian(row_i, col_j, block_jacobian1, block_hessian_ii, mult_eq_part);
                    for (int i = 0; i < vert_i_dim_unfixed; ++i)
                    {
                        for (int j = 0; j < i + 1; ++j)
                        {
                            values[nnz_idx] += block_hessian_ii.coeffRef(i, j);
                            ++nnz_idx;
                        }
                        // nnz_idx += vert_i_dim_unfixed * (vert_i_dim_unfixed + 1) / 2;
                    }
                }
                else
                {
                    Eigen::Map<Eigen::MatrixXd> block_hessian(values.data() + nnz_idx, vert_i_dim_unfixed, vert_j_dim_unfixed);
                    edge->computeEqualityHessianInc(row_i, col_j, block_jacobian1, block_hessian, mult_eq_part);
                    nnz_idx += vert_i_dim_unfixed * vert_j_dim_unfixed;
                }
            }
        }
    }
}

int HyperGraphOptimizationProblemEdgeBased::computeSparseHessianInequalitiesNNZ(bool lower_part_only)
{
    int nnz = 0;

    assert(_graph.hasEdgeSet());
    if (!_graph_precomputed) precomputeGraphQuantities();

    if (getInequalityDimension() == 0) return nnz;

    OptimizationEdgeSet::Ptr edges = _graph.getEdgeSet();

    // Iterate inequality edges
    for (BaseEdge::Ptr& edge : edges->getInequalityEdgesRef())
    {
        if (edge->isLinear()) continue;

        for (int row_i = 0; row_i < edge->getNumVertices(); ++row_i)
        {
            const VertexInterface* vertex1 = edge->getVertexRaw(row_i);

            int vert_i_dim_unfixed = vertex1->getDimensionUnfixed();
            if (vert_i_dim_unfixed == 0) continue;

            int col_end = lower_part_only ? row_i + 1 : edge->getNumVertices();
            for (int col_j = 0; col_j < col_end; ++col_j)  // hessian is symmetric
            {
                const VertexInterface* vertex2 = edge->getVertexRaw(col_j);

                int vert_j_dim_unfixed = vertex2->getDimensionUnfixed();
                if (vert_j_dim_unfixed == 0) continue;

                if (lower_part_only && vertex1 == vertex2)
                {
                    int n = vertex1->getDimensionUnfixed();
                    nnz += n * (n + 1) / 2;  // lower triangular including diagonal
                }
                else
                    nnz += vertex1->getDimensionUnfixed() * vertex2->getDimensionUnfixed();  // block size
            }
        }
    }
    // Iterate mixed edges
    for (BaseMixedEdge::Ptr& edge : edges->getMixedEdgesRef())
    {
        if (getInequalityDimension() == 0 || edge->isInequalityLinear()) continue;

        for (int row_i = 0; row_i < edge->getNumVertices(); ++row_i)
        {
            const VertexInterface* vertex1 = edge->getVertexRaw(row_i);

            int vert_i_dim_unfixed = vertex1->getDimensionUnfixed();
            if (vert_i_dim_unfixed == 0) continue;

            int col_end = lower_part_only ? row_i + 1 : edge->getNumVertices();
            for (int col_j = 0; col_j < col_end; ++col_j)  // hessian is symmetric
            {
                const VertexInterface* vertex2 = edge->getVertexRaw(col_j);

                int vert_j_dim_unfixed = vertex2->getDimensionUnfixed();
                if (vert_j_dim_unfixed == 0) continue;

                if (lower_part_only && vertex1 == vertex2)
                {
                    int n = vertex1->getDimensionUnfixed();
                    nnz += n * (n + 1) / 2;  // lower triangular including diagonal
                }
                else
                    nnz += vertex1->getDimensionUnfixed() * vertex2->getDimensionUnfixed();  // block size
            }
        }
    }

    return nnz;
}
void HyperGraphOptimizationProblemEdgeBased::computeSparseHessianInequalitiesStructure(Eigen::Ref<Eigen::VectorXi> i_row,
                                                                                       Eigen::Ref<Eigen::VectorXi> j_col, bool lower_part_only)
{
    assert(_graph.hasEdgeSet());
    if (!_graph_precomputed) precomputeGraphQuantities();

    assert(i_row.size() == computeSparseHessianInequalitiesNNZ(lower_part_only));
    assert(i_row.size() == j_col.size());

    OptimizationEdgeSet::Ptr edges = _graph.getEdgeSet();

    int nnz_idx = 0;

    // Iterate inequality edges
    for (BaseEdge::Ptr& edge : edges->getInequalityEdgesRef())
    {
        if (edge->isLinear()) continue;

        for (int row_i = 0; row_i < edge->getNumVertices(); ++row_i)
        {
            const VertexInterface* vertex1 = edge->getVertexRaw(row_i);

            int vert_i_dim_unfixed = vertex1->getDimensionUnfixed();
            if (vert_i_dim_unfixed == 0) continue;

            int col_end = lower_part_only ? row_i + 1 : edge->getNumVertices();
            for (int col_j = 0; col_j < col_end; ++col_j)  // hessian is symmetric
            {
                const VertexInterface* vertex2 = edge->getVertexRaw(col_j);

                int vert_j_dim_unfixed = vertex2->getDimensionUnfixed();
                if (vert_j_dim_unfixed == 0) continue;

                if (lower_part_only && vertex1 == vertex2)
                {
                    for (int v1_idx = 0; v1_idx < vert_i_dim_unfixed; ++v1_idx)
                    {
                        for (int v2_idx = 0; v2_idx < v1_idx + 1; ++v2_idx)
                        {
                            i_row[nnz_idx] = vertex1->getVertexIdx() + v1_idx;
                            j_col[nnz_idx] = vertex2->getVertexIdx() + v2_idx;
                            ++nnz_idx;
                        }
                    }
                }
                else
                {
                    for (int v1_idx = 0; v1_idx < vert_i_dim_unfixed; ++v1_idx)
                    {
                        for (int v2_idx = 0; v2_idx < vert_j_dim_unfixed; ++v2_idx)
                        {
                            i_row[nnz_idx] = vertex1->getVertexIdx() + v1_idx;
                            j_col[nnz_idx] = vertex2->getVertexIdx() + v2_idx;
                            ++nnz_idx;
                        }
                    }
                }
            }
        }
    }

    // Iterate mixed edges
    for (BaseMixedEdge::Ptr& edge : edges->getMixedEdgesRef())
    {
        if (getInequalityDimension() == 0 || edge->isInequalityLinear()) continue;

        for (int row_i = 0; row_i < edge->getNumVertices(); ++row_i)
        {
            const VertexInterface* vertex1 = edge->getVertexRaw(row_i);

            int vert_i_dim_unfixed = vertex1->getDimensionUnfixed();
            if (vert_i_dim_unfixed == 0) continue;

            int col_end = lower_part_only ? row_i + 1 : edge->getNumVertices();
            for (int col_j = 0; col_j < col_end; ++col_j)  // hessian is symmetric
            {
                const VertexInterface* vertex2 = edge->getVertexRaw(col_j);

                int vert_j_dim_unfixed = vertex2->getDimensionUnfixed();
                if (vert_j_dim_unfixed == 0) continue;

                if (lower_part_only && vertex1 == vertex2)
                {
                    for (int v1_idx = 0; v1_idx < vert_i_dim_unfixed; ++v1_idx)
                    {
                        for (int v2_idx = 0; v2_idx < v1_idx + 1; ++v2_idx)
                        {
                            i_row[nnz_idx] = vertex1->getVertexIdx() + v1_idx;
                            j_col[nnz_idx] = vertex2->getVertexIdx() + v2_idx;
                            ++nnz_idx;
                        }
                    }
                }
                else
                {
                    for (int v1_idx = 0; v1_idx < vert_i_dim_unfixed; ++v1_idx)
                    {
                        for (int v2_idx = 0; v2_idx < vert_j_dim_unfixed; ++v2_idx)
                        {
                            i_row[nnz_idx] = vertex1->getVertexIdx() + v1_idx;
                            j_col[nnz_idx] = vertex2->getVertexIdx() + v2_idx;
                            ++nnz_idx;
                        }
                    }
                }
            }
        }
    }
}
void HyperGraphOptimizationProblemEdgeBased::computeSparseHessianInequalitiesValues(Eigen::Ref<Eigen::VectorXd> values, const double* multipliers,
                                                                                    bool lower_part_only)
{
    assert(values.size() == computeSparseHessianInequalitiesNNZ(lower_part_only));

    values.setZero();

    if (getInequalityDimension() == 0) return;

    int nnz_idx = 0;

    OptimizationEdgeSet::Ptr edges = _graph.getEdgeSet();

    // Iterate inequality edges
    for (BaseEdge::Ptr& edge : edges->getInequalityEdgesRef())
    {
        if (edge->isLinear()) continue;

        const double* mult_ineq_part = multipliers ? multipliers + edge->getEdgeIdx() : nullptr;

        for (int row_i = 0; row_i < edge->getNumVertices(); ++row_i)
        {
            const VertexInterface* vertex1 = edge->getVertexRaw(row_i);

            int vert_i_dim_unfixed = vertex1->getDimensionUnfixed();
            if (vert_i_dim_unfixed == 0) continue;

            // compute jacobian to speed up hessian computation
            Eigen::MatrixXd block_jacobian1(edge->getDimension(), vertex1->getDimensionUnfixed());
            edge->computeJacobian(row_i, block_jacobian1);

            int col_end = lower_part_only ? row_i + 1 : edge->getNumVertices();
            for (int col_j = 0; col_j < col_end; ++col_j)  // hessian is symmetric
            {
                const VertexInterface* vertex2 = edge->getVertexRaw(col_j);

                int vert_j_dim_unfixed = vertex2->getDimensionUnfixed();
                if (vert_j_dim_unfixed == 0) continue;

                // warning: we require the hessian to be initialized with zeros
                // and that edge->computeHessian only adds (+=) values!!
                if (lower_part_only && vertex1 == vertex2)
                {
                    // we only need the lower triangular part, so cache Hessian
                    Eigen::MatrixXd block_hessian_ii(vert_i_dim_unfixed, vert_j_dim_unfixed);
                    edge->computeHessian(row_i, col_j, block_jacobian1, block_hessian_ii, mult_ineq_part);
                    for (int i = 0; i < vert_i_dim_unfixed; ++i)
                    {
                        for (int j = 0; j < i + 1; ++j)
                        {
                            values[nnz_idx] += block_hessian_ii.coeffRef(i, j);
                            ++nnz_idx;
                        }
                        // nnz_idx += vert_i_dim_unfixed * (vert_i_dim_unfixed + 1) / 2;
                    }
                }
                else
                {
                    Eigen::Map<Eigen::MatrixXd> block_hessian(values.data() + nnz_idx, vert_i_dim_unfixed, vert_j_dim_unfixed);
                    edge->computeHessianInc(row_i, col_j, block_jacobian1, block_hessian, mult_ineq_part);
                    nnz_idx += vert_i_dim_unfixed * vert_j_dim_unfixed;
                }
            }
        }
    }

    // Iterate mixed edges
    for (BaseMixedEdge::Ptr& edge : edges->getMixedEdgesRef())
    {
        if (edge->getInequalityDimension() == 0 || edge->isInequalityLinear()) continue;

        const double* mult_ineq_part = multipliers ? multipliers + edge->getEdgeInequalityIdx() : nullptr;

        for (int row_i = 0; row_i < edge->getNumVertices(); ++row_i)
        {
            const VertexInterface* vertex1 = edge->getVertexRaw(row_i);

            int vert_i_dim_unfixed = vertex1->getDimensionUnfixed();
            if (vert_i_dim_unfixed == 0) continue;

            // compute jacobian to speed up hessian computation
            Eigen::MatrixXd block_jacobian1(edge->getInequalityDimension(), vertex1->getDimensionUnfixed());
            edge->computeInequalityJacobian(row_i, block_jacobian1);

            int col_end = lower_part_only ? row_i + 1 : edge->getNumVertices();
            for (int col_j = 0; col_j < col_end; ++col_j)  // hessian is symmetric
            {
                const VertexInterface* vertex2 = edge->getVertexRaw(col_j);

                int vert_j_dim_unfixed = vertex2->getDimensionUnfixed();
                if (vert_j_dim_unfixed == 0) continue;

                // warning: we require the hessian to be initialized with zeros
                // and that edge->computeHessian only adds (+=) values!!
                if (lower_part_only && vertex1 == vertex2)
                {
                    // we only need the lower triangular part, so cache Hessian
                    Eigen::MatrixXd block_hessian_ii(vert_i_dim_unfixed, vert_j_dim_unfixed);
                    edge->computeInequalityHessian(row_i, col_j, block_jacobian1, block_hessian_ii, mult_ineq_part);
                    for (int i = 0; i < vert_i_dim_unfixed; ++i)
                    {
                        for (int j = 0; j < i + 1; ++j)
                        {
                            values[nnz_idx] += block_hessian_ii.coeffRef(i, j);
                            ++nnz_idx;
                        }
                        // nnz_idx += vert_i_dim_unfixed * (vert_i_dim_unfixed + 1) / 2;
                    }
                }
                else
                {
                    Eigen::Map<Eigen::MatrixXd> block_hessian(values.data() + nnz_idx, vert_i_dim_unfixed, vert_j_dim_unfixed);
                    edge->computeInequalityHessianInc(row_i, col_j, block_jacobian1, block_hessian, mult_ineq_part);
                    nnz_idx += vert_i_dim_unfixed * vert_j_dim_unfixed;
                }
            }
        }
    }
}

void HyperGraphOptimizationProblemEdgeBased::computeSparseHessiansNNZ(int& nnz_obj, int& nnz_eq, int& nnz_ineq, bool lower_part_only)
{
    // we should not forget to check for isObjectiveLinear/isEqualityLinear/isInequalityLinear to match number of NNZ
    nnz_obj  = computeSparseHessianObjectiveNNZ(lower_part_only);
    nnz_eq   = computeSparseHessianEqualitiesNNZ(lower_part_only);
    nnz_ineq = computeSparseHessianInequalitiesNNZ(lower_part_only);
}

void HyperGraphOptimizationProblemEdgeBased::computeSparseHessiansStructure(
    Eigen::Ref<Eigen::VectorXi> i_row_obj, Eigen::Ref<Eigen::VectorXi> j_col_obj, Eigen::Ref<Eigen::VectorXi> i_row_eq,
    Eigen::Ref<Eigen::VectorXi> j_col_eq, Eigen::Ref<Eigen::VectorXi> i_row_ineq, Eigen::Ref<Eigen::VectorXi> j_col_ineq, bool lower_part_only)
{
    computeSparseHessianObjectiveStructure(i_row_obj, j_col_obj, lower_part_only);
    computeSparseHessianEqualitiesStructure(i_row_eq, j_col_eq, lower_part_only);
    computeSparseHessianInequalitiesStructure(i_row_ineq, j_col_ineq, lower_part_only);
}
void HyperGraphOptimizationProblemEdgeBased::computeSparseHessiansValues(Eigen::Ref<Eigen::VectorXd> values_obj,
                                                                         Eigen::Ref<Eigen::VectorXd> values_eq,
                                                                         Eigen::Ref<Eigen::VectorXd> values_ineq, double multiplier_obj,
                                                                         const double* multipliers_eq, const double* multipliers_ineq,
                                                                         bool lower_part_only)
{
    assert(values_obj.size() == computeSparseHessianObjectiveNNZ(lower_part_only));
    assert(values_eq.size() == computeSparseHessianEqualitiesNNZ(lower_part_only));
    assert(values_ineq.size() == computeSparseHessianInequalitiesNNZ(lower_part_only));

    values_obj.setZero();
    values_eq.setZero();
    values_ineq.setZero();

    int obj_nnz_idx  = 0;
    int eq_nnz_idx   = 0;
    int ineq_nnz_idx = 0;

    OptimizationEdgeSet::Ptr edges = _graph.getEdgeSet();

    // Iterate plain objective edges
    for (BaseEdge::Ptr& edge : edges->getObjectiveEdgesRef())
    {
        if (edge->isLinear()) continue;

        for (int row_i = 0; row_i < edge->getNumVertices(); ++row_i)
        {
            const VertexInterface* vertex1 = edge->getVertexRaw(row_i);

            int vert_i_dim_unfixed = vertex1->getDimensionUnfixed();
            if (vert_i_dim_unfixed == 0) continue;

            // compute jacobian to speed up hessian computation
            Eigen::MatrixXd block_jacobian1(edge->getDimension(), vertex1->getDimensionUnfixed());
            edge->computeJacobian(row_i, block_jacobian1);

            int col_end = lower_part_only ? row_i + 1 : edge->getNumVertices();
            for (int col_j = 0; col_j < col_end; ++col_j)  // hessian is symmetric
            {
                const VertexInterface* vertex2 = edge->getVertexRaw(col_j);

                int vert_j_dim_unfixed = vertex2->getDimensionUnfixed();
                if (vert_j_dim_unfixed == 0) continue;

                // warning: we require the hessian to be initialized with zeros
                // and that edge->computeHessian only adds (+=) values!!
                if (lower_part_only && vertex1 == vertex2)
                {
                    // we only need the lower triangular part, so cache Hessian
                    Eigen::MatrixXd block_hessian_ii(vert_i_dim_unfixed, vert_j_dim_unfixed);
                    edge->computeHessian(row_i, col_j, block_jacobian1, block_hessian_ii, nullptr, multiplier_obj);
                    for (int i = 0; i < vert_i_dim_unfixed; ++i)
                    {
                        for (int j = 0; j < i + 1; ++j)
                        {
                            values_obj[obj_nnz_idx] += block_hessian_ii.coeffRef(i, j);
                            ++obj_nnz_idx;
                        }
                        // nnz_idx += vert_i_dim_unfixed * (vert_i_dim_unfixed + 1) / 2;
                    }
                }
                else
                {
                    Eigen::Map<Eigen::MatrixXd> block_hessian(values_obj.data() + obj_nnz_idx, vert_i_dim_unfixed, vert_j_dim_unfixed);
                    edge->computeHessianInc(row_i, col_j, block_jacobian1, block_hessian, nullptr, multiplier_obj);
                    obj_nnz_idx += vert_i_dim_unfixed * vert_j_dim_unfixed;
                }
            }
        }
    }

    // Iterate lsq objective edges
    for (BaseEdge::Ptr& edge : edges->getLsqObjectiveEdgesRef())
    {
        // if (edge->isLinear()) continue; // we must not check this one, since Linear-squared is not linear ;-)

        for (int row_i = 0; row_i < edge->getNumVertices(); ++row_i)
        {
            const VertexInterface* vertex1 = edge->getVertexRaw(row_i);

            int vert_i_dim_unfixed = vertex1->getDimensionUnfixed();
            if (vert_i_dim_unfixed == 0) continue;

            // compute jacobian to speed up hessian computation
            Eigen::MatrixXd block_jacobian1(edge->getDimension(), vertex1->getDimensionUnfixed());
            edge->computeJacobian(row_i, block_jacobian1);

            int col_end = lower_part_only ? row_i + 1 : edge->getNumVertices();
            for (int col_j = 0; col_j < col_end; ++col_j)  // hessian is symmetric
            {
                const VertexInterface* vertex2 = edge->getVertexRaw(col_j);

                int vert_j_dim_unfixed = vertex2->getDimensionUnfixed();
                if (vert_j_dim_unfixed == 0) continue;

                // we need the second jacobian for applying the chain rule
                Eigen::MatrixXd block_jacobian2(edge->getDimension(), vertex2->getDimensionUnfixed());
                edge->computeJacobian(col_j, block_jacobian2);

                // approximate H = 2*J^T*J
                if (lower_part_only && vertex1 == vertex2)
                {
                    // we only need the lower triangular part, so cache Hessian
                    Eigen::MatrixXd block_hessian_ii(vert_i_dim_unfixed, vert_j_dim_unfixed);
                    block_hessian_ii.noalias() = 2.0 * multiplier_obj * block_jacobian1.transpose() * block_jacobian2;
                    for (int i = 0; i < vert_i_dim_unfixed; ++i)
                    {
                        for (int j = 0; j < i + 1; ++j)
                        {
                            values_obj[obj_nnz_idx] += block_hessian_ii.coeffRef(i, j);
                            ++obj_nnz_idx;
                        }
                        // nnz_idx += vert_i_dim_unfixed * (vert_i_dim_unfixed + 1) / 2;
                    }
                }
                else
                {
                    Eigen::Map<Eigen::MatrixXd> block_hessian(values_obj.data() + obj_nnz_idx, vert_i_dim_unfixed, vert_j_dim_unfixed);
                    block_hessian += 2.0 * multiplier_obj * block_jacobian1.transpose() * block_jacobian2;
                    obj_nnz_idx += vert_i_dim_unfixed * vert_j_dim_unfixed;
                }
            }
        }
    }

    // Iterate equality edges
    for (BaseEdge::Ptr& edge : edges->getEqualityEdgesRef())
    {
        if (edge->isLinear()) continue;

        const double* mult_eq_part = multipliers_eq ? multipliers_eq + edge->getEdgeIdx() : nullptr;

        for (int row_i = 0; row_i < edge->getNumVertices(); ++row_i)
        {
            const VertexInterface* vertex1 = edge->getVertexRaw(row_i);

            int vert_i_dim_unfixed = vertex1->getDimensionUnfixed();
            if (vert_i_dim_unfixed == 0) continue;

            // compute jacobian to speed up hessian computation
            Eigen::MatrixXd block_jacobian1(edge->getDimension(), vertex1->getDimensionUnfixed());
            edge->computeJacobian(row_i, block_jacobian1);

            int col_end = lower_part_only ? row_i + 1 : edge->getNumVertices();
            for (int col_j = 0; col_j < col_end; ++col_j)  // hessian is symmetric
            {
                const VertexInterface* vertex2 = edge->getVertexRaw(col_j);

                int vert_j_dim_unfixed = vertex2->getDimensionUnfixed();
                if (vert_j_dim_unfixed == 0) continue;

                // warning: we require the hessian to be initialized with zeros
                // and that edge->computeHessian only adds (+=) values!!
                if (lower_part_only && vertex1 == vertex2)
                {
                    // we only need the lower triangular part, so cache Hessian
                    Eigen::MatrixXd block_hessian_ii(vert_i_dim_unfixed, vert_j_dim_unfixed);
                    edge->computeHessian(row_i, col_j, block_jacobian1, block_hessian_ii, mult_eq_part);
                    for (int i = 0; i < vert_i_dim_unfixed; ++i)
                    {
                        for (int j = 0; j < i + 1; ++j)
                        {
                            values_eq[eq_nnz_idx] += block_hessian_ii.coeffRef(i, j);
                            ++eq_nnz_idx;
                        }
                        // nnz_idx += vert_i_dim_unfixed * (vert_i_dim_unfixed + 1) / 2;
                    }
                }
                else
                {
                    Eigen::Map<Eigen::MatrixXd> block_hessian(values_eq.data() + eq_nnz_idx, vert_i_dim_unfixed, vert_j_dim_unfixed);
                    edge->computeHessianInc(row_i, col_j, block_jacobian1, block_hessian, mult_eq_part);
                    eq_nnz_idx += vert_i_dim_unfixed * vert_j_dim_unfixed;
                }
            }
        }
    }

    // Iterate inequality edges
    for (BaseEdge::Ptr& edge : edges->getInequalityEdgesRef())
    {
        if (edge->isLinear()) continue;

        const double* mult_ineq_part = multipliers_ineq ? multipliers_ineq + edge->getEdgeIdx() : nullptr;

        for (int row_i = 0; row_i < edge->getNumVertices(); ++row_i)
        {
            const VertexInterface* vertex1 = edge->getVertexRaw(row_i);

            int vert_i_dim_unfixed = vertex1->getDimensionUnfixed();
            if (vert_i_dim_unfixed == 0) continue;

            // compute jacobian to speed up hessian computation
            Eigen::MatrixXd block_jacobian1(edge->getDimension(), vertex1->getDimensionUnfixed());
            edge->computeJacobian(row_i, block_jacobian1);

            int col_end = lower_part_only ? row_i + 1 : edge->getNumVertices();
            for (int col_j = 0; col_j < col_end; ++col_j)  // hessian is symmetric
            {
                const VertexInterface* vertex2 = edge->getVertexRaw(col_j);

                int vert_j_dim_unfixed = vertex2->getDimensionUnfixed();
                if (vert_j_dim_unfixed == 0) continue;

                // warning: we require the hessian to be initialized with zeros
                // and that edge->computeHessian only adds (+=) values!!
                if (lower_part_only && vertex1 == vertex2)
                {
                    // we only need the lower triangular part, so cache Hessian
                    Eigen::MatrixXd block_hessian_ii(vert_i_dim_unfixed, vert_j_dim_unfixed);
                    edge->computeHessian(row_i, col_j, block_jacobian1, block_hessian_ii, mult_ineq_part);
                    for (int i = 0; i < vert_i_dim_unfixed; ++i)
                    {
                        for (int j = 0; j < i + 1; ++j)
                        {
                            values_ineq[ineq_nnz_idx] += block_hessian_ii.coeffRef(i, j);
                            ++ineq_nnz_idx;
                        }
                        // nnz_idx += vert_i_dim_unfixed * (vert_i_dim_unfixed + 1) / 2;
                    }
                }
                else
                {
                    Eigen::Map<Eigen::MatrixXd> block_hessian(values_ineq.data() + ineq_nnz_idx, vert_i_dim_unfixed, vert_j_dim_unfixed);
                    edge->computeHessianInc(row_i, col_j, block_jacobian1, block_hessian, mult_ineq_part);
                    ineq_nnz_idx += vert_i_dim_unfixed * vert_j_dim_unfixed;
                }
            }
        }
    }

    // Iterate mixed edges
    for (BaseMixedEdge::Ptr& edge : edges->getMixedEdgesRef())
    {
        // if (edge->getObjectiveDimension() == 0 || edge->isObjectiveLinear()) continue;

        const double* mult_eq_part   = multipliers_eq ? multipliers_eq + edge->getEdgeEqualityIdx() : nullptr;
        const double* mult_ineq_part = multipliers_ineq ? multipliers_ineq + edge->getEdgeInequalityIdx() : nullptr;

        for (int row_i = 0; row_i < edge->getNumVertices(); ++row_i)
        {
            const VertexInterface* vertex1 = edge->getVertexRaw(row_i);

            int vert_i_dim_unfixed = vertex1->getDimensionUnfixed();
            if (vert_i_dim_unfixed == 0) continue;

            // compute jacobian to speed up hessian computation
            Eigen::MatrixXd obj_block_jacobian1(edge->getObjectiveDimension(), vertex1->getDimensionUnfixed());
            Eigen::MatrixXd eq_block_jacobian1(edge->getEqualityDimension(), vertex1->getDimensionUnfixed());
            Eigen::MatrixXd ineq_block_jacobian1(edge->getInequalityDimension(), vertex1->getDimensionUnfixed());
            edge->computeJacobians(row_i, obj_block_jacobian1, eq_block_jacobian1, ineq_block_jacobian1);

            int col_end = lower_part_only ? row_i + 1 : edge->getNumVertices();
            for (int col_j = 0; col_j < col_end; ++col_j)  // hessian is symmetric
            {
                const VertexInterface* vertex2 = edge->getVertexRaw(col_j);

                int vert_j_dim_unfixed = vertex2->getDimensionUnfixed();
                if (vert_j_dim_unfixed == 0) continue;

                if (edge->isObjectiveLeastSquaresForm())
                {
                    if (!edge->isObjectiveLinear() && edge->getObjectiveDimension() != 0)
                    {
                        // we need the second jacobian for applying the chain rule
                        Eigen::MatrixXd obj_block_jacobian2(edge->getObjectiveDimension(), vertex2->getDimensionUnfixed());
                        edge->computeObjectiveJacobian(col_j, obj_block_jacobian2);

                        // approximate H = 2*J^T*J
                        if (lower_part_only && vertex1 == vertex2)
                        {
                            // we only need the lower triangular part, so cache Hessian
                            Eigen::MatrixXd block_hessian_ii(vert_i_dim_unfixed, vert_j_dim_unfixed);
                            block_hessian_ii.noalias() = 2.0 * multiplier_obj * obj_block_jacobian1.transpose() * obj_block_jacobian2;
                            for (int i = 0; i < vert_i_dim_unfixed; ++i)
                            {
                                for (int j = 0; j < i + 1; ++j)
                                {
                                    values_obj[obj_nnz_idx] += block_hessian_ii.coeffRef(i, j);
                                    ++obj_nnz_idx;
                                }
                                // nnz_idx += vert_i_dim_unfixed * (vert_i_dim_unfixed + 1) / 2;
                            }
                        }
                        else
                        {
                            Eigen::Map<Eigen::MatrixXd> obj_block_hessian(values_obj.data() + obj_nnz_idx, vert_i_dim_unfixed, vert_j_dim_unfixed);
                            obj_block_hessian += 2.0 * multiplier_obj * obj_block_jacobian1.transpose() * obj_block_jacobian2;
                            obj_nnz_idx += vert_i_dim_unfixed * vert_j_dim_unfixed;
                        }
                    }

                    if (!edge->isEqualityLinear() && !edge->isInequalityLinear())
                    {
                        if (lower_part_only && vertex1 == vertex2)
                        {
                            // we only need the lower triangular part, so cache Hessian
                            Eigen::MatrixXd eq_block_hessian_ii(vert_i_dim_unfixed, vert_j_dim_unfixed);
                            Eigen::MatrixXd ineq_block_hessian_ii(vert_i_dim_unfixed, vert_j_dim_unfixed);
                            edge->computeConstraintHessians(row_i, col_j, eq_block_jacobian1, ineq_block_jacobian1, eq_block_hessian_ii,
                                                            ineq_block_hessian_ii, mult_eq_part, mult_ineq_part);
                            for (int i = 0; i < vert_i_dim_unfixed; ++i)
                            {
                                for (int j = 0; j < i + 1; ++j)
                                {
                                    values_eq[eq_nnz_idx] += eq_block_hessian_ii.coeffRef(i, j);
                                    ++eq_nnz_idx;
                                    values_ineq[ineq_nnz_idx] += ineq_block_hessian_ii.coeffRef(i, j);
                                    ++ineq_nnz_idx;
                                }
                                // nnz_idx += vert_i_dim_unfixed * (vert_i_dim_unfixed + 1) / 2;
                            }
                        }
                        else
                        {
                            Eigen::Map<Eigen::MatrixXd> eq_block_hessian(values_eq.data() + eq_nnz_idx, vert_i_dim_unfixed, vert_j_dim_unfixed);
                            Eigen::Map<Eigen::MatrixXd> ineq_block_hessian(values_ineq.data() + ineq_nnz_idx, vert_i_dim_unfixed, vert_j_dim_unfixed);
                            edge->computeConstraintHessiansInc(row_i, col_j, eq_block_jacobian1, ineq_block_jacobian1, eq_block_hessian,
                                                               ineq_block_hessian, mult_eq_part, mult_ineq_part);
                            eq_nnz_idx += vert_i_dim_unfixed * vert_j_dim_unfixed;
                            ineq_nnz_idx += vert_i_dim_unfixed * vert_j_dim_unfixed;
                        }
                    }
                    else
                    {
                        if (!edge->isEqualityLinear())
                        {
                            if (lower_part_only && vertex1 == vertex2)
                            {
                                // we only need the lower triangular part, so cache Hessian
                                Eigen::MatrixXd block_hessian_ii(vert_i_dim_unfixed, vert_j_dim_unfixed);
                                edge->computeEqualityHessian(row_i, col_j, eq_block_jacobian1, block_hessian_ii, mult_eq_part);
                                for (int i = 0; i < vert_i_dim_unfixed; ++i)
                                {
                                    for (int j = 0; j < i + 1; ++j)
                                    {
                                        values_eq[eq_nnz_idx] += block_hessian_ii.coeffRef(i, j);
                                        ++eq_nnz_idx;
                                    }
                                    // nnz_idx += vert_i_dim_unfixed * (vert_i_dim_unfixed + 1) / 2;
                                }
                            }
                            else
                            {
                                Eigen::Map<Eigen::MatrixXd> eq_block_hessian(values_eq.data() + eq_nnz_idx, vert_i_dim_unfixed, vert_j_dim_unfixed);
                                edge->computeEqualityHessianInc(row_i, col_j, eq_block_jacobian1, eq_block_hessian, mult_eq_part);
                                eq_nnz_idx += vert_i_dim_unfixed * vert_j_dim_unfixed;
                            }
                        }
                        else if (!edge->isInequalityLinear())
                        {
                            if (lower_part_only && vertex1 == vertex2)
                            {
                                // we only need the lower triangular part, so cache Hessian
                                Eigen::MatrixXd block_hessian_ii(vert_i_dim_unfixed, vert_j_dim_unfixed);
                                edge->computeInequalityHessian(row_i, col_j, ineq_block_jacobian1, block_hessian_ii, mult_ineq_part);
                                for (int i = 0; i < vert_i_dim_unfixed; ++i)
                                {
                                    for (int j = 0; j < i + 1; ++j)
                                    {
                                        values_ineq[ineq_nnz_idx] += block_hessian_ii.coeffRef(i, j);
                                        ++ineq_nnz_idx;
                                    }
                                    // nnz_idx += vert_i_dim_unfixed * (vert_i_dim_unfixed + 1) / 2;
                                }
                            }
                            else
                            {
                                Eigen::Map<Eigen::MatrixXd> ineq_block_hessian(values_ineq.data() + ineq_nnz_idx, vert_i_dim_unfixed,
                                                                               vert_j_dim_unfixed);
                                edge->computeInequalityHessianInc(row_i, col_j, ineq_block_jacobian1, ineq_block_hessian, mult_ineq_part);
                                ineq_nnz_idx += vert_i_dim_unfixed * vert_j_dim_unfixed;
                            }
                        }
                    }
                }
                else  // TODO(roesmann): check all reasonable cases:
                {
                    if (!edge->isObjectiveLinear() && !edge->isEqualityLinear() && !edge->isInequalityLinear() &&
                        edge->getObjectiveDimension() != 0 && edge->getEqualityDimension() != 0 && edge->getInequalityDimension() != 0)
                    {
                        if (lower_part_only && vertex1 == vertex2)
                        {
                            // we only need the lower triangular part, so cache Hessian
                            Eigen::MatrixXd obj_block_hessian_ii(vert_i_dim_unfixed, vert_j_dim_unfixed);
                            Eigen::MatrixXd eq_block_hessian_ii(vert_i_dim_unfixed, vert_j_dim_unfixed);
                            Eigen::MatrixXd ineq_block_hessian_ii(vert_i_dim_unfixed, vert_j_dim_unfixed);
                            edge->computeHessians(row_i, col_j, obj_block_jacobian1, eq_block_jacobian1, ineq_block_jacobian1, obj_block_hessian_ii,
                                                  eq_block_hessian_ii, ineq_block_hessian_ii, mult_eq_part, mult_ineq_part, multiplier_obj);
                            for (int i = 0; i < vert_i_dim_unfixed; ++i)
                            {
                                for (int j = 0; j < i + 1; ++j)
                                {
                                    values_obj[obj_nnz_idx] += obj_block_hessian_ii.coeffRef(i, j);
                                    ++obj_nnz_idx;
                                    values_eq[eq_nnz_idx] += eq_block_hessian_ii.coeffRef(i, j);
                                    ++eq_nnz_idx;
                                    values_ineq[ineq_nnz_idx] += ineq_block_hessian_ii.coeffRef(i, j);
                                    ++ineq_nnz_idx;
                                }
                                // nnz_idx += vert_i_dim_unfixed * (vert_i_dim_unfixed + 1) / 2;
                            }
                        }
                        else
                        {
                            Eigen::Map<Eigen::MatrixXd> obj_block_hessian(values_obj.data() + obj_nnz_idx, vert_i_dim_unfixed, vert_j_dim_unfixed);
                            Eigen::Map<Eigen::MatrixXd> eq_block_hessian(values_eq.data() + eq_nnz_idx, vert_i_dim_unfixed, vert_j_dim_unfixed);
                            Eigen::Map<Eigen::MatrixXd> ineq_block_hessian(values_ineq.data() + ineq_nnz_idx, vert_i_dim_unfixed, vert_j_dim_unfixed);
                            edge->computeHessiansInc(row_i, col_j, obj_block_jacobian1, eq_block_jacobian1, ineq_block_jacobian1, obj_block_hessian,
                                                     eq_block_hessian, ineq_block_hessian, mult_eq_part, mult_ineq_part, multiplier_obj);
                            obj_nnz_idx += vert_i_dim_unfixed * vert_j_dim_unfixed;
                            eq_nnz_idx += vert_i_dim_unfixed * vert_j_dim_unfixed;
                            ineq_nnz_idx += vert_i_dim_unfixed * vert_j_dim_unfixed;
                        }
                    }
                    else if (!edge->isEqualityLinear() && !edge->isInequalityLinear() && edge->getEqualityDimension() != 0 &&
                             edge->getInequalityDimension() != 0)
                    {
                        if (lower_part_only && vertex1 == vertex2)
                        {
                            // we only need the lower triangular part, so cache Hessian
                            Eigen::MatrixXd eq_block_hessian_ii(vert_i_dim_unfixed, vert_j_dim_unfixed);
                            Eigen::MatrixXd ineq_block_hessian_ii(vert_i_dim_unfixed, vert_j_dim_unfixed);
                            edge->computeConstraintHessians(row_i, col_j, eq_block_jacobian1, ineq_block_jacobian1, eq_block_hessian_ii,
                                                            ineq_block_hessian_ii, mult_eq_part, mult_ineq_part);
                            for (int i = 0; i < vert_i_dim_unfixed; ++i)
                            {
                                for (int j = 0; j < i + 1; ++j)
                                {
                                    values_eq[eq_nnz_idx] += eq_block_hessian_ii.coeffRef(i, j);
                                    ++eq_nnz_idx;
                                    values_ineq[ineq_nnz_idx] += ineq_block_hessian_ii.coeffRef(i, j);
                                    ++ineq_nnz_idx;
                                }
                                // nnz_idx += vert_i_dim_unfixed * (vert_i_dim_unfixed + 1) / 2;
                            }
                        }
                        else
                        {
                            Eigen::Map<Eigen::MatrixXd> eq_block_hessian(values_eq.data() + eq_nnz_idx, vert_i_dim_unfixed, vert_j_dim_unfixed);
                            Eigen::Map<Eigen::MatrixXd> ineq_block_hessian(values_ineq.data() + ineq_nnz_idx, vert_i_dim_unfixed, vert_j_dim_unfixed);
                            edge->computeConstraintHessiansInc(row_i, col_j, eq_block_jacobian1, ineq_block_jacobian1, eq_block_hessian,
                                                               ineq_block_hessian, mult_eq_part, mult_ineq_part);
                            eq_nnz_idx += vert_i_dim_unfixed * vert_j_dim_unfixed;
                            ineq_nnz_idx += vert_i_dim_unfixed * vert_j_dim_unfixed;
                        }
                    }
                    else
                    {
                        if (!edge->isObjectiveLinear() && edge->getObjectiveDimension() != 0)
                        {
                            if (lower_part_only && vertex1 == vertex2)
                            {
                                // we only need the lower triangular part, so cache Hessian
                                Eigen::MatrixXd block_hessian_ii(vert_i_dim_unfixed, vert_j_dim_unfixed);
                                edge->computeObjectiveHessian(row_i, col_j, obj_block_jacobian1, block_hessian_ii, nullptr, multiplier_obj);
                                for (int i = 0; i < vert_i_dim_unfixed; ++i)
                                {
                                    for (int j = 0; j < i + 1; ++j)
                                    {
                                        values_obj[obj_nnz_idx] += block_hessian_ii.coeffRef(i, j);
                                        ++obj_nnz_idx;
                                    }
                                    // nnz_idx += vert_i_dim_unfixed * (vert_i_dim_unfixed + 1) / 2;
                                }
                            }
                            else
                            {
                                Eigen::Map<Eigen::MatrixXd> obj_block_hessian(values_obj.data() + obj_nnz_idx, vert_i_dim_unfixed,
                                                                              vert_j_dim_unfixed);
                                edge->computeObjectiveHessianInc(row_i, col_j, obj_block_jacobian1, obj_block_hessian, nullptr, multiplier_obj);
                                obj_nnz_idx += vert_i_dim_unfixed * vert_j_dim_unfixed;
                            }
                        }

                        if (!edge->isEqualityLinear() && edge->getEqualityDimension() != 0)
                        {
                            if (lower_part_only && vertex1 == vertex2)
                            {
                                // we only need the lower triangular part, so cache Hessian
                                Eigen::MatrixXd block_hessian_ii(vert_i_dim_unfixed, vert_j_dim_unfixed);
                                edge->computeEqualityHessian(row_i, col_j, eq_block_jacobian1, block_hessian_ii, mult_eq_part);
                                for (int i = 0; i < vert_i_dim_unfixed; ++i)
                                {
                                    for (int j = 0; j < i + 1; ++j)
                                    {
                                        values_eq[eq_nnz_idx] += block_hessian_ii.coeffRef(i, j);
                                        ++eq_nnz_idx;
                                    }
                                    // nnz_idx += vert_i_dim_unfixed * (vert_i_dim_unfixed + 1) / 2;
                                }
                            }
                            else
                            {
                                Eigen::Map<Eigen::MatrixXd> eq_block_hessian(values_eq.data() + eq_nnz_idx, vert_i_dim_unfixed, vert_j_dim_unfixed);
                                edge->computeEqualityHessianInc(row_i, col_j, eq_block_jacobian1, eq_block_hessian, mult_eq_part);
                                eq_nnz_idx += vert_i_dim_unfixed * vert_j_dim_unfixed;
                            }
                        }
                        else if (!edge->isInequalityLinear() && edge->getInequalityDimension() != 0)
                        {
                            if (lower_part_only && vertex1 == vertex2)
                            {
                                // we only need the lower triangular part, so cache Hessian
                                Eigen::MatrixXd block_hessian_ii(vert_i_dim_unfixed, vert_j_dim_unfixed);
                                edge->computeInequalityHessian(row_i, col_j, ineq_block_jacobian1, block_hessian_ii, mult_ineq_part);
                                for (int i = 0; i < vert_i_dim_unfixed; ++i)
                                {
                                    for (int j = 0; j < i + 1; ++j)
                                    {
                                        values_ineq[ineq_nnz_idx] += block_hessian_ii.coeffRef(i, j);
                                        ++ineq_nnz_idx;
                                    }
                                    // nnz_idx += vert_i_dim_unfixed * (vert_i_dim_unfixed + 1) / 2;
                                }
                            }
                            else
                            {
                                Eigen::Map<Eigen::MatrixXd> ineq_block_hessian(values_ineq.data() + ineq_nnz_idx, vert_i_dim_unfixed,
                                                                               vert_j_dim_unfixed);
                                edge->computeInequalityHessianInc(row_i, col_j, ineq_block_jacobian1, ineq_block_hessian, mult_ineq_part);
                                ineq_nnz_idx += vert_i_dim_unfixed * vert_j_dim_unfixed;
                            }
                        }
                    }
                }
            }
        }
    }
}

void HyperGraphOptimizationProblemEdgeBased::computeSparseHessianLagrangian(Eigen::SparseMatrix<double, Eigen::ColMajor, long long>& H,
                                                                            const double* multipliers_eq, const double* multipliers_ineq,
                                                                            const Eigen::VectorXi* col_nnz, bool upper_part_only)
{
    H.setZero();
    // TODO(roesmann): what is most efficient, if we just want to udpate the Hessian?
    if (col_nnz) H.reserve(*col_nnz);

    OptimizationEdgeSet::Ptr edges = _graph.getEdgeSet();

    // Iterate plain objective edges
    for (BaseEdge::Ptr& edge : edges->getObjectiveEdgesRef())
    {
        if (edge->isLinear()) continue;

        for (int row_i = 0; row_i < edge->getNumVertices(); ++row_i)
        {
            const VertexInterface* vertex1 = edge->getVertexRaw(row_i);

            int vert_i_dim_unfixed = vertex1->getDimensionUnfixed();
            if (vert_i_dim_unfixed == 0) continue;

            // compute jacobian to speed up hessian computation
            Eigen::MatrixXd block_jacobian1(edge->getDimension(), vertex1->getDimensionUnfixed());
            edge->computeJacobian(row_i, block_jacobian1);

            int col_start = upper_part_only ? row_i : 0;
            for (int col_j = col_start; col_j < edge->getNumVertices(); ++col_j)  // hessian is symmetric
            {
                const VertexInterface* vertex2 = edge->getVertexRaw(col_j);

                int vert_j_dim_unfixed = vertex2->getDimensionUnfixed();
                if (vert_j_dim_unfixed == 0) continue;

                // warning: we require the hessian to be initialized with zeros
                // and that edge->computeHessian only adds (+=) values!!
                Eigen::MatrixXd block_hessian(vert_i_dim_unfixed, vert_j_dim_unfixed);
                edge->computeHessian(row_i, col_j, block_jacobian1, block_hessian);

                // insert values
                if (upper_part_only && vertex1 == vertex2)
                {
                    // block on the main diagonal, so use upper part only
                    for (int i = 0; i < block_hessian.rows(); ++i)
                    {
                        for (int j = i; j < block_hessian.cols(); ++j)
                        {
                            H.coeffRef(vertex1->getVertexIdx() + i, vertex2->getVertexIdx() + j) += block_hessian.coeffRef(i, j);
                        }
                    }
                }
                else
                {
                    for (int i = 0; i < block_hessian.rows(); ++i)
                    {
                        for (int j = 0; j < block_hessian.cols(); ++j)
                        {
                            H.coeffRef(vertex1->getVertexIdx() + i, vertex2->getVertexIdx() + j) += block_hessian.coeffRef(i, j);
                        }
                    }
                }
            }
        }
    }

    // Iterate lsq objective edges
    for (BaseEdge::Ptr& edge : edges->getLsqObjectiveEdgesRef())
    {
        // if (edge->isLinear()) continue; // we must not check this one, since Linear-squared is not linear ;-)

        for (int row_i = 0; row_i < edge->getNumVertices(); ++row_i)
        {
            const VertexInterface* vertex1 = edge->getVertexRaw(row_i);

            int vert_i_dim_unfixed = vertex1->getDimensionUnfixed();
            if (vert_i_dim_unfixed == 0) continue;

            // compute jacobian to speed up hessian computation
            Eigen::MatrixXd block_jacobian1(edge->getDimension(), vertex1->getDimensionUnfixed());
            edge->computeJacobian(row_i, block_jacobian1);

            int col_start = upper_part_only ? row_i : 0;
            for (int col_j = col_start; col_j < edge->getNumVertices(); ++col_j)  // hessian is symmetric
            {
                const VertexInterface* vertex2 = edge->getVertexRaw(col_j);

                int vert_j_dim_unfixed = vertex2->getDimensionUnfixed();
                if (vert_j_dim_unfixed == 0) continue;

                // we need the second jacobian for applying the chain rule
                Eigen::MatrixXd block_jacobian2(edge->getDimension(), vertex2->getDimensionUnfixed());
                edge->computeJacobian(col_j, block_jacobian2);

                // approximate H = 2*J^T*J
                Eigen::MatrixXd block_hessian(vert_i_dim_unfixed, vert_j_dim_unfixed);
                block_hessian = 2.0 * block_jacobian1.transpose() * block_jacobian2;

                // insert values
                if (upper_part_only && vertex1 == vertex2)
                {
                    // block on the main diagonal, so use upper part only
                    for (int i = 0; i < block_hessian.rows(); ++i)
                    {
                        for (int j = i; j < block_hessian.cols(); ++j)
                        {
                            H.coeffRef(vertex1->getVertexIdx() + i, vertex2->getVertexIdx() + j) += block_hessian.coeffRef(i, j);
                        }
                    }
                }
                else
                {
                    for (int i = 0; i < block_hessian.rows(); ++i)
                    {
                        for (int j = 0; j < block_hessian.cols(); ++j)
                        {
                            H.coeffRef(vertex1->getVertexIdx() + i, vertex2->getVertexIdx() + j) += block_hessian.coeffRef(i, j);
                        }
                    }
                }
            }
        }
    }

    // Iterate equality edges
    for (BaseEdge::Ptr& edge : edges->getEqualityEdgesRef())
    {
        if (edge->isLinear()) continue;

        const double* mult_eq_part = multipliers_eq ? multipliers_eq + edge->getEdgeIdx() : nullptr;

        for (int row_i = 0; row_i < edge->getNumVertices(); ++row_i)
        {
            const VertexInterface* vertex1 = edge->getVertexRaw(row_i);

            int vert_i_dim_unfixed = vertex1->getDimensionUnfixed();
            if (vert_i_dim_unfixed == 0) continue;

            // compute jacobian to speed up hessian computation
            Eigen::MatrixXd block_jacobian1(edge->getDimension(), vertex1->getDimensionUnfixed());
            edge->computeJacobian(row_i, block_jacobian1);

            int col_start = upper_part_only ? row_i : 0;
            for (int col_j = col_start; col_j < edge->getNumVertices(); ++col_j)  // hessian is symmetric
            {
                const VertexInterface* vertex2 = edge->getVertexRaw(col_j);

                int vert_j_dim_unfixed = vertex2->getDimensionUnfixed();
                if (vert_j_dim_unfixed == 0) continue;

                // warning: we require the hessian to be initialized with zeros
                // and that edge->computeHessian only adds (+=) values!!
                Eigen::MatrixXd block_hessian(vert_i_dim_unfixed, vert_j_dim_unfixed);
                edge->computeHessian(row_i, col_j, block_jacobian1, block_hessian, mult_eq_part);

                // insert values
                if (upper_part_only && vertex1 == vertex2)
                {
                    // block on the main diagonal, so use upper part only
                    for (int i = 0; i < block_hessian.rows(); ++i)
                    {
                        for (int j = i; j < block_hessian.cols(); ++j)
                        {
                            H.coeffRef(vertex1->getVertexIdx() + i, vertex2->getVertexIdx() + j) += block_hessian.coeffRef(i, j);
                        }
                    }
                }
                else
                {
                    for (int i = 0; i < block_hessian.rows(); ++i)
                    {
                        for (int j = 0; j < block_hessian.cols(); ++j)
                        {
                            H.coeffRef(vertex1->getVertexIdx() + i, vertex2->getVertexIdx() + j) += block_hessian.coeffRef(i, j);
                        }
                    }
                }
            }
        }
    }

    // Iterate inequality edges
    for (BaseEdge::Ptr& edge : edges->getInequalityEdgesRef())
    {
        if (edge->isLinear()) continue;

        const double* mult_ineq_part = multipliers_ineq ? multipliers_ineq + edge->getEdgeIdx() : nullptr;

        for (int row_i = 0; row_i < edge->getNumVertices(); ++row_i)
        {
            const VertexInterface* vertex1 = edge->getVertexRaw(row_i);

            int vert_i_dim_unfixed = vertex1->getDimensionUnfixed();
            if (vert_i_dim_unfixed == 0) continue;

            // compute jacobian to speed up hessian computation
            Eigen::MatrixXd block_jacobian1(edge->getDimension(), vertex1->getDimensionUnfixed());
            edge->computeJacobian(row_i, block_jacobian1);

            int col_start = upper_part_only ? row_i : 0;
            for (int col_j = col_start; col_j < edge->getNumVertices(); ++col_j)  // hessian is symmetric
            {
                const VertexInterface* vertex2 = edge->getVertexRaw(col_j);

                int vert_j_dim_unfixed = vertex2->getDimensionUnfixed();
                if (vert_j_dim_unfixed == 0) continue;

                // warning: we require the hessian to be initialized with zeros
                // and that edge->computeHessian only adds (+=) values!!
                Eigen::MatrixXd block_hessian(vert_i_dim_unfixed, vert_j_dim_unfixed);
                edge->computeHessian(row_i, col_j, block_jacobian1, block_hessian, mult_ineq_part);

                // insert values
                if (upper_part_only && vertex1 == vertex2)
                {
                    // block on the main diagonal, so use upper part only
                    for (int i = 0; i < block_hessian.rows(); ++i)
                    {
                        for (int j = i; j < block_hessian.cols(); ++j)
                        {
                            H.coeffRef(vertex1->getVertexIdx() + i, vertex2->getVertexIdx() + j) += block_hessian.coeffRef(i, j);
                        }
                    }
                }
                else
                {
                    for (int i = 0; i < block_hessian.rows(); ++i)
                    {
                        for (int j = 0; j < block_hessian.cols(); ++j)
                        {
                            H.coeffRef(vertex1->getVertexIdx() + i, vertex2->getVertexIdx() + j) += block_hessian.coeffRef(i, j);
                        }
                    }
                }
            }
        }
    }

    // Iterate mixed edges
    for (BaseMixedEdge::Ptr& edge : edges->getMixedEdgesRef())
    {
        // if (edge->getObjectiveDimension() == 0 || edge->isObjectiveLinear()) continue;

        const double* mult_eq_part   = multipliers_eq ? multipliers_eq + edge->getEdgeEqualityIdx() : nullptr;
        const double* mult_ineq_part = multipliers_ineq ? multipliers_ineq + edge->getEdgeInequalityIdx() : nullptr;

        for (int row_i = 0; row_i < edge->getNumVertices(); ++row_i)
        {
            const VertexInterface* vertex1 = edge->getVertexRaw(row_i);

            int vert_i_dim_unfixed = vertex1->getDimensionUnfixed();
            if (vert_i_dim_unfixed == 0) continue;

            // compute jacobian to speed up hessian computation
            Eigen::MatrixXd obj_block_jacobian1(edge->getObjectiveDimension(), vertex1->getDimensionUnfixed());
            Eigen::MatrixXd eq_block_jacobian1(edge->getEqualityDimension(), vertex1->getDimensionUnfixed());
            Eigen::MatrixXd ineq_block_jacobian1(edge->getInequalityDimension(), vertex1->getDimensionUnfixed());
            edge->computeJacobians(row_i, obj_block_jacobian1, eq_block_jacobian1, ineq_block_jacobian1);

            int col_start = upper_part_only ? row_i : 0;
            for (int col_j = col_start; col_j < edge->getNumVertices(); ++col_j)  // hessian is symmetric
            {
                const VertexInterface* vertex2 = edge->getVertexRaw(col_j);

                int vert_j_dim_unfixed = vertex2->getDimensionUnfixed();
                if (vert_j_dim_unfixed == 0) continue;

                Eigen::MatrixXd block_hessian(vert_i_dim_unfixed, vert_j_dim_unfixed);
                block_hessian.setZero();

                if (edge->isObjectiveLeastSquaresForm())
                {
                    if (!edge->isObjectiveLinear() && edge->getObjectiveDimension() != 0)
                    {
                        // we need the second jacobian for applying the chain rule
                        Eigen::MatrixXd obj_block_jacobian2(edge->getObjectiveDimension(), vertex2->getDimensionUnfixed());
                        edge->computeObjectiveJacobian(col_j, obj_block_jacobian2);

                        // approximate H = 2*J^T*J
                        block_hessian += 2.0 * obj_block_jacobian1.transpose() * obj_block_jacobian2;
                    }

                    if (!edge->isEqualityLinear() && !edge->isInequalityLinear())
                    {
                        edge->computeConstraintHessiansInc(row_i, col_j, eq_block_jacobian1, ineq_block_jacobian1, block_hessian, block_hessian,
                                                           mult_eq_part, mult_ineq_part);
                    }
                    else
                    {
                        if (!edge->isEqualityLinear())
                        {
                            edge->computeEqualityHessianInc(row_i, col_j, eq_block_jacobian1, block_hessian, mult_eq_part);
                        }
                        else if (!edge->isInequalityLinear())
                        {
                            edge->computeInequalityHessianInc(row_i, col_j, ineq_block_jacobian1, block_hessian, mult_ineq_part);
                        }
                    }
                }
                else  // TODO(roesmann): check all reasonable cases:
                {
                    if (!edge->isObjectiveLinear() && !edge->isEqualityLinear() && !edge->isInequalityLinear() &&
                        edge->getObjectiveDimension() != 0 && edge->getEqualityDimension() != 0 && edge->getInequalityDimension() != 0)
                    {
                        edge->computeHessiansInc(row_i, col_j, obj_block_jacobian1, eq_block_jacobian1, ineq_block_jacobian1, block_hessian,
                                                 block_hessian, block_hessian, mult_eq_part, mult_ineq_part, 1.0);
                    }
                    else if (!edge->isEqualityLinear() && !edge->isInequalityLinear() && edge->getEqualityDimension() != 0 &&
                             edge->getInequalityDimension() != 0)
                    {
                        edge->computeConstraintHessiansInc(row_i, col_j, eq_block_jacobian1, ineq_block_jacobian1, block_hessian, block_hessian,
                                                           mult_eq_part, mult_ineq_part);
                    }
                    else
                    {
                        if (!edge->isObjectiveLinear() && edge->getObjectiveDimension() != 0)
                        {
                            edge->computeObjectiveHessianInc(row_i, col_j, obj_block_jacobian1, block_hessian, nullptr, 1.0);
                        }

                        if (!edge->isEqualityLinear() && edge->getEqualityDimension() != 0)
                        {
                            edge->computeEqualityHessianInc(row_i, col_j, eq_block_jacobian1, block_hessian, mult_eq_part);
                        }
                        else if (!edge->isInequalityLinear() && edge->getInequalityDimension() != 0)
                        {
                            edge->computeInequalityHessianInc(row_i, col_j, ineq_block_jacobian1, block_hessian, mult_ineq_part);
                        }
                    }
                }
                // insert values
                if (upper_part_only && vertex1 == vertex2)
                {
                    // block on the main diagonal, so use upper part only
                    for (int i = 0; i < block_hessian.rows(); ++i)
                    {
                        for (int j = i; j < block_hessian.cols(); ++j)
                        {
                            H.coeffRef(vertex1->getVertexIdx() + i, vertex2->getVertexIdx() + j) += block_hessian.coeffRef(i, j);
                        }
                    }
                }
                else
                {
                    for (int i = 0; i < block_hessian.rows(); ++i)
                    {
                        for (int j = 0; j < block_hessian.cols(); ++j)
                        {
                            H.coeffRef(vertex1->getVertexIdx() + i, vertex2->getVertexIdx() + j) += block_hessian.coeffRef(i, j);
                        }
                    }
                }
            }
        }
    }
}

void HyperGraphOptimizationProblemEdgeBased::computeSparseHessianLagrangianNNZperCol(Eigen::Ref<Eigen::VectorXi> col_nnz, bool upper_part_only)
{
    assert(col_nnz.size() == getParameterDimension());

    col_nnz.setZero();
    // TODO(roesmann): what is most efficient, if we just want to udpate the Hessian?

    // TODO(roesmann): with the following strategy, we overestimate the number of nnz, since vertices are shared among multiple edges...
    // however, we do not want to create the complete sparsity pattern here by creating a matrix of zeros and ones (tagging).

    OptimizationEdgeSet::Ptr edges = _graph.getEdgeSet();

    // Iterate plain objective edges
    for (BaseEdge::Ptr& edge : edges->getObjectiveEdgesRef())
    {
        if (edge->isLinear()) continue;

        for (int row_i = 0; row_i < edge->getNumVertices(); ++row_i)
        {
            const VertexInterface* vertex1 = edge->getVertexRaw(row_i);

            int vert_i_dim_unfixed = vertex1->getDimensionUnfixed();
            if (vert_i_dim_unfixed == 0) continue;

            int col_start = upper_part_only ? row_i : 0;
            for (int col_j = col_start; col_j < edge->getNumVertices(); ++col_j)  // hessian is symmetric
            {
                const VertexInterface* vertex2 = edge->getVertexRaw(col_j);

                int vert_j_dim_unfixed = vertex2->getDimensionUnfixed();
                if (vert_j_dim_unfixed == 0) continue;

                if (upper_part_only && vertex1 == vertex2)
                {
                    // just the upper triangular part of the bock matrix including diagonal, so N*(N+1)/2 elements
                    for (int l = 0; l < vert_j_dim_unfixed; ++l)
                    {
                        col_nnz(vertex2->getVertexIdx() + l) += l + 1;  // first column 0, second 1, ...., vert_i_dim_unfixed
                    }
                }
                else
                {
                    col_nnz.segment(vertex2->getVertexIdx(), vert_j_dim_unfixed).array() += vert_i_dim_unfixed;
                }
            }
        }
    }

    // Iterate lsq objective edges
    for (BaseEdge::Ptr& edge : edges->getLsqObjectiveEdgesRef())
    {
        // if (edge->isLinear()) continue; // we must not check this one, since Linear-squared is not linear ;-)

        for (int row_i = 0; row_i < edge->getNumVertices(); ++row_i)
        {
            const VertexInterface* vertex1 = edge->getVertexRaw(row_i);

            int vert_i_dim_unfixed = vertex1->getDimensionUnfixed();
            if (vert_i_dim_unfixed == 0) continue;

            int col_start = upper_part_only ? row_i : 0;
            for (int col_j = col_start; col_j < edge->getNumVertices(); ++col_j)  // hessian is symmetric
            {
                const VertexInterface* vertex2 = edge->getVertexRaw(col_j);

                int vert_j_dim_unfixed = vertex2->getDimensionUnfixed();
                if (vert_j_dim_unfixed == 0) continue;

                if (upper_part_only && vertex1 == vertex2)
                {
                    // just the upper triangular part of the bock matrix including diagonal, so N*(N+1)/2 elements
                    for (int l = 0; l < vert_j_dim_unfixed; ++l)
                    {
                        col_nnz(vertex2->getVertexIdx() + l) += l + 1;  // first column 0, second 1, ...., vert_i_dim_unfixed
                    }
                }
                else
                {
                    col_nnz.segment(vertex2->getVertexIdx(), vert_j_dim_unfixed).array() += vert_i_dim_unfixed;
                }
            }
        }
    }

    // Iterate equality edges
    for (BaseEdge::Ptr& edge : edges->getEqualityEdgesRef())
    {
        if (edge->isLinear()) continue;

        for (int row_i = 0; row_i < edge->getNumVertices(); ++row_i)
        {
            const VertexInterface* vertex1 = edge->getVertexRaw(row_i);

            int vert_i_dim_unfixed = vertex1->getDimensionUnfixed();
            if (vert_i_dim_unfixed == 0) continue;

            int col_start = upper_part_only ? row_i : 0;
            for (int col_j = col_start; col_j < edge->getNumVertices(); ++col_j)  // hessian is symmetric
            {
                const VertexInterface* vertex2 = edge->getVertexRaw(col_j);

                int vert_j_dim_unfixed = vertex2->getDimensionUnfixed();
                if (vert_j_dim_unfixed == 0) continue;

                if (upper_part_only && vertex1 == vertex2)
                {
                    // just the upper triangular part of the bock matrix including diagonal, so N*(N+1)/2 elements
                    for (int l = 0; l < vert_j_dim_unfixed; ++l)
                    {
                        col_nnz(vertex2->getVertexIdx() + l) += l + 1;  // first column 0, second 1, ...., vert_i_dim_unfixed
                    }
                }
                else
                {
                    col_nnz.segment(vertex2->getVertexIdx(), vert_j_dim_unfixed).array() += vert_i_dim_unfixed;
                }
            }
        }
    }

    // Iterate inequality edges
    for (BaseEdge::Ptr& edge : edges->getInequalityEdgesRef())
    {
        if (edge->isLinear()) continue;

        for (int row_i = 0; row_i < edge->getNumVertices(); ++row_i)
        {
            const VertexInterface* vertex1 = edge->getVertexRaw(row_i);

            int vert_i_dim_unfixed = vertex1->getDimensionUnfixed();
            if (vert_i_dim_unfixed == 0) continue;

            int col_start = upper_part_only ? row_i : 0;
            for (int col_j = col_start; col_j < edge->getNumVertices(); ++col_j)  // hessian is symmetric
            {
                const VertexInterface* vertex2 = edge->getVertexRaw(col_j);

                int vert_j_dim_unfixed = vertex2->getDimensionUnfixed();
                if (vert_j_dim_unfixed == 0) continue;

                if (upper_part_only && vertex1 == vertex2)
                {
                    // just the upper triangular part of the bock matrix including diagonal, so N*(N+1)/2 elements
                    for (int l = 0; l < vert_j_dim_unfixed; ++l)
                    {
                        col_nnz(vertex2->getVertexIdx() + l) += l + 1;  // first column 0, second 1, ...., vert_i_dim_unfixed
                    }
                }
                else
                {
                    col_nnz.segment(vertex2->getVertexIdx(), vert_j_dim_unfixed).array() += vert_i_dim_unfixed;
                }
            }
        }
    }

    // Iterate mixed edges
    for (BaseMixedEdge::Ptr& edge : edges->getMixedEdgesRef())
    {
        // if (edge->getObjectiveDimension() == 0 || edge->isObjectiveLinear()) continue;

        for (int row_i = 0; row_i < edge->getNumVertices(); ++row_i)
        {
            const VertexInterface* vertex1 = edge->getVertexRaw(row_i);

            int vert_i_dim_unfixed = vertex1->getDimensionUnfixed();
            if (vert_i_dim_unfixed == 0) continue;

            int col_start = upper_part_only ? row_i : 0;
            for (int col_j = col_start; col_j < edge->getNumVertices(); ++col_j)  // hessian is symmetric
            {
                const VertexInterface* vertex2 = edge->getVertexRaw(col_j);

                int vert_j_dim_unfixed = vertex2->getDimensionUnfixed();
                if (vert_j_dim_unfixed == 0) continue;

                if (upper_part_only && vertex1 == vertex2)
                {
                    // just the upper triangular part of the bock matrix including diagonal, so N*(N+1)/2 elements
                    for (int l = 0; l < vert_j_dim_unfixed; ++l)
                    {
                        col_nnz(vertex2->getVertexIdx() + l) += l + 1;  // first column 0, second 1, ...., vert_i_dim_unfixed
                    }
                }
                else
                {
                    col_nnz.segment(vertex2->getVertexIdx(), vert_j_dim_unfixed).array() += vert_i_dim_unfixed;
                }
            }
        }
    }
}

void HyperGraphOptimizationProblemEdgeBased::computeSparseJacobianTwoSideBoundedLinearForm(Eigen::SparseMatrix<double, Eigen::ColMajor, long long>& A,
                                                                                           bool include_finite_bounds, const Eigen::VectorXi* col_nnz)
{
    OptimizationEdgeSet::Ptr edges = _graph.getEdgeSet();

    A.setZero();  // we need to reset since we call insert below
    if (col_nnz) A.reserve(*col_nnz);

    int inequality_row_start = getEqualityDimension();
    int bounds_row_start     = inequality_row_start + getInequalityDimension();

    // Iterate equality edges
    for (BaseEdge::Ptr& edge : edges->getEqualityEdgesRef())
    {
        for (int vert_idx = 0; vert_idx < edge->getNumVertices(); ++vert_idx)
        {
            const VertexInterface* vertex = edge->getVertexRaw(vert_idx);

            int vert_dim_unfixed = vertex->getDimensionUnfixed();
            if (vert_dim_unfixed == 0) continue;

            // TODO(roesmann): Check if we can avoid the temporary here...
            Eigen::MatrixXd block_jacobian(edge->getDimension(), vert_dim_unfixed);
            edge->computeJacobian(vert_idx, block_jacobian);

            // Now insert values
            for (int i = 0; i < block_jacobian.rows(); ++i)
            {
                for (int j = 0; j < block_jacobian.cols(); ++j)
                {
                    A.insert(edge->getEdgeIdx() + i, vertex->getVertexIdx() + j) = block_jacobian.coeffRef(i, j);
                }
            }
        }
    }

    // Iterate inequality edges
    for (BaseEdge::Ptr& edge : edges->getInequalityEdgesRef())
    {
        int A_edge_idx = inequality_row_start + edge->getEdgeIdx();
        for (int vert_idx = 0; vert_idx < edge->getNumVertices(); ++vert_idx)
        {
            const VertexInterface* vertex = edge->getVertexRaw(vert_idx);

            int vert_dim_unfixed = vertex->getDimensionUnfixed();
            if (vert_dim_unfixed == 0) continue;

            // TODO(roesmann): Check if we can avoid the temporary here...
            Eigen::MatrixXd block_jacobian(edge->getDimension(), vert_dim_unfixed);
            edge->computeJacobian(vert_idx, block_jacobian);

            // Now insert values
            for (int i = 0; i < block_jacobian.rows(); ++i)
            {
                for (int j = 0; j < block_jacobian.cols(); ++j)
                {
                    A.insert(A_edge_idx + i, vertex->getVertexIdx() + j) = block_jacobian.coeffRef(i, j);
                }
            }
        }
    }

    // Iterate mixed edges
    for (BaseMixedEdge::Ptr& edge : edges->getMixedEdgesRef())
    {
        int A_edge_eq_idx   = edge->getEdgeEqualityIdx();
        int A_edge_ineq_idx = inequality_row_start + edge->getEdgeInequalityIdx();

        for (int vert_idx = 0; vert_idx < edge->getNumVertices(); ++vert_idx)
        {
            const VertexInterface* vertex = edge->getVertexRaw(vert_idx);

            int vert_dim_unfixed = vertex->getDimensionUnfixed();
            if (vert_dim_unfixed == 0) continue;

            // TODO(roesmann): Check if we can avoid the temporary here...
            Eigen::MatrixXd block_jacobian_eq(edge->getEqualityDimension(), vert_dim_unfixed);
            Eigen::MatrixXd block_jacobian_ineq(edge->getInequalityDimension(), vert_dim_unfixed);
            edge->computeConstraintJacobians(vert_idx, block_jacobian_eq, block_jacobian_ineq);

            // Now insert values
            for (int i = 0; i < block_jacobian_eq.rows(); ++i)
            {
                for (int j = 0; j < block_jacobian_eq.cols(); ++j)
                {
                    A.insert(A_edge_eq_idx + i, vertex->getVertexIdx() + j) = block_jacobian_eq.coeffRef(i, j);
                }
            }
            for (int i = 0; i < block_jacobian_ineq.rows(); ++i)
            {
                for (int j = 0; j < block_jacobian_ineq.cols(); ++j)
                {
                    A.insert(A_edge_ineq_idx + i, vertex->getVertexIdx() + j) = block_jacobian_ineq.coeffRef(i, j);
                }
            }
        }
    }

    if (include_finite_bounds)
    {
        int row_idx = bounds_row_start;
        for (const VertexInterface* vertex : _graph.getVertexSet()->getActiveVertices())
        {
            int vert_idx = vertex->getVertexIdx();
            int free_idx = 0;
            for (int i = 0; i < vertex->getDimension(); ++i)
            {
                if (vertex->isFixedComponent(i)) continue;

                if (vertex->hasFiniteLowerBound(i) || vertex->hasFiniteUpperBound(i))
                {
                    A.insert(row_idx, vert_idx + free_idx) = 1;
                    ++row_idx;
                }
                ++free_idx;
            }
        }
    }
}

void HyperGraphOptimizationProblemEdgeBased::computeSparseJacobianTwoSideBoundedLinearFormNNZPerColumn(Eigen::Ref<Eigen::VectorXi> col_nnz,
                                                                                                       bool include_finite_bounds)
{
    assert(col_nnz.size() == getParameterDimension());

    col_nnz.setZero();

    OptimizationEdgeSet::Ptr edges = _graph.getEdgeSet();

    // Iterate equality edges
    for (BaseEdge::Ptr& edge : edges->getEqualityEdgesRef())
    {
        for (int vert_idx = 0; vert_idx < edge->getNumVertices(); ++vert_idx)
        {
            const VertexInterface* vertex = edge->getVertexRaw(vert_idx);
            int vert_dim_unfixed          = vertex->getDimensionUnfixed();
            if (vert_dim_unfixed == 0) continue;

            col_nnz.segment(vertex->getVertexIdx(), vert_dim_unfixed).array() += edge->getDimension();
        }
    }

    // Iterate equality edges
    for (BaseEdge::Ptr& edge : edges->getInequalityEdgesRef())
    {
        for (int vert_idx = 0; vert_idx < edge->getNumVertices(); ++vert_idx)
        {
            const VertexInterface* vertex = edge->getVertexRaw(vert_idx);
            int vert_dim_unfixed          = vertex->getDimensionUnfixed();
            if (vert_dim_unfixed == 0) continue;

            col_nnz.segment(vertex->getVertexIdx(), vert_dim_unfixed).array() += edge->getDimension();
        }
    }

    // Iterate mixed edges
    for (BaseMixedEdge::Ptr& edge : edges->getMixedEdgesRef())
    {
        for (int vert_idx = 0; vert_idx < edge->getNumVertices(); ++vert_idx)
        {
            const VertexInterface* vertex = edge->getVertexRaw(vert_idx);

            int vert_dim_unfixed = vertex->getDimensionUnfixed();
            if (vert_dim_unfixed == 0) continue;

            col_nnz.segment(vertex->getVertexIdx(), vert_dim_unfixed).array() += edge->getEqualityDimension() + edge->getInequalityDimension();
        }
    }

    if (include_finite_bounds)
    {
        for (const VertexInterface* vertex : _graph.getVertexSet()->getActiveVertices())
        {
            if (!vertex->hasFixedComponents() || !vertex->hasFiniteBounds()) continue;

            int vert_idx = vertex->getVertexIdx();
            int free_idx = 0;
            for (int i = 0; i < vertex->getDimension(); ++i)
            {
                if (vertex->isFixedComponent(i)) continue;

                if (vertex->hasFiniteLowerBound(i) || vertex->hasFiniteUpperBound(i))
                {
                    col_nnz[vert_idx + free_idx] += 1;
                }
                ++free_idx;
            }
        }
    }
}

int HyperGraphOptimizationProblemEdgeBased::computeSparseJacobianTwoSideBoundedLinearFormNNZ(bool include_finite_bounds)
{
    int nnz = 0;
    nnz += computeSparseJacobianEqualitiesNNZ();
    nnz += computeSparseJacobianInequalitiesNNZ();
    if (include_finite_bounds) nnz += finiteCombinedBoundsDimension();
    return nnz;
}

void HyperGraphOptimizationProblemEdgeBased::computeSparseJacobianTwoSideBoundedLinearFormStructure(Eigen::Ref<Eigen::VectorXi> i_row,
                                                                                                    Eigen::Ref<Eigen::VectorXi> j_col,
                                                                                                    bool include_finite_bounds)
{
    assert(i_row.size() == computeCombinedSparseJacobiansNNZ(include_finite_bounds));
    assert(j_col.size() == i_row.size());

    OptimizationEdgeSet::Ptr edges = _graph.getEdgeSet();

    int nnz_idx = 0;

    int inequality_row_start = getEqualityDimension();
    int bounds_row_start     = inequality_row_start + getInequalityDimension();

    // Iterate equality edges
    for (BaseEdge::Ptr& edge : edges->getEqualityEdgesRef())
    {
        for (int vert_idx = 0; vert_idx < edge->getNumVertices(); ++vert_idx)
        {
            const VertexInterface* vertex = edge->getVertexRaw(vert_idx);
            // Iterate all free variables
            int idx_free = 0;
            for (int i = 0; i < vertex->getDimension(); ++i)
            {
                if (!vertex->isFixedComponent(i))
                {
                    // Iterate inner edge dimension
                    for (int j = 0; j < edge->getDimension(); ++j)
                    {
                        i_row[nnz_idx] = edge->getEdgeIdx() + j;
                        j_col[nnz_idx] = vertex->getVertexIdx() + idx_free;
                        ++nnz_idx;
                    }
                    ++idx_free;
                }
            }
        }
    }

    // Iterate equality edges
    for (BaseEdge::Ptr& edge : edges->getInequalityEdgesRef())
    {
        for (int vert_idx = 0; vert_idx < edge->getNumVertices(); ++vert_idx)
        {
            const VertexInterface* vertex = edge->getVertexRaw(vert_idx);
            // Iterate all free variables
            int idx_free = 0;
            for (int i = 0; i < vertex->getDimension(); ++i)
            {
                if (!vertex->isFixedComponent(i))
                {
                    // Iterate inner edge dimension
                    for (int j = 0; j < edge->getDimension(); ++j)
                    {
                        i_row[nnz_idx] = inequality_row_start + edge->getEdgeIdx() + j;
                        j_col[nnz_idx] = vertex->getVertexIdx() + idx_free;
                        ++nnz_idx;
                    }
                    ++idx_free;
                }
            }
        }
    }

    // Iterate mixed edges
    for (BaseMixedEdge::Ptr& edge : edges->getMixedEdgesRef())
    {
        for (int vert_idx = 0; vert_idx < edge->getNumVertices(); ++vert_idx)
        {
            const VertexInterface* vertex = edge->getVertexRaw(vert_idx);

            // we need to do the following in the same order as in computeSparseJacobianTwoSideBoundedLinearFormValues()
            // Iterate all free variables
            int idx_free = 0;
            for (int i = 0; i < vertex->getDimension(); ++i)
            {
                if (!vertex->isFixedComponent(i))
                {
                    for (int j = 0; j < edge->getEqualityDimension(); ++j)
                    {
                        i_row[nnz_idx] = edge->getEdgeEqualityIdx() + j;
                        j_col[nnz_idx] = vertex->getVertexIdx() + idx_free;
                        ++nnz_idx;
                    }
                    ++idx_free;
                }
            }

            // Iterate all free variables
            idx_free = 0;
            for (int i = 0; i < vertex->getDimension(); ++i)
            {
                if (!vertex->isFixedComponent(i))
                {
                    for (int j = 0; j < edge->getInequalityDimension(); ++j)
                    {
                        i_row[nnz_idx] = inequality_row_start + edge->getEdgeInequalityIdx() + j;
                        j_col[nnz_idx] = vertex->getVertexIdx() + idx_free;
                        ++nnz_idx;
                    }
                    ++idx_free;
                }
            }
        }
    }

    if (include_finite_bounds)
    {
        int row_idx = bounds_row_start;
        for (const VertexInterface* vertex : _graph.getVertexSet()->getActiveVertices())
        {
            int vert_idx = vertex->getVertexIdx();
            int free_idx = 0;
            for (int i = 0; i < vertex->getDimension(); ++i)
            {
                if (vertex->isFixedComponent(i)) continue;

                if (vertex->hasFiniteLowerBound(i) || vertex->hasFiniteUpperBound(i))
                {
                    i_row[nnz_idx] = row_idx;
                    j_col[nnz_idx] = vert_idx + free_idx;

                    ++row_idx;
                    ++nnz_idx;
                }
                ++free_idx;
            }
        }
    }

    assert(nnz_idx == i_row.size());
}

void HyperGraphOptimizationProblemEdgeBased::computeSparseJacobianTwoSideBoundedLinearFormValues(Eigen::Ref<Eigen::VectorXd> values,
                                                                                                 bool include_finite_bounds)
{
    assert(values.size() == computeSparseJacobianTwoSideBoundedLinearFormNNZ(include_finite_bounds));

    int nnz_idx = 0;

    OptimizationEdgeSet::Ptr edges = _graph.getEdgeSet();

    // Iterate equality edges
    for (BaseEdge::Ptr& edge : edges->getEqualityEdgesRef())
    {
        for (int vert_idx = 0; vert_idx < edge->getNumVertices(); ++vert_idx)
        {
            int vert_dim_unfixed = edge->getVertexRaw(vert_idx)->getDimensionUnfixed();
            if (vert_dim_unfixed == 0) continue;

            Eigen::Map<Eigen::MatrixXd> block_jacobian(values.data() + nnz_idx, edge->getDimension(), vert_dim_unfixed);
            edge->computeJacobian(vert_idx, block_jacobian);
            nnz_idx += block_jacobian.rows() * block_jacobian.cols();
        }
    }

    // Iterate inequality edges
    for (BaseEdge::Ptr& edge : edges->getInequalityEdgesRef())
    {
        for (int vert_idx = 0; vert_idx < edge->getNumVertices(); ++vert_idx)
        {
            int vert_dim_unfixed = edge->getVertexRaw(vert_idx)->getDimensionUnfixed();
            if (vert_dim_unfixed == 0) continue;

            Eigen::Map<Eigen::MatrixXd> block_jacobian(values.data() + nnz_idx, edge->getDimension(), vert_dim_unfixed);
            edge->computeJacobian(vert_idx, block_jacobian);
            nnz_idx += block_jacobian.rows() * block_jacobian.cols();
        }
    }

    // Iterate mixed edges
    for (BaseMixedEdge::Ptr& edge : edges->getMixedEdgesRef())
    {
        if (edge->getEqualityDimension() == 0 && edge->getInequalityDimension() == 0) continue;
        // TODO(roesmann) implement the other cases as well

        for (int vert_idx = 0; vert_idx < edge->getNumVertices(); ++vert_idx)
        {
            int vert_dim_unfixed = edge->getVertexRaw(vert_idx)->getDimensionUnfixed();
            if (vert_dim_unfixed == 0) continue;

            Eigen::Map<Eigen::MatrixXd> block_eq(values.data() + nnz_idx, edge->getEqualityDimension(), vert_dim_unfixed);
            nnz_idx += block_eq.rows() * block_eq.cols();
            Eigen::Map<Eigen::MatrixXd> block_ineq(values.data() + nnz_idx, edge->getInequalityDimension(), vert_dim_unfixed);
            nnz_idx += block_ineq.rows() * block_ineq.cols();

            edge->computeConstraintJacobians(vert_idx, block_eq, block_ineq);
        }
    }

    if (include_finite_bounds)
    {
        values.tail(finiteCombinedBoundsDimension()).setOnes();
    }

    assert(nnz_idx + finiteCombinedBoundsDimension() - 1 == values.size());
}

void HyperGraphOptimizationProblemEdgeBased::computeSparseJacobianTwoSideBoundedLinearFormAndHessianLagrangian(
    Eigen::SparseMatrix<double, Eigen::ColMajor, long long>& H, const double* multipliers_eq, const double* multipliers_ineq,
    Eigen::SparseMatrix<double, Eigen::ColMajor, long long>& A, bool include_finite_bounds, const Eigen::VectorXi* col_nnz_H,
    const Eigen::VectorXi* col_nnz_A, bool upper_part_only_H)
{
    H.setZero();
    // TODO(roesmann): what is most efficient, if we just want to udpate the Hessian?
    if (col_nnz_H) H.reserve(*col_nnz_H);

    A.setZero();  // we need to reset since we call insert below
    if (col_nnz_A) A.reserve(*col_nnz_A);

    int inequality_row_start = getEqualityDimension();
    int bounds_row_start     = inequality_row_start + getInequalityDimension();

    OptimizationEdgeSet::Ptr edges = _graph.getEdgeSet();

    // Iterate plain objective edges
    for (BaseEdge::Ptr& edge : edges->getObjectiveEdgesRef())
    {
        for (int row_i = 0; row_i < edge->getNumVertices(); ++row_i)
        {
            const VertexInterface* vertex1 = edge->getVertexRaw(row_i);

            int vert_i_dim_unfixed = vertex1->getDimensionUnfixed();
            if (vert_i_dim_unfixed == 0) continue;

            // compute jacobian to speed up hessian computation
            Eigen::MatrixXd block_jacobian1(edge->getDimension(), vertex1->getDimensionUnfixed());
            edge->computeJacobian(row_i, block_jacobian1);

            if (edge->isLinear()) continue;

            int col_start = upper_part_only_H ? row_i : 0;
            for (int col_j = col_start; col_j < edge->getNumVertices(); ++col_j)  // hessian is symmetric
            {
                const VertexInterface* vertex2 = edge->getVertexRaw(col_j);

                int vert_j_dim_unfixed = vertex2->getDimensionUnfixed();
                if (vert_j_dim_unfixed == 0) continue;

                // warning: we require the hessian to be initialized with zeros
                // and that edge->computeHessian only adds (+=) values!!
                Eigen::MatrixXd block_hessian(vert_i_dim_unfixed, vert_j_dim_unfixed);
                edge->computeHessian(row_i, col_j, block_jacobian1, block_hessian);

                // insert values
                if (upper_part_only_H && vertex1 == vertex2)
                {
                    // block on the main diagonal, so use upper part only
                    for (int i = 0; i < block_hessian.rows(); ++i)
                    {
                        for (int j = i; j < block_hessian.cols(); ++j)
                        {
                            H.coeffRef(vertex1->getVertexIdx() + i, vertex2->getVertexIdx() + j) += block_hessian.coeffRef(i, j);
                        }
                    }
                }
                else
                {
                    for (int i = 0; i < block_hessian.rows(); ++i)
                    {
                        for (int j = 0; j < block_hessian.cols(); ++j)
                        {
                            H.coeffRef(vertex1->getVertexIdx() + i, vertex2->getVertexIdx() + j) += block_hessian.coeffRef(i, j);
                        }
                    }
                }
            }
        }
    }

    // Iterate lsq objective edges
    for (BaseEdge::Ptr& edge : edges->getLsqObjectiveEdgesRef())
    {
        // if (edge->isLinear()) continue; // we must not check this one, since Linear-squared is not linear ;-)

        for (int row_i = 0; row_i < edge->getNumVertices(); ++row_i)
        {
            const VertexInterface* vertex1 = edge->getVertexRaw(row_i);

            int vert_i_dim_unfixed = vertex1->getDimensionUnfixed();
            if (vert_i_dim_unfixed == 0) continue;

            // compute jacobian to speed up hessian computation
            Eigen::MatrixXd block_jacobian1(edge->getDimension(), vertex1->getDimensionUnfixed());
            edge->computeJacobian(row_i, block_jacobian1);

            if (edge->isLinear()) continue;

            int col_start = upper_part_only_H ? row_i : 0;
            for (int col_j = col_start; col_j < edge->getNumVertices(); ++col_j)  // hessian is symmetric
            {
                const VertexInterface* vertex2 = edge->getVertexRaw(col_j);

                int vert_j_dim_unfixed = vertex2->getDimensionUnfixed();
                if (vert_j_dim_unfixed == 0) continue;

                // we need the second jacobian for applying the chain rule
                Eigen::MatrixXd block_jacobian2(edge->getDimension(), vertex2->getDimensionUnfixed());
                edge->computeJacobian(col_j, block_jacobian2);

                // approximate H = 2*J^T*J
                Eigen::MatrixXd block_hessian(vert_i_dim_unfixed, vert_j_dim_unfixed);
                block_hessian = 2.0 * block_jacobian1.transpose() * block_jacobian2;

                // insert values
                if (upper_part_only_H && vertex1 == vertex2)
                {
                    // block on the main diagonal, so use upper part only
                    for (int i = 0; i < block_hessian.rows(); ++i)
                    {
                        for (int j = i; j < block_hessian.cols(); ++j)
                        {
                            H.coeffRef(vertex1->getVertexIdx() + i, vertex2->getVertexIdx() + j) += block_hessian.coeffRef(i, j);
                        }
                    }
                }
                else
                {
                    for (int i = 0; i < block_hessian.rows(); ++i)
                    {
                        for (int j = 0; j < block_hessian.cols(); ++j)
                        {
                            H.coeffRef(vertex1->getVertexIdx() + i, vertex2->getVertexIdx() + j) += block_hessian.coeffRef(i, j);
                        }
                    }
                }
            }
        }
    }

    // Iterate equality edges
    for (BaseEdge::Ptr& edge : edges->getEqualityEdgesRef())
    {
        const double* mult_eq_part = multipliers_eq ? multipliers_eq + edge->getEdgeIdx() : nullptr;

        for (int row_i = 0; row_i < edge->getNumVertices(); ++row_i)
        {
            const VertexInterface* vertex1 = edge->getVertexRaw(row_i);

            int vert_i_dim_unfixed = vertex1->getDimensionUnfixed();
            if (vert_i_dim_unfixed == 0) continue;

            // compute jacobian to speed up hessian computation
            Eigen::MatrixXd block_jacobian1(edge->getDimension(), vertex1->getDimensionUnfixed());
            edge->computeJacobian(row_i, block_jacobian1);

            // Now insert values to Jacobian matrix A
            for (int i = 0; i < block_jacobian1.rows(); ++i)
            {
                for (int j = 0; j < block_jacobian1.cols(); ++j)
                {
                    A.insert(edge->getEdgeIdx() + i, vertex1->getVertexIdx() + j) = block_jacobian1.coeffRef(i, j);
                }
            }

            if (edge->isLinear()) continue;

            int col_start = upper_part_only_H ? row_i : 0;
            for (int col_j = col_start; col_j < edge->getNumVertices(); ++col_j)  // hessian is symmetric
            {
                const VertexInterface* vertex2 = edge->getVertexRaw(col_j);

                int vert_j_dim_unfixed = vertex2->getDimensionUnfixed();
                if (vert_j_dim_unfixed == 0) continue;

                // warning: we require the hessian to be initialized with zeros
                // and that edge->computeHessian only adds (+=) values!!
                Eigen::MatrixXd block_hessian(vert_i_dim_unfixed, vert_j_dim_unfixed);
                edge->computeHessian(row_i, col_j, block_jacobian1, block_hessian, mult_eq_part);

                // insert values
                if (upper_part_only_H && vertex1 == vertex2)
                {
                    // block on the main diagonal, so use upper part only
                    for (int i = 0; i < block_hessian.rows(); ++i)
                    {
                        for (int j = i; j < block_hessian.cols(); ++j)
                        {
                            H.coeffRef(vertex1->getVertexIdx() + i, vertex2->getVertexIdx() + j) += block_hessian.coeffRef(i, j);
                        }
                    }
                }
                else
                {
                    for (int i = 0; i < block_hessian.rows(); ++i)
                    {
                        for (int j = 0; j < block_hessian.cols(); ++j)
                        {
                            H.coeffRef(vertex1->getVertexIdx() + i, vertex2->getVertexIdx() + j) += block_hessian.coeffRef(i, j);
                        }
                    }
                }
            }
        }
    }

    // Iterate inequality edges
    for (BaseEdge::Ptr& edge : edges->getInequalityEdgesRef())
    {
        const double* mult_ineq_part = multipliers_ineq ? multipliers_ineq + edge->getEdgeIdx() : nullptr;
        int A_edge_idx               = inequality_row_start + edge->getEdgeIdx();

        for (int row_i = 0; row_i < edge->getNumVertices(); ++row_i)
        {
            const VertexInterface* vertex1 = edge->getVertexRaw(row_i);

            int vert_i_dim_unfixed = vertex1->getDimensionUnfixed();
            if (vert_i_dim_unfixed == 0) continue;

            // compute jacobian to speed up hessian computation
            Eigen::MatrixXd block_jacobian1(edge->getDimension(), vertex1->getDimensionUnfixed());
            edge->computeJacobian(row_i, block_jacobian1);

            // Now insert values to Jacobian A
            for (int i = 0; i < block_jacobian1.rows(); ++i)
            {
                for (int j = 0; j < block_jacobian1.cols(); ++j)
                {
                    A.insert(A_edge_idx + i, vertex1->getVertexIdx() + j) = block_jacobian1.coeffRef(i, j);
                }
            }

            if (edge->isLinear()) continue;

            int col_start = upper_part_only_H ? row_i : 0;
            for (int col_j = col_start; col_j < edge->getNumVertices(); ++col_j)  // hessian is symmetric
            {
                const VertexInterface* vertex2 = edge->getVertexRaw(col_j);

                int vert_j_dim_unfixed = vertex2->getDimensionUnfixed();
                if (vert_j_dim_unfixed == 0) continue;

                // warning: we require the hessian to be initialized with zeros
                // and that edge->computeHessian only adds (+=) values!!
                Eigen::MatrixXd block_hessian(vert_i_dim_unfixed, vert_j_dim_unfixed);
                edge->computeHessian(row_i, col_j, block_jacobian1, block_hessian, mult_ineq_part);

                // insert values
                if (upper_part_only_H && vertex1 == vertex2)
                {
                    // block on the main diagonal, so use upper part only
                    for (int i = 0; i < block_hessian.rows(); ++i)
                    {
                        for (int j = i; j < block_hessian.cols(); ++j)
                        {
                            H.coeffRef(vertex1->getVertexIdx() + i, vertex2->getVertexIdx() + j) += block_hessian.coeffRef(i, j);
                        }
                    }
                }
                else
                {
                    for (int i = 0; i < block_hessian.rows(); ++i)
                    {
                        for (int j = 0; j < block_hessian.cols(); ++j)
                        {
                            H.coeffRef(vertex1->getVertexIdx() + i, vertex2->getVertexIdx() + j) += block_hessian.coeffRef(i, j);
                        }
                    }
                }
            }
        }
    }

    // Iterate mixed edges
    for (BaseMixedEdge::Ptr& edge : edges->getMixedEdgesRef())
    {
        // if (edge->getObjectiveDimension() == 0 || edge->isObjectiveLinear()) continue;

        const double* mult_eq_part   = multipliers_eq ? multipliers_eq + edge->getEdgeEqualityIdx() : nullptr;
        const double* mult_ineq_part = multipliers_ineq ? multipliers_ineq + edge->getEdgeInequalityIdx() : nullptr;

        int A_edge_eq_idx   = edge->getEdgeEqualityIdx();
        int A_edge_ineq_idx = inequality_row_start + edge->getEdgeInequalityIdx();

        for (int row_i = 0; row_i < edge->getNumVertices(); ++row_i)
        {
            const VertexInterface* vertex1 = edge->getVertexRaw(row_i);

            int vert_i_dim_unfixed = vertex1->getDimensionUnfixed();
            if (vert_i_dim_unfixed == 0) continue;

            // compute jacobian to speed up hessian computation
            Eigen::MatrixXd obj_block_jacobian1(edge->getObjectiveDimension(), vertex1->getDimensionUnfixed());
            Eigen::MatrixXd eq_block_jacobian1(edge->getEqualityDimension(), vertex1->getDimensionUnfixed());
            Eigen::MatrixXd ineq_block_jacobian1(edge->getInequalityDimension(), vertex1->getDimensionUnfixed());
            edge->computeJacobians(row_i, obj_block_jacobian1, eq_block_jacobian1, ineq_block_jacobian1);

            // Now insert values to Jacobian A
            for (int i = 0; i < eq_block_jacobian1.rows(); ++i)
            {
                for (int j = 0; j < eq_block_jacobian1.cols(); ++j)
                {
                    A.insert(A_edge_eq_idx + i, vertex1->getVertexIdx() + j) = eq_block_jacobian1.coeffRef(i, j);
                }
            }
            for (int i = 0; i < ineq_block_jacobian1.rows(); ++i)
            {
                for (int j = 0; j < ineq_block_jacobian1.cols(); ++j)
                {
                    A.insert(A_edge_ineq_idx + i, vertex1->getVertexIdx() + j) = ineq_block_jacobian1.coeffRef(i, j);
                }
            }

            int col_start = upper_part_only_H ? row_i : 0;
            for (int col_j = col_start; col_j < edge->getNumVertices(); ++col_j)  // hessian is symmetric
            {
                const VertexInterface* vertex2 = edge->getVertexRaw(col_j);

                int vert_j_dim_unfixed = vertex2->getDimensionUnfixed();
                if (vert_j_dim_unfixed == 0) continue;

                Eigen::MatrixXd block_hessian(vert_i_dim_unfixed, vert_j_dim_unfixed);
                block_hessian.setZero();

                if (edge->isObjectiveLeastSquaresForm())
                {
                    if (!edge->isObjectiveLinear() && edge->getObjectiveDimension() != 0)
                    {
                        // we need the second jacobian for applying the chain rule
                        Eigen::MatrixXd obj_block_jacobian2(edge->getObjectiveDimension(), vertex2->getDimensionUnfixed());
                        edge->computeObjectiveJacobian(col_j, obj_block_jacobian2);

                        // approximate H = 2*J^T*J
                        block_hessian += 2.0 * obj_block_jacobian1.transpose() * obj_block_jacobian2;
                    }

                    if (!edge->isEqualityLinear() && !edge->isInequalityLinear())
                    {
                        edge->computeConstraintHessiansInc(row_i, col_j, eq_block_jacobian1, ineq_block_jacobian1, block_hessian, block_hessian,
                                                           mult_eq_part, mult_ineq_part);
                    }
                    else
                    {
                        if (!edge->isEqualityLinear())
                        {
                            edge->computeEqualityHessianInc(row_i, col_j, eq_block_jacobian1, block_hessian, mult_eq_part);
                        }
                        else if (!edge->isInequalityLinear())
                        {
                            edge->computeInequalityHessianInc(row_i, col_j, ineq_block_jacobian1, block_hessian, mult_ineq_part);
                        }
                    }
                }
                else  // TODO(roesmann): check all reasonable cases:
                {
                    if (!edge->isObjectiveLinear() && !edge->isEqualityLinear() && !edge->isInequalityLinear() &&
                        edge->getObjectiveDimension() != 0 && edge->getEqualityDimension() != 0 && edge->getInequalityDimension() != 0)
                    {
                        edge->computeHessiansInc(row_i, col_j, obj_block_jacobian1, eq_block_jacobian1, ineq_block_jacobian1, block_hessian,
                                                 block_hessian, block_hessian, mult_eq_part, mult_ineq_part, 1.0);
                    }
                    else if (!edge->isEqualityLinear() && !edge->isInequalityLinear() && edge->getEqualityDimension() != 0 &&
                             edge->getInequalityDimension() != 0)
                    {
                        edge->computeConstraintHessiansInc(row_i, col_j, eq_block_jacobian1, ineq_block_jacobian1, block_hessian, block_hessian,
                                                           mult_eq_part, mult_ineq_part);
                    }
                    else
                    {
                        if (!edge->isObjectiveLinear() && edge->getObjectiveDimension() != 0)
                        {
                            edge->computeObjectiveHessianInc(row_i, col_j, obj_block_jacobian1, block_hessian, nullptr, 1.0);
                        }

                        if (!edge->isEqualityLinear() && edge->getEqualityDimension() != 0)
                        {
                            edge->computeEqualityHessianInc(row_i, col_j, eq_block_jacobian1, block_hessian, mult_eq_part);
                        }
                        else if (!edge->isInequalityLinear() && edge->getInequalityDimension() != 0)
                        {
                            edge->computeInequalityHessianInc(row_i, col_j, ineq_block_jacobian1, block_hessian, mult_ineq_part);
                        }
                    }
                }
                // insert values
                if (upper_part_only_H && vertex1 == vertex2)
                {
                    // block on the main diagonal, so use upper part only
                    for (int i = 0; i < block_hessian.rows(); ++i)
                    {
                        for (int j = i; j < block_hessian.cols(); ++j)
                        {
                            H.coeffRef(vertex1->getVertexIdx() + i, vertex2->getVertexIdx() + j) += block_hessian.coeffRef(i, j);
                        }
                    }
                }
                else
                {
                    for (int i = 0; i < block_hessian.rows(); ++i)
                    {
                        for (int j = 0; j < block_hessian.cols(); ++j)
                        {
                            H.coeffRef(vertex1->getVertexIdx() + i, vertex2->getVertexIdx() + j) += block_hessian.coeffRef(i, j);
                        }
                    }
                }
            }
        }
    }

    if (include_finite_bounds)
    {
        int row_idx = bounds_row_start;
        for (const VertexInterface* vertex : _graph.getVertexSet()->getActiveVertices())
        {
            int vert_idx = vertex->getVertexIdx();
            int free_idx = 0;
            for (int i = 0; i < vertex->getDimension(); ++i)
            {
                if (vertex->isFixedComponent(i)) continue;

                if (vertex->hasFiniteLowerBound(i) || vertex->hasFiniteUpperBound(i))
                {
                    A.insert(row_idx, vert_idx + free_idx) = 1;
                    ++row_idx;
                }
                ++free_idx;
            }
        }
    }
}

void HyperGraphOptimizationProblemEdgeBased::computeSparseJacobianTwoSideBoundedLinearFormAndHessianObjective(
    Eigen::SparseMatrix<double, Eigen::ColMajor, long long>& H, Eigen::SparseMatrix<double, Eigen::ColMajor, long long>& A,
    bool include_finite_bounds, const Eigen::VectorXi* col_nnz_H, const Eigen::VectorXi* col_nnz_A, bool upper_part_only_H)
{
    H.setZero();
    // TODO(roesmann): what is most efficient, if we just want to udpate the Hessian?
    if (col_nnz_H) H.reserve(*col_nnz_H);

    A.setZero();  // we need to reset since we call insert below
    if (col_nnz_A) A.reserve(*col_nnz_A);

    int inequality_row_start = getEqualityDimension();
    int bounds_row_start     = inequality_row_start + getInequalityDimension();

    OptimizationEdgeSet::Ptr edges = _graph.getEdgeSet();

    // Iterate plain objective edges
    for (BaseEdge::Ptr& edge : edges->getObjectiveEdgesRef())
    {
        for (int row_i = 0; row_i < edge->getNumVertices(); ++row_i)
        {
            const VertexInterface* vertex1 = edge->getVertexRaw(row_i);

            int vert_i_dim_unfixed = vertex1->getDimensionUnfixed();
            if (vert_i_dim_unfixed == 0) continue;

            // compute jacobian to speed up hessian computation
            Eigen::MatrixXd block_jacobian1(edge->getDimension(), vertex1->getDimensionUnfixed());
            edge->computeJacobian(row_i, block_jacobian1);

            if (edge->isLinear()) continue;

            int col_start = upper_part_only_H ? row_i : 0;
            for (int col_j = col_start; col_j < edge->getNumVertices(); ++col_j)  // hessian is symmetric
            {
                const VertexInterface* vertex2 = edge->getVertexRaw(col_j);

                int vert_j_dim_unfixed = vertex2->getDimensionUnfixed();
                if (vert_j_dim_unfixed == 0) continue;

                // warning: we require the hessian to be initialized with zeros
                // and that edge->computeHessian only adds (+=) values!!
                Eigen::MatrixXd block_hessian(vert_i_dim_unfixed, vert_j_dim_unfixed);
                edge->computeHessian(row_i, col_j, block_jacobian1, block_hessian);

                // insert values
                if (upper_part_only_H && vertex1 == vertex2)
                {
                    // block on the main diagonal, so use upper part only
                    for (int i = 0; i < block_hessian.rows(); ++i)
                    {
                        for (int j = i; j < block_hessian.cols(); ++j)
                        {
                            H.coeffRef(vertex1->getVertexIdx() + i, vertex2->getVertexIdx() + j) += block_hessian.coeffRef(i, j);
                        }
                    }
                }
                else
                {
                    for (int i = 0; i < block_hessian.rows(); ++i)
                    {
                        for (int j = 0; j < block_hessian.cols(); ++j)
                        {
                            H.coeffRef(vertex1->getVertexIdx() + i, vertex2->getVertexIdx() + j) += block_hessian.coeffRef(i, j);
                        }
                    }
                }
            }
        }
    }

    // Iterate lsq objective edges
    for (BaseEdge::Ptr& edge : edges->getLsqObjectiveEdgesRef())
    {
        // if (edge->isLinear()) continue; // we must not check this one, since Linear-squared is not linear ;-)

        for (int row_i = 0; row_i < edge->getNumVertices(); ++row_i)
        {
            const VertexInterface* vertex1 = edge->getVertexRaw(row_i);

            int vert_i_dim_unfixed = vertex1->getDimensionUnfixed();
            if (vert_i_dim_unfixed == 0) continue;

            // compute jacobian to speed up hessian computation
            Eigen::MatrixXd block_jacobian1(edge->getDimension(), vertex1->getDimensionUnfixed());
            edge->computeJacobian(row_i, block_jacobian1);

            if (edge->isLinear()) continue;

            int col_start = upper_part_only_H ? row_i : 0;
            for (int col_j = col_start; col_j < edge->getNumVertices(); ++col_j)  // hessian is symmetric
            {
                const VertexInterface* vertex2 = edge->getVertexRaw(col_j);

                int vert_j_dim_unfixed = vertex2->getDimensionUnfixed();
                if (vert_j_dim_unfixed == 0) continue;

                // we need the second jacobian for applying the chain rule
                Eigen::MatrixXd block_jacobian2(edge->getDimension(), vertex2->getDimensionUnfixed());
                edge->computeJacobian(col_j, block_jacobian2);

                // approximate H = 2*J^T*J
                Eigen::MatrixXd block_hessian(vert_i_dim_unfixed, vert_j_dim_unfixed);
                block_hessian = 2.0 * block_jacobian1.transpose() * block_jacobian2;

                // insert values
                if (upper_part_only_H && vertex1 == vertex2)
                {
                    // block on the main diagonal, so use upper part only
                    for (int i = 0; i < block_hessian.rows(); ++i)
                    {
                        for (int j = i; j < block_hessian.cols(); ++j)
                        {
                            H.coeffRef(vertex1->getVertexIdx() + i, vertex2->getVertexIdx() + j) += block_hessian.coeffRef(i, j);
                        }
                    }
                }
                else
                {
                    for (int i = 0; i < block_hessian.rows(); ++i)
                    {
                        for (int j = 0; j < block_hessian.cols(); ++j)
                        {
                            H.coeffRef(vertex1->getVertexIdx() + i, vertex2->getVertexIdx() + j) += block_hessian.coeffRef(i, j);
                        }
                    }
                }
            }
        }
    }

    // Iterate equality edges
    for (BaseEdge::Ptr& edge : edges->getEqualityEdgesRef())
    {
        for (int row_i = 0; row_i < edge->getNumVertices(); ++row_i)
        {
            const VertexInterface* vertex1 = edge->getVertexRaw(row_i);

            int vert_i_dim_unfixed = vertex1->getDimensionUnfixed();
            if (vert_i_dim_unfixed == 0) continue;

            // compute jacobian to speed up hessian computation
            Eigen::MatrixXd block_jacobian1(edge->getDimension(), vertex1->getDimensionUnfixed());
            edge->computeJacobian(row_i, block_jacobian1);

            // Now insert values to Jacobian matrix A
            for (int i = 0; i < block_jacobian1.rows(); ++i)
            {
                for (int j = 0; j < block_jacobian1.cols(); ++j)
                {
                    A.insert(edge->getEdgeIdx() + i, vertex1->getVertexIdx() + j) = block_jacobian1.coeffRef(i, j);
                }
            }
        }
    }

    // Iterate inequality edges
    for (BaseEdge::Ptr& edge : edges->getInequalityEdgesRef())
    {
        int A_edge_idx = inequality_row_start + edge->getEdgeIdx();

        for (int row_i = 0; row_i < edge->getNumVertices(); ++row_i)
        {
            const VertexInterface* vertex1 = edge->getVertexRaw(row_i);

            int vert_i_dim_unfixed = vertex1->getDimensionUnfixed();
            if (vert_i_dim_unfixed == 0) continue;

            // compute jacobian to speed up hessian computation
            Eigen::MatrixXd block_jacobian1(edge->getDimension(), vertex1->getDimensionUnfixed());
            edge->computeJacobian(row_i, block_jacobian1);

            // Now insert values to Jacobian A
            for (int i = 0; i < block_jacobian1.rows(); ++i)
            {
                for (int j = 0; j < block_jacobian1.cols(); ++j)
                {
                    A.insert(A_edge_idx + i, vertex1->getVertexIdx() + j) = block_jacobian1.coeffRef(i, j);
                }
            }
        }
    }

    // Iterate mixed edges
    for (BaseMixedEdge::Ptr& edge : edges->getMixedEdgesRef())
    {
        // if (edge->getObjectiveDimension() == 0 || edge->isObjectiveLinear()) continue;
        int A_edge_eq_idx   = edge->getEdgeEqualityIdx();
        int A_edge_ineq_idx = inequality_row_start + edge->getEdgeInequalityIdx();

        for (int row_i = 0; row_i < edge->getNumVertices(); ++row_i)
        {
            const VertexInterface* vertex1 = edge->getVertexRaw(row_i);

            int vert_i_dim_unfixed = vertex1->getDimensionUnfixed();
            if (vert_i_dim_unfixed == 0) continue;

            // compute jacobian to speed up hessian computation
            Eigen::MatrixXd obj_block_jacobian1(edge->getObjectiveDimension(), vertex1->getDimensionUnfixed());
            Eigen::MatrixXd eq_block_jacobian1(edge->getEqualityDimension(), vertex1->getDimensionUnfixed());
            Eigen::MatrixXd ineq_block_jacobian1(edge->getInequalityDimension(), vertex1->getDimensionUnfixed());
            edge->computeJacobians(row_i, obj_block_jacobian1, eq_block_jacobian1, ineq_block_jacobian1);

            // Now insert values to Jacobian A
            for (int i = 0; i < eq_block_jacobian1.rows(); ++i)
            {
                for (int j = 0; j < eq_block_jacobian1.cols(); ++j)
                {
                    A.insert(A_edge_eq_idx + i, vertex1->getVertexIdx() + j) = eq_block_jacobian1.coeffRef(i, j);
                }
            }
            for (int i = 0; i < ineq_block_jacobian1.rows(); ++i)
            {
                for (int j = 0; j < ineq_block_jacobian1.cols(); ++j)
                {
                    A.insert(A_edge_ineq_idx + i, vertex1->getVertexIdx() + j) = ineq_block_jacobian1.coeffRef(i, j);
                }
            }

            if (edge->isObjectiveLinear() || edge->getObjectiveDimension() == 0) continue;

            int col_start = upper_part_only_H ? row_i : 0;
            for (int col_j = col_start; col_j < edge->getNumVertices(); ++col_j)  // hessian is symmetric
            {
                const VertexInterface* vertex2 = edge->getVertexRaw(col_j);

                int vert_j_dim_unfixed = vertex2->getDimensionUnfixed();
                if (vert_j_dim_unfixed == 0) continue;

                Eigen::MatrixXd block_hessian(vert_i_dim_unfixed, vert_j_dim_unfixed);
                block_hessian.setZero();

                if (edge->isObjectiveLeastSquaresForm())
                {
                    if (!edge->isObjectiveLinear() && edge->getObjectiveDimension() != 0)
                    {
                        // we need the second jacobian for applying the chain rule
                        Eigen::MatrixXd obj_block_jacobian2(edge->getObjectiveDimension(), vertex2->getDimensionUnfixed());
                        edge->computeObjectiveJacobian(col_j, obj_block_jacobian2);

                        // approximate H = 2*J^T*J
                        block_hessian += 2.0 * obj_block_jacobian1.transpose() * obj_block_jacobian2;
                    }
                }
                else  // TODO(roesmann): check all reasonable cases:
                {
                    edge->computeObjectiveHessianInc(row_i, col_j, obj_block_jacobian1, block_hessian, nullptr, 1.0);
                }

                // insert values
                if (upper_part_only_H && vertex1 == vertex2)
                {
                    // block on the main diagonal, so use upper part only
                    for (int i = 0; i < block_hessian.rows(); ++i)
                    {
                        for (int j = i; j < block_hessian.cols(); ++j)
                        {
                            H.coeffRef(vertex1->getVertexIdx() + i, vertex2->getVertexIdx() + j) += block_hessian.coeffRef(i, j);
                        }
                    }
                }
                else
                {
                    for (int i = 0; i < block_hessian.rows(); ++i)
                    {
                        for (int j = 0; j < block_hessian.cols(); ++j)
                        {
                            H.coeffRef(vertex1->getVertexIdx() + i, vertex2->getVertexIdx() + j) += block_hessian.coeffRef(i, j);
                        }
                    }
                }
            }
        }
    }

    if (include_finite_bounds)
    {
        int row_idx = bounds_row_start;
        for (const VertexInterface* vertex : _graph.getVertexSet()->getActiveVertices())
        {
            int vert_idx = vertex->getVertexIdx();
            int free_idx = 0;
            for (int i = 0; i < vertex->getDimension(); ++i)
            {
                if (vertex->isFixedComponent(i)) continue;

                if (vertex->hasFiniteLowerBound(i) || vertex->hasFiniteUpperBound(i))
                {
                    A.insert(row_idx, vert_idx + free_idx) = 1;
                    ++row_idx;
                }
                ++free_idx;
            }
        }
    }
}

}  // namespace corbo
