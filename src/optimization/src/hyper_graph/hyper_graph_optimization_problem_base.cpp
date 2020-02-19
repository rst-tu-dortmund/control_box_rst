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

#include <corbo-optimization/hyper_graph/hyper_graph_optimization_problem_base.h>

namespace corbo {

void BaseHyperGraphOptimizationProblem::precomputeVertexQuantities()
{
    assert(_graph.hasVertexSet());

    // Vertex set
    VertexSetInterface::Ptr vertices = _graph.getVertexSet();
    if (vertices->isModified())
    {
        vertices->computeActiveVertices();
        vertices->setModified(false);  // avoid recomputation of active vertices
        _dim_param = vertices->getParameterDimension();
        vertices->computeVertexIndices();
    }
}

void BaseHyperGraphOptimizationProblem::precomputeEdgeQuantities()
{
    assert(_graph.hasEdgeSet());

    OptimizationEdgeSet::Ptr edges = _graph.getEdgeSet();
    if (edges->isModified())
    {
        edges->getDimensions(_dim_non_lsq_obj, _dim_lsq_obj, _dim_eq, _dim_ineq);
        edges->computeEdgeIndices();
        edges->setModified(false);
    }
}

void BaseHyperGraphOptimizationProblem::precomputeGraphQuantities()
{
    assert(_graph.hasEdgeSet());
    assert(_graph.hasVertexSet());

    precomputeVertexQuantities();
    precomputeEdgeQuantities();

    _graph_precomputed = true;
}

void BaseHyperGraphOptimizationProblem::clear()
{
    if (_graph.hasEdgeSet()) _graph.getEdgeSet()->clear();
    _graph_precomputed = false;

    _dim_param       = 0;
    _dim_non_lsq_obj = 0;
    _dim_lsq_obj     = 0;
    _dim_eq          = 0;
    _dim_ineq        = 0;
}

// implmements interface method
double BaseHyperGraphOptimizationProblem::computeValueNonLsqObjective()
{
    assert(_graph.hasEdgeSet());
    if (!_graph_precomputed) precomputeGraphQuantities();

    double value = 0;

    OptimizationEdgeSet::Ptr edges = _graph.getEdgeSet();
    // Iterate plain objective edges
    for (BaseEdge::Ptr& edge : edges->getObjectiveEdgesRef())
    {
        value += edge->computeSumOfValues();
    }
    // Iterate mixed edges
    for (BaseMixedEdge::Ptr& edge : edges->getMixedEdgesRef())
    {
        // we consider only non-least squares objectives here
        if (edge->isObjectiveLeastSquaresForm() || edge->getObjectiveDimension() == 0) continue;

        edge->precompute();
        value += edge->computeSumOfObjectiveValues();
    }
    return value;
}

void BaseHyperGraphOptimizationProblem::computeValuesLsqObjective(Eigen::Ref<Eigen::VectorXd> values)
{
    assert(_graph.hasEdgeSet());
    if (!_graph_precomputed) precomputeGraphQuantities();
    OptimizationEdgeSet::Ptr edges = _graph.getEdgeSet();
    // Iterate plain objective edges
    for (BaseEdge::Ptr& edge : edges->getLsqObjectiveEdgesRef())
    {
        edge->computeValues(values.segment(edge->getEdgeIdx(), edge->getDimension()));
    }
    // Iterate mixed edges
    for (BaseMixedEdge::Ptr& edge : edges->getMixedEdgesRef())
    {
        // we consider only least squares objectives here
        if (!edge->isObjectiveLeastSquaresForm() || edge->getObjectiveDimension() == 0) continue;

        edge->precompute();
        edge->computeObjectiveValues(values.segment(edge->getEdgeObjectiveIdx(), edge->getObjectiveDimension()));
    }
}

double BaseHyperGraphOptimizationProblem::computeValueObjective()
{
    assert(_graph.hasEdgeSet());
    if (!_graph_precomputed) precomputeGraphQuantities();

    double value = 0;

    OptimizationEdgeSet::Ptr edges = _graph.getEdgeSet();
    // Iterate plain objective edges
    for (BaseEdge::Ptr& edge : edges->getObjectiveEdgesRef())
    {
        PRINT_DEBUG_COND_ONCE(edge->isLeastSquaresForm(), "BaseHyperGraphOptimizationProblem::computeValueObjective(): "
                                                              << "least-squares edge found in non-lsq container");
        value += edge->computeSumOfValues();
    }
    for (BaseEdge::Ptr& edge : edges->getLsqObjectiveEdgesRef())
    {
        PRINT_DEBUG_COND_ONCE(!edge->isLeastSquaresForm(), "BaseHyperGraphOptimizationProblem::computeValueObjective(): "
                                                               << "non-least-squares edge found in lsq container");
        value += edge->computeSquaredNormOfValues();
    }
    // Iterate mixed edges
    for (BaseMixedEdge::Ptr& edge : edges->getMixedEdgesRef())
    {
        if (edge->getObjectiveDimension() == 0) continue;

        edge->precompute();
        if (edge->isObjectiveLeastSquaresForm())
            value += edge->computeSquaredNormOfObjectiveValues();
        else
            value += edge->computeSumOfObjectiveValues();
    }
    return value;
}

void BaseHyperGraphOptimizationProblem::computeValuesEquality(Eigen::Ref<Eigen::VectorXd> values)
{
    assert(_graph.hasEdgeSet());
    if (!_graph_precomputed) precomputeGraphQuantities();
    OptimizationEdgeSet::Ptr edges = _graph.getEdgeSet();
    // Iterate plain equality edges
    for (BaseEdge::Ptr& edge : edges->getEqualityEdgesRef())
    {
        edge->computeValues(values.segment(edge->getEdgeIdx(), edge->getDimension()));
    }
    // Iterate mixed edges
    for (BaseMixedEdge::Ptr& edge : edges->getMixedEdgesRef())
    {
        if (edge->getEqualityDimension() == 0) continue;

        edge->precompute();
        edge->computeEqualityValues(values.segment(edge->getEdgeEqualityIdx(), edge->getEqualityDimension()));
    }
}

void BaseHyperGraphOptimizationProblem::computeValuesInequality(Eigen::Ref<Eigen::VectorXd> values)
{
    assert(_graph.hasEdgeSet());
    if (!_graph_precomputed) precomputeGraphQuantities();
    OptimizationEdgeSet::Ptr edges = _graph.getEdgeSet();
    // Iterate plain inequality edges
    for (BaseEdge::Ptr& edge : edges->getInequalityEdgesRef())
    {
        edge->computeValues(values.segment(edge->getEdgeIdx(), edge->getDimension()));
    }
    // Iterate mixed edges
    for (BaseMixedEdge::Ptr& edge : edges->getMixedEdgesRef())
    {
        if (edge->getInequalityDimension() == 0) continue;

        edge->precompute();
        edge->computeInequalityValues(values.segment(edge->getEdgeInequalityIdx(), edge->getInequalityDimension()));
    }
}

void BaseHyperGraphOptimizationProblem::computeValues(double& non_lsq_obj_value, Eigen::Ref<Eigen::VectorXd> lsq_obj_values,
                                                      Eigen::Ref<Eigen::VectorXd> eq_values, Eigen::Ref<Eigen::VectorXd> ineq_values)
{
    assert(_graph.hasEdgeSet());
    if (!_graph_precomputed) precomputeGraphQuantities();

    OptimizationEdgeSet::Ptr edges = _graph.getEdgeSet();
    // Iterate plain objective edges
    non_lsq_obj_value = 0;
    for (BaseEdge::Ptr& edge : edges->getObjectiveEdgesRef())
    {
        non_lsq_obj_value += edge->computeSumOfValues();
    }
    // Iterate lsq objective edges
    for (BaseEdge::Ptr& edge : edges->getLsqObjectiveEdgesRef())
    {
        edge->computeValues(lsq_obj_values.segment(edge->getEdgeIdx(), edge->getDimension()));
    }
    // Iterate plain equality edges
    for (BaseEdge::Ptr& edge : edges->getEqualityEdgesRef())
    {
        edge->computeValues(eq_values.segment(edge->getEdgeIdx(), edge->getDimension()));
    }
    // Iterate plain inequality edges
    for (BaseEdge::Ptr& edge : edges->getInequalityEdgesRef())
    {
        edge->computeValues(ineq_values.segment(edge->getEdgeIdx(), edge->getDimension()));
    }
    // Iterate mixed edges
    for (BaseMixedEdge::Ptr& edge : edges->getMixedEdgesRef())
    {
        edge->precompute();

        if (edge->isObjectiveLeastSquaresForm())
        {
            edge->computeObjectiveValues(lsq_obj_values.segment(edge->getEdgeObjectiveIdx(), edge->getObjectiveDimension()));
        }
        else
        {
            non_lsq_obj_value += edge->computeSumOfObjectiveValues();
        }

        edge->computeEqualityValues(eq_values.segment(edge->getEdgeEqualityIdx(), edge->getEqualityDimension()));
        edge->computeInequalityValues(ineq_values.segment(edge->getEdgeInequalityIdx(), edge->getInequalityDimension()));
    }
}

int BaseHyperGraphOptimizationProblem::finiteCombinedBoundsDimension()
{
    assert(_graph.hasVertexSet());
    if (!_graph_precomputed) precomputeVertexQuantities();

    int dim = 0;
    // Iterate vertices
    for (const VertexInterface* vertex : _graph.getVertexSet()->getActiveVertices())
    {
        dim += vertex->getNumberFiniteBounds(true);  // only unfixed elements
    }
    return dim;
}

int BaseHyperGraphOptimizationProblem::finiteBoundsDimension()
{
    assert(_graph.hasVertexSet());
    if (!_graph_precomputed) precomputeVertexQuantities();

    int dim = 0;
    // Iterate vertices
    for (const VertexInterface* vertex : _graph.getVertexSet()->getActiveVertices())
    {
        dim += vertex->getNumberFiniteLowerBounds(true);  // only unfixed elements
        dim += vertex->getNumberFiniteUpperBounds(true);  // only unfixed elements
    }
    return dim;
}

void BaseHyperGraphOptimizationProblem::computeValuesActiveInequality(Eigen::Ref<Eigen::VectorXd> values, double weight)
{
    // we don't have any better strategy than default currently
    computeValuesInequality(values);
    for (int i = 0; i < values.size(); ++i)
    {
        if (values[i] < 0)
            values[i] = 0;
        else
            values[i] *= weight;
    }
}

void BaseHyperGraphOptimizationProblem::computeDistanceFiniteCombinedBounds(Eigen::Ref<Eigen::VectorXd> values)
{
    assert(_graph.hasVertexSet());
    if (!_graph_precomputed) precomputeVertexQuantities();

    int idx = 0;
    for (const VertexInterface* vertex : _graph.getVertexSet()->getActiveVertices())
    {
        for (int i = 0; i < vertex->getDimension(); ++i)
        {
            if (vertex->isFixedComponent(i)) continue;

            if (vertex->hasFiniteLowerBound(i) || vertex->hasFiniteUpperBound(i))
            {
                if (vertex->getData()[i] < vertex->getLowerBounds()[i])
                    values[idx] = vertex->getLowerBounds()[i] - vertex->getData()[i];
                else if (vertex->getData()[i] > vertex->getUpperBounds()[i])
                    values[idx] = vertex->getData()[i] - vertex->getUpperBounds()[i];
                else
                    values[idx] = 0;
                ++idx;
            }
        }
    }
}

void BaseHyperGraphOptimizationProblem::computeLowerAndUpperBoundDiff(Eigen::Ref<Eigen::VectorXd> lb_minus_x, Eigen::Ref<Eigen::VectorXd> ub_minus_x)
{
    assert(_graph.hasVertexSet());
    if (!_graph_precomputed) precomputeVertexQuantities();

    int idx = 0;
    for (const VertexInterface* vertex : _graph.getVertexSet()->getActiveVertices())
    {
        for (int i = 0; i < vertex->getDimension(); ++i)
        {
            if (vertex->isFixedComponent(i)) continue;

            if (vertex->hasFiniteLowerBound(i) || vertex->hasFiniteUpperBound(i))
            {
                lb_minus_x[idx] = vertex->getLowerBounds()[i] - vertex->getData()[i];
                ub_minus_x[idx] = vertex->getUpperBounds()[i] - vertex->getData()[i];
                ++idx;
            }
        }
    }
}

void BaseHyperGraphOptimizationProblem::getParametersAndBoundsFinite(Eigen::Ref<Eigen::VectorXd> lb_finite_bounds,
                                                                     Eigen::Ref<Eigen::VectorXd> ub_finite_bounds,
                                                                     Eigen::Ref<Eigen::VectorXd> x_finite_bounds)
{
    assert(_graph.hasVertexSet());
    if (!_graph_precomputed) precomputeVertexQuantities();

    int idx = 0;
    for (const VertexInterface* vertex : _graph.getVertexSet()->getActiveVertices())
    {
        for (int i = 0; i < vertex->getDimension(); ++i)
        {
            if (vertex->isFixedComponent(i)) continue;

            if (vertex->hasFiniteLowerBound(i) || vertex->hasFiniteUpperBound(i))
            {
                lb_finite_bounds[idx] = vertex->getLowerBounds()[i];
                ub_finite_bounds[idx] = vertex->getUpperBounds()[i];
                x_finite_bounds[idx]  = vertex->getData()[i];
                ++idx;
            }
        }
    }
}

void BaseHyperGraphOptimizationProblem::computeDenseJacobianFiniteCombinedBounds(Eigen::Ref<Eigen::MatrixXd> jacobian, double weight)
{
    assert(_graph.hasVertexSet());
    if (!_graph_precomputed) precomputeVertexQuantities();
    assert(jacobian.rows() == finiteCombinedBoundsDimension());
    assert(jacobian.cols() == getParameterDimension());

    jacobian.setZero();

    int row_idx = 0;
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
                    jacobian(row_idx, vert_idx + free_idx) = -weight;
                }
                else if (vertex->getData()[i] > vertex->getUpperBounds()[i])
                {
                    jacobian(row_idx, vert_idx + free_idx) = weight;
                }
                ++row_idx;
            }
            ++free_idx;
        }
    }
}

int BaseHyperGraphOptimizationProblem::computeSparseJacobianFiniteCombinedBoundsNNZ() { return finiteCombinedBoundsDimension(); }

void BaseHyperGraphOptimizationProblem::computeSparseJacobianFiniteCombinedBoundsStructure(Eigen::Ref<Eigen::VectorXi> i_row,
                                                                                           Eigen::Ref<Eigen::VectorXi> j_col)
{
    assert(i_row.size() == computeSparseJacobianFiniteCombinedBoundsNNZ());
    assert(j_col.size() == i_row.size());

    int row_idx = 0;
    for (const VertexInterface* vertex : _graph.getVertexSet()->getActiveVertices())
    {
        int vert_idx = vertex->getVertexIdx();
        int free_idx = 0;
        for (int i = 0; i < vertex->getDimension(); ++i)
        {
            if (vertex->isFixedComponent(i)) continue;

            if (vertex->hasFiniteLowerBound(i) || vertex->hasFiniteUpperBound(i))
            {
                i_row[row_idx] = row_idx;
                j_col[row_idx] = vert_idx + free_idx;
                ++row_idx;
            }
            ++free_idx;
        }
    }
}

void BaseHyperGraphOptimizationProblem::computeSparseJacobianFiniteCombinedBoundsValues(Eigen::Ref<Eigen::VectorXd> values, double weight)
{
    assert(values.size() == computeSparseJacobianFiniteCombinedBoundsNNZ());

    int row_idx = 0;
    for (const VertexInterface* vertex : _graph.getVertexSet()->getActiveVertices())
    {
        for (int i = 0; i < vertex->getDimension(); ++i)
        {
            if (vertex->isFixedComponent(i)) continue;

            if (vertex->hasFiniteLowerBound(i) || vertex->hasFiniteUpperBound(i))
            {
                if (vertex->getData()[i] < vertex->getLowerBounds()[i])
                {
                    values(row_idx) = -weight;
                }
                else if (vertex->getData()[i] > vertex->getUpperBounds()[i])
                {
                    values(row_idx) = weight;
                }
                else
                {
                    values(row_idx) = 0.0;
                }
                ++row_idx;
            }
        }
    }
}

void BaseHyperGraphOptimizationProblem::computeDenseJacobianFiniteCombinedBoundsIdentity(Eigen::Ref<Eigen::MatrixXd> jacobian)
{
    assert(_graph.hasVertexSet());
    if (!_graph_precomputed) precomputeVertexQuantities();
    assert(jacobian.rows() == finiteCombinedBoundsDimension());
    assert(jacobian.cols() == getParameterDimension());

    jacobian.setZero();

    int row_idx = 0;
    for (const VertexInterface* vertex : _graph.getVertexSet()->getActiveVertices())
    {
        int vert_idx = vertex->getVertexIdx();
        int free_idx = 0;
        for (int i = 0; i < vertex->getDimension(); ++i)
        {
            if (vertex->isFixedComponent(i)) continue;

            if (vertex->hasFiniteLowerBound(i) || vertex->hasFiniteUpperBound(i))
            {
                jacobian(row_idx, vert_idx + free_idx) = 1;
                ++row_idx;
            }
            ++free_idx;
        }
    }
}

bool BaseHyperGraphOptimizationProblem::checkIfAllUnfixedParam(std::function<bool(double, int)> fun)
{
    assert(_graph.hasVertexSet());
    if (!_graph_precomputed) precomputeVertexQuantities();

    int idx = 0;
    for (const VertexInterface* vertex : _graph.getVertexSet()->getActiveVertices())
    {
        for (int i = 0; i < vertex->getDimension(); ++i)
        {
            if (vertex->isFixedComponent(i)) continue;

            if (vertex->hasFiniteLowerBound(i) || vertex->hasFiniteUpperBound(i))
            {
                if (!fun(vertex->getData()[i], idx)) return false;
                ++idx;
            }
        }
    }
    return true;
}

}  // namespace corbo
