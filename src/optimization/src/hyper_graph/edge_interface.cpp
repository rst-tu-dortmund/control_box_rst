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

#include <corbo-optimization/hyper_graph/edge_interface.h>
#include <corbo-optimization/hyper_graph/vector_vertex.h>

#include <corbo-core/console.h>

namespace corbo {

constexpr const double HESSIAN_DELTA = 1e-2;

int EdgeInterface::getNumFiniteVerticesLowerBounds() const
{
    int num_lb = 0;
    for (int i = 0; i < getNumVertices(); ++i)
    {
        assert(getVertex(i));
        num_lb += getVertex(i)->getNumberFiniteLowerBounds(true);
    }
    return num_lb;
}
int EdgeInterface::getNumFiniteVerticesUpperBounds() const
{
    int num_ub = 0;
    for (int i = 0; i < getNumVertices(); ++i)
    {
        assert(getVertex(i));
        num_ub += getVertex(i)->getNumberFiniteUpperBounds(true);
    }
    return num_ub;
}

void BaseEdge::computeJacobian(int vtx_idx, Eigen::Ref<Eigen::MatrixXd> block_jacobian, const double* multipliers)
{
    if (getNumVertices() == 0) return;

    assert(vtx_idx >= 0 && vtx_idx < getNumVertices());
    assert(block_jacobian.rows() == getDimension());
    // TODO(roesmann) fixed vertices
    assert(block_jacobian.cols() == getVertexRaw(vtx_idx)->getDimensionUnfixed());

    constexpr const double delta     = 1e-9;
    constexpr const double neg2delta = -2 * delta;
    constexpr const double scalar    = 1.0 / (2 * delta);

    VertexInterface* vertex = getVertexRaw(vtx_idx);

    Eigen::VectorXd values1(getDimension());
    Eigen::VectorXd values2(getDimension());

    int col_idx = 0;
    for (int i = 0; i < vertex->getDimension(); ++i)
    {
        if (vertex->isFixedComponent(i)) continue;

        vertex->plus(i, delta);  // increment current value by delta
        computeValues(values2);

        vertex->plus(i, neg2delta);  // subtract 2*delta
        computeValues(values1);
        block_jacobian.col(col_idx) = scalar * (values2 - values1);

        vertex->plus(i, delta);  // revert offset

        // increase column index
        ++col_idx;
    }

    if (multipliers)
    {
        // Eigen::Map<const Eigen::VectorXd> mult_vec(multipliers, getDimension());
        for (int i = 0; i < getDimension(); ++i) block_jacobian.row(i) *= multipliers[i];
    }
}

void BaseEdge::computeHessianInc(int vtx_idx_i, int vtx_idx_j, Eigen::Ref<Eigen::MatrixXd> block_hessian_ij, const double* multipliers, double weight)
{
    if (getNumVertices() == 0) return;

    assert(vtx_idx_i >= 0 && vtx_idx_i < getNumVertices());
    assert(vtx_idx_j >= 0 && vtx_idx_j < getNumVertices());
    assert(block_hessian_ij.rows() == getVertexRaw(vtx_idx_i)->getDimensionUnfixed());
    assert(block_hessian_ij.cols() == getVertexRaw(vtx_idx_j)->getDimensionUnfixed());

    if (isLinear()) return;  // we do not need to set to zero, since we assume that it is already initialized with zeros!

    // parameter
    // TODO(roesmann): improve hessian computation: delta=1e-2 is not much,
    //                 but the results are numerically not stable
    constexpr const double delta = HESSIAN_DELTA;

    double scalar = 1.0 / delta;
    if (weight != 1.0) scalar *= weight;

    VertexInterface* vertex_i = getVertexRaw(vtx_idx_i);
    VertexInterface* vertex_j = getVertexRaw(vtx_idx_j);

    Eigen::MatrixXd jacobian1(getDimension(), vertex_i->getDimensionUnfixed());
    Eigen::MatrixXd jacobian2(getDimension(), vertex_i->getDimensionUnfixed());

    computeJacobian(vtx_idx_i, jacobian1, multipliers);

    int param_idx_j = 0;
    for (int j = 0; j < vertex_j->getDimension(); ++j)
    {
        if (vertex_j->isFixedComponent(j)) continue;

        vertex_j->plus(j, delta);  // increment current value by delta

        computeJacobian(vtx_idx_i, jacobian2);  // always use jacobian of vertex1

        if (multipliers)
        {
            for (int val_idx = 0; val_idx < getDimension(); ++val_idx)
                block_hessian_ij.col(param_idx_j) += scalar * multipliers[val_idx] * (jacobian2.row(val_idx) - jacobian1.row(val_idx)).transpose();
        }
        else
        {
            for (int val_idx = 0; val_idx < getDimension(); ++val_idx)
                block_hessian_ij.col(param_idx_j) += scalar * (jacobian2.row(val_idx) - jacobian1.row(val_idx)).transpose();
        }
        vertex_j->plus(j, -delta);  // revert offset

        // increase vertex value index
        ++param_idx_j;
    }
}

void BaseEdge::computeHessianInc(int vtx_idx_i, int vtx_idx_j, const Eigen::Ref<const Eigen::MatrixXd>& block_jacobian_i,
                                 Eigen::Ref<Eigen::MatrixXd> block_hessian_ij, const double* multipliers, double weight)
{
    if (getNumVertices() == 0) return;

    assert(vtx_idx_i >= 0 && vtx_idx_i < getNumVertices());
    assert(vtx_idx_j >= 0 && vtx_idx_j < getNumVertices());
    assert(block_jacobian_i.rows() == getDimension());
    assert(block_jacobian_i.cols() == getVertexRaw(vtx_idx_i)->getDimensionUnfixed());
    assert(block_hessian_ij.rows() == getVertexRaw(vtx_idx_i)->getDimensionUnfixed());
    assert(block_hessian_ij.cols() == getVertexRaw(vtx_idx_j)->getDimensionUnfixed());

    if (isLinear()) return;  // we do not need to set to zero, since we assume that it is already initialized with zeros!

    // parameter
    // TODO(roesmann): improve hessian computation: delta=1e-2 is not much,
    //                 but the results are numerically not stable
    constexpr const double delta = HESSIAN_DELTA;

    double scalar = 1.0 / delta;
    if (weight != 1.0) scalar *= weight;

    VertexInterface* vertex_i = getVertexRaw(vtx_idx_i);
    VertexInterface* vertex_j = getVertexRaw(vtx_idx_j);

    Eigen::MatrixXd jacobian2(getDimension(), vertex_i->getDimensionUnfixed());

    int param_idx_j = 0;
    for (int j = 0; j < vertex_j->getDimension(); ++j)
    {
        if (vertex_j->isFixedComponent(j)) continue;

        vertex_j->plus(j, delta);  // increment current value by delta

        computeJacobian(vtx_idx_i, jacobian2);  // always use jacobian of vertex1

        if (multipliers)
        {
            for (int val_idx = 0; val_idx < getDimension(); ++val_idx)
                block_hessian_ij.col(param_idx_j) +=
                    scalar * multipliers[val_idx] * (jacobian2.row(val_idx) - block_jacobian_i.row(val_idx)).transpose();
        }
        else
        {
            for (int val_idx = 0; val_idx < getDimension(); ++val_idx)
                block_hessian_ij.col(param_idx_j) += scalar * (jacobian2.row(val_idx) - block_jacobian_i.row(val_idx)).transpose();
        }
        vertex_j->plus(j, -delta);  // revert offset

        // increase vertex value index
        ++param_idx_j;
    }
}

void BaseEdge::computeHessian(int vtx_idx_i, int vtx_idx_j, const Eigen::Ref<const Eigen::MatrixXd>& block_jacobian_i,
                              Eigen::Ref<Eigen::MatrixXd> block_hessian_ij, const double* multipliers, double weight)
{
    if (getNumVertices() == 0) return;

    assert(vtx_idx_i >= 0 && vtx_idx_i < getNumVertices());
    assert(vtx_idx_j >= 0 && vtx_idx_j < getNumVertices());
    assert(block_jacobian_i.rows() == getDimension());
    assert(block_jacobian_i.cols() == getVertexRaw(vtx_idx_i)->getDimensionUnfixed());
    assert(block_hessian_ij.rows() == getVertexRaw(vtx_idx_i)->getDimensionUnfixed());
    assert(block_hessian_ij.cols() == getVertexRaw(vtx_idx_j)->getDimensionUnfixed());

    if (isLinear()) return;  // we do not need to set to zero, since we assume that it is already initialized with zeros!

    // parameter
    // TODO(roesmann): improve hessian computation: delta=1e-2 is not much,
    //                 but the results are numerically not stable
    constexpr const double delta = HESSIAN_DELTA;

    double scalar = 1.0 / delta;
    if (weight != 1.0) scalar *= weight;

    VertexInterface* vertex_i = getVertexRaw(vtx_idx_i);
    VertexInterface* vertex_j = getVertexRaw(vtx_idx_j);

    Eigen::MatrixXd jacobian2(getDimension(), vertex_i->getDimensionUnfixed());

    int param_idx_j = 0;
    for (int j = 0; j < vertex_j->getDimension(); ++j)
    {
        if (vertex_j->isFixedComponent(j)) continue;

        vertex_j->plus(j, delta);  // increment current value by delta

        computeJacobian(vtx_idx_i, jacobian2);  // always use jacobian of vertex1

        if (multipliers)
        {
            block_hessian_ij.col(param_idx_j) = scalar * multipliers[0] * (jacobian2.row(0) - block_jacobian_i.row(0)).transpose();
            for (int val_idx = 1; val_idx < getDimension(); ++val_idx)
                block_hessian_ij.col(param_idx_j) +=
                    scalar * multipliers[val_idx] * (jacobian2.row(val_idx) - block_jacobian_i.row(val_idx)).transpose();
        }
        else
        {
            block_hessian_ij.col(param_idx_j) = scalar * (jacobian2.row(0) - block_jacobian_i.row(0)).transpose();
            for (int val_idx = 1; val_idx < getDimension(); ++val_idx)
                block_hessian_ij.col(param_idx_j) += scalar * (jacobian2.row(val_idx) - block_jacobian_i.row(val_idx)).transpose();
        }
        vertex_j->plus(j, -delta);  // revert offset

        // increase vertex value index
        ++param_idx_j;
    }
}

void BaseMixedEdge::computeObjectiveJacobian(int vtx_idx, Eigen::Ref<Eigen::MatrixXd> block_jacobian, const double* multipliers)
{
    if (getNumVertices() == 0 || getObjectiveDimension() < 1) return;

    assert(vtx_idx >= 0 && vtx_idx < getNumVertices());
    assert(block_jacobian.rows() == getObjectiveDimension());
    // TODO(roesmann) fixed vertices
    assert(block_jacobian.cols() == getVertexRaw(vtx_idx)->getDimensionUnfixed());

    constexpr const double delta     = 1e-9;
    constexpr const double neg2delta = -2 * delta;
    constexpr const double scalar    = 1.0 / (2 * delta);

    VertexInterface* vertex = getVertexRaw(vtx_idx);

    Eigen::VectorXd values1(getObjectiveDimension());
    Eigen::VectorXd values2(getObjectiveDimension());

    int col_idx = 0;
    for (int i = 0; i < vertex->getDimension(); ++i)
    {
        if (vertex->isFixedComponent(i)) continue;

        vertex->plus(i, delta);  // increment current value by delta
        precompute();
        computeObjectiveValues(values2);

        vertex->plus(i, neg2delta);  // subtract 2*delta
        precompute();
        computeObjectiveValues(values1);
        block_jacobian.col(col_idx) = scalar * (values2 - values1);

        vertex->plus(i, delta);  // revert offset

        // increase column index
        ++col_idx;
    }

    if (multipliers)
    {
        // Eigen::Map<const Eigen::VectorXd> mult_vec(multipliers, getObjectiveDimension());
        for (int i = 0; i < getObjectiveDimension(); ++i) block_jacobian.row(i) *= multipliers[i];
    }
}
void BaseMixedEdge::computeEqualityJacobian(int vtx_idx, Eigen::Ref<Eigen::MatrixXd> block_jacobian, const double* multipliers)
{
    if (getNumVertices() == 0 || getEqualityDimension() < 1) return;

    assert(vtx_idx >= 0 && vtx_idx < getNumVertices());
    assert(block_jacobian.rows() == getEqualityDimension());
    // TODO(roesmann) fixed vertices
    assert(block_jacobian.cols() == getVertexRaw(vtx_idx)->getDimensionUnfixed());

    constexpr const double delta     = 1e-9;
    constexpr const double neg2delta = -2 * delta;
    constexpr const double scalar    = 1.0 / (2 * delta);

    VertexInterface* vertex = getVertexRaw(vtx_idx);

    Eigen::VectorXd values1(getEqualityDimension());
    Eigen::VectorXd values2(getEqualityDimension());

    int col_idx = 0;
    for (int i = 0; i < vertex->getDimension(); ++i)
    {
        if (vertex->isFixedComponent(i)) continue;

        vertex->plus(i, delta);  // increment current value by delta
        precompute();
        computeEqualityValues(values2);

        vertex->plus(i, neg2delta);  // subtract 2*delta
        precompute();
        computeEqualityValues(values1);
        block_jacobian.col(col_idx) = scalar * (values2 - values1);

        vertex->plus(i, delta);  // revert offset

        // increase column index
        ++col_idx;
    }

    if (multipliers)
    {
        // Eigen::Map<const Eigen::VectorXd> mult_vec(multipliers, getEqualityDimension());
        for (int i = 0; i < getEqualityDimension(); ++i) block_jacobian.row(i) *= multipliers[i];
    }
}
void BaseMixedEdge::computeInequalityJacobian(int vtx_idx, Eigen::Ref<Eigen::MatrixXd> block_jacobian, const double* multipliers)
{
    if (getNumVertices() == 0 || getInequalityDimension() < 1) return;

    assert(vtx_idx >= 0 && vtx_idx < getNumVertices());
    assert(block_jacobian.rows() == getInequalityDimension());
    // TODO(roesmann) fixed vertices
    assert(block_jacobian.cols() == getVertexRaw(vtx_idx)->getDimensionUnfixed());

    constexpr const double delta     = 1e-9;
    constexpr const double neg2delta = -2 * delta;
    constexpr const double scalar    = 1.0 / (2 * delta);

    VertexInterface* vertex = getVertexRaw(vtx_idx);

    Eigen::VectorXd values1(getInequalityDimension());
    Eigen::VectorXd values2(getInequalityDimension());

    int col_idx = 0;
    for (int i = 0; i < vertex->getDimension(); ++i)
    {
        if (vertex->isFixedComponent(i)) continue;

        vertex->plus(i, delta);  // increment current value by delta
        precompute();
        computeInequalityValues(values2);

        vertex->plus(i, neg2delta);  // subtract 2*delta
        precompute();
        computeInequalityValues(values1);
        block_jacobian.col(col_idx) = scalar * (values2 - values1);

        vertex->plus(i, delta);  // revert offset

        // increase column index
        ++col_idx;
    }

    if (multipliers)
    {
        // Eigen::Map<const Eigen::VectorXd> mult_vec(multipliers, getInequalityDimension());
        for (int i = 0; i < getInequalityDimension(); ++i) block_jacobian.row(i) *= multipliers[i];
    }
}

void BaseMixedEdge::computeJacobians(int vtx_idx, Eigen::Ref<Eigen::MatrixXd> obj_jacobian, Eigen::Ref<Eigen::MatrixXd> eq_jacobian,
                                     Eigen::Ref<Eigen::MatrixXd> ineq_jacobian, const double* obj_multipliers, const double* eq_multipliers,
                                     const double* ineq_multipliers)
{
    if (getNumVertices() == 0) return;

    int obj_dim  = getObjectiveDimension();
    int eq_dim   = getEqualityDimension();
    int ineq_dim = getInequalityDimension();

    assert(vtx_idx >= 0 && vtx_idx < getNumVertices());
    assert(obj_jacobian.rows() == obj_dim);
    assert(obj_jacobian.cols() == getVertexRaw(vtx_idx)->getDimensionUnfixed());
    assert(eq_jacobian.rows() == eq_dim);
    assert(eq_jacobian.cols() == getVertexRaw(vtx_idx)->getDimensionUnfixed());
    assert(ineq_jacobian.rows() == ineq_dim);
    assert(ineq_jacobian.cols() == getVertexRaw(vtx_idx)->getDimensionUnfixed());

    constexpr const double delta     = 1e-9;
    constexpr const double neg2delta = -2 * delta;
    constexpr const double scalar    = 1.0 / (2 * delta);

    VertexInterface* vertex = getVertexRaw(vtx_idx);

    Eigen::VectorXd obj_values1(obj_dim);
    Eigen::VectorXd obj_values2(obj_dim);
    Eigen::VectorXd eq_values1(eq_dim);
    Eigen::VectorXd eq_values2(eq_dim);
    Eigen::VectorXd ineq_values1(ineq_dim);
    Eigen::VectorXd ineq_values2(ineq_dim);

    int col_idx = 0;
    for (int i = 0; i < vertex->getDimension(); ++i)
    {
        if (vertex->isFixedComponent(i)) continue;

        vertex->plus(i, delta);  // increment current value by delta
        precompute();
        if (obj_dim > 0) computeObjectiveValues(obj_values2);
        if (eq_dim > 0) computeEqualityValues(eq_values2);
        if (ineq_dim > 0) computeInequalityValues(ineq_values2);

        vertex->plus(i, neg2delta);  // subtract 2*delta
        precompute();
        if (obj_dim > 0) computeObjectiveValues(obj_values1);
        if (eq_dim > 0) computeEqualityValues(eq_values1);
        if (ineq_dim > 0) computeInequalityValues(ineq_values1);

        if (obj_dim > 0) obj_jacobian.col(col_idx)   = scalar * (obj_values2 - obj_values1);
        if (eq_dim > 0) eq_jacobian.col(col_idx)     = scalar * (eq_values2 - eq_values1);
        if (ineq_dim > 0) ineq_jacobian.col(col_idx) = scalar * (ineq_values2 - ineq_values1);

        vertex->plus(i, delta);  // revert offset

        // increase column index
        ++col_idx;
    }

    if (obj_multipliers && obj_dim > 0)
    {
        for (int i = 0; i < obj_dim; ++i) obj_jacobian.row(i) *= obj_multipliers[i];
    }
    if (eq_multipliers && eq_dim > 0)
    {
        for (int i = 0; i < eq_dim; ++i) eq_jacobian.row(i) *= eq_multipliers[i];
    }
    if (ineq_multipliers && ineq_dim > 0)
    {
        for (int i = 0; i < ineq_dim; ++i) ineq_jacobian.row(i) *= ineq_multipliers[i];
    }
}

void BaseMixedEdge::computeConstraintJacobians(int vtx_idx, Eigen::Ref<Eigen::MatrixXd> eq_jacobian, Eigen::Ref<Eigen::MatrixXd> ineq_jacobian,
                                               const double* eq_multipliers, const double* ineq_multipliers)
{
    if (getNumVertices() == 0) return;

    int eq_dim   = getEqualityDimension();
    int ineq_dim = getInequalityDimension();

    assert(vtx_idx >= 0 && vtx_idx < getNumVertices());
    assert(eq_jacobian.rows() == eq_dim);
    assert(eq_jacobian.cols() == getVertexRaw(vtx_idx)->getDimensionUnfixed());
    assert(ineq_jacobian.rows() == ineq_dim);
    assert(ineq_jacobian.cols() == getVertexRaw(vtx_idx)->getDimensionUnfixed());

    constexpr const double delta     = 1e-9;
    constexpr const double neg2delta = -2 * delta;
    constexpr const double scalar    = 1.0 / (2 * delta);

    VertexInterface* vertex = getVertexRaw(vtx_idx);

    Eigen::VectorXd eq_values1(eq_dim);
    Eigen::VectorXd eq_values2(eq_dim);
    Eigen::VectorXd ineq_values1(ineq_dim);
    Eigen::VectorXd ineq_values2(ineq_dim);

    int col_idx = 0;
    for (int i = 0; i < vertex->getDimension(); ++i)
    {
        if (vertex->isFixedComponent(i)) continue;

        vertex->plus(i, delta);  // increment current value by delta
        precompute();
        if (eq_dim > 0) computeEqualityValues(eq_values2);
        if (ineq_dim > 0) computeInequalityValues(ineq_values2);

        vertex->plus(i, neg2delta);  // subtract 2*delta
        precompute();
        if (eq_dim > 0) computeEqualityValues(eq_values1);
        if (ineq_dim > 0) computeInequalityValues(ineq_values1);

        if (eq_dim > 0) eq_jacobian.col(col_idx)     = scalar * (eq_values2 - eq_values1);
        if (ineq_dim > 0) ineq_jacobian.col(col_idx) = scalar * (ineq_values2 - ineq_values1);

        vertex->plus(i, delta);  // revert offset

        // increase column index
        ++col_idx;
    }

    if (eq_multipliers && eq_dim > 0)
    {
        for (int i = 0; i < eq_dim; ++i) eq_jacobian.row(i) *= eq_multipliers[i];
    }
    if (ineq_multipliers && ineq_dim > 0)
    {
        for (int i = 0; i < ineq_dim; ++i) ineq_jacobian.row(i) *= ineq_multipliers[i];
    }
}

void BaseMixedEdge::computeObjectiveHessian(int vtx_idx_i, int vtx_idx_j, const Eigen::Ref<const Eigen::MatrixXd>& block_jacobian_i,
                                            Eigen::Ref<Eigen::MatrixXd> block_hessian_ij, const double* multipliers, double weight)
{
    if (getNumVertices() == 0) return;

    assert(vtx_idx_i >= 0 && vtx_idx_i < getNumVertices());
    assert(vtx_idx_j >= 0 && vtx_idx_j < getNumVertices());
    assert(block_jacobian_i.rows() == getObjectiveDimension());
    assert(block_jacobian_i.cols() == getVertexRaw(vtx_idx_i)->getDimensionUnfixed());
    assert(block_hessian_ij.rows() == getVertexRaw(vtx_idx_i)->getDimensionUnfixed());
    assert(block_hessian_ij.cols() == getVertexRaw(vtx_idx_j)->getDimensionUnfixed());

    if (isObjectiveLinear()) return;  // we do not need to set to zero, since we assume that it is already initialized with zeros!

    // parameter
    // TODO(roesmann): improve hessian computation: delta=1e-2 is not much,
    //                 but the results are numerically not stable
    constexpr const double delta = HESSIAN_DELTA;

    double scalar = 1.0 / delta;
    if (weight != 1.0) scalar *= weight;

    VertexInterface* vertex_i = getVertexRaw(vtx_idx_i);
    VertexInterface* vertex_j = getVertexRaw(vtx_idx_j);

    Eigen::MatrixXd jacobian2(getObjectiveDimension(), vertex_i->getDimensionUnfixed());

    int param_idx_j = 0;
    for (int j = 0; j < vertex_j->getDimension(); ++j)
    {
        if (vertex_j->isFixedComponent(j)) continue;

        vertex_j->plus(j, delta);  // increment current value by delta

        computeObjectiveJacobian(vtx_idx_i, jacobian2);  // always use jacobian of vertex1

        if (multipliers)
        {
            block_hessian_ij.col(param_idx_j) = scalar * multipliers[0] * (jacobian2.row(0) - block_jacobian_i.row(0)).transpose();
            for (int val_idx = 1; val_idx < getObjectiveDimension(); ++val_idx)
                block_hessian_ij.col(param_idx_j) +=
                    scalar * multipliers[val_idx] * (jacobian2.row(val_idx) - block_jacobian_i.row(val_idx)).transpose();
        }
        else
        {
            block_hessian_ij.col(param_idx_j) = scalar * (jacobian2.row(0) - block_jacobian_i.row(0)).transpose();
            for (int val_idx = 1; val_idx < getObjectiveDimension(); ++val_idx)
                block_hessian_ij.col(param_idx_j) += scalar * (jacobian2.row(val_idx) - block_jacobian_i.row(val_idx)).transpose();
        }
        vertex_j->plus(j, -delta);  // revert offset

        // increase vertex value index
        ++param_idx_j;
    }
}
void BaseMixedEdge::computeEqualityHessian(int vtx_idx_i, int vtx_idx_j, const Eigen::Ref<const Eigen::MatrixXd>& block_jacobian_i,
                                           Eigen::Ref<Eigen::MatrixXd> block_hessian_ij, const double* multipliers, double weight)
{
    if (getNumVertices() == 0) return;

    assert(vtx_idx_i >= 0 && vtx_idx_i < getNumVertices());
    assert(vtx_idx_j >= 0 && vtx_idx_j < getNumVertices());
    assert(block_jacobian_i.rows() == getEqualityDimension());
    assert(block_jacobian_i.cols() == getVertexRaw(vtx_idx_i)->getDimensionUnfixed());
    assert(block_hessian_ij.rows() == getVertexRaw(vtx_idx_i)->getDimensionUnfixed());
    assert(block_hessian_ij.cols() == getVertexRaw(vtx_idx_j)->getDimensionUnfixed());

    if (isEqualityLinear()) return;  // we do not need to set to zero, since we assume that it is already initialized with zeros!

    // parameter
    // TODO(roesmann): improve hessian computation: delta=1e-2 is not much,
    //                 but the results are numerically not stable
    constexpr const double delta = HESSIAN_DELTA;

    double scalar = 1.0 / delta;
    if (weight != 1.0) scalar *= weight;

    VertexInterface* vertex_i = getVertexRaw(vtx_idx_i);
    VertexInterface* vertex_j = getVertexRaw(vtx_idx_j);

    Eigen::MatrixXd jacobian2(getEqualityDimension(), vertex_i->getDimensionUnfixed());

    int param_idx_j = 0;
    for (int j = 0; j < vertex_j->getDimension(); ++j)
    {
        if (vertex_j->isFixedComponent(j)) continue;

        vertex_j->plus(j, delta);  // increment current value by delta

        computeEqualityJacobian(vtx_idx_i, jacobian2);  // always use jacobian of vertex1

        if (multipliers)
        {
            block_hessian_ij.col(param_idx_j) = scalar * multipliers[0] * (jacobian2.row(0) - block_jacobian_i.row(0)).transpose();
            for (int val_idx = 1; val_idx < getEqualityDimension(); ++val_idx)
                block_hessian_ij.col(param_idx_j) +=
                    scalar * multipliers[val_idx] * (jacobian2.row(val_idx) - block_jacobian_i.row(val_idx)).transpose();
        }
        else
        {
            block_hessian_ij.col(param_idx_j) = scalar * (jacobian2.row(0) - block_jacobian_i.row(0)).transpose();
            for (int val_idx = 1; val_idx < getEqualityDimension(); ++val_idx)
                block_hessian_ij.col(param_idx_j) += scalar * (jacobian2.row(val_idx) - block_jacobian_i.row(val_idx)).transpose();
        }
        vertex_j->plus(j, -delta);  // revert offset

        // increase vertex value index
        ++param_idx_j;
    }
}
void BaseMixedEdge::computeInequalityHessian(int vtx_idx_i, int vtx_idx_j, const Eigen::Ref<const Eigen::MatrixXd>& block_jacobian_i,
                                             Eigen::Ref<Eigen::MatrixXd> block_hessian_ij, const double* multipliers, double weight)
{
    if (getNumVertices() == 0) return;

    assert(vtx_idx_i >= 0 && vtx_idx_i < getNumVertices());
    assert(vtx_idx_j >= 0 && vtx_idx_j < getNumVertices());
    assert(block_jacobian_i.rows() == getInequalityDimension());
    assert(block_jacobian_i.cols() == getVertexRaw(vtx_idx_i)->getDimensionUnfixed());
    assert(block_hessian_ij.rows() == getVertexRaw(vtx_idx_i)->getDimensionUnfixed());
    assert(block_hessian_ij.cols() == getVertexRaw(vtx_idx_j)->getDimensionUnfixed());

    if (isInequalityLinear()) return;  // we do not need to set to zero, since we assume that it is already initialized with zeros!

    // parameter
    // TODO(roesmann): improve hessian computation: delta=1e-2 is not much,
    //                 but the results are numerically not stable
    constexpr const double delta = HESSIAN_DELTA;

    double scalar = 1.0 / delta;
    if (weight != 1.0) scalar *= weight;

    VertexInterface* vertex_i = getVertexRaw(vtx_idx_i);
    VertexInterface* vertex_j = getVertexRaw(vtx_idx_j);

    Eigen::MatrixXd jacobian2(getInequalityDimension(), vertex_i->getDimensionUnfixed());

    int param_idx_j = 0;
    for (int j = 0; j < vertex_j->getDimension(); ++j)
    {
        if (vertex_j->isFixedComponent(j)) continue;

        vertex_j->plus(j, delta);  // increment current value by delta

        computeInequalityJacobian(vtx_idx_i, jacobian2);  // always use jacobian of vertex1

        if (multipliers)
        {
            block_hessian_ij.col(param_idx_j) = scalar * multipliers[0] * (jacobian2.row(0) - block_jacobian_i.row(0)).transpose();
            for (int val_idx = 1; val_idx < getInequalityDimension(); ++val_idx)
                block_hessian_ij.col(param_idx_j) +=
                    scalar * multipliers[val_idx] * (jacobian2.row(val_idx) - block_jacobian_i.row(val_idx)).transpose();
        }
        else
        {
            block_hessian_ij.col(param_idx_j) = scalar * (jacobian2.row(0) - block_jacobian_i.row(0)).transpose();
            for (int val_idx = 1; val_idx < getInequalityDimension(); ++val_idx)
                block_hessian_ij.col(param_idx_j) += scalar * (jacobian2.row(val_idx) - block_jacobian_i.row(val_idx)).transpose();
        }
        vertex_j->plus(j, -delta);  // revert offset

        // increase vertex value index
        ++param_idx_j;
    }
}
void BaseMixedEdge::computeHessians(int vtx_idx_i, int vtx_idx_j, const Eigen::Ref<const Eigen::MatrixXd>& obj_jacobian_i,
                                    const Eigen::Ref<const Eigen::MatrixXd>& eq_jacobian_i, const Eigen::Ref<const Eigen::MatrixXd>& ineq_jacobian_i,
                                    Eigen::Ref<Eigen::MatrixXd> obj_hessian_ij, Eigen::Ref<Eigen::MatrixXd> eq_hessian_ij,
                                    Eigen::Ref<Eigen::MatrixXd> ineq_hessian_ij, const double* multipliers_eq, const double* multipliers_ineq,
                                    double weight_obj)
{
    if (getNumVertices() == 0) return;

    assert(vtx_idx_i >= 0 && vtx_idx_i < getNumVertices());
    assert(vtx_idx_j >= 0 && vtx_idx_j < getNumVertices());
    assert(obj_jacobian_i.rows() == getObjectiveDimension());
    assert(obj_jacobian_i.cols() == getVertexRaw(vtx_idx_i)->getDimensionUnfixed());
    assert(eq_jacobian_i.rows() == getEqualityDimension());
    assert(eq_jacobian_i.cols() == getVertexRaw(vtx_idx_i)->getDimensionUnfixed());
    assert(ineq_jacobian_i.rows() == getInequalityDimension());
    assert(ineq_jacobian_i.cols() == getVertexRaw(vtx_idx_i)->getDimensionUnfixed());
    assert(obj_hessian_ij.rows() == getVertexRaw(vtx_idx_i)->getDimensionUnfixed());
    assert(obj_hessian_ij.cols() == getVertexRaw(vtx_idx_j)->getDimensionUnfixed());
    assert(eq_hessian_ij.rows() == getVertexRaw(vtx_idx_i)->getDimensionUnfixed());
    assert(eq_hessian_ij.cols() == getVertexRaw(vtx_idx_j)->getDimensionUnfixed());
    assert(ineq_hessian_ij.rows() == getVertexRaw(vtx_idx_i)->getDimensionUnfixed());
    assert(ineq_hessian_ij.cols() == getVertexRaw(vtx_idx_j)->getDimensionUnfixed());

    // parameter
    // TODO(roesmann): improve hessian computation: delta=1e-2 is not much,
    //                 but the results are numerically not stable
    constexpr const double delta  = HESSIAN_DELTA;
    constexpr const double scalar = 1.0 / delta;

    VertexInterface* vertex_i = getVertexRaw(vtx_idx_i);
    VertexInterface* vertex_j = getVertexRaw(vtx_idx_j);

    Eigen::MatrixXd obj_jacobian_j(getObjectiveDimension(), vertex_i->getDimensionUnfixed());
    Eigen::MatrixXd eq_jacobian_j(getEqualityDimension(), vertex_i->getDimensionUnfixed());
    Eigen::MatrixXd ineq_jacobian_j(getInequalityDimension(), vertex_i->getDimensionUnfixed());

    int param_idx_j = 0;
    for (int j = 0; j < vertex_j->getDimension(); ++j)
    {
        if (vertex_j->isFixedComponent(j)) continue;

        vertex_j->plus(j, delta);  // increment current value by delta

        computeJacobians(vtx_idx_i, obj_jacobian_j, eq_jacobian_j, ineq_jacobian_j);

        if (!isObjectiveLinear() && getObjectiveDimension() > 0)
        {
            obj_hessian_ij.col(param_idx_j) = weight_obj * scalar * (obj_jacobian_j.row(0) - obj_jacobian_i.row(0)).transpose();
            for (int val_idx = 1; val_idx < getObjectiveDimension(); ++val_idx)
                obj_hessian_ij.col(param_idx_j) += weight_obj * scalar * (obj_jacobian_j.row(val_idx) - obj_jacobian_i.row(val_idx)).transpose();
        }
        if (!isEqualityLinear() && getEqualityDimension() > 0)
        {
            if (multipliers_eq)
            {
                eq_hessian_ij.col(param_idx_j) = scalar * multipliers_eq[0] * (eq_jacobian_j.row(0) - eq_jacobian_i.row(0)).transpose();
                for (int val_idx = 1; val_idx < getEqualityDimension(); ++val_idx)
                    eq_hessian_ij.col(param_idx_j) +=
                        scalar * multipliers_eq[val_idx] * (eq_jacobian_j.row(val_idx) - eq_jacobian_i.row(val_idx)).transpose();
            }
            else
            {
                eq_hessian_ij.col(param_idx_j) = scalar * (eq_jacobian_j.row(0) - eq_jacobian_i.row(0)).transpose();
                for (int val_idx = 1; val_idx < getEqualityDimension(); ++val_idx)
                    eq_hessian_ij.col(param_idx_j) += scalar * (eq_jacobian_j.row(val_idx) - eq_jacobian_i.row(val_idx)).transpose();
            }
        }
        if (!isInequalityLinear() && getInequalityDimension() > 0)
        {
            if (multipliers_ineq)
            {
                ineq_hessian_ij.col(param_idx_j) = scalar * multipliers_ineq[0] * (ineq_jacobian_j.row(0) - ineq_jacobian_i.row(0)).transpose();
                for (int val_idx = 1; val_idx < getInequalityDimension(); ++val_idx)
                    ineq_hessian_ij.col(param_idx_j) +=
                        scalar * multipliers_ineq[val_idx] * (ineq_jacobian_j.row(val_idx) - ineq_jacobian_i.row(val_idx)).transpose();
            }
            else
            {
                ineq_hessian_ij.col(param_idx_j) = scalar * (ineq_jacobian_j.row(0) - ineq_jacobian_i.row(0)).transpose();
                for (int val_idx = 1; val_idx < getInequalityDimension(); ++val_idx)
                    ineq_hessian_ij.col(param_idx_j) += scalar * (ineq_jacobian_j.row(val_idx) - ineq_jacobian_i.row(val_idx)).transpose();
            }
        }

        vertex_j->plus(j, -delta);  // revert offset

        // increase vertex value index
        ++param_idx_j;
    }
}

void BaseMixedEdge::computeConstraintHessians(int vtx_idx_i, int vtx_idx_j, const Eigen::Ref<const Eigen::MatrixXd>& eq_jacobian_i,
                                              const Eigen::Ref<const Eigen::MatrixXd>& ineq_jacobian_i, Eigen::Ref<Eigen::MatrixXd> eq_hessian_ij,
                                              Eigen::Ref<Eigen::MatrixXd> ineq_hessian_ij, const double* multipliers_eq,
                                              const double* multipliers_ineq)
{
    if (getNumVertices() == 0) return;

    assert(vtx_idx_i >= 0 && vtx_idx_i < getNumVertices());
    assert(vtx_idx_j >= 0 && vtx_idx_j < getNumVertices());
    assert(eq_jacobian_i.rows() == getEqualityDimension());
    assert(eq_jacobian_i.cols() == getVertexRaw(vtx_idx_i)->getDimensionUnfixed());
    assert(ineq_jacobian_i.rows() == getInequalityDimension());
    assert(ineq_jacobian_i.cols() == getVertexRaw(vtx_idx_i)->getDimensionUnfixed());
    assert(eq_hessian_ij.rows() == getVertexRaw(vtx_idx_i)->getDimensionUnfixed());
    assert(eq_hessian_ij.cols() == getVertexRaw(vtx_idx_j)->getDimensionUnfixed());
    assert(ineq_hessian_ij.rows() == getVertexRaw(vtx_idx_i)->getDimensionUnfixed());
    assert(ineq_hessian_ij.cols() == getVertexRaw(vtx_idx_j)->getDimensionUnfixed());

    // parameter
    // TODO(roesmann): improve hessian computation: delta=1e-2 is not much,
    //                 but the results are numerically not stable
    constexpr const double delta  = HESSIAN_DELTA;
    constexpr const double scalar = 1.0 / delta;

    VertexInterface* vertex_i = getVertexRaw(vtx_idx_i);
    VertexInterface* vertex_j = getVertexRaw(vtx_idx_j);

    Eigen::MatrixXd eq_jacobian_j(getEqualityDimension(), vertex_i->getDimensionUnfixed());
    Eigen::MatrixXd ineq_jacobian_j(getInequalityDimension(), vertex_i->getDimensionUnfixed());

    int param_idx_j = 0;
    for (int j = 0; j < vertex_j->getDimension(); ++j)
    {
        if (vertex_j->isFixedComponent(j)) continue;

        vertex_j->plus(j, delta);  // increment current value by delta

        computeConstraintJacobians(vtx_idx_i, eq_jacobian_j, ineq_jacobian_j);

        if (!isEqualityLinear() && getEqualityDimension() > 0)
        {
            if (multipliers_eq)
            {
                eq_hessian_ij.col(param_idx_j) = scalar * multipliers_eq[0] * (eq_jacobian_j.row(0) - eq_jacobian_i.row(0)).transpose();
                for (int val_idx = 1; val_idx < getEqualityDimension(); ++val_idx)
                    eq_hessian_ij.col(param_idx_j) +=
                        scalar * multipliers_eq[val_idx] * (eq_jacobian_j.row(val_idx) - eq_jacobian_i.row(val_idx)).transpose();
            }
            else
            {
                eq_hessian_ij.col(param_idx_j) = scalar * (eq_jacobian_j.row(0) - eq_jacobian_i.row(0)).transpose();
                for (int val_idx = 1; val_idx < getEqualityDimension(); ++val_idx)
                    eq_hessian_ij.col(param_idx_j) += scalar * (eq_jacobian_j.row(val_idx) - eq_jacobian_i.row(val_idx)).transpose();
            }
        }
        if (!isInequalityLinear() && getInequalityDimension() > 0)
        {
            if (multipliers_ineq)
            {
                ineq_hessian_ij.col(param_idx_j) = scalar * multipliers_ineq[0] * (ineq_jacobian_j.row(0) - ineq_jacobian_i.row(0)).transpose();
                for (int val_idx = 1; val_idx < getInequalityDimension(); ++val_idx)
                    ineq_hessian_ij.col(param_idx_j) +=
                        scalar * multipliers_ineq[val_idx] * (ineq_jacobian_j.row(val_idx) - ineq_jacobian_i.row(val_idx)).transpose();
            }
            else
            {
                ineq_hessian_ij.col(param_idx_j) = scalar * (ineq_jacobian_j.row(0) - ineq_jacobian_i.row(0)).transpose();
                for (int val_idx = 1; val_idx < getInequalityDimension(); ++val_idx)
                    ineq_hessian_ij.col(param_idx_j) += scalar * (ineq_jacobian_j.row(val_idx) - ineq_jacobian_i.row(val_idx)).transpose();
            }
        }

        vertex_j->plus(j, -delta);  // revert offset

        // increase vertex value index
        ++param_idx_j;
    }
}

void BaseMixedEdge::computeObjectiveHessianInc(int vtx_idx_i, int vtx_idx_j, const Eigen::Ref<const Eigen::MatrixXd>& block_jacobian_i,
                                               Eigen::Ref<Eigen::MatrixXd> block_hessian_ij, const double* multipliers, double weight)
{
    if (getNumVertices() == 0) return;

    assert(vtx_idx_i >= 0 && vtx_idx_i < getNumVertices());
    assert(vtx_idx_j >= 0 && vtx_idx_j < getNumVertices());
    assert(block_jacobian_i.rows() == getObjectiveDimension());
    assert(block_jacobian_i.cols() == getVertexRaw(vtx_idx_i)->getDimensionUnfixed());
    assert(block_hessian_ij.rows() == getVertexRaw(vtx_idx_i)->getDimensionUnfixed());
    assert(block_hessian_ij.cols() == getVertexRaw(vtx_idx_j)->getDimensionUnfixed());

    if (isObjectiveLinear()) return;  // we do not need to set to zero, since we assume that it is already initialized with zeros!

    // parameter
    // TODO(roesmann): improve hessian computation: delta=1e-2 is not much,
    //                 but the results are numerically not stable
    constexpr const double delta = HESSIAN_DELTA;

    double scalar = 1.0 / delta;
    if (weight != 1.0) scalar *= weight;

    VertexInterface* vertex_i = getVertexRaw(vtx_idx_i);
    VertexInterface* vertex_j = getVertexRaw(vtx_idx_j);

    Eigen::MatrixXd jacobian2(getObjectiveDimension(), vertex_i->getDimensionUnfixed());

    int param_idx_j = 0;
    for (int j = 0; j < vertex_j->getDimension(); ++j)
    {
        if (vertex_j->isFixedComponent(j)) continue;

        vertex_j->plus(j, delta);  // increment current value by delta

        computeObjectiveJacobian(vtx_idx_i, jacobian2);  // always use jacobian of vertex1

        if (multipliers)
        {
            for (int val_idx = 0; val_idx < getObjectiveDimension(); ++val_idx)
                block_hessian_ij.col(param_idx_j) +=
                    scalar * multipliers[val_idx] * (jacobian2.row(val_idx) - block_jacobian_i.row(val_idx)).transpose();
        }
        else
        {
            for (int val_idx = 0; val_idx < getObjectiveDimension(); ++val_idx)
                block_hessian_ij.col(param_idx_j) += scalar * (jacobian2.row(val_idx) - block_jacobian_i.row(val_idx)).transpose();
        }
        vertex_j->plus(j, -delta);  // revert offset

        // increase vertex value index
        ++param_idx_j;
    }
}
void BaseMixedEdge::computeEqualityHessianInc(int vtx_idx_i, int vtx_idx_j, const Eigen::Ref<const Eigen::MatrixXd>& block_jacobian_i,
                                              Eigen::Ref<Eigen::MatrixXd> block_hessian_ij, const double* multipliers, double weight)
{
    if (getNumVertices() == 0) return;

    assert(vtx_idx_i >= 0 && vtx_idx_i < getNumVertices());
    assert(vtx_idx_j >= 0 && vtx_idx_j < getNumVertices());
    assert(block_jacobian_i.rows() == getEqualityDimension());
    assert(block_jacobian_i.cols() == getVertexRaw(vtx_idx_i)->getDimensionUnfixed());
    assert(block_hessian_ij.rows() == getVertexRaw(vtx_idx_i)->getDimensionUnfixed());
    assert(block_hessian_ij.cols() == getVertexRaw(vtx_idx_j)->getDimensionUnfixed());

    if (isEqualityLinear()) return;  // we do not need to set to zero, since we assume that it is already initialized with zeros!

    // parameter
    // TODO(roesmann): improve hessian computation: delta=1e-2 is not much,
    //                 but the results are numerically not stable
    constexpr const double delta = HESSIAN_DELTA;

    double scalar = 1.0 / delta;
    if (weight != 1.0) scalar *= weight;

    VertexInterface* vertex_i = getVertexRaw(vtx_idx_i);
    VertexInterface* vertex_j = getVertexRaw(vtx_idx_j);

    Eigen::MatrixXd jacobian2(getEqualityDimension(), vertex_i->getDimensionUnfixed());

    int param_idx_j = 0;
    for (int j = 0; j < vertex_j->getDimension(); ++j)
    {
        if (vertex_j->isFixedComponent(j)) continue;

        vertex_j->plus(j, delta);  // increment current value by delta

        computeEqualityJacobian(vtx_idx_i, jacobian2);  // always use jacobian of vertex1

        if (multipliers)
        {
            for (int val_idx = 0; val_idx < getEqualityDimension(); ++val_idx)
                block_hessian_ij.col(param_idx_j) +=
                    scalar * multipliers[val_idx] * (jacobian2.row(val_idx) - block_jacobian_i.row(val_idx)).transpose();
        }
        else
        {
            for (int val_idx = 0; val_idx < getEqualityDimension(); ++val_idx)
                block_hessian_ij.col(param_idx_j) += scalar * (jacobian2.row(val_idx) - block_jacobian_i.row(val_idx)).transpose();
        }
        vertex_j->plus(j, -delta);  // revert offset

        // increase vertex value index
        ++param_idx_j;
    }
}
void BaseMixedEdge::computeInequalityHessianInc(int vtx_idx_i, int vtx_idx_j, const Eigen::Ref<const Eigen::MatrixXd>& block_jacobian_i,
                                                Eigen::Ref<Eigen::MatrixXd> block_hessian_ij, const double* multipliers, double weight)
{
    if (getNumVertices() == 0) return;

    assert(vtx_idx_i >= 0 && vtx_idx_i < getNumVertices());
    assert(vtx_idx_j >= 0 && vtx_idx_j < getNumVertices());
    assert(block_jacobian_i.rows() == getInequalityDimension());
    assert(block_jacobian_i.cols() == getVertexRaw(vtx_idx_i)->getDimensionUnfixed());
    assert(block_hessian_ij.rows() == getVertexRaw(vtx_idx_i)->getDimensionUnfixed());
    assert(block_hessian_ij.cols() == getVertexRaw(vtx_idx_j)->getDimensionUnfixed());

    if (isInequalityLinear()) return;  // we do not need to set to zero, since we assume that it is already initialized with zeros!

    // parameter
    // TODO(roesmann): improve hessian computation: delta=1e-2 is not much,
    //                 but the results are numerically not stable
    constexpr const double delta = HESSIAN_DELTA;

    double scalar = 1.0 / delta;
    if (weight != 1.0) scalar *= weight;

    VertexInterface* vertex_i = getVertexRaw(vtx_idx_i);
    VertexInterface* vertex_j = getVertexRaw(vtx_idx_j);

    Eigen::MatrixXd jacobian2(getInequalityDimension(), vertex_i->getDimensionUnfixed());

    int param_idx_j = 0;
    for (int j = 0; j < vertex_j->getDimension(); ++j)
    {
        if (vertex_j->isFixedComponent(j)) continue;

        vertex_j->plus(j, delta);  // increment current value by delta

        computeInequalityJacobian(vtx_idx_i, jacobian2);  // always use jacobian of vertex1

        if (multipliers)
        {
            for (int val_idx = 0; val_idx < getInequalityDimension(); ++val_idx)
                block_hessian_ij.col(param_idx_j) +=
                    scalar * multipliers[val_idx] * (jacobian2.row(val_idx) - block_jacobian_i.row(val_idx)).transpose();
        }
        else
        {
            for (int val_idx = 0; val_idx < getInequalityDimension(); ++val_idx)
                block_hessian_ij.col(param_idx_j) += scalar * (jacobian2.row(val_idx) - block_jacobian_i.row(val_idx)).transpose();
        }
        vertex_j->plus(j, -delta);  // revert offset

        // increase vertex value index
        ++param_idx_j;
    }
}
void BaseMixedEdge::computeHessiansInc(int vtx_idx_i, int vtx_idx_j, const Eigen::Ref<const Eigen::MatrixXd>& obj_jacobian_i,
                                       const Eigen::Ref<const Eigen::MatrixXd>& eq_jacobian_i,
                                       const Eigen::Ref<const Eigen::MatrixXd>& ineq_jacobian_i, Eigen::Ref<Eigen::MatrixXd> obj_hessian_ij,
                                       Eigen::Ref<Eigen::MatrixXd> eq_hessian_ij, Eigen::Ref<Eigen::MatrixXd> ineq_hessian_ij,
                                       const double* multipliers_eq, const double* multipliers_ineq, double weight_obj)
{
    if (getNumVertices() == 0) return;

    assert(vtx_idx_i >= 0 && vtx_idx_i < getNumVertices());
    assert(vtx_idx_j >= 0 && vtx_idx_j < getNumVertices());
    assert(obj_jacobian_i.rows() == getObjectiveDimension());
    assert(obj_jacobian_i.cols() == getVertexRaw(vtx_idx_i)->getDimensionUnfixed());
    assert(eq_jacobian_i.rows() == getEqualityDimension());
    assert(eq_jacobian_i.cols() == getVertexRaw(vtx_idx_i)->getDimensionUnfixed());
    assert(ineq_jacobian_i.rows() == getInequalityDimension());
    assert(ineq_jacobian_i.cols() == getVertexRaw(vtx_idx_i)->getDimensionUnfixed());
    assert(obj_hessian_ij.rows() == getVertexRaw(vtx_idx_i)->getDimensionUnfixed());
    assert(obj_hessian_ij.cols() == getVertexRaw(vtx_idx_j)->getDimensionUnfixed());
    assert(eq_hessian_ij.rows() == getVertexRaw(vtx_idx_i)->getDimensionUnfixed());
    assert(eq_hessian_ij.cols() == getVertexRaw(vtx_idx_j)->getDimensionUnfixed());
    assert(ineq_hessian_ij.rows() == getVertexRaw(vtx_idx_i)->getDimensionUnfixed());
    assert(ineq_hessian_ij.cols() == getVertexRaw(vtx_idx_j)->getDimensionUnfixed());

    // parameter
    // TODO(roesmann): improve hessian computation: delta=1e-2 is not much,
    //                 but the results are numerically not stable
    constexpr const double delta  = HESSIAN_DELTA;
    constexpr const double scalar = 1.0 / delta;

    VertexInterface* vertex_i = getVertexRaw(vtx_idx_i);
    VertexInterface* vertex_j = getVertexRaw(vtx_idx_j);

    Eigen::MatrixXd obj_jacobian_j(getObjectiveDimension(), vertex_i->getDimensionUnfixed());
    Eigen::MatrixXd eq_jacobian_j(getEqualityDimension(), vertex_i->getDimensionUnfixed());
    Eigen::MatrixXd ineq_jacobian_j(getInequalityDimension(), vertex_i->getDimensionUnfixed());

    int param_idx_j = 0;
    for (int j = 0; j < vertex_j->getDimension(); ++j)
    {
        if (vertex_j->isFixedComponent(j)) continue;

        vertex_j->plus(j, delta);  // increment current value by delta

        computeJacobians(vtx_idx_i, obj_jacobian_j, eq_jacobian_j, ineq_jacobian_j);

        if (!isObjectiveLinear())
        {
            for (int val_idx = 0; val_idx < getObjectiveDimension(); ++val_idx)
                obj_hessian_ij.col(param_idx_j) += weight_obj * scalar * (obj_jacobian_j.row(val_idx) - obj_jacobian_i.row(val_idx)).transpose();
        }
        if (!isEqualityLinear())
        {
            if (multipliers_eq)
            {
                for (int val_idx = 0; val_idx < getEqualityDimension(); ++val_idx)
                    eq_hessian_ij.col(param_idx_j) +=
                        scalar * multipliers_eq[val_idx] * (eq_jacobian_j.row(val_idx) - eq_jacobian_i.row(val_idx)).transpose();
            }
            else
            {
                for (int val_idx = 0; val_idx < getEqualityDimension(); ++val_idx)
                    eq_hessian_ij.col(param_idx_j) += scalar * (eq_jacobian_j.row(val_idx) - eq_jacobian_i.row(val_idx)).transpose();
            }
        }
        if (!isInequalityLinear())
        {
            if (multipliers_ineq)
            {
                for (int val_idx = 0; val_idx < getInequalityDimension(); ++val_idx)
                    ineq_hessian_ij.col(param_idx_j) +=
                        scalar * multipliers_ineq[val_idx] * (ineq_jacobian_j.row(val_idx) - ineq_jacobian_i.row(val_idx)).transpose();
            }
            else
            {
                for (int val_idx = 0; val_idx < getInequalityDimension(); ++val_idx)
                    ineq_hessian_ij.col(param_idx_j) += scalar * (ineq_jacobian_j.row(val_idx) - ineq_jacobian_i.row(val_idx)).transpose();
            }
        }

        vertex_j->plus(j, -delta);  // revert offset

        // increase vertex value index
        ++param_idx_j;
    }
}

void BaseMixedEdge::computeConstraintHessiansInc(int vtx_idx_i, int vtx_idx_j, const Eigen::Ref<const Eigen::MatrixXd>& eq_jacobian_i,
                                                 const Eigen::Ref<const Eigen::MatrixXd>& ineq_jacobian_i, Eigen::Ref<Eigen::MatrixXd> eq_hessian_ij,
                                                 Eigen::Ref<Eigen::MatrixXd> ineq_hessian_ij, const double* multipliers_eq,
                                                 const double* multipliers_ineq)
{
    if (getNumVertices() == 0) return;

    assert(vtx_idx_i >= 0 && vtx_idx_i < getNumVertices());
    assert(vtx_idx_j >= 0 && vtx_idx_j < getNumVertices());
    assert(eq_jacobian_i.rows() == getEqualityDimension());
    assert(eq_jacobian_i.cols() == getVertexRaw(vtx_idx_i)->getDimensionUnfixed());
    assert(ineq_jacobian_i.rows() == getInequalityDimension());
    assert(ineq_jacobian_i.cols() == getVertexRaw(vtx_idx_i)->getDimensionUnfixed());
    assert(eq_hessian_ij.rows() == getVertexRaw(vtx_idx_i)->getDimensionUnfixed());
    assert(eq_hessian_ij.cols() == getVertexRaw(vtx_idx_j)->getDimensionUnfixed());
    assert(ineq_hessian_ij.rows() == getVertexRaw(vtx_idx_i)->getDimensionUnfixed());
    assert(ineq_hessian_ij.cols() == getVertexRaw(vtx_idx_j)->getDimensionUnfixed());

    // parameter
    // TODO(roesmann): improve hessian computation: delta=1e-2 is not much,
    //                 but the results are numerically not stable
    constexpr const double delta  = HESSIAN_DELTA;
    constexpr const double scalar = 1.0 / delta;

    VertexInterface* vertex_i = getVertexRaw(vtx_idx_i);
    VertexInterface* vertex_j = getVertexRaw(vtx_idx_j);

    Eigen::MatrixXd eq_jacobian_j(getEqualityDimension(), vertex_i->getDimensionUnfixed());
    Eigen::MatrixXd ineq_jacobian_j(getInequalityDimension(), vertex_i->getDimensionUnfixed());

    int param_idx_j = 0;
    for (int j = 0; j < vertex_j->getDimension(); ++j)
    {
        if (vertex_j->isFixedComponent(j)) continue;

        vertex_j->plus(j, delta);  // increment current value by delta

        computeConstraintJacobians(vtx_idx_i, eq_jacobian_j, ineq_jacobian_j);

        if (!isEqualityLinear())
        {
            if (multipliers_eq)
            {
                for (int val_idx = 0; val_idx < getEqualityDimension(); ++val_idx)
                    eq_hessian_ij.col(param_idx_j) +=
                        scalar * multipliers_eq[val_idx] * (eq_jacobian_j.row(val_idx) - eq_jacobian_i.row(val_idx)).transpose();
            }
            else
            {
                for (int val_idx = 0; val_idx < getEqualityDimension(); ++val_idx)
                    eq_hessian_ij.col(param_idx_j) += scalar * (eq_jacobian_j.row(val_idx) - eq_jacobian_i.row(val_idx)).transpose();
            }
        }
        if (!isInequalityLinear())
        {
            if (multipliers_ineq)
            {
                for (int val_idx = 0; val_idx < getInequalityDimension(); ++val_idx)
                    ineq_hessian_ij.col(param_idx_j) +=
                        scalar * multipliers_ineq[val_idx] * (ineq_jacobian_j.row(val_idx) - ineq_jacobian_i.row(val_idx)).transpose();
            }
            else
            {
                for (int val_idx = 0; val_idx < getInequalityDimension(); ++val_idx)
                    ineq_hessian_ij.col(param_idx_j) += scalar * (ineq_jacobian_j.row(val_idx) - ineq_jacobian_i.row(val_idx)).transpose();
            }
        }

        vertex_j->plus(j, -delta);  // revert offset

        // increase vertex value index
        ++param_idx_j;
    }
}

}  // namespace corbo
