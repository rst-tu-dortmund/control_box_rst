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

#include <corbo-optimization/hyper_graph/hyper_graph.h>

#include <corbo-core/console.h>
#include <corbo-core/types.h>

#include <algorithm>

namespace corbo {

// generic lambdas are c++14, hence we create a local function here
template <typename EdgeContainer>
bool check_edge_vertices(const std::vector<VertexInterface*>& vertices, EdgeContainer& edges)
{
    for (int i = 0; i < edges.size(); ++i)
    {
        EdgeInterface* edge = edges[i].get();
        for (int i = 0; i < edge->getNumVertices(); ++i)
        {
            auto it = std::find(vertices.begin(), vertices.end(), edge->getVertexRaw(i));
            if (it == vertices.end()) return false;
        }
    }
    return true;
}

bool HyperGraph::checkGraphConsistency()
{
    if (!hasVertexSet() || !hasEdgeSet()) return false;

    std::vector<VertexInterface*> vertices;
    getVertexSet()->getVertices(vertices);

    if (!check_edge_vertices(vertices, getEdgeSet()->getObjectiveEdgesRef())) return false;
    if (!check_edge_vertices(vertices, getEdgeSet()->getEqualityEdgesRef())) return false;
    if (!check_edge_vertices(vertices, getEdgeSet()->getInequalityEdgesRef())) return false;
    if (!check_edge_vertices(vertices, getEdgeSet()->getMixedEdgesRef())) return false;
    return true;
}

// void HyperGraph2::addVertex(VertexInterface* vtx)
//{
//    _vertices.push_back(vtx);
//    _vertex_idx_precomputed    = false;
//    _finite_bounds_precomputed = false;
//}
// void HyperGraph2::addObjectiveEdge(EdgeInterface::UPtr edge)
//{
//    // register edge at adjacent vertices
//    for (int i = 0; i < edge->getNumVertices(); ++i)
//    {
//        edge->getVertexRaw(i)->registerObjectiveEdge(edge.get());
//    }
//    _edges_objective.push_back(std::move(edge));
//    _edge_obj_idx_precomputed = false;
//}
// void HyperGraph2::addEqualityConstraintEdge(EdgeInterface::UPtr edge)
//{
//    // register edge at adjacent vertices
//    for (int i = 0; i < edge->getNumVertices(); ++i) edge->getVertexRaw(i)->registerEqualityEdge(edge.get());
//    _edges_equalities.push_back(std::move(edge));
//    _edge_eq_idx_precomputed = false;
//}
// void HyperGraph2::addInequalityConstraintEdge(EdgeInterface::UPtr edge)
//{
//    // register edge at adjacent vertices
//    for (int i = 0; i < edge->getNumVertices(); ++i) edge->getVertexRaw(i)->registerInequalityEdge(edge.get());
//    _edges_inequalities.push_back(std::move(edge));
//    _edge_ineq_idx_precomputed = false;
//}

// void HyperGraph2::computeDenseJacobianActiveConstraints(VertexContainer& vertices, EdgeContainer& edges, Eigen::Ref<Eigen::MatrixXd> jacobian,
//                                                        double weight)
//{
//    if (edges.empty() || vertices.empty()) return;

//    constexpr const double delta     = 1e-9;
//    constexpr const double neg2delta = -2 * delta;
//    constexpr const double scalar    = 1.0 / (2 * delta);

//    // just iterate edges
//    for (EdgeInterface::UPtr& edge : edges)
//    {
//        assert(!edge->isLeastSquaresForm() && "computeDenseJacobianActive() currently does not support least-squares forms");

//        // compute values and check which values are active
//        edge->computeValues();
//        Eigen::Array<bool, -1, 1> active = Eigen::Map<Eigen::ArrayXd>(edge->getValuesRaw(), edge->getDimension()) > 0.0;

//        if (!active.any()) continue;  // all parts are inactive -> zero jacobian

//        for (int i = 0; i < edge->getNumVertices(); ++i)
//        {
//            if (edge->getVertexRaw(i)->getDimensionUnfixed() == 0) continue;

//            edge->computeJacobian(i, jacobian.block(edge->getEdgeIdx(), edge->getVertexRaw(i)->getVertexIdx(), edge->getDimension(),
//                                                    edge->getVertexRaw(i)->getDimensionUnfixed()),
//                                  nullptr);

//            // now set values to zero if inactive or multiply with weight otherwise
//            for (int j = 0; j < edge->getDimension(); ++j)
//            {
//                if (active[j])
//                    jacobian.block(edge->getEdgeIdx() + j, edge->getVertexRaw(i)->getVertexIdx(), 1, edge->getVertexRaw(i)->getDimensionUnfixed())
//                    *=
//                        weight;
//                else
//                    jacobian.block(edge->getEdgeIdx() + j, edge->getVertexRaw(i)->getVertexIdx(), 1, edge->getVertexRaw(i)->getDimensionUnfixed())
//                        .setZero();
//            }
//        }
//    }
//}

// void HyperGraph2::computeValuesObjective(Eigen::Ref<Eigen::VectorXd> values)
//{
//    if (_edges_objective.empty()) return;

//    assert(values.size() == objectiveDimension());

//    // we require that edge indices are valid
//    if (!_edge_obj_idx_precomputed) computeObjectiveEdgeIndices();

//    for (EdgeInterface::UPtr& edge : _edges_objective)
//    {
//        // TODO(roesmann) avoid copy by changing edge interface to directly accept external value vector in computeValues
//        edge->computeValues();

//        if (!_retain_lsq_form && edge->isLeastSquaresForm())
//        {
//            values(edge->getEdgeIdx()) = Eigen::Map<const Eigen::VectorXd>(edge->getValues(), edge->getDimension()).squaredNorm();
//        }
//        else
//        {
//            values.segment(edge->getEdgeIdx(), edge->getDimension()) = Eigen::Map<const Eigen::VectorXd>(edge->getValues(), edge->getDimension());
//        }
//    }
//}

// void HyperGraph2::computeValuesEquality(Eigen::Ref<Eigen::VectorXd> values)
//{
//    if (_edges_equalities.empty()) return;

//    assert(values.size() == equalityDimension());

//    // we require that edge indices are valid
//    if (!_edge_eq_idx_precomputed) computeEqualityEdgeIndices();

//    for (EdgeInterface::UPtr& edge : _edges_equalities)
//    {
//        // TODO(roesmann) avoid copy by changing edge interface to directly accept external value vector in computeValues
//        edge->computeValues();

//        if (!_retain_lsq_form && edge->isLeastSquaresForm())
//        {
//            values(edge->getEdgeIdx()) = Eigen::Map<const Eigen::VectorXd>(edge->getValues(), edge->getDimension()).squaredNorm();
//        }
//        else
//        {
//            values.segment(edge->getEdgeIdx(), edge->getDimension()) = Eigen::Map<const Eigen::VectorXd>(edge->getValues(), edge->getDimension());
//        }
//    }
//}

// void HyperGraph2::computeValuesInequality(Eigen::Ref<Eigen::VectorXd> values)
//{
//    if (_edges_inequalities.empty()) return;

//    assert(values.size() == inequalityDimension());

//    // we require that edge indices are valid
//    if (!_edge_ineq_idx_precomputed) computeInequalityEdgeIndices();

//    for (EdgeInterface::UPtr& edge : _edges_inequalities)
//    {
//        // TODO(roesmann) avoid copy by changing edge interface to directly accept external value vector in computeValues
//        edge->computeValues();

//        if (!_retain_lsq_form && edge->isLeastSquaresForm())
//        {
//            values(edge->getEdgeIdx()) = Eigen::Map<const Eigen::VectorXd>(edge->getValues(), edge->getDimension()).squaredNorm();
//        }
//        else
//        {
//            values.segment(edge->getEdgeIdx(), edge->getDimension()) = Eigen::Map<const Eigen::VectorXd>(edge->getValues(), edge->getDimension());
//        }
//    }
//}

// void HyperGraph2::computeValuesActiveInequality(Eigen::Ref<Eigen::VectorXd> values, double weight)
//{
//    if (_edges_inequalities.empty()) return;

//    assert(values.size() == inequalityDimension());

//    // we require that edge indices are valid
//    if (!_edge_ineq_idx_precomputed) computeInequalityEdgeIndices();

//    for (EdgeInterface::UPtr& edge : _edges_inequalities)
//    {
//        // TODO(roesmann) avoid copy by changing edge interface to directly accept external value vector in computeValues
//        edge->computeValues();

//        if (!_retain_lsq_form && edge->isLeastSquaresForm())
//        {
//            values(edge->getEdgeIdx()) = Eigen::Map<const Eigen::VectorXd>(edge->getValues(), edge->getDimension()).squaredNorm();
//            if (values(edge->getEdgeIdx()) < 0)
//                values(edge->getEdgeIdx()) = 0;
//            else
//                values(edge->getEdgeIdx()) *= weight;
//        }
//        else
//        {
//            Eigen::Map<const Eigen::VectorXd> edge_val(edge->getValues(), edge->getDimension());
//            values.segment(edge->getEdgeIdx(), edge->getDimension()) = (edge_val.array() < 0).select(0.0, edge_val * weight);
//        }
//    }
//}

// void HyperGraph2::computeValuesBounds(Eigen::Ref<Eigen::VectorXd> lb, Eigen::Ref<Eigen::VectorXd> ub)
//{
//    if (_vertices.empty()) return;

//    assert(lb.size() == parameterDimension());
//    assert(ub.size() == parameterDimension());

//    if (!_vertex_idx_precomputed) computeVertexIndices();

//    for (const VertexInterface* vertex : _vertices)
//    {
//        if (vertex->getDimensionUnfixed() == vertex->getDimension())
//        {
//            lb.segment(vertex->getVertexIdx(), vertex->getDimension()) = vertex->getLowerBoundsMap();
//            ub.segment(vertex->getVertexIdx(), vertex->getDimension()) = vertex->getUpperBoundsMap();
//        }
//        else
//        {
//            int param_idx = 0;
//            for (int i = 0; i < vertex->getDimension(); ++i)
//            {
//                if (!vertex->isFixedComponent(i))
//                {
//                    lb(vertex->getVertexIdx() + param_idx) = vertex->getLowerBounds()[i];
//                    ub(vertex->getVertexIdx() + param_idx) = vertex->getUpperBounds()[i];
//                    ++param_idx;
//                }
//            }
//        }
//    }
//}

// void HyperGraph2::computeValuesLowerBounds(Eigen::Ref<Eigen::VectorXd> values)
//{
//    if (_vertices.empty()) return;

//    assert(values.size() == parameterDimension());

//    if (!_vertex_idx_precomputed) computeVertexIndices();

//    for (const VertexInterface* vertex : _vertices)
//    {
//        if (vertex->getDimensionUnfixed() == vertex->getDimension())
//        {
//            values.segment(vertex->getVertexIdx(), vertex->getDimension()) = vertex->getLowerBoundsMap();
//        }
//        else
//        {
//            int param_idx = 0;
//            for (int i = 0; i < vertex->getDimension(); ++i)
//            {
//                if (!vertex->isFixedComponent(i))
//                {
//                    values(vertex->getVertexIdx() + param_idx) = vertex->getLowerBounds()[i];
//                    ++param_idx;
//                }
//            }
//        }
//    }
//}

// void HyperGraph2::computeValuesUpperBounds(Eigen::Ref<Eigen::VectorXd> values)
//{
//    if (_vertices.empty()) return;

//    assert(values.size() == parameterDimension());

//    if (!_vertex_idx_precomputed) computeVertexIndices();

//    for (const VertexInterface* vertex : _vertices)
//    {
//        if (vertex->getDimensionUnfixed() == vertex->getDimension())
//        {
//            values.segment(vertex->getVertexIdx(), vertex->getDimension()) = vertex->getUpperBoundsMap();
//        }

//        else
//        {
//            int param_idx = 0;
//            for (int i = 0; i < vertex->getDimension(); ++i)
//            {
//                if (!vertex->isFixedComponent(i))
//                {
//                    values(vertex->getVertexIdx() + param_idx) = vertex->getUpperBounds()[i];
//                    ++param_idx;
//                }
//            }
//        }
//    }
//}

// void HyperGraph2::computeDistanceFiniteCombinedBounds(Eigen::Ref<Eigen::VectorXd> values)
//{
//    assert(values.size() == finiteCombinedBoundsDimension());

//    int idx = 0;
//    for (const VertexInterface* vertex : _vertices)
//    {
//        if (!vertex->hasFiniteBounds()) continue;

//        for (int i = 0; i < vertex->getDimension(); ++i)
//        {
//            if (vertex->isFixedComponent(i)) continue;

//            if (vertex->getLowerBounds()[i] > -CORBO_INF_DBL || vertex->getUpperBounds()[i] < CORBO_INF_DBL)
//            {
//                if (vertex->getData()[i] < vertex->getLowerBounds()[i])
//                    values(idx) = vertex->getLowerBounds()[i] - vertex->getData()[i];
//                else if (vertex->getData()[i] > vertex->getUpperBounds()[i])
//                    values(idx) = vertex->getData()[i] - vertex->getUpperBounds()[i];
//                else
//                    values(idx) = 0;

//                ++idx;  // only count finite bounds
//            }
//        }
//    }
//}

// void HyperGraph2::computeDenseJacobianObjective(Eigen::Ref<Eigen::MatrixXd> jacobian, const double* multipliers)
//{
//    // we require that edge indices are valid
//    if (_edges_objective.empty() || _vertices.empty()) return;

//    assert(jacobian.rows() == objectiveDimension());
//    assert(jacobian.cols() == parameterDimension());

//    // we require that edge indices are valid
//    if (!_edge_obj_idx_precomputed) computeObjectiveEdgeIndices();
//    // if custom jacobians are specified, we also require valid vertex indices
//    if ((_num_custom_jacobians_obj > 0 || _edge_based) && !_vertex_idx_precomputed) computeVertexIndices();

//    computeDenseJacobian<EdgeType::Objective>(_vertices, _edges_objective, jacobian, _num_custom_jacobians_obj, _edge_based, multipliers);
//}

// void HyperGraph2::computeDenseJacobianEqualities(Eigen::Ref<Eigen::MatrixXd> jacobian, const double* multipliers)
//{
//    if (_edges_equalities.empty() || _vertices.empty()) return;

//    assert(jacobian.rows() == equalityDimension());
//    assert(jacobian.cols() == parameterDimension());

//    // we require that edge indices are valid
//    if (!_edge_eq_idx_precomputed) computeEqualityEdgeIndices();
//    // if custom jacobians are specified, we also require valid vertex indices
//    if ((_num_custom_jacobians_eq > 0 || _edge_based) && !_vertex_idx_precomputed) computeVertexIndices();

//    computeDenseJacobian<EdgeType::Equality>(_vertices, _edges_equalities, jacobian, _num_custom_jacobians_eq, _edge_based, multipliers);
//}

// void HyperGraph2::computeDenseJacobianInequalities(Eigen::Ref<Eigen::MatrixXd> jacobian, const double* multipliers)
//{
//    if (_edges_inequalities.empty() || _vertices.empty()) return;

//    assert(jacobian.rows() == inequalityDimension());
//    assert(jacobian.cols() == parameterDimension());

//    // we require that edge indices are valid
//    if (!_edge_ineq_idx_precomputed) computeInequalityEdgeIndices();
//    // if custom jacobians are specified, we also require valid vertex indices
//    if ((_num_custom_jacobians_ineq > 0 || _edge_based) && !_vertex_idx_precomputed) computeVertexIndices();

//    computeDenseJacobian<EdgeType::Inequality>(_vertices, _edges_inequalities, jacobian, _num_custom_jacobians_ineq, _edge_based, multipliers);
//}

//// only active constraints
// void HyperGraph2::computeDenseJacobianActiveInequalities(Eigen::Ref<Eigen::MatrixXd> jacobian, double weight)
//{
//    if (_edges_inequalities.empty() || _vertices.empty()) return;

//    assert(jacobian.rows() == inequalityDimension());
//    assert(jacobian.cols() == parameterDimension());

//    // we require that edge indices are valid
//    if (!_edge_ineq_idx_precomputed) computeInequalityEdgeIndices();
//    // if custom jacobians are specified, we also require valid vertex indices
//    if (!_vertex_idx_precomputed) computeVertexIndices();

//    computeDenseJacobianActiveConstraints(_vertices, _edges_inequalities, jacobian, weight);
//}

// void HyperGraph2::computeDenseJacobianFiniteCombinedBounds(Eigen::Ref<Eigen::MatrixXd> jacobian, double weight)
//{
//    if (_vertices.empty()) return;

//    assert(jacobian.rows() == finiteCombinedBoundsDimension());
//    assert(jacobian.cols() == parameterDimension());

//    // if custom jacobians are specified, we also require valid vertex indices
//    if (!_vertex_idx_precomputed) computeVertexIndices();

//    int idx = 0;
//    for (const VertexInterface* vertex : _vertices)
//    {
//        if (!vertex->hasFiniteBounds()) continue;

//        for (int i = 0; i < vertex->getDimension(); ++i)
//        {
//            if (vertex->isFixedComponent(i)) continue;

//            if (vertex->getLowerBounds()[i] > -CORBO_INF_DBL || vertex->getUpperBounds()[i] < CORBO_INF_DBL)
//            {
//                if (vertex->getData()[i] < vertex->getLowerBounds()[i])
//                    jacobian(idx, vertex->getVertexIdx() + i) = -weight;  // -1 * weight
//                else if (vertex->getData()[i] > vertex->getUpperBounds()[i])
//                    jacobian(idx, vertex->getVertexIdx() + i) = weight;  // +1 * weight
//                // else
//                //    jacobian(idx, vertex->vertexIdx() + i) = 0.0;   // jacobian is already initialized with zeros!

//                ++idx;  // only count finite bounds
//            }
//        }
//    }
//}

// void HyperGraph2::computeDenseHessian(VertexContainer& vertices, EdgeContainer& edges, Eigen::Ref<Eigen::MatrixXd> hessian,
//                                      const Eigen::Ref<const Eigen::MatrixXd> jacobian, const double* multipliers, bool jacob_scaled)
//{
//    if (edges.empty() || vertices.empty()) return;

//    for (EdgeInterface::UPtr& edge : edges)
//    {
//        // Hessian for linear equations is zero!
//        if (edge->isLinear()) continue;

//        assert((!edge->isLeastSquaresForm() || _retain_lsq_form) &&
//               "Hessian computation for automatic least-square form resolution not yet implemented.");

//        for (int i = 0; i < edge->getNumVertices(); ++i)
//        {
//            VertexInterface* vertex_i = edge->getVertexRaw(i);

//            if (vertex_i->getDimensionUnfixed() == 0) continue;

//            for (int j = i; j < edge->getNumVertices(); ++j)
//            {
//                VertexInterface* vertex_j = edge->getVertexRaw(j);

//                if (vertex_j->getDimensionUnfixed() == 0) continue;

//                if (multipliers)
//                {
//                    edge->computeHessian(
//                        i, j, jacobian.block(edge->getEdgeIdx(), vertex_i->getVertexIdx(), edge->getDimension(), vertex_i->getDimensionUnfixed()),
//                        hessian.block(vertex_i->getVertexIdx(), vertex_j->getVertexIdx(), vertex_i->getDimensionUnfixed(),
//                                      vertex_j->getDimensionUnfixed()),
//                        multipliers + edge->getEdgeIdx(), jacob_scaled);
//                }
//                else
//                {
//                    edge->computeHessian(
//                        i, j, jacobian.block(edge->getEdgeIdx(), vertex_i->getVertexIdx(), edge->getDimension(), vertex_i->getDimensionUnfixed()),
//                        hessian.block(vertex_i->getVertexIdx(), vertex_j->getVertexIdx(), vertex_i->getDimensionUnfixed(),
//                                      vertex_j->getDimensionUnfixed()),
//                        nullptr, false);
//                }

//                if (i != j)
//                {
//                    // assume symmetry:
//                    hessian
//                        .block(vertex_j->getVertexIdx(), vertex_i->getVertexIdx(), vertex_j->getDimensionUnfixed(), vertex_i->getDimensionUnfixed())
//                        .noalias() += hessian
//                                          .block(vertex_i->getVertexIdx(), vertex_j->getVertexIdx(), vertex_i->getDimensionUnfixed(),
//                                                 vertex_j->getDimensionUnfixed())
//                                          .transpose();
//                }
//            }
//        }
//    }
//}

// void HyperGraph2::computeDenseHessianObjective(const Eigen::Ref<const Eigen::MatrixXd>& jacobian, Eigen::Ref<Eigen::MatrixXd> hessian,
//                                               const double* multiplier, bool jacob_scaled)
//{
//    if (_edges_objective.empty() || _vertices.empty()) return;

//    assert(hessian.rows() == parameterDimension());
//    assert(hessian.cols() == parameterDimension());

//    // we require that edge indices are valid
//    if (!_edge_obj_idx_precomputed) computeObjectiveEdgeIndices();
//    // if custom jacobians are specified, we also require valid vertex indices
//    // if ((_num_custom_jacobians_obj > 0 || edge_based) && !_vertex_idx_precomputed) computeVertexIndices();
//    if (!_vertex_idx_precomputed) computeVertexIndices();

//    computeDenseHessian(_vertices, _edges_objective, hessian, jacobian, multiplier, jacob_scaled);
//}

// void HyperGraph2::computeEdgeIndices(EdgeContainer& edges, int* dim, int* num_custom_jacobians)
//{
//    int aux_num_custom_jacobians = 0;
//    if (edges.empty())
//    {
//        if (dim) *dim                                   = 0;
//        if (num_custom_jacobians) *num_custom_jacobians = aux_num_custom_jacobians;
//        return;
//    }

//    edges.front()->_edge_jacobi_idx = 0;
//    if (edges.front()->providesJacobian()) ++aux_num_custom_jacobians;  // useful byproduct
//    for (int i = 1; i < (int)edges.size(); ++i)
//    {
//        if (!_retain_lsq_form && edges[i]->isLeastSquaresForm())
//        {
//            edges[i]->_edge_jacobi_idx = edges[i - 1]->_edge_jacobi_idx + 1;  // resulting dimension is always 1 -> f(x)^T * f(x)
//        }
//        else
//        {
//            edges[i]->_edge_jacobi_idx = edges[i - 1]->_edge_jacobi_idx + edges[i - 1]->getDimension();
//        }
//        if (edges[i]->providesJacobian()) ++aux_num_custom_jacobians;  // useful byproduct
//    }

//    if (dim) *dim = edges.back()->_edge_jacobi_idx + (!_retain_lsq_form && edges.back()->isLeastSquaresForm() ? 1 : edges.back()->getDimension());
//    if (num_custom_jacobians) *num_custom_jacobians = aux_num_custom_jacobians;
//}

// void HyperGraph2::computeVertexIndices()
//{
//    if (_vertices.empty()) return;

//    _vertices.front()->_vertex_idx = 0;
//    for (int i = 1; i < (int)_vertices.size(); ++i)
//    {
//        _vertices[i]->_vertex_idx = _vertices[i - 1]->_vertex_idx + _vertices[i - 1]->getDimensionUnfixed();
//    }

//    _vertex_idx_precomputed = true;
//}

// void HyperGraph2::computeFiniteBoundDimensions()
//{
//    _dim_finite_bounds          = 0;
//    _dim_finite_combined_bounds = 0;
//    for (const VertexInterface* vertex : _vertices)
//    {
//        if (!vertex->hasFiniteBounds()) continue;
//        if (vertex->getDimensionUnfixed() == 0) continue;

//        for (int i = 0; i < vertex->getDimension(); ++i)
//        {
//            if (vertex->isFixedComponent(i)) continue;

//            if (vertex->getLowerBounds()[i] > -CORBO_INF_DBL)
//            {
//                ++_dim_finite_bounds;
//                ++_dim_finite_combined_bounds;
//                if (vertex->getUpperBounds()[i] < CORBO_INF_DBL) ++_dim_finite_bounds;
//            }
//            else if (vertex->getUpperBounds()[i] < CORBO_INF_DBL)
//            {
//                ++_dim_finite_bounds;
//                ++_dim_finite_combined_bounds;
//            }
//        }
//    }
//    _finite_bounds_precomputed = true;
//}

// int HyperGraph2::finiteBoundsDimension()
//{
//    if (_finite_bounds_precomputed) return _dim_finite_bounds;
//    computeFiniteBoundDimensions();
//    return _dim_finite_bounds;
//}

// int HyperGraph2::finiteCombinedBoundsDimension()
//{
//    if (_finite_bounds_precomputed) return _dim_finite_combined_bounds;
//    computeFiniteBoundDimensions();
//    return _dim_finite_combined_bounds;
//}

// bool HyperGraph2::isLeastSquaresProblem() const
//{
//    for (const EdgeInterface::UPtr& edge : _edges_objective)
//    {
//        if (!edge->isLeastSquaresForm()) return false;
//    }
//    return _retain_lsq_form;
//}

// void HyperGraph2::backupParameters()
//{
//    for (VertexInterface* vertex : _vertices) vertex->push();
//}

// void HyperGraph2::restoreBackupParameters(bool keep_backup)
//{
//    if (keep_backup)
//        for (VertexInterface* vertex : _vertices) vertex->top();
//    else
//        for (VertexInterface* vertex : _vertices) vertex->pop();
//}

// void HyperGraph2::discardBackupParameters()
//{
//    for (VertexInterface* vertex : _vertices) vertex->discardTop();
//}

// void HyperGraph2::applyIncrement(const Eigen::Ref<const Eigen::VectorXd>& increment)
//{
//    assert(increment.size() == parameterDimension());

//    if (!_vertex_idx_precomputed) computeVertexIndices();

//    for (VertexInterface* vertex : _vertices)
//    {
//        if (vertex->getDimensionUnfixed() != 0) vertex->plusUnfixed(increment.segment(vertex->getVertexIdx(),
//        vertex->getDimensionUnfixed()).data());
//    }
//}

// double HyperGraph2::getLowerBound(int idx)
//{
//    if (_vertices.empty()) return -CORBO_INF_DBL;

//    if (!_vertex_idx_precomputed) computeVertexIndices();

//    for (const VertexInterface* vertex : _vertices)
//    {
//        int vtx_idx = vertex->getVertexIdx();
//        if (vertex->getDimensionUnfixed() == vertex->getDimension())
//        {
//            if (vtx_idx + vertex->getDimension() > idx)
//            {
//                return vertex->getLowerBounds()[idx - vtx_idx];
//            }
//        }
//        else
//        {
//            int param_idx = 0;
//            for (int i = 0; i < vertex->getDimension(); ++i)
//            {
//                if (!vertex->isFixedComponent(i))
//                {
//                    if (vtx_idx + param_idx == idx) return vertex->getLowerBounds()[i];
//                    ++param_idx;
//                }
//            }
//        }
//    }
//    return -CORBO_INF_DBL;
//}

// double HyperGraph2::getUpperBound(int idx)
//{
//    if (_vertices.empty()) return CORBO_INF_DBL;

//    if (!_vertex_idx_precomputed) computeVertexIndices();

//    for (const VertexInterface* vertex : _vertices)
//    {
//        int vtx_idx = vertex->getVertexIdx();
//        if (vertex->getDimensionUnfixed() == vertex->getDimension())
//        {
//            if (vtx_idx + vertex->getDimension() > idx)
//            {
//                return vertex->getUpperBounds()[idx - vtx_idx];
//            }
//        }
//        else
//        {
//            int param_idx = 0;
//            for (int i = 0; i < vertex->getDimension(); ++i)
//            {
//                if (!vertex->isFixedComponent(i))
//                {
//                    if (vtx_idx + param_idx == idx) return vertex->getUpperBounds()[i];
//                    ++param_idx;
//                }
//            }
//        }
//    }
//    return CORBO_INF_DBL;
//}
// void HyperGraph2::setLowerBound(int idx, double lb)
//{
//    if (_vertices.empty()) return;

//    if (!_vertex_idx_precomputed) computeVertexIndices();

//    for (VertexInterface* vertex : _vertices)
//    {
//        int vtx_idx = vertex->getVertexIdx();
//        if (vertex->getDimensionUnfixed() == vertex->getDimension())
//        {
//            if (vtx_idx + vertex->getDimension() > idx)
//            {
//                vertex->setLowerBound(idx - vtx_idx, lb);
//                return;
//            }
//        }
//        else
//        {
//            int param_idx = 0;
//            for (int i = 0; i < vertex->getDimension(); ++i)
//            {
//                if (!vertex->isFixedComponent(i))
//                {
//                    if (vtx_idx + param_idx == idx)
//                    {
//                        vertex->setLowerBound(i, lb);
//                        return;
//                    }
//                    ++param_idx;
//                }
//            }
//        }
//    }
//}

// void HyperGraph2::setUpperBound(int idx, double ub)
//{
//    if (_vertices.empty()) return;

//    if (!_vertex_idx_precomputed) computeVertexIndices();

//    for (VertexInterface* vertex : _vertices)
//    {
//        int vtx_idx = vertex->getVertexIdx();
//        if (vertex->getDimensionUnfixed() == vertex->getDimension())
//        {
//            if (vtx_idx + vertex->getDimension() > idx)
//            {
//                vertex->setUpperBound(idx - vtx_idx, ub);
//                return;
//            }
//        }
//        else
//        {
//            int param_idx = 0;
//            for (int i = 0; i < vertex->getDimension(); ++i)
//            {
//                if (!vertex->isFixedComponent(i))
//                {
//                    if (vtx_idx + param_idx == idx)
//                    {
//                        vertex->setUpperBound(i, ub);
//                        return;
//                    }
//                    ++param_idx;
//                }
//            }
//        }
//    }
//}

// double HyperGraph2::getParameterValue(int idx)
//{
//    if (_vertices.empty()) return CORBO_MAX_DBL;

//    if (!_vertex_idx_precomputed) computeVertexIndices();

//    for (const VertexInterface* vertex : _vertices)
//    {
//        int vtx_idx = vertex->getVertexIdx();
//        if (vertex->getDimensionUnfixed() == vertex->getDimension())
//        {
//            if (vtx_idx + vertex->getDimension() > idx)
//            {
//                return vertex->getData()[idx - vtx_idx];
//            }
//        }
//        else
//        {
//            int param_idx = 0;
//            for (int i = 0; i < vertex->getDimension(); ++i)
//            {
//                if (!vertex->isFixedComponent(i))
//                {
//                    if (vtx_idx + param_idx == idx) return vertex->getData()[i];
//                    ++param_idx;
//                }
//            }
//        }
//    }
//    return CORBO_MAX_DBL;
//}

// void HyperGraph2::setParameterValue(int idx, double x)
//{
//    if (_vertices.empty()) return;

//    if (!_vertex_idx_precomputed) computeVertexIndices();

//    for (VertexInterface* vertex : _vertices)
//    {
//        int vtx_idx = vertex->getVertexIdx();
//        if (vertex->getDimensionUnfixed() == vertex->getDimension())
//        {
//            if (vtx_idx + vertex->getDimension() > idx)
//            {
//                vertex->setData(idx - vtx_idx, x);
//                return;
//            }
//        }
//        else
//        {
//            int param_idx = 0;
//            for (int i = 0; i < vertex->getDimension(); ++i)
//            {
//                if (!vertex->isFixedComponent(i))
//                {
//                    if (vtx_idx + param_idx == idx)
//                    {
//                        vertex->setData(i, x);
//                        return;
//                    }
//                    ++param_idx;
//                }
//            }
//        }
//    }
//}

// void HyperGraph2::resetStructure()
//{
//    _edge_obj_idx_precomputed  = false;
//    _edge_eq_idx_precomputed   = false;
//    _edge_ineq_idx_precomputed = false;
//    _vertex_idx_precomputed    = false;
//}

// void HyperGraph2::clearConnectedEdges()
//{
//    for (VertexInterface* vertex : _vertices) vertex->clearConnectedEdges();
//}

// void HyperGraph2::clear()
//{
//    clearConnectedEdges();
//    _edges_objective.clear();
//    _edges_equalities.clear();
//    _edges_inequalities.clear();
//    _vertices.clear();
//    resetStructure();
//}

}  // namespace corbo
