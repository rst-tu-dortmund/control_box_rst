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

#include <corbo-core/macros.h>

namespace corbo {

template <HyperGraph2::EdgeType type>
void HyperGraph2::computeDenseJacobian(VertexContainer& vertices, EdgeContainer& edges, Eigen::Ref<Eigen::MatrixXd> jacobian,
                                       int num_custom_jacobians, bool edge_based, const double* multipliers)
{
    if (edges.empty() || vertices.empty()) return;

    constexpr const double delta     = 1e-9;
    constexpr const double neg2delta = -2 * delta;
    constexpr const double scalar    = 1.0 / (2 * delta);

    // first store analytic jacobians if any
    if (num_custom_jacobians > 0 || edge_based)
    {
        for (EdgeInterface::UPtr& edge : edges)
        {
            if (edge->providesJacobian() || edge_based)
            {
                for (int i = 0; i < edge->getNumVertices(); ++i)
                {
                    if (edge->getVertexRaw(i)->getDimensionUnfixed() == 0) continue;

                    // in least squares form, we must compute the jacobian of f(x)^2, hence -> 2 * f(x)^T * Jacobian(f(x)).
                    if (!_retain_lsq_form && edge->isLeastSquaresForm())
                    {
                        Eigen::MatrixXd jacobf(edge->getDimension(), edge->getVertexRaw(i)->getDimensionUnfixed());
                        // we require a temporary matrix since we need to apply the chain rule in order to get a
                        // dimension of 1
                        edge->computeJacobian(i, jacobf, nullptr);
                        edge->computeValues();
                        const Eigen::Map<const Eigen::VectorXd> values(edge->getValues(), edge->getDimension());
                        if (multipliers)
                        {
                            jacobian.block(edge->getEdgeIdx(), edge->getVertexRaw(i)->getVertexIdx(), 1, edge->getVertexRaw(i)->getDimensionUnfixed())
                                .noalias() = 2.0 * multipliers[edge->getEdgeIdx()] * values.transpose() * jacobf;
                        }
                        else
                        {
                            jacobian.block(edge->getEdgeIdx(), edge->getVertexRaw(i)->getVertexIdx(), 1, edge->getVertexRaw(i)->getDimensionUnfixed())
                                .noalias() = 2.0 * values.transpose() * jacobf;
                        }
                    }
                    else
                    {
                        if (multipliers)
                            edge->computeJacobian(i, jacobian.block(edge->getEdgeIdx(), edge->getVertexRaw(i)->getVertexIdx(), edge->getDimension(),
                                                                    edge->getVertexRaw(i)->getDimensionUnfixed()),
                                                  multipliers + edge->getEdgeIdx());
                        else
                            edge->computeJacobian(i, jacobian.block(edge->getEdgeIdx(), edge->getVertexRaw(i)->getVertexIdx(), edge->getDimension(),
                                                                    edge->getVertexRaw(i)->getDimensionUnfixed()),
                                                  nullptr);
                    }
                }
            }
        }
        if (num_custom_jacobians == edges.size() || edge_based) return;  // we have everything completed already
    }

    // now let's iterate vertices
    int col_idx = 0;
    for (VertexInterface* vertex : vertices)
    {
        if (vertex->getDimensionUnfixed() == 0) continue;  // without increasing col_idx

        int num_connected_edges;
        static_if(type == EdgeType::Objective)
        {
            num_connected_edges = (int)vertex->getConnectedObjectiveEdgesRef().size();
            if (vertex->getNumObjectiveEdgesWithCustomJacobian() == num_connected_edges)
            {
                ++col_idx;
                continue;  // skip what we already have...
            }
        }
        else static_if(type == EdgeType::Equality)
        {
            num_connected_edges = (int)vertex->getConnectedEqualityEdgesRef().size();
            if (vertex->getNumEqualityEdgesWithCustomJacobian() == num_connected_edges)
            {
                ++col_idx;
                continue;  // skip what we already have...
            }
        }
        else
        {
            num_connected_edges = (int)vertex->getConnectedInequalityEdgesRef().size();
            if (vertex->getNumInequalityEdgesWithCustomJacobian() == num_connected_edges)
            {
                ++col_idx;
                continue;  // skip what we already have...
            }
        }

        std::vector<Eigen::VectorXd> values2(num_connected_edges);  // if edges provide custom jacobians, we reserve too much.

        for (int i = 0; i < vertex->getDimension(); ++i)
        {
            if (vertex->isFixedComponent(i)) continue;  // continue without increasing col_idx

            std::set<EdgeInterface*>* connected_edges;
            static_if(type == EdgeType::Objective) connected_edges     = &vertex->getConnectedObjectiveEdgesRef();
            else static_if(type == EdgeType::Equality) connected_edges = &vertex->getConnectedEqualityEdgesRef();
            else connected_edges                                       = &vertex->getConnectedInequalityEdgesRef();

            vertex->plus(i, delta);  // increment current value by delta
            // compute values for all connected edges
            int edge_idx = 0;
            for (EdgeInterface* edge : *connected_edges)
            {
                if (!edge->providesJacobian())
                {
                    edge->computeValues();
                    values2[edge_idx] = Eigen::Map<const Eigen::VectorXd>(edge->getValues(), edge->getDimension());
                }
                ++edge_idx;
            }
            vertex->plus(i, neg2delta);  // subtract 2*delta
            edge_idx = 0;
            for (EdgeInterface* edge : *connected_edges)
            {
                if (!edge->providesJacobian())
                {
                    edge->computeValues();
                    const Eigen::Map<const Eigen::VectorXd> values(edge->getValues(), edge->getDimension());
                    // in least squares form, we must compute the jacobian of f(x)^2, hence -> 2 * f(x)^T * Jacobian(f(x)).
                    if (!_retain_lsq_form && edge->isLeastSquaresForm())
                    {
                        // get [1 x 1] jacobi block in the jacobi matrix
                        if (multipliers)
                            jacobian(edge->getEdgeIdx(), col_idx) =
                                2.0 * values.transpose() * scalar * multipliers[edge->getEdgeIdx()] * (values2[edge_idx] - values);
                        else
                            jacobian(edge->getEdgeIdx(), col_idx) = 2.0 * values.transpose() * scalar * (values2[edge_idx] - values);
                    }
                    else
                    {
                        // get [edge_dim x 1] jacobi block in the jacobi matrix
                        if (multipliers)
                        {
                            const Eigen::Map<const Eigen::VectorXd> values_multp(multipliers + edge->getEdgeIdx(), edge->getDimension());
                            jacobian.block(edge->getEdgeIdx(), col_idx, edge->getDimension(), 1).noalias() =
                                scalar * (values2[edge_idx] - values).cwiseProduct(values_multp);
                        }
                        else
                        {
                            jacobian.block(edge->getEdgeIdx(), col_idx, edge->getDimension(), 1).noalias() = scalar * (values2[edge_idx] - values);
                        }
                    }
                }
                ++edge_idx;
            }
            vertex->plus(i, delta);  // revert offset

            // increase column index
            ++col_idx;
        }
    }
}

}  // namespace corbo
