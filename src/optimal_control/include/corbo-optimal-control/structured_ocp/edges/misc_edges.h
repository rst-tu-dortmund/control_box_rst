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

#ifndef SRC_OPTIMAL_CONTROL_INCLUDE_CORBO_OPTIMAL_CONTROL_STRUCTURED_OCP_EDGES_MISC_EDGES_H_
#define SRC_OPTIMAL_CONTROL_INCLUDE_CORBO_OPTIMAL_CONTROL_STRUCTURED_OCP_EDGES_MISC_EDGES_H_

#include <corbo-optimization/hyper_graph/edge.h>

#include <corbo-optimal-control/functions/stage_functions.h>
#include <corbo-optimization/hyper_graph/scalar_vertex.h>
#include <corbo-optimization/hyper_graph/vector_vertex.h>
#include <corbo-systems/system_dynamics_interface.h>

#include <functional>
#include <memory>

namespace corbo {

class TwoScalarEqualEdge : public Edge<ScalarVertex, ScalarVertex>
{
 public:
    using Ptr  = std::shared_ptr<TwoScalarEqualEdge>;
    using UPtr = std::unique_ptr<TwoScalarEqualEdge>;

    explicit TwoScalarEqualEdge(ScalarVertex& s1, ScalarVertex& s2) : Edge<ScalarVertex, ScalarVertex>(s1, s2) {}

    // implements interface method
    int getDimension() const override { return 1; }
    // implement in child class:
    bool isLinear() const override { return true; }

    // bool providesJacobian() const override { return false; }

    bool isLeastSquaresForm() const override { return false; }

    void computeValues(Eigen::Ref<Eigen::VectorXd> values) override
    {
        const ScalarVertex* s1 = static_cast<const ScalarVertex*>(_vertices[0]);
        const ScalarVertex* s2 = static_cast<const ScalarVertex*>(_vertices[1]);

        values.coeffRef(0) = s2->value() - s1->value();
    }

 public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

}  // namespace corbo

#endif  // SRC_OPTIMAL_CONTROL_INCLUDE_CORBO_OPTIMAL_CONTROL_STRUCTURED_OCP_EDGES_MISC_EDGES_H_
