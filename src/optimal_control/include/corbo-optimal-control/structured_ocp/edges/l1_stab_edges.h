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

#ifndef SRC_OPTIMAL_CONTROL_INCLUDE_CORBO_OPTIMAL_CONTROL_STRUCTURED_OCP_EDGES_L1_STAB_EDGES_H_
#define SRC_OPTIMAL_CONTROL_INCLUDE_CORBO_OPTIMAL_CONTROL_STRUCTURED_OCP_EDGES_L1_STAB_EDGES_H_

#include <corbo-optimization/hyper_graph/edge.h>

#include <corbo-numerics/finite_differences_collocation.h>
#include <corbo-optimal-control/functions/stage_functions.h>
#include <corbo-optimization/hyper_graph/scalar_vertex.h>
#include <corbo-optimization/hyper_graph/vector_vertex.h>
#include <corbo-systems/system_dynamics_interface.h>

#include <functional>
#include <memory>

namespace corbo {

class L1StabCostEdge : public Edge<VectorVertex>
{
 public:
    using Ptr  = std::shared_ptr<L1StabCostEdge>;
    using UPtr = std::unique_ptr<L1StabCostEdge>;

    explicit L1StabCostEdge(VectorVertex& s, double delta_pow_k) : Edge<VectorVertex>(s)
    {
        _s           = static_cast<const VectorVertex*>(_vertices[0]);
        _delta_pow_k = delta_pow_k;
    }

    // implements interface method
    int getDimension() const override { return 1; }
    // implement in child class:
    bool isLinear() const override { return false; }

    // bool providesJacobian() const override { return false; }

    bool isLeastSquaresForm() const override { return false; }

    void computeValues(Eigen::Ref<Eigen::VectorXd> values) override
    {
        assert(_s);
        values[0] = _delta_pow_k * _s->values().sum();
    }

 private:
    const VectorVertex* _s = nullptr;

    double _delta_pow_k = 1.6;

 public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

class L1StabInequalityEdge : public Edge<VectorVertex, VectorVertex>
{
 public:
    using Ptr  = std::shared_ptr<L1StabInequalityEdge>;
    using UPtr = std::unique_ptr<L1StabInequalityEdge>;

    explicit L1StabInequalityEdge(VectorVertex& x, VectorVertex& s, const Eigen::VectorXd* xref) : Edge<VectorVertex, VectorVertex>(x, s)
    {
        _x    = static_cast<const VectorVertex*>(_vertices[0]);
        _s    = static_cast<const VectorVertex*>(_vertices[1]);
        _xref = xref;
        assert(_x);
        _dim_x = _x->getDimension();
    }

    // implements interface method
    int getDimension() const override
    {
        return 2 * _dim_x;  // TODO(roesmann) or unfixed?
    }
    // implement in child class:
    bool isLinear() const override { return false; }

    // bool providesJacobian() const override { return false; }

    bool isLeastSquaresForm() const override { return false; }

    void computeValues(Eigen::Ref<Eigen::VectorXd> values) override
    {
        assert(_s && _x);
        assert(_s->getDimension() == _dim_x);
        if (_xref)
        {
            assert(_xref->size() == _dim_x);
            values.head(_dim_x).noalias() = _x->values() - *_xref - _s->values();
            values.tail(_dim_x).noalias() = *_xref - _x->values() - _s->values();
        }
        else
        {
            values.head(_dim_x).noalias() = _x->values() - _s->values();
            values.tail(_dim_x).noalias() = -_x->values() - _s->values();
        }
    }

 private:
    const VectorVertex* _x       = nullptr;
    const VectorVertex* _s       = nullptr;
    const Eigen::VectorXd* _xref = nullptr;

    int _dim_x = 0;

 public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

}  // namespace corbo

#endif  // SRC_OPTIMAL_CONTROL_INCLUDE_CORBO_OPTIMAL_CONTROL_STRUCTURED_OCP_EDGES_L1_STAB_EDGES_H_
