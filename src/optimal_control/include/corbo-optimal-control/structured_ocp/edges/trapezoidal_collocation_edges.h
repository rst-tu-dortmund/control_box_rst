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

#ifndef SRC_OPTIMAL_CONTROL_INCLUDE_CORBO_OPTIMAL_CONTROL_STRUCTURED_OCP_EDGES_TRAPEZOIDAL_COLLOCATION_EDGES_H_
#define SRC_OPTIMAL_CONTROL_INCLUDE_CORBO_OPTIMAL_CONTROL_STRUCTURED_OCP_EDGES_TRAPEZOIDAL_COLLOCATION_EDGES_H_

#include <corbo-optimization/hyper_graph/edge.h>

#include <corbo-optimal-control/functions/stage_functions.h>
#include <corbo-optimization/hyper_graph/scalar_vertex.h>
#include <corbo-optimization/hyper_graph/vector_vertex.h>
#include <corbo-systems/system_dynamics_interface.h>

#include <functional>
#include <memory>

namespace corbo {

class TrapezoidalCollocationDynamicsOnlyEdge : public Edge<VectorVertex, VectorVertex, VectorVertex, VectorVertex, ScalarVertex>
{
 public:
    using Ptr  = std::shared_ptr<TrapezoidalCollocationDynamicsOnlyEdge>;
    using UPtr = std::unique_ptr<TrapezoidalCollocationDynamicsOnlyEdge>;

    explicit TrapezoidalCollocationDynamicsOnlyEdge(SystemDynamicsInterface::Ptr dynamics, VectorVertex& x1, VectorVertex& u1, VectorVertex& u2,
                                                    VectorVertex& x2, ScalarVertex& dt)
        : Edge<VectorVertex, VectorVertex, VectorVertex, VectorVertex, ScalarVertex>(x1, u1, u2, x2, dt), _dynamics(dynamics)
    {
        _f1.resize(_dynamics->getStateDimension());
        _f2.resize(_dynamics->getStateDimension());
    }

    virtual ~TrapezoidalCollocationDynamicsOnlyEdge() = default;

    // implements interface method
    int getDimension() const override { return _dynamics ? _dynamics->getStateDimension() : 0; }
    // implement in child class:
    bool isLinear() const override { return false; }

    // bool providesJacobian() const override { return false; }

    bool isLeastSquaresForm() const override { return false; }

    void computeValues(Eigen::Ref<Eigen::VectorXd> values) override
    {
        assert(_dynamics);
        assert(_f1.size() == _dynamics->getStateDimension());
        assert(_f2.size() == _dynamics->getStateDimension());

        const VectorVertex* x1 = static_cast<const VectorVertex*>(_vertices[0]);
        const VectorVertex* u1 = static_cast<const VectorVertex*>(_vertices[1]);
        const VectorVertex* u2 = static_cast<const VectorVertex*>(_vertices[2]);
        const VectorVertex* x2 = static_cast<const VectorVertex*>(_vertices[3]);
        const ScalarVertex* dt = static_cast<const ScalarVertex*>(_vertices[4]);

        _dynamics->dynamics(x1->values(), u1->values(), _f1);
        _dynamics->dynamics(x2->values(), u2->values(), _f2);

        values.noalias() = 0.5 * dt->value() * (_f1 + _f2) - (x2->values() - x1->values());
    }

 protected:
 private:
    SystemDynamicsInterface::Ptr _dynamics;

    Eigen::VectorXd _f1;
    Eigen::VectorXd _f2;

 public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

// TODO(roesmann): maybe a mixed edge for the uncompressed HS with linear control to compute uc just once.
class TrapezoidalCollocationIntegralCostEdge : public Edge<VectorVertex, VectorVertex, VectorVertex, VectorVertex, ScalarVertex>
{
 public:
    using Ptr  = std::shared_ptr<TrapezoidalCollocationIntegralCostEdge>;
    using UPtr = std::unique_ptr<TrapezoidalCollocationIntegralCostEdge>;

    explicit TrapezoidalCollocationIntegralCostEdge(VectorVertex& x1, VectorVertex& u1, VectorVertex& u2, VectorVertex& x2, ScalarVertex& dt,
                                                    StageCost::Ptr stage_cost, int k)
        : Edge<VectorVertex, VectorVertex, VectorVertex, VectorVertex, ScalarVertex>(x1, u1, u2, x2, dt), _stage_cost(stage_cost), _k(k)
    {
        assert(stage_cost);
        assert(stage_cost->getIntegralStateControlTermDimension(_k) == 1);

        _cost1.resize(1);
        _cost2.resize(1);
    }

    virtual ~TrapezoidalCollocationIntegralCostEdge() = default;

    // implements interface method
    int getDimension() const override { return 1; }

    // implement in child class:
    bool isLinear() const override { return false; }

    // bool providesJacobian() const override { return false; }
    bool isLeastSquaresForm() const override { return false; }

    void computeValues(Eigen::Ref<Eigen::VectorXd> values) override
    {
        const VectorVertex* x1 = static_cast<const VectorVertex*>(_vertices[0]);
        const VectorVertex* u1 = static_cast<const VectorVertex*>(_vertices[1]);
        const VectorVertex* u2 = static_cast<const VectorVertex*>(_vertices[2]);
        const VectorVertex* x2 = static_cast<const VectorVertex*>(_vertices[3]);
        const ScalarVertex* dt = static_cast<const ScalarVertex*>(_vertices[4]);

        _stage_cost->computeIntegralStateControlTerm(_k, x1->values(), u1->values(), _cost1);
        _stage_cost->computeIntegralStateControlTerm(_k, x2->values(), u2->values(), _cost2);
        values[0] = 0.5 * dt->value() * (_cost1[0] + _cost2[0]);
    }

 private:
    StageCost::Ptr _stage_cost;

    int _k = 0;

    Eigen::VectorXd _cost1;
    Eigen::VectorXd _cost2;

 public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

class TrapezoidalCollocationIntegralEqualityDynamicsEdge : public Edge<VectorVertex, VectorVertex, VectorVertex, VectorVertex, ScalarVertex>
{
 public:
    using Ptr  = std::shared_ptr<TrapezoidalCollocationIntegralEqualityDynamicsEdge>;
    using UPtr = std::unique_ptr<TrapezoidalCollocationIntegralEqualityDynamicsEdge>;

    explicit TrapezoidalCollocationIntegralEqualityDynamicsEdge(SystemDynamicsInterface::Ptr dynamics, VectorVertex& x1, VectorVertex& u1,
                                                                VectorVertex& u2, VectorVertex& x2, ScalarVertex& dt,
                                                                StageEqualityConstraint::Ptr stage_eq, int k)
        : Edge<VectorVertex, VectorVertex, VectorVertex, VectorVertex, ScalarVertex>(x1, u1, u2, x2, dt),
          _dynamics(dynamics),
          _stage_eq(stage_eq),
          _k(k)
    {
        _dim_dyn    = _dynamics->getStateDimension();
        _dim_int_eq = _stage_eq ? _stage_eq->getIntegralStateControlTermDimension(_k) : 0;
        _dim_eq     = _dim_dyn + _dim_int_eq;

        _eq1.resize(_dim_int_eq);
        _eq2.resize(_dim_int_eq);

        _f1.resize(_dim_dyn);
        _f2.resize(_dim_dyn);
    }

    virtual ~TrapezoidalCollocationIntegralEqualityDynamicsEdge() = default;

    // implements interface method
    int getDimension() const override { return _dim_eq; }

    // implement in child class:
    bool isLinear() const override { return false; }

    // bool providesJacobian() const override { return false; }
    bool isLeastSquaresForm() const override { return false; }

    void computeValues(Eigen::Ref<Eigen::VectorXd> values) override
    {
        const VectorVertex* x1 = static_cast<const VectorVertex*>(_vertices[0]);
        const VectorVertex* u1 = static_cast<const VectorVertex*>(_vertices[1]);
        const VectorVertex* u2 = static_cast<const VectorVertex*>(_vertices[2]);
        const VectorVertex* x2 = static_cast<const VectorVertex*>(_vertices[3]);
        const ScalarVertex* dt = static_cast<const ScalarVertex*>(_vertices[4]);

        _dynamics->dynamics(x1->values(), u1->values(), _f1);
        _dynamics->dynamics(x2->values(), u2->values(), _f2);

        values.head(_dim_dyn).noalias() = 0.5 * dt->value() * (_f1 + _f2) - (x2->values() - x1->values());

        if (_dim_int_eq > 0)
        {
            _stage_eq->computeIntegralStateControlTerm(_k, x1->values(), u1->values(), _eq1);
            _stage_eq->computeIntegralStateControlTerm(_k, x2->values(), u2->values(), _eq2);
            values.tail(_dim_int_eq).noalias() = 0.5 * dt->value() * (_eq1 + _eq2);
        }
    }

 private:
    SystemDynamicsInterface::Ptr _dynamics;

    StageEqualityConstraint::Ptr _stage_eq;

    int _dim_eq     = 0;
    int _dim_int_eq = 0;
    int _dim_dyn    = 0;

    Eigen::VectorXd _eq1;
    Eigen::VectorXd _eq2;

    int _k = 0;

    Eigen::VectorXd _f1;
    Eigen::VectorXd _f2;

 public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

class TrapezoidalCollocationIntegralEqualityEdge : public Edge<VectorVertex, VectorVertex, VectorVertex, VectorVertex, ScalarVertex>
{
 public:
    using Ptr  = std::shared_ptr<TrapezoidalCollocationIntegralEqualityEdge>;
    using UPtr = std::unique_ptr<TrapezoidalCollocationIntegralEqualityEdge>;

    explicit TrapezoidalCollocationIntegralEqualityEdge(VectorVertex& x1, VectorVertex& u1, VectorVertex& u2, VectorVertex& x2, ScalarVertex& dt,
                                                        StageEqualityConstraint::Ptr stage_eq, int k)
        : Edge<VectorVertex, VectorVertex, VectorVertex, VectorVertex, ScalarVertex>(x1, u1, u2, x2, dt), _stage_eq(stage_eq), _k(k)
    {
        _dim_eq = _stage_eq ? _stage_eq->getIntegralStateControlTermDimension(_k) : 0;

        _eq1.resize(_dim_eq);
        _eq2.resize(_dim_eq);
    }

    virtual ~TrapezoidalCollocationIntegralEqualityEdge() = default;

    // implements interface method
    int getDimension() const override { return _dim_eq; }

    // implement in child class:
    bool isLinear() const override { return false; }

    // bool providesJacobian() const override { return false; }
    bool isLeastSquaresForm() const override { return false; }

    void computeValues(Eigen::Ref<Eigen::VectorXd> values) override
    {
        const VectorVertex* x1 = static_cast<const VectorVertex*>(_vertices[0]);
        const VectorVertex* u1 = static_cast<const VectorVertex*>(_vertices[1]);
        const VectorVertex* u2 = static_cast<const VectorVertex*>(_vertices[2]);
        const VectorVertex* x2 = static_cast<const VectorVertex*>(_vertices[3]);
        const ScalarVertex* dt = static_cast<const ScalarVertex*>(_vertices[4]);

        _stage_eq->computeIntegralStateControlTerm(_k, x1->values(), u1->values(), _eq1);
        _stage_eq->computeIntegralStateControlTerm(_k, x2->values(), u2->values(), _eq2);
        values.noalias() = 0.5 * dt->value() * (_eq1 + _eq2);
    }

 private:
    StageEqualityConstraint::Ptr _stage_eq;

    int _dim_eq = 0;

    Eigen::VectorXd _eq1;
    Eigen::VectorXd _eq2;

    int _k = 0;

 public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

class TrapezoidalCollocationIntegralInequalityEdge : public Edge<VectorVertex, VectorVertex, VectorVertex, VectorVertex, ScalarVertex>
{
 public:
    using Ptr  = std::shared_ptr<TrapezoidalCollocationIntegralInequalityEdge>;
    using UPtr = std::unique_ptr<TrapezoidalCollocationIntegralInequalityEdge>;

    explicit TrapezoidalCollocationIntegralInequalityEdge(VectorVertex& x1, VectorVertex& u1, VectorVertex& u2, VectorVertex& x2, ScalarVertex& dt,
                                                          StageInequalityConstraint::Ptr stage_ineq, int k)
        : Edge<VectorVertex, VectorVertex, VectorVertex, VectorVertex, ScalarVertex>(x1, u1, u2, x2, dt), _stage_ineq(stage_ineq), _k(k)
    {
        _dim_ineq = _stage_ineq ? _stage_ineq->getIntegralStateControlTermDimension(_k) : 0;

        _ineq1.resize(_dim_ineq);
        _ineq2.resize(_dim_ineq);
    }

    virtual ~TrapezoidalCollocationIntegralInequalityEdge() = default;

    // implements interface method
    int getDimension() const override { return _dim_ineq; }

    // implement in child class:
    bool isLinear() const override { return false; }

    // bool providesJacobian() const override { return false; }

    bool isLeastSquaresForm() const override { return false; }

    void computeValues(Eigen::Ref<Eigen::VectorXd> values) override
    {
        const VectorVertex* x1 = static_cast<const VectorVertex*>(_vertices[0]);
        const VectorVertex* u1 = static_cast<const VectorVertex*>(_vertices[1]);
        const VectorVertex* u2 = static_cast<const VectorVertex*>(_vertices[2]);
        const VectorVertex* x2 = static_cast<const VectorVertex*>(_vertices[3]);
        const ScalarVertex* dt = static_cast<const ScalarVertex*>(_vertices[4]);

        _stage_ineq->computeIntegralStateControlTerm(_k, x1->values(), u1->values(), _ineq1);
        _stage_ineq->computeIntegralStateControlTerm(_k, x2->values(), u2->values(), _ineq2);
        values.noalias() = 0.5 * dt->value() * (_ineq1 + _ineq2);
    }

 private:
    StageInequalityConstraint::Ptr _stage_ineq;

    int _dim_ineq = 0;

    Eigen::VectorXd _ineq1;
    Eigen::VectorXd _ineq2;

    int _k = 0;

 public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

}  // namespace corbo

#endif  // SRC_OPTIMAL_CONTROL_INCLUDE_CORBO_OPTIMAL_CONTROL_STRUCTURED_OCP_EDGES_TRAPEZOIDAL_COLLOCATION_EDGES_H_
