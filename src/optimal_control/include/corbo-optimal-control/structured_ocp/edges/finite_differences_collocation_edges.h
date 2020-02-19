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

#ifndef SRC_OPTIMAL_CONTROL_INCLUDE_CORBO_OPTIMAL_CONTROL_STRUCTURED_OCP_EDGES_FINITE_DIFFERENCES_COLLOCATION_EDGES_H_
#define SRC_OPTIMAL_CONTROL_INCLUDE_CORBO_OPTIMAL_CONTROL_STRUCTURED_OCP_EDGES_FINITE_DIFFERENCES_COLLOCATION_EDGES_H_

#include <corbo-optimization/hyper_graph/edge.h>

#include <corbo-numerics/finite_differences_collocation.h>
#include <corbo-optimal-control/functions/stage_functions.h>
#include <corbo-optimization/hyper_graph/scalar_vertex.h>
#include <corbo-optimization/hyper_graph/vector_vertex.h>
#include <corbo-systems/system_dynamics_interface.h>

#include <functional>
#include <memory>

namespace corbo {

// TODO(roesmann): if we combine StageCost / StageEqualtiy ... we can merge some of the edges below

class FDCollocationEdge : public Edge<VectorVertex, VectorVertex, VectorVertex, ScalarVertex>
{
 public:
    using Ptr  = std::shared_ptr<FDCollocationEdge>;
    using UPtr = std::unique_ptr<FDCollocationEdge>;

    explicit FDCollocationEdge(SystemDynamicsInterface::Ptr dynamics, VectorVertex& x1, VectorVertex& u1, VectorVertex& x2, ScalarVertex& dt)
        : Edge<VectorVertex, VectorVertex, VectorVertex, ScalarVertex>(x1, u1, x2, dt), _dynamics(dynamics)
    {
        _x1 = static_cast<const VectorVertex*>(_vertices[0]);
        _u1 = static_cast<const VectorVertex*>(_vertices[1]);
        _x2 = static_cast<const VectorVertex*>(_vertices[2]);
        _dt = static_cast<const ScalarVertex*>(_vertices[3]);
    }

    // implements interface method
    int getDimension() const override
    {
        assert(_dynamics);
        return _dynamics->getStateDimension();
    }
    // implement in child class:
    bool isLinear() const override { return false; }

    // bool providesJacobian() const override { return false; }

    bool isLeastSquaresForm() const override { return false; }

    void computeValues(Eigen::Ref<Eigen::VectorXd> values) override
    {
        assert(_fd_eval);
        assert(_dynamics);
        assert(_x1->getDimension() == _dynamics->getStateDimension());
        assert(_u1->getDimension() == _dynamics->getInputDimension());
        assert(_x2->getDimension() == _dynamics->getStateDimension());

        _fd_eval->computeEqualityConstraint(_x1->values(), _u1->values(), _x2->values(), _dt->value(), *_dynamics, values);
    }

    void setFiniteDifferencesCollocationMethod(FiniteDifferencesCollocationInterface::Ptr fd_eval) { _fd_eval = fd_eval; }

 private:
    SystemDynamicsInterface::Ptr _dynamics;
    FiniteDifferencesCollocationInterface::Ptr _fd_eval;

    const VectorVertex* _x1 = nullptr;
    const VectorVertex* _u1 = nullptr;
    const VectorVertex* _x2 = nullptr;
    const ScalarVertex* _dt = nullptr;

 public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

// TODO(roesmann): maybe a mixed edge for the uncompressed HS with linear control to compute uc just once.
class TrapezoidalIntegralCostEdge : public Edge<VectorVertex, VectorVertex, VectorVertex, ScalarVertex>
{
 public:
    using Ptr  = std::shared_ptr<TrapezoidalIntegralCostEdge>;
    using UPtr = std::unique_ptr<TrapezoidalIntegralCostEdge>;

    explicit TrapezoidalIntegralCostEdge(VectorVertex& x1, VectorVertex& u1, VectorVertex& x2, ScalarVertex& dt, StageCost::Ptr stage_cost, int k)
        : Edge<VectorVertex, VectorVertex, VectorVertex, ScalarVertex>(x1, u1, x2, dt), _stage_cost(stage_cost), _k(k)
    {
        assert(stage_cost);
        assert(stage_cost->getIntegralStateControlTermDimension(_k) == 1);

        _cost1.resize(1);
        _cost2.resize(1);
    }

    virtual ~TrapezoidalIntegralCostEdge() = default;

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
        const VectorVertex* x2 = static_cast<const VectorVertex*>(_vertices[2]);
        const ScalarVertex* dt = static_cast<const ScalarVertex*>(_vertices[3]);

        _stage_cost->computeIntegralStateControlTerm(_k, x1->values(), u1->values(), _cost1);
        _stage_cost->computeIntegralStateControlTerm(_k, x2->values(), u1->values(), _cost2);
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

class TrapezoidalIntegralEqualityDynamicsEdge : public Edge<VectorVertex, VectorVertex, VectorVertex, ScalarVertex>
{
 public:
    using Ptr  = std::shared_ptr<TrapezoidalIntegralEqualityDynamicsEdge>;
    using UPtr = std::unique_ptr<TrapezoidalIntegralEqualityDynamicsEdge>;

    explicit TrapezoidalIntegralEqualityDynamicsEdge(SystemDynamicsInterface::Ptr dynamics, VectorVertex& x1, VectorVertex& u1, VectorVertex& x2,
                                                     ScalarVertex& dt, StageEqualityConstraint::Ptr stage_eq, int k)
        : Edge<VectorVertex, VectorVertex, VectorVertex, ScalarVertex>(x1, u1, x2, dt), _dynamics(dynamics), _stage_eq(stage_eq), _k(k)
    {
        _dim_dyn    = _dynamics->getStateDimension();
        _dim_int_eq = _stage_eq ? _stage_eq->getIntegralStateControlTermDimension(_k) : 0;
        _dim_eq     = _dim_dyn + _dim_int_eq;

        _eq1.resize(_dim_int_eq);
        _eq2.resize(_dim_int_eq);
    }

    virtual ~TrapezoidalIntegralEqualityDynamicsEdge() = default;

    // implements interface method
    int getDimension() const override { return _dim_eq; }

    // implement in child class:
    bool isLinear() const override { return false; }

    // bool providesJacobian() const override { return false; }
    bool isLeastSquaresForm() const override { return false; }

    void computeValues(Eigen::Ref<Eigen::VectorXd> values) override
    {
        assert(_dynamics);
        assert(_fd_eval);
        const VectorVertex* x1 = static_cast<const VectorVertex*>(_vertices[0]);
        const VectorVertex* u1 = static_cast<const VectorVertex*>(_vertices[1]);
        const VectorVertex* x2 = static_cast<const VectorVertex*>(_vertices[2]);
        const ScalarVertex* dt = static_cast<const ScalarVertex*>(_vertices[3]);

        _fd_eval->computeEqualityConstraint(x1->values(), u1->values(), x2->values(), dt->value(), *_dynamics, values.head(_dim_dyn));

        if (_dim_int_eq > 0)
        {
            _stage_eq->computeIntegralStateControlTerm(_k, x1->values(), u1->values(), _eq1);
            _stage_eq->computeIntegralStateControlTerm(_k, x2->values(), u1->values(), _eq2);
            values.tail(_dim_int_eq).noalias() = 0.5 * dt->value() * (_eq1 + _eq2);
        }
    }

    void setFiniteDifferencesCollocationMethod(FiniteDifferencesCollocationInterface::Ptr fd_eval) { _fd_eval = fd_eval; }

 private:
    SystemDynamicsInterface::Ptr _dynamics;
    FiniteDifferencesCollocationInterface::Ptr _fd_eval;

    StageEqualityConstraint::Ptr _stage_eq;

    int _dim_eq     = 0;
    int _dim_int_eq = 0;
    int _dim_dyn    = 0;

    Eigen::VectorXd _eq1;
    Eigen::VectorXd _eq2;

    int _k = 0;

 public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

class TrapezoidalIntegralEqualityEdge : public Edge<VectorVertex, VectorVertex, VectorVertex, ScalarVertex>
{
 public:
    using Ptr  = std::shared_ptr<TrapezoidalIntegralEqualityEdge>;
    using UPtr = std::unique_ptr<TrapezoidalIntegralEqualityEdge>;

    explicit TrapezoidalIntegralEqualityEdge(VectorVertex& x1, VectorVertex& u1, VectorVertex& x2, ScalarVertex& dt,
                                             StageEqualityConstraint::Ptr stage_eq, int k)
        : Edge<VectorVertex, VectorVertex, VectorVertex, ScalarVertex>(x1, u1, x2, dt), _stage_eq(stage_eq), _k(k)
    {
        _dim_eq = _stage_eq ? _stage_eq->getIntegralStateControlTermDimension(_k) : 0;

        _eq1.resize(_dim_eq);
        _eq2.resize(_dim_eq);
    }

    virtual ~TrapezoidalIntegralEqualityEdge() = default;

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
        const VectorVertex* x2 = static_cast<const VectorVertex*>(_vertices[2]);
        const ScalarVertex* dt = static_cast<const ScalarVertex*>(_vertices[3]);

        _stage_eq->computeIntegralStateControlTerm(_k, x1->values(), u1->values(), _eq1);
        _stage_eq->computeIntegralStateControlTerm(_k, x2->values(), u1->values(), _eq2);
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

class TrapezoidalIntegralInequalityEdge : public Edge<VectorVertex, VectorVertex, VectorVertex, ScalarVertex>
{
 public:
    using Ptr  = std::shared_ptr<TrapezoidalIntegralInequalityEdge>;
    using UPtr = std::unique_ptr<TrapezoidalIntegralInequalityEdge>;

    explicit TrapezoidalIntegralInequalityEdge(VectorVertex& x1, VectorVertex& u1, VectorVertex& x2, ScalarVertex& dt,
                                               StageInequalityConstraint::Ptr stage_ineq, int k)
        : Edge<VectorVertex, VectorVertex, VectorVertex, ScalarVertex>(x1, u1, x2, dt), _stage_ineq(stage_ineq), _k(k)
    {
        _dim_ineq = _stage_ineq ? _stage_ineq->getIntegralStateControlTermDimension(_k) : 0;

        _ineq1.resize(_dim_ineq);
        _ineq2.resize(_dim_ineq);
    }

    virtual ~TrapezoidalIntegralInequalityEdge() = default;

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
        const VectorVertex* x2 = static_cast<const VectorVertex*>(_vertices[2]);
        const ScalarVertex* dt = static_cast<const ScalarVertex*>(_vertices[3]);

        _stage_ineq->computeIntegralStateControlTerm(_k, x1->values(), u1->values(), _ineq1);
        _stage_ineq->computeIntegralStateControlTerm(_k, x2->values(), u1->values(), _ineq2);
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

class LeftSumCostEdge : public Edge<VectorVertex, VectorVertex, ScalarVertex>
{
 public:
    using Ptr  = std::shared_ptr<LeftSumCostEdge>;
    using UPtr = std::unique_ptr<LeftSumCostEdge>;

    explicit LeftSumCostEdge(VectorVertex& x1, VectorVertex& u1, ScalarVertex& dt, StageCost::Ptr stage_cost, int k)
        : Edge<VectorVertex, VectorVertex, ScalarVertex>(x1, u1, dt), _stage_cost(stage_cost), _k(k)
    {
        assert(stage_cost);
        assert(stage_cost->getIntegralStateControlTermDimension(_k) == 1);
    }

    virtual ~LeftSumCostEdge() = default;

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
        const ScalarVertex* dt = static_cast<const ScalarVertex*>(_vertices[2]);

        _stage_cost->computeIntegralStateControlTerm(_k, x1->values(), u1->values(), values);
        values *= dt->value();
    }

 private:
    StageCost::Ptr _stage_cost;

    int _k = 0;

 public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

class LeftSumEqualityEdge : public Edge<VectorVertex, VectorVertex, ScalarVertex>
{
 public:
    using Ptr  = std::shared_ptr<LeftSumEqualityEdge>;
    using UPtr = std::unique_ptr<LeftSumEqualityEdge>;

    explicit LeftSumEqualityEdge(VectorVertex& x1, VectorVertex& u1, ScalarVertex& dt, StageEqualityConstraint::Ptr stage_eq, int k)
        : Edge<VectorVertex, VectorVertex, ScalarVertex>(x1, u1, dt), _stage_eq(stage_eq), _k(k)
    {
        _dim_eq = _stage_eq ? _stage_eq->getIntegralStateControlTermDimension(_k) : 0;
    }

    virtual ~LeftSumEqualityEdge() = default;

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
        const ScalarVertex* dt = static_cast<const ScalarVertex*>(_vertices[2]);

        _stage_eq->computeIntegralStateControlTerm(_k, x1->values(), u1->values(), values);
        values *= dt->value();
    }

 private:
    StageEqualityConstraint::Ptr _stage_eq;

    int _dim_eq = 0;

    int _k = 0;

 public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

class LeftSumInequalityEdge : public Edge<VectorVertex, VectorVertex, ScalarVertex>
{
 public:
    using Ptr  = std::shared_ptr<LeftSumInequalityEdge>;
    using UPtr = std::unique_ptr<LeftSumInequalityEdge>;

    explicit LeftSumInequalityEdge(VectorVertex& x1, VectorVertex& u1, ScalarVertex& dt, StageInequalityConstraint::Ptr stage_ineq, int k)
        : Edge<VectorVertex, VectorVertex, ScalarVertex>(x1, u1, dt), _stage_ineq(stage_ineq), _k(k)
    {
        _dim_ineq = _stage_ineq ? _stage_ineq->getIntegralStateControlTermDimension(_k) : 0;
    }

    virtual ~LeftSumInequalityEdge() = default;

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
        const ScalarVertex* dt = static_cast<const ScalarVertex*>(_vertices[2]);

        _stage_ineq->computeIntegralStateControlTerm(_k, x1->values(), u1->values(), values);
        values *= dt->value();
    }

 private:
    StageInequalityConstraint::Ptr _stage_ineq;

    int _dim_ineq = 0;

    int _k = 0;

 public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

}  // namespace corbo

#endif  // SRC_OPTIMAL_CONTROL_INCLUDE_CORBO_OPTIMAL_CONTROL_STRUCTURED_OCP_EDGES_FINITE_DIFFERENCES_COLLOCATION_EDGES_H_
