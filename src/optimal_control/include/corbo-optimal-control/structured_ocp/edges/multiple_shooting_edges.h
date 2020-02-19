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

#ifndef SRC_OPTIMAL_CONTROL_INCLUDE_CORBO_OPTIMAL_CONTROL_STRUCTURED_OCP_EDGES_MULTIPLE_SHOOTING_EDGES_H_
#define SRC_OPTIMAL_CONTROL_INCLUDE_CORBO_OPTIMAL_CONTROL_STRUCTURED_OCP_EDGES_MULTIPLE_SHOOTING_EDGES_H_

#include <corbo-optimization/hyper_graph/edge.h>

#include <corbo-numerics/integrator_interface.h>
#include <corbo-optimal-control/functions/stage_functions.h>
#include <corbo-optimization/hyper_graph/scalar_vertex.h>
#include <corbo-optimization/hyper_graph/vector_vertex.h>
#include <corbo-systems/system_dynamics_interface.h>

#include <functional>
#include <memory>

namespace corbo {

// Refer to MSVariableDynamicsOnlyEdge for the variable dt case.
// we could always use the variabledt edge, since dt is fixed in that case.
// However, we need some benchmarks how much time it requires to always iterate the fixed vertex.. (probably not much?)
class MSDynamicsOnlyEdge : public Edge<VectorVertex, VectorVertex, VectorVertex>
{
 public:
    using Ptr  = std::shared_ptr<MSDynamicsOnlyEdge>;
    using UPtr = std::unique_ptr<MSDynamicsOnlyEdge>;

    explicit MSDynamicsOnlyEdge(SystemDynamicsInterface::Ptr dynamics, VectorVertex& x1, VectorVertex& u1, VectorVertex& x2, double dt)
        : Edge<VectorVertex, VectorVertex, VectorVertex>(x1, u1, x2), _dt(dt), _dynamics(dynamics)
    {
        _x1 = static_cast<const VectorVertex*>(_vertices[0]);
        _u1 = static_cast<const VectorVertex*>(_vertices[1]);
        _x2 = static_cast<const VectorVertex*>(_vertices[2]);
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
        assert(_integrator);
        assert(_dynamics);
        assert(_x1->getDimension() == _dynamics->getStateDimension());
        assert(_u1->getDimension() == _dynamics->getInputDimension());
        assert(_x2->getDimension() == _dynamics->getStateDimension());

        _integrator->computeEqualityConstraint(_x1->values(), _u1->values(), _x2->values(), _dt, *_dynamics, values);
    }

    void setIntegrator(NumericalIntegratorExplicitInterface::Ptr integrator) { _integrator = integrator; }

 private:
    double _dt;
    SystemDynamicsInterface::Ptr _dynamics;
    NumericalIntegratorExplicitInterface::Ptr _integrator;

    const VectorVertex* _x1 = nullptr;
    const VectorVertex* _u1 = nullptr;
    const VectorVertex* _x2 = nullptr;

 public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

class MSVariableDynamicsOnlyEdge : public Edge<VectorVertex, VectorVertex, VectorVertex, ScalarVertex>
{
 public:
    using Ptr  = std::shared_ptr<MSVariableDynamicsOnlyEdge>;
    using UPtr = std::unique_ptr<MSVariableDynamicsOnlyEdge>;

    explicit MSVariableDynamicsOnlyEdge(SystemDynamicsInterface::Ptr dynamics, VectorVertex& x1, VectorVertex& u1, VectorVertex& x2, ScalarVertex& dt)
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
        assert(_integrator);
        assert(_dynamics);
        assert(_x1->getDimension() == _dynamics->getStateDimension());
        assert(_u1->getDimension() == _dynamics->getInputDimension());
        assert(_x2->getDimension() == _dynamics->getStateDimension());

        _integrator->computeEqualityConstraint(_x1->values(), _u1->values(), _x2->values(), _dt->value(), *_dynamics, values);
    }

    void setIntegrator(NumericalIntegratorExplicitInterface::Ptr integrator) { _integrator = integrator; }

 private:
    SystemDynamicsInterface::Ptr _dynamics;
    NumericalIntegratorExplicitInterface::Ptr _integrator;

    const VectorVertex* _x1 = nullptr;
    const VectorVertex* _u1 = nullptr;
    const VectorVertex* _x2 = nullptr;
    const ScalarVertex* _dt = nullptr;

 public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

class MultipleShootingEdgeSingleControl : public MixedEdge<VectorVertex, VectorVertex, ScalarVertex, VectorVertex>
{
 public:
    using Ptr  = std::shared_ptr<MultipleShootingEdgeSingleControl>;
    using UPtr = std::unique_ptr<MultipleShootingEdgeSingleControl>;

    explicit MultipleShootingEdgeSingleControl(SystemDynamicsInterface::Ptr dynamics, StageCost::ConstPtr stage_cost,
                                               StageEqualityConstraint::ConstPtr stage_eq, StageInequalityConstraint::ConstPtr stage_ineq, int k,
                                               VectorVertex& x_k, VectorVertex& u_k, ScalarVertex& dt_k, VectorVertex& x_kp1)
        : MixedEdge<VectorVertex, VectorVertex, ScalarVertex, VectorVertex>(x_k, u_k, dt_k, x_kp1),
          _stage_cost(stage_cost),
          _stage_eq(stage_eq),
          _stage_ineq(stage_ineq),
          _dynamics(dynamics),
          _k(k)
    {
        assert(_dynamics);

        _obj_dim      = 1;  // we have only a single obejtive value
        _dyn_dim      = _dynamics->getStateDimension();
        _other_eq_dim = _stage_eq ? _stage_eq->getIntegralStateControlTermDimension(_k) : 0;
        _ineq_dim     = _stage_ineq ? _stage_ineq->getIntegralStateControlTermDimension(_k) : 0;

        _current.resize(_obj_dim + _dyn_dim + _other_eq_dim + _ineq_dim);

        // TODO(roesmann) see notes at member declaration;
        _values.resize(_current.size());

        _xkvert    = static_cast<const VectorVertex*>(_vertices[0]);
        _ukvert    = static_cast<const VectorVertex*>(_vertices[1]);
        _dtvert    = static_cast<const ScalarVertex*>(_vertices[2]);
        _xnextvert = static_cast<const VectorVertex*>(_vertices[3]);

        configureIntegrand();
    }

    // implements interface method
    int getObjectiveDimension() const override { return _obj_dim; }
    // implements interface method
    int getEqualityDimension() const override { return _dyn_dim + _other_eq_dim; }
    // implements interface method
    int getInequalityDimension() const override { return _ineq_dim; }

    // implement in child class:
    // TODO(roesmann) implement if appropriate
    bool isObjectiveLinear() const override { return false; }
    bool isEqualityLinear() const override { return false; }
    bool isInequalityLinear() const override { return false; }

    // bool providesJacobian() const override { return false; }

    // TODO(roesmann): integral of squared functions is not always the same as the squared integral of functions
    bool isObjectiveLeastSquaresForm() const override { return false; }

    void precompute() override
    {
        assert(_integrator);
        assert(_dynamics);

        assert(_xkvert->getDimension() == _dynamics->getStateDimension());
        assert(_ukvert->getDimension() == _dynamics->getInputDimension());

        // get current values (we must avoid alias, hence we need to copy)
        if (_obj_dim > 0) _current[0] = 0;
        _current.segment(_obj_dim, _dyn_dim) = _xkvert->values();
        int tail = _other_eq_dim + _ineq_dim;
        if (tail > 0) _current.tail(tail).setZero();

        // TODO(roesmann) check if we can ensure to avoid alias in the integrator part to speed up simple integrators
        // first simple tests with replacing _values by _current led to different (wrong) results
        // _integrator->integrate(_current, _dtvert->value(), _integrand, _values);
        _integrator->solveIVP(_current, _dtvert->value(), _integrand, _values);
    }
    void computeObjectiveValues(Eigen::Ref<Eigen::VectorXd> obj_values) override
    {
        if (_obj_dim > 0) obj_values[0] = _values[0];
    }
    void computeEqualityValues(Eigen::Ref<Eigen::VectorXd> eq_values) override
    {
        // system dynamics result must comply with subsequent state
        assert(_xnextvert->values().size() == _dyn_dim);

        // _values.segment(_obj_dim, _dyn_dim) -= x_next->values();
        eq_values.head(_dyn_dim).noalias() = _values.segment(_obj_dim, _dyn_dim) - _xnextvert->values();
        eq_values.tail(_other_eq_dim)      = _values.segment(_obj_dim + _dyn_dim, _other_eq_dim);
    }
    void computeInequalityValues(Eigen::Ref<Eigen::VectorXd> ineq_values) override
    {
        if (_ineq_dim > 0) ineq_values = _values.tail(_ineq_dim);
    }

    void setIntegrator(NumericalIntegratorExplicitInterface::Ptr integrator) { _integrator = integrator; }

 protected:
    void configureIntegrand()
    {
        _integrand = [&](const Eigen::VectorXd& current, Eigen::Ref<Eigen::VectorXd> result) {
            assert(_ukvert);
            assert(_dtvert);

            int idx = 0;
            _dynamics->dynamics(current.segment(_obj_dim, _dyn_dim), _ukvert->values(), result.segment(_obj_dim, _dyn_dim));
            if (_obj_dim > 0)
            {
                assert(_obj_dim == 1);
                _stage_cost->computeIntegralStateControlTerm(_k, current.segment(_obj_dim, _dyn_dim), _ukvert->values(), result.head(1));
                idx += 1;
            }
            idx += _dyn_dim;
            if (_other_eq_dim > 0)
            {
                _stage_eq->computeIntegralStateControlTerm(_k, current.segment(_obj_dim, _dyn_dim), _ukvert->values(),
                                                           result.segment(idx, _other_eq_dim));
                idx += _other_eq_dim;
            }
            if (_ineq_dim)
            {
                _stage_ineq->computeIntegralStateControlTerm(_k, current.segment(_obj_dim, _dyn_dim), _ukvert->values(),
                                                             result.segment(idx, _ineq_dim));
                idx += _ineq_dim;
            }
            assert(idx == result.size());
        };
    }

 private:
    const VectorVertex* _xkvert    = nullptr;
    const VectorVertex* _ukvert    = nullptr;
    const ScalarVertex* _dtvert    = nullptr;
    const VectorVertex* _xnextvert = nullptr;

    int _dyn_dim      = 0;
    int _obj_dim      = 0;
    int _other_eq_dim = 0;
    int _ineq_dim     = 0;

    Eigen::VectorXd _current;
    Eigen::VectorXd _values;  // TODO(roesmann): if integrators are alias safe we can remove _values in favor of just _current

    std::function<void(const Eigen::VectorXd&, Eigen::Ref<Eigen::VectorXd>)> _integrand;

    StageCost::ConstPtr _stage_cost;
    StageEqualityConstraint::ConstPtr _stage_eq;
    StageInequalityConstraint::ConstPtr _stage_ineq;

    SystemDynamicsInterface::Ptr _dynamics;
    NumericalIntegratorExplicitInterface::Ptr _integrator;

    int _k = 0;

 public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

class MultipleShootingEdgeSingleControlTimeScaling : public MixedEdge<VectorVertex, VectorVertex, ScalarVertex, VectorVertex>
{
 public:
    using Ptr  = std::shared_ptr<MultipleShootingEdgeSingleControlTimeScaling>;
    using UPtr = std::unique_ptr<MultipleShootingEdgeSingleControlTimeScaling>;

    explicit MultipleShootingEdgeSingleControlTimeScaling(SystemDynamicsInterface::Ptr dynamics, StageCost::ConstPtr stage_cost,
                                                          StageEqualityConstraint::ConstPtr stage_eq, StageInequalityConstraint::ConstPtr stage_ineq,
                                                          int k, VectorVertex& x_k, VectorVertex& u_k, ScalarVertex& time, VectorVertex& x_kp1,
                                                          double dt)
        : MixedEdge<VectorVertex, VectorVertex, ScalarVertex, VectorVertex>(x_k, u_k, time, x_kp1),
          _stage_cost(stage_cost),
          _stage_eq(stage_eq),
          _stage_ineq(stage_ineq),
          _dynamics(dynamics),
          _dt(dt),
          _k(k)
    {
        assert(_dynamics);

        _obj_dim      = 1;  // we have only a single obejtive value
        _dyn_dim      = _dynamics->getStateDimension();
        _other_eq_dim = _stage_eq ? _stage_eq->getIntegralStateControlTermDimension(_k) : 0;
        _ineq_dim     = _stage_ineq ? _stage_ineq->getIntegralStateControlTermDimension(_k) : 0;

        _current.resize(_obj_dim + _dyn_dim + _other_eq_dim + _ineq_dim);

        // TODO(roesmann) see notes at member declaration;
        _values.resize(_current.size());

        _xkvert    = static_cast<const VectorVertex*>(_vertices[0]);
        _ukvert    = static_cast<const VectorVertex*>(_vertices[1]);
        _timevert  = static_cast<const ScalarVertex*>(_vertices[2]);
        _xnextvert = static_cast<const VectorVertex*>(_vertices[3]);

        configureIntegrand();
    }

    // implements interface method
    int getObjectiveDimension() const override { return _obj_dim; }
    // implements interface method
    int getEqualityDimension() const override { return _dyn_dim + _other_eq_dim; }
    // implements interface method
    int getInequalityDimension() const override { return _ineq_dim; }

    // implement in child class:
    // TODO(roesmann) implement if appropriate
    bool isObjectiveLinear() const override { return false; }
    bool isEqualityLinear() const override { return false; }
    bool isInequalityLinear() const override { return false; }

    // bool providesJacobian() const override { return false; }

    // TODO(roesmann): integral of squared functions is not always the same as the squared integral of functions
    bool isObjectiveLeastSquaresForm() const override { return false; }

    void precompute() override
    {
        assert(_integrator);
        assert(_dynamics);

        assert(_xkvert->getDimension() == _dynamics->getStateDimension());
        assert(_ukvert->getDimension() == _dynamics->getInputDimension());

        // get current values (we must avoid alias, hence we need to copy)
        if (_obj_dim > 0) _current[0] = 0;
        _current.segment(_obj_dim, _dyn_dim) = _xkvert->values();
        int tail = _other_eq_dim + _ineq_dim;
        if (tail > 0) _current.tail(tail).setZero();

        // TODO(roesmann) check if we can ensure to avoid alias in the integrator part to speed up simple integrators
        // first simple tests with replacing _values by _current led to different (wrong) results
        // _integrator->integrate(_current, _dtvert->value(), _integrand, _values);
        _integrator->solveIVP(_current, _dt, _integrand, _values);
    }
    void computeObjectiveValues(Eigen::Ref<Eigen::VectorXd> obj_values) override
    {
        if (_obj_dim > 0) obj_values[0] = _values[0];
    }
    void computeEqualityValues(Eigen::Ref<Eigen::VectorXd> eq_values) override
    {
        // system dynamics result must comply with subsequent state
        assert(_xnextvert->values().size() == _dyn_dim);

        // _values.segment(_obj_dim, _dyn_dim) -= x_next->values();
        eq_values.head(_dyn_dim).noalias() = _values.segment(_obj_dim, _dyn_dim) - _xnextvert->values();
        eq_values.tail(_other_eq_dim)      = _values.segment(_obj_dim + _dyn_dim, _other_eq_dim);
    }
    void computeInequalityValues(Eigen::Ref<Eigen::VectorXd> ineq_values) override
    {
        if (_ineq_dim > 0) ineq_values = _values.tail(_ineq_dim);
    }

    void setIntegrator(NumericalIntegratorExplicitInterface::Ptr integrator) { _integrator = integrator; }

 protected:
    void configureIntegrand()
    {
        _integrand = [&](const Eigen::VectorXd& current, Eigen::Ref<Eigen::VectorXd> result) {
            assert(_ukvert);
            assert(_timevert);

            int idx = 0;
            _dynamics->dynamics(current.segment(_obj_dim, _dyn_dim), _ukvert->values(), result.segment(_obj_dim, _dyn_dim));
            // time scaling:
            result.segment(_obj_dim, _dyn_dim) *= _timevert->value();

            if (_obj_dim > 0)
            {
                assert(_obj_dim == 1);
                _stage_cost->computeIntegralStateControlTerm(_k, current.segment(_obj_dim, _dyn_dim), _ukvert->values(), result.head(1));
                idx += 1;
            }
            idx += _dyn_dim;
            if (_other_eq_dim > 0)
            {
                _stage_eq->computeIntegralStateControlTerm(_k, current.segment(_obj_dim, _dyn_dim), _ukvert->values(),
                                                           result.segment(idx, _other_eq_dim));
                idx += _other_eq_dim;
            }
            if (_ineq_dim)
            {
                _stage_ineq->computeIntegralStateControlTerm(_k, current.segment(_obj_dim, _dyn_dim), _ukvert->values(),
                                                             result.segment(idx, _ineq_dim));
                idx += _ineq_dim;
            }
            assert(idx == result.size());
        };
    }

 private:
    const VectorVertex* _xkvert    = nullptr;
    const VectorVertex* _ukvert    = nullptr;
    const ScalarVertex* _timevert  = nullptr;
    const VectorVertex* _xnextvert = nullptr;

    int _dyn_dim      = 0;
    int _obj_dim      = 0;
    int _other_eq_dim = 0;
    int _ineq_dim     = 0;

    Eigen::VectorXd _current;
    Eigen::VectorXd _values;  // TODO(roesmann): if integrators are alias safe we can remove _values in favor of just _current

    std::function<void(const Eigen::VectorXd&, Eigen::Ref<Eigen::VectorXd>)> _integrand;

    StageCost::ConstPtr _stage_cost;
    StageEqualityConstraint::ConstPtr _stage_eq;
    StageInequalityConstraint::ConstPtr _stage_ineq;

    SystemDynamicsInterface::Ptr _dynamics;
    NumericalIntegratorExplicitInterface::Ptr _integrator;

    double _dt = 0.1;

    int _k = 0;

 public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

// vertex order:
// vertex: sk
// vertex: sk+1
// vertex: u[0]
// vertex: u[1]
// vertex: u[...]
// vertex: u[nc-1]
// vertex: dt[0]
// vertex: dt[1]
// vertex: dt[...]
// vertex: dt[nt-1]
class MultipleShootingEdge : public MixedEdge<>
{
 public:
    using Ptr  = std::shared_ptr<MultipleShootingEdge>;
    using UPtr = std::unique_ptr<MultipleShootingEdge>;

    explicit MultipleShootingEdge(SystemDynamicsInterface::Ptr dynamics, StageCost::ConstPtr stage_cost, StageEqualityConstraint::ConstPtr stage_eq,
                                  StageInequalityConstraint::ConstPtr stage_ineq, int k, int n_interval, int nc_interval, int nt_interval,
                                  bool eval_intermediate_constr)
        : MixedEdge<>(2 + nc_interval + nt_interval),
          _stage_cost(stage_cost),
          _stage_eq(stage_eq),
          _stage_ineq(stage_ineq),
          _n_interval(n_interval),
          _nc_interval(nc_interval),
          _nt_interval(nt_interval),
          _eval_intermediate_constr(eval_intermediate_constr),
          _dynamics(dynamics),
          _k(k)
    {
    }

    void finalize()
    {
        assert(_dynamics);
        assert((int)_vertices.size() == 2 + _nc_interval + _nt_interval);

        _shooting_node = static_cast<const VectorVertex*>(_vertices[0]);

        _int_obj_dim = (_stage_cost && _stage_cost->getIntegralStateControlTermDimension(_k) > 0) ? 1 : 0;  // we have only a single obejtive value
        _dyn_dim     = _dynamics->getStateDimension();

        _other_int_eq_dim = _stage_eq ? _stage_eq->getIntegralStateControlTermDimension(_k) : 0;
        _int_ineq_dim     = _stage_ineq ? _stage_ineq->getIntegralStateControlTermDimension(_k) : 0;

        _current.resize(_int_obj_dim + _dyn_dim + _other_int_eq_dim + _int_ineq_dim);

        // TODO(roesmann) see notes at member declaration;
        _values.resize(_current.size());

        // chache total inequality dimension
        _total_dim_obj  = _int_obj_dim;
        _total_dim_eq   = _dyn_dim + _other_int_eq_dim;
        _total_dim_ineq = _int_ineq_dim;

        activateIntermediateStateCostAndConstraints();

        configureIntegrand();
    }

    // implements interface method
    int getObjectiveDimension() const override { return _total_dim_obj; }
    // implements interface method
    int getEqualityDimension() const override { return _total_dim_eq; }
    // implements interface method
    int getInequalityDimension() const override { return _total_dim_ineq; }

    // implement in child class:
    bool isObjectiveLinear() const override { return false; }
    bool isEqualityLinear() const override { return false; }
    bool isInequalityLinear() const override { return false; }

    // bool providesJacobian() const override { return false; }

    // TODO(roesmann): integral of squared functions is not always the same as the squared integral of functions
    bool isObjectiveLeastSquaresForm() const override { return false; }

    void precompute() override
    {
        assert(_integrator);
        assert(_dynamics);
        assert((int)_vertices.size() == 2 + _nc_interval + _nt_interval);

        assert(_shooting_node->getDimension() == _dynamics->getStateDimension());

        // get current values (we must avoid alias, hence we need to copy)
        if (_int_obj_dim > 0) _current[0] = 0;
        _current.segment(_int_obj_dim, _dyn_dim) = _shooting_node->values();
        int tail = _other_int_eq_dim + _int_ineq_dim;
        if (tail > 0) _current.tail(tail).setZero();

        _current_idx = 0;
        for (; _current_idx < _n_interval - 1; ++_current_idx)
        {
            if (_current_idx < _nc_interval) _ukvert = static_cast<VectorVertex*>(_vertices[_current_idx + 2]);
            if (_current_idx < _nt_interval) _dtvert = static_cast<ScalarVertex*>(_vertices[_current_idx + _nc_interval + 2]);

            assert(_ukvert->getDimension() == _dynamics->getInputDimension());

            // compute integral stuff
            // _integrator->integrate(_current, _dtvert->value(), _integrand, _values);
            _integrator->solveIVP(_current, _dtvert->value(), _integrand, _values);

            // TODO(roesmann) check if we can ensure to avoid alias in the integrator part to speed up simple integrators
            // first simple tests with replacing _values by _current led to different (wrong) results
            if (_current_idx < _n_interval - 2)
            {
                _current = _values;

                // also compute intermediate state dependend stuff
                if (_eval_intermediate_constr)
                {
                    const int k                          = _k + _current_idx;
                    Eigen::Ref<const Eigen::VectorXd> xk = _values.segment(_int_obj_dim, _dyn_dim);

                    if (_nonint_obj_dim > 0)
                    {
                        Eigen::Ref<Eigen::VectorXd> cost_segment = _nonint_obj_values.segment(_current_idx * _nonint_obj_dim, _nonint_obj_dim);
                        _stage_cost->computeConcatenatedNonIntegralStateTerms(k, xk, _ukvert->values(), _dtvert->value(), cost_segment, true);
                    }

                    if (_nonint_eq_dim > 0)
                    {
                        Eigen::Ref<Eigen::VectorXd> eq_segment = _nonint_eq_values.segment(_current_idx * _nonint_eq_dim, _nonint_eq_dim);
                        _stage_eq->computeConcatenatedNonIntegralStateTerms(k, xk, _ukvert->values(), _dtvert->value(), eq_segment, false);
                    }

                    if (_nonint_ineq_dim > 0)
                    {
                        Eigen::Ref<Eigen::VectorXd> ineq_segment = _nonint_eq_values.segment(_current_idx * _nonint_ineq_dim, _nonint_ineq_dim);
                        _stage_ineq->computeConcatenatedNonIntegralStateTerms(k, xk, _ukvert->values(), _dtvert->value(), ineq_segment, false);
                    }

                    // lower bounds
                    if (_dim_lb_x > 0)
                    {
                        int lb_idx = 0;
                        for (int j = 0; j < _shooting_node->getDimension(); ++j)
                        {
                            if (_shooting_node->hasFiniteLowerBound(j))
                            {
                                _lb_ineq_values(_current_idx * _dim_lb_x + lb_idx) = _shooting_node->lowerBound()[j] - xk[j];  // x_min - x <= 0
                                ++lb_idx;
                            }
                        }
                    }

                    // upper bounds
                    if (_dim_ub_x > 0)
                    {
                        int ub_idx = 0;
                        for (int j = 0; j < _shooting_node->getDimension(); ++j)
                        {
                            if (_shooting_node->hasFiniteUpperBound(j))
                            {
                                _ub_ineq_values(_current_idx * _dim_ub_x + ub_idx) = xk[j] - _shooting_node->upperBound()[j];  // x - x_max <= 0
                                ++ub_idx;
                            }
                        }
                    }
                }
            }
        }
    }

    void computeObjectiveValues(Eigen::Ref<Eigen::VectorXd> obj_values) override
    {
        if (_int_obj_dim > 0 && _nonint_obj_dim > 0)
        {
            obj_values[0] = _values[0] + _nonint_obj_values.sum();
        }
        else if (_int_obj_dim > 0)
        {
            obj_values[0] = _values[0];
        }
        else if (_nonint_obj_dim > 0)
        {
            obj_values[0] = _nonint_obj_values.sum();
        }
    }
    void computeEqualityValues(Eigen::Ref<Eigen::VectorXd> eq_values) override
    {
        // system dynamics result must comply with subsequent state
        const VectorVertex* x_next = static_cast<const VectorVertex*>(_vertices[1]);
        assert(x_next->values().size() == _dyn_dim);

        eq_values.head(_dyn_dim).noalias() = _values.segment(_int_obj_dim, _dyn_dim) - x_next->values();

        if (_other_int_eq_dim > 0) eq_values.segment(_dyn_dim, _other_int_eq_dim) = _values.segment(_int_obj_dim + _dyn_dim, _other_int_eq_dim);

        if (_eval_intermediate_constr)
        {
            assert(_num_intermediate_states > 0);
            if (_nonint_eq_dim > 0) eq_values.segment(_dyn_dim + _other_int_eq_dim, _num_intermediate_states * _nonint_eq_dim) = _nonint_eq_values;
        }

        assert(eq_values.size() == _dyn_dim + _other_int_eq_dim + _num_intermediate_states * _nonint_eq_dim);
    }
    void computeInequalityValues(Eigen::Ref<Eigen::VectorXd> ineq_values) override
    {
        if (_int_ineq_dim > 0) ineq_values.head(_int_ineq_dim) = _values.tail(_int_ineq_dim);

        if (_eval_intermediate_constr)
        {
            assert(_num_intermediate_states > 0);
            int idx = _int_ineq_dim;
            if (_nonint_ineq_dim > 0)
            {
                ineq_values.segment(idx, _num_intermediate_states * _nonint_ineq_dim) = _nonint_ineq_values;
                idx += _num_intermediate_states * _nonint_ineq_dim;
            }

            // bounds
            if (_dim_lb_x > 0)
            {
                ineq_values.segment(idx, _num_intermediate_states * _dim_lb_x) = _lb_ineq_values;
                // idx += _num_intermediate_states * _nonint_ineq_dim;
            }
            if (_dim_ub_x > 0)
            {
                ineq_values.tail(_num_intermediate_states * _dim_ub_x) = _ub_ineq_values;
                // idx += _num_intermediate_states * _nonint_ineq_dim;
            }
        }

        assert(ineq_values.size() == _int_ineq_dim + _num_intermediate_states * (_nonint_ineq_dim + _dim_lb_x + _dim_ub_x));
    }

    void setIntegrator(NumericalIntegratorExplicitInterface::Ptr integrator) { _integrator = integrator; }

 protected:
    void activateIntermediateStateCostAndConstraints()
    {
        if (!_eval_intermediate_constr) return;
        _num_intermediate_states = _n_interval - 2;  // no start and state
        if (_num_intermediate_states < 1)
        {
            _eval_intermediate_constr = false;
            return;
        }

        // TODO(roesmann): we assume for now that all terms have the same dimension _k!!!
        _nonint_obj_dim  = _stage_cost ? _stage_cost->getConcatenatedNonIntegralStateTermDimension(_k, true) : 0;
        _nonint_eq_dim   = _stage_eq ? _stage_eq->getConcatenatedNonIntegralStateTermDimension(_k, false) : 0;
        _nonint_ineq_dim = _stage_ineq ? _stage_ineq->getConcatenatedNonIntegralStateTermDimension(_k, false) : 0;

        // consider number of finite bounds for all vertices
        // we apply the bounds of vertex x1
        const VectorVertex* x1 = static_cast<const VectorVertex*>(_vertices[0]);
        _dim_lb_x              = x1->getNumberFiniteLowerBounds(false);
        _dim_ub_x              = x1->getNumberFiniteUpperBounds(false);

        // update overall dimensions
        _total_dim_obj  = (_int_obj_dim > 0 || _nonint_obj_dim > 0) ? 1 : 0;
        _total_dim_eq   = _other_int_eq_dim + _dyn_dim + _num_intermediate_states * _nonint_eq_dim;
        _total_dim_ineq = _int_ineq_dim + _num_intermediate_states * (_nonint_ineq_dim + _dim_lb_x + _dim_ub_x);

        // allocate caches
        _nonint_obj_values.resize(_num_intermediate_states * _nonint_obj_dim);
        _nonint_eq_values.resize(_num_intermediate_states * _nonint_eq_dim);
        _nonint_ineq_values.resize(_num_intermediate_states * _nonint_ineq_dim);
        _lb_ineq_values.resize(_num_intermediate_states * _dim_lb_x);
        _ub_ineq_values.resize(_num_intermediate_states * _dim_ub_x);
    }

    void configureIntegrand()
    {
        _integrand = [&](const Eigen::VectorXd& current, Eigen::Ref<Eigen::VectorXd> result) {
            assert(_ukvert);
            assert(_dtvert);

            const int k_integrand = _k + _current_idx;

            int idx = 0;
            _dynamics->dynamics(current.segment(_int_obj_dim, _dyn_dim), _ukvert->values(), result.segment(_int_obj_dim, _dyn_dim));
            if (_int_obj_dim > 0)
            {
                assert(_int_obj_dim == 1);
                _stage_cost->computeIntegralStateControlTerm(k_integrand, current.segment(_int_obj_dim, _dyn_dim), _ukvert->values(), result.head(1));
                idx += 1;
            }
            idx += _dyn_dim;
            if (_other_int_eq_dim > 0)
            {
                _stage_eq->computeIntegralStateControlTerm(k_integrand, current.segment(_int_obj_dim, _dyn_dim), _ukvert->values(),
                                                           result.segment(idx, _other_int_eq_dim));
                idx += _other_int_eq_dim;
            }
            if (_int_ineq_dim)
            {
                _stage_ineq->computeIntegralStateControlTerm(k_integrand, current.segment(_int_obj_dim, _dyn_dim), _ukvert->values(),
                                                             result.segment(idx, _int_ineq_dim));
                idx += _int_ineq_dim;
            }
            assert(idx == result.size());
        };
    }

 private:
    const VectorVertex* _shooting_node = nullptr;
    VectorVertex* _ukvert              = nullptr;
    ScalarVertex* _dtvert              = nullptr;

    int _dyn_dim = 0;

    int _int_obj_dim      = 0;
    int _other_int_eq_dim = 0;
    int _int_ineq_dim     = 0;

    int _nonint_obj_dim  = 0;
    int _nonint_eq_dim   = 0;
    int _nonint_ineq_dim = 0;

    int _dim_lb_x = 0;
    int _dim_ub_x = 0;

    int _total_dim_obj  = 0;
    int _total_dim_eq   = 0;
    int _total_dim_ineq = 0;

    int _n_interval              = 0;
    int _nc_interval             = 0;
    int _nt_interval             = 0;
    int _num_intermediate_states = 0;
    int _current_idx             = 0;

    bool _eval_intermediate_constr = false;

    Eigen::VectorXd _current;
    Eigen::VectorXd _values;  // TODO(roesmann): if integrators are alias safe we can remove _values in favor of just _current

    Eigen::VectorXd _nonint_obj_values;
    Eigen::VectorXd _nonint_ineq_values;
    Eigen::VectorXd _nonint_eq_values;
    Eigen::VectorXd _lb_ineq_values;
    Eigen::VectorXd _ub_ineq_values;

    std::function<void(const Eigen::VectorXd&, Eigen::Ref<Eigen::VectorXd>)> _integrand;

    StageCost::ConstPtr _stage_cost;
    StageEqualityConstraint::ConstPtr _stage_eq;
    StageInequalityConstraint::ConstPtr _stage_ineq;

    SystemDynamicsInterface::Ptr _dynamics;
    NumericalIntegratorExplicitInterface::Ptr _integrator;

    int _k = 0;

 public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

class MSDynamicsOnlyMultiControlsEdge : public Edge<>
{
 public:
    using Ptr  = std::shared_ptr<MSDynamicsOnlyMultiControlsEdge>;
    using UPtr = std::unique_ptr<MSDynamicsOnlyMultiControlsEdge>;

    explicit MSDynamicsOnlyMultiControlsEdge(SystemDynamicsInterface::Ptr dynamics, int num_controls) : Edge<>(3 + num_controls), _dynamics(dynamics)
    {
    }

    void finalize()
    {
        _x1 = static_cast<const VectorVertex*>(_vertices[0]);
        _x2 = static_cast<const VectorVertex*>(_vertices[1]);
        _dt = static_cast<const ScalarVertex*>(_vertices[2]);

        _f.resize(_dynamics->getStateDimension());
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
        assert(_integrator);
        assert(_dynamics);
        assert(_x1->getDimension() == _dynamics->getStateDimension());
        assert(_x2->getDimension() == _dynamics->getStateDimension());
        assert(_vertices.size() > 4);

        const VectorVertex* u = static_cast<const VectorVertex*>(_vertices[3]);
        _integrator->solveIVP(_x1->values(), u->values(), _dt->value(), *_dynamics, _f);

        for (int i = 4; i < (int)_vertices.size(); ++i)
        {
            u = static_cast<const VectorVertex*>(_vertices[i]);
            _integrator->solveIVP(_f, u->values(), _dt->value(), *_dynamics, values);
            _f = values;  // TODO(roesmann): avoid this by checking the indices maybe?!
        }
        values -= _x2->values();
    }

    void setIntegrator(NumericalIntegratorExplicitInterface::Ptr integrator) { _integrator = integrator; }

 private:
    SystemDynamicsInterface::Ptr _dynamics;
    NumericalIntegratorExplicitInterface::Ptr _integrator;

    Eigen::VectorXd _f;

    const VectorVertex* _x1 = nullptr;
    const VectorVertex* _x2 = nullptr;
    const ScalarVertex* _dt = nullptr;

 public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

class MSDynamicsOnlyMultiControlsMultiDtsEdge : public Edge<>
{
 public:
    using Ptr  = std::shared_ptr<MSDynamicsOnlyMultiControlsMultiDtsEdge>;
    using UPtr = std::unique_ptr<MSDynamicsOnlyMultiControlsMultiDtsEdge>;

    explicit MSDynamicsOnlyMultiControlsMultiDtsEdge(SystemDynamicsInterface::Ptr dynamics, VectorVertex& x1, VectorVertex& x2, int num_controls)
        : Edge<>(2 + 2 * num_controls), _dynamics(dynamics)
    {
        setVertex(0, x1);
        setVertex(1, x2);
        _x1 = static_cast<const VectorVertex*>(_vertices[0]);
        _x2 = static_cast<const VectorVertex*>(_vertices[1]);

        _f.resize(dynamics->getStateDimension());
    }

    void finalize() {}

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
        assert(_integrator);
        assert(_dynamics);
        assert(_x1->getDimension() == _dynamics->getStateDimension());
        assert(_x2->getDimension() == _dynamics->getStateDimension());
        assert(_vertices.size() > 5);

        const VectorVertex* u  = static_cast<const VectorVertex*>(_vertices[2]);
        const ScalarVertex* dt = static_cast<const ScalarVertex*>(_vertices[3]);
        _integrator->solveIVP(_x1->values(), u->values(), dt->value(), *_dynamics, values);

        for (int i = 4; i < (int)_vertices.size() - 2; i = i + 2)
        {
            u  = static_cast<const VectorVertex*>(_vertices[i]);
            dt = static_cast<const ScalarVertex*>(_vertices[i + 1]);
            _integrator->solveIVP(values, u->values(), dt->value(), *_dynamics, _f);
        }
        u  = static_cast<const VectorVertex*>(_vertices[(int)_vertices.size() - 2]);
        dt = static_cast<const ScalarVertex*>(_vertices.back());
        _integrator->solveIVP(_f, u->values(), dt->value(), *_dynamics, values);
        values -= _x2->values();
    }

    void setIntegrator(NumericalIntegratorExplicitInterface::Ptr integrator) { _integrator = integrator; }

 private:
    SystemDynamicsInterface::Ptr _dynamics;
    NumericalIntegratorExplicitInterface::Ptr _integrator;

    Eigen::VectorXd _f;

    const VectorVertex* _x1 = nullptr;
    const VectorVertex* _x2 = nullptr;

 public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

}  // namespace corbo

#endif  // SRC_OPTIMAL_CONTROL_INCLUDE_CORBO_OPTIMAL_CONTROL_STRUCTURED_OCP_EDGES_MULTIPLE_SHOOTING_EDGES_H_
