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

#ifndef SRC_OPTIMAL_CONTROL_INCLUDE_CORBO_OPTIMAL_CONTROL_STRUCTURED_OCP_EDGES_COLLOCATION_EDGES_H_
#define SRC_OPTIMAL_CONTROL_INCLUDE_CORBO_OPTIMAL_CONTROL_STRUCTURED_OCP_EDGES_COLLOCATION_EDGES_H_

#include <corbo-optimization/hyper_graph/edge.h>

#include <corbo-numerics/dynamics_eval_interface.h>
#include <corbo-numerics/quadrature_interface.h>
#include <corbo-optimal-control/functions/stage_functions.h>
#include <corbo-optimization/hyper_graph/scalar_vertex.h>
#include <corbo-optimization/hyper_graph/vector_vertex.h>
#include <corbo-systems/system_dynamics_interface.h>

#include <functional>
#include <memory>

namespace corbo {

class QuadratureCollocationDynamicsOnly : public Edge<>
{
 public:
    using Ptr  = std::shared_ptr<QuadratureCollocationDynamicsOnly>;
    using UPtr = std::unique_ptr<QuadratureCollocationDynamicsOnly>;

    explicit QuadratureCollocationDynamicsOnly(SystemDynamicsInterface::Ptr dynamics, QuadratureCollocationInterface::Ptr quadrature)
        : Edge<>(5 + quadrature->getNumIntermediateStates() + quadrature->getNumIntermediateControls()),  // x1, u1, x2, u2, dt + intermediate
          _dynamics(dynamics),
          _collocation(quadrature),
          _num_intermediate_x(quadrature->getNumIntermediateStates()),
          _num_intermediate_u(quadrature->getNumIntermediateControls())
    {
    }

    virtual ~QuadratureCollocationDynamicsOnly() = default;

    void finalize()
    {
        _dim_x = _dynamics->getStateDimension();

        assert(_vertices.size() == 5 + _num_intermediate_x + _num_intermediate_u);

        _x1_points.clear();
        _u1_points.clear();

        _x2 = &static_cast<const VectorVertex*>(_vertices[0])->values();
        _u2 = &static_cast<const VectorVertex*>(_vertices[1])->values();
        _dt = static_cast<const ScalarVertex*>(_vertices[2]);

        int idx = 3;
        for (int i = 0; i < _num_intermediate_x + 1; ++i)
        {
            _x1_points.push_back(&static_cast<const VectorVertex*>(_vertices[idx++])->values());
        }
        for (int i = 0; i < _num_intermediate_u + 1; ++i)
        {
            _u1_points.push_back(&static_cast<const VectorVertex*>(_vertices[idx++])->values());
        }
    }

    // implements interface method
    int getDimension() const override { return _dim_x; }
    // implement in child class:
    bool isLinear() const override { return false; }

    // bool providesJacobian() const override { return false; }

    bool isLeastSquaresForm() const override { return false; }

    void computeValues(Eigen::Ref<Eigen::VectorXd> values) override
    {
        assert(_collocation);
        assert(_dynamics);
        assert(!_x1_points.empty());
        assert(!_u1_points.empty());
        assert(_x2 && _u2 && _dt);
        assert(_x1_points.front()->size() == _dynamics->getStateDimension());
        assert(_u1_points.front()->size() == _dynamics->getInputDimension());
        assert(_x2->size() == _dynamics->getStateDimension());
        assert(_u2->size() == _dynamics->getInputDimension());

        // perform quadrature with given control mid points:
        //        _collocation->quadrature(_x1->values(), _u1->values(), _x2->values(), _u2->values(), _dt->value(), *_dynamics, _intermediate_u,
        //        values,
        //                                 _intermediate_x);
        _collocation->quadrature(_x1_points, _u1_points, *_u2, *_x2, _dt->value(), *_dynamics, values);

        values.noalias() -= (*_x2 - *_x1_points.front());
    }

 protected:
 private:
    SystemDynamicsInterface::Ptr _dynamics;
    QuadratureCollocationInterface::Ptr _collocation;

    int _num_intermediate_x = 0;
    int _num_intermediate_u = 0;
    int _dim_x              = 0;

    std::vector<const Eigen::VectorXd*> _x1_points;
    std::vector<const Eigen::VectorXd*> _u1_points;
    const Eigen::VectorXd* _x2 = nullptr;
    const Eigen::VectorXd* _u2 = nullptr;

    const ScalarVertex* _dt = nullptr;

 public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

class QuadratureCollocationEdge : public MixedEdge<>
{
 public:
    using Ptr  = std::shared_ptr<QuadratureCollocationEdge>;
    using UPtr = std::unique_ptr<QuadratureCollocationEdge>;

    explicit QuadratureCollocationEdge(SystemDynamicsInterface::Ptr dynamics, QuadratureCollocationInterface::Ptr quadrature,
                                       StageCost::Ptr stage_cost, StageEqualityConstraint::Ptr stage_eq, StageInequalityConstraint::Ptr stage_ineq,
                                       bool eval_intermediate_constr, int k)
        : MixedEdge<>(5 + quadrature->getNumIntermediateStates() + quadrature->getNumIntermediateControls()),
          _dynamics(dynamics),
          _collocation(quadrature),
          _stage_cost(stage_cost),
          _stage_eq(stage_eq),
          _stage_ineq(stage_ineq),
          _k(k),
          _eval_intermediate_constr(eval_intermediate_constr),
          _num_intermediate_x(quadrature->getNumIntermediateStates()),
          _num_intermediate_u(quadrature->getNumIntermediateControls())
    {
        assert(quadrature->isIntermediateControlSubjectToOptim());
    }

    void finalize()
    {
        int x_dim = _dynamics->getStateDimension();

        _dim_obj      = (_stage_cost && _stage_cost->getIntegralStateControlTermDimension(_k) > 0) ? 1 : 0;  // we have only a single obejtive value
        _dim_int_eq   = _stage_eq ? _stage_eq->getIntegralStateControlTermDimension(_k) : 0;
        _dim_eq       = x_dim + _dim_int_eq;
        _dim_int_ineq = _stage_ineq ? _stage_ineq->getIntegralStateControlTermDimension(_k) : 0;
        _dim_ineq     = _dim_int_ineq;
        _dim_dyn      = _dynamics->getStateDimension();

        assert(_vertices.size() == 5 + _num_intermediate_x + _num_intermediate_u);

        _x1_points.clear();
        _u1_points.clear();

        _x2 = &static_cast<const VectorVertex*>(_vertices[0])->values();
        _u2 = &static_cast<const VectorVertex*>(_vertices[1])->values();
        _dt = static_cast<const ScalarVertex*>(_vertices[2]);

        int idx = 3;
        _x1     = static_cast<const VectorVertex*>(_vertices[idx]);
        for (int i = 0; i < _num_intermediate_x + 1; ++i)
        {
            _x1_points.push_back(&static_cast<const VectorVertex*>(_vertices[idx++])->values());
        }
        for (int i = 0; i < _num_intermediate_u + 1; ++i)
        {
            _u1_points.push_back(&static_cast<const VectorVertex*>(_vertices[idx++])->values());
        }

        _quadrature_values.resize(_dim_obj + _dim_dyn + _dim_int_eq + _dim_int_ineq);

        configureIntegrand();
        PRINT_WARNING_COND_NAMED(_eval_intermediate_constr, "_eval_intermediate_constr not yet implemented");
        // activateIntermediateConstraints();
    }

    virtual ~QuadratureCollocationEdge() = default;

    // implements interface method
    int getObjectiveDimension() const override { return _dim_obj; }
    // implements interface method
    int getEqualityDimension() const override { return _dim_eq; }
    int getInequalityDimension() const override { return _dim_ineq; }

    // implement in child class:
    bool isObjectiveLinear() const override { return false; }
    bool isEqualityLinear() const override { return false; }
    bool isInequalityLinear() const override { return false; }

    // bool providesJacobian() const override { return false; }

    bool isObjectiveLeastSquaresForm() const override { return false; }

    void precompute() override
    {
        assert(_collocation);
        assert(_dynamics);
        assert(!_x1_points.empty());
        assert(!_u1_points.empty());
        assert(_x2 && _u2 && _dt);
        assert(_x1_points.front()->size() == _dynamics->getStateDimension());
        assert(_u1_points.front()->size() == _dynamics->getInputDimension());
        assert(_x2->size() == _dynamics->getStateDimension());
        assert(_u2->size() == _dynamics->getInputDimension());

        _collocation->quadrature(_x1_points, _u1_points, *_u2, *_x2, _dt->value(), _integrand, _quadrature_values);
    }

    void computeObjectiveValues(Eigen::Ref<Eigen::VectorXd> obj_values) override { obj_values = _quadrature_values.head(1); }

    void computeEqualityValues(Eigen::Ref<Eigen::VectorXd> eq_values) override
    {
        eq_values.head(_dim_dyn).noalias() = *_x2 - _x1->values() - _quadrature_values.segment(_dim_obj, _dim_dyn);
        if (_dim_int_eq > 0)
        {
            eq_values.segment(_dim_dyn, _dim_int_eq).noalias() = _quadrature_values.segment(_dim_obj + _dim_dyn, _dim_int_eq);
        }

        // in case we consider intermediate constraints
        //        if (_dim_nonint_eq > 0)
        //        {
        //            int cur_idx = _dim_dyn + _dim_int_eq;
        //            for (int i = 0; i < _intermediate_x; ++i)
        //            {
        //                _stage_eq->computeConcatenatedNonIntegralStateTerms(_k, *_intermediate_x[i], *_intermediate_u[i], _dt->value(),
        //                                                                    eq_values.segment(cur_idx, _dim_nonint_eq));
        //                cur_idx += _dim_nonint_eq;
        //            }
        //            assert(cur_idx == getEqualityDimension());
        //        }
    }
    void computeInequalityValues(Eigen::Ref<Eigen::VectorXd> ineq_values) override
    {
        if (_dim_int_ineq > 0)
        {
            ineq_values = _quadrature_values.tail(_dim_int_ineq);
        }

        // in case we consider intermediate constraints, the following dimensions are > 0
        //        int cur_idx = _dim_int_eq;
        //        if (_dim_nonint_ineq > 0)
        //        {
        //            for (int i = 0; i < (int)_intermediate_x.size(); ++i)
        //            {
        //                _stage_ineq->computeConcatenatedNonIntegralStateTerms(_k, *_intermediate_x[i], *_intermediate_u[i], _dt->value(),
        //                                                                      ineq_values.segment(cur_idx, _dim_nonint_ineq));
        //                cur_idx += _dim_nonint_eq;
        //            }
        //        }
        //
        //        // lower bounds
        //        if (_dim_lb_x > 0)
        //        {
        //            for (int i = 0; i < (int)_intermediate_x.size(); ++i)
        //            {
        //                for (int j = 0; j < _intermediate_x[i]->size(); ++j)
        //                {
        //                    if (_x1->hasFiniteLowerBound(j))
        //                    {
        //                        ineq_values(cur_idx++) = _x1->lowerBound()[j] - (*_intermediate_x[i])[j];  // x_min - x <= 0
        //                    }
        //                }
        //            }
        //        }
        //
        //        // upper bounds
        //        if (_dim_ub_x > 0)
        //        {
        //            for (int i = 0; i < (int)_intermediate_x.size(); ++i)
        //            {
        //                for (int j = 0; j < _intermediate_x[i]->size(); ++j)
        //                {
        //                    if (_x1->hasFiniteUpperBound(j))
        //                    {
        //                        ineq_values(cur_idx++) = (*_intermediate_x[i])[j] - _x1->upperBound()[j];  // x - x_max <= 0
        //                    }
        //                }
        //            }
        //        }
        //        assert(cur_idx == getInequalityDimension());
    }

 protected:
    //    void activateIntermediateConstraints()
    //    {
    //        // get number inequalities in case _eval_intermediate_constr is turned on
    //        // TODO(roesmann) currently we only evaluate intermediate states, since u is usually just a linear function
    //        if (_eval_intermediate_constr)
    //        {
    //            _dim_nonint_eq        = _stage_eq ? _stage_eq->getConcatenatedNonIntegralStateTermDimension(_k) : 0;
    //            _dim_nonint_ineq      = _stage_ineq ? _stage_ineq->getConcatenatedNonIntegralStateTermDimension(_k) : 0;
    //
    //            // consider number of finite bounds for all vertices
    //            // we apply the bounds of vertex x1
    //            _dim_lb_x = _x1_points[0]->getNumberFiniteLowerBounds(false);
    //            _dim_ub_x = _x1_points[0]->getNumberFiniteUpperBounds(false);
    //
    //            // update overall dimensions
    //            _dim_eq   = _dim_int_eq + _dim_dyn + _num_intermediate_x * _dim_nonint_eq;
    //            _dim_ineq = _dim_int_ineq + _num_intermediate_x * (_dim_nonint_ineq + _dim_lb_x + _dim_ub_x);
    //        }
    //    }

    void configureIntegrand()
    {
        _integrand = [&](const Eigen::VectorXd& x, const Eigen::VectorXd& u, Eigen::Ref<Eigen::VectorXd> result) {
            int idx = 0;
            _dynamics->dynamics(x, u, result.segment(_dim_obj, _dim_dyn));

            if (_dim_obj > 0)
            {
                assert(_dim_obj == 1);
                _stage_cost->computeIntegralStateControlTerm(_k, x, u, result.head(1));
                idx += 1;
            }
            idx += _dim_dyn;

            if (_dim_int_eq > 0)
            {
                _stage_eq->computeIntegralStateControlTerm(_k, x, u, result.segment(idx, _dim_int_eq));
                idx += _dim_int_eq;
            }

            if (_dim_int_ineq > 0)
            {
                _stage_ineq->computeIntegralStateControlTerm(_k, x, u, result.segment(idx, _dim_int_eq));
                idx += _dim_int_eq;
            }
            assert(idx == result.size());
        };
    }

 private:
    SystemDynamicsInterface::Ptr _dynamics;
    QuadratureCollocationInterface::Ptr _collocation;

    StageCost::Ptr _stage_cost;
    StageEqualityConstraint::Ptr _stage_eq;
    StageInequalityConstraint::Ptr _stage_ineq;

    int _dim_obj      = 1;
    int _dim_eq       = 0;
    int _dim_int_eq   = 0;
    int _dim_ineq     = 0;
    int _dim_int_ineq = 0;
    int _dim_dyn      = 0;

    int _dim_nonint_eq   = 0;
    int _dim_nonint_ineq = 0;

    int _dim_lb_x = 0;
    int _dim_ub_x = 0;

    int _k = 0;

    bool _eval_intermediate_constr = false;

    int _num_intermediate_x;
    int _num_intermediate_u;

    Eigen::VectorXd _quadrature_values;

    std::function<void(const Eigen::VectorXd&, const Eigen::VectorXd&, Eigen::Ref<Eigen::VectorXd>)> _integrand;

    const VectorVertex* _x1 = nullptr;
    std::vector<const Eigen::VectorXd*> _x1_points;
    std::vector<const Eigen::VectorXd*> _u1_points;
    const Eigen::VectorXd* _x2 = nullptr;
    const Eigen::VectorXd* _u2 = nullptr;

    const ScalarVertex* _dt = nullptr;

 public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

class CompressedCollocationEdge : public Edge<VectorVertex, VectorVertex, ScalarVertex, VectorVertex, VectorVertex>
{
 public:
    using Ptr  = std::shared_ptr<CompressedCollocationEdge>;
    using UPtr = std::unique_ptr<CompressedCollocationEdge>;

    explicit CompressedCollocationEdge(SystemDynamicsInterface::Ptr dynamics, VectorVertex& x1, VectorVertex& u1, ScalarVertex& dt, VectorVertex& u2,
                                       VectorVertex& x2)
        : Edge<VectorVertex, VectorVertex, ScalarVertex, VectorVertex, VectorVertex>(x1, u1, dt, u2, x2), _dynamics(dynamics)
    {
        _x1 = static_cast<const VectorVertex*>(_vertices[0]);
        _u1 = static_cast<const VectorVertex*>(_vertices[1]);
        _dt = static_cast<const ScalarVertex*>(_vertices[2]);
        _u2 = static_cast<const VectorVertex*>(_vertices[3]);
        _x2 = static_cast<const VectorVertex*>(_vertices[4]);
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
        assert(_collocation);
        assert(_dynamics);
        assert(_x1->getDimension() == _dynamics->getStateDimension());
        assert(_u1->getDimension() == _dynamics->getInputDimension());
        assert(_x2->getDimension() == _dynamics->getStateDimension());
        assert(_u2->getDimension() == _dynamics->getInputDimension());

        _collocation->computeEqualityConstraint(_x1->values(), _u1->values(), _x2->values(), _u2->values(), _dt->value(), *_dynamics, values);
    }

    void setCollocationMethod(QuadratureCollocationInterface::Ptr collocation)
    {
        assert(!collocation->isIntermediateControlSubjectToOptim() && collocation->isSupportingCompressedStatesMode());
        _collocation = collocation;
    }

 private:
    SystemDynamicsInterface::Ptr _dynamics;
    QuadratureCollocationInterface::Ptr _collocation;

    const VectorVertex* _x1 = nullptr;
    const VectorVertex* _u1 = nullptr;
    const VectorVertex* _x2 = nullptr;
    const VectorVertex* _u2 = nullptr;
    const ScalarVertex* _dt = nullptr;

 public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

class CompressedCollocationMultipleControlsEdge : public Edge<>
{
 public:
    using Ptr  = std::shared_ptr<CompressedCollocationMultipleControlsEdge>;
    using UPtr = std::unique_ptr<CompressedCollocationMultipleControlsEdge>;

    explicit CompressedCollocationMultipleControlsEdge(SystemDynamicsInterface::Ptr dynamics, QuadratureCollocationInterface::Ptr quadrature)
        : Edge<>(5 + quadrature->getNumIntermediatePoints()),
          _dynamics(dynamics),
          _collocation(quadrature),
          _num_intermediate_points(quadrature->getNumIntermediatePoints())
    {
        assert(quadrature->isIntermediateControlSubjectToOptim());
    }

    virtual ~CompressedCollocationMultipleControlsEdge() { clearInternalBuffer(); }

    void finalize()
    {
        _dim_x = _dynamics->getStateDimension();

        _x1 = static_cast<const VectorVertex*>(_vertices[0]);
        _u1 = static_cast<const VectorVertex*>(_vertices[1]);
        _dt = static_cast<const ScalarVertex*>(_vertices[2]);
        _u2 = static_cast<const VectorVertex*>(_vertices[3]);
        _x2 = static_cast<const VectorVertex*>(_vertices[4]);

        int idx = 5;
        _intermediate_u.clear();
        for (int i = 0; i < _num_intermediate_points; ++i)
        {
            VectorVertex* interm_u = static_cast<VectorVertex*>(_vertices[idx + i]);
            _intermediate_u.push_back(&interm_u->values());
        }

        clearInternalBuffer();

        for (int i = 0; i < _num_intermediate_points; ++i)
        {
            _intermediate_x.push_back(new Eigen::VectorXd(_dim_x));
        }
    }

    // implements interface method
    int getDimension() const override { return _dim_x; }
    // implement in child class:
    bool isLinear() const override { return false; }

    // bool providesJacobian() const override { return false; }

    bool isLeastSquaresForm() const override { return false; }

    void computeValues(Eigen::Ref<Eigen::VectorXd> values) override
    {
        assert(_collocation);
        assert(_dynamics);
        assert(_x1->getDimension() == _dynamics->getStateDimension());
        assert(_u1->getDimension() == _dynamics->getInputDimension());
        assert(_x2->getDimension() == _dynamics->getStateDimension());
        assert(_u2->getDimension() == _dynamics->getInputDimension());

        // perform quadrature with given control mid points:
        _collocation->quadrature(_x1->values(), _u1->values(), _x2->values(), _u2->values(), _dt->value(), *_dynamics, _intermediate_u, values,
                                 _intermediate_x);

        values.noalias() -= (_x2->values() - _x1->values());
    }

 protected:
    void clearInternalBuffer()
    {
        for (int i = 0; i < _intermediate_x.size(); ++i) delete _intermediate_x[i];
        _intermediate_x.clear();
    }

 private:
    SystemDynamicsInterface::Ptr _dynamics;
    QuadratureCollocationInterface::Ptr _collocation;

    int _num_intermediate_points;
    int _dim_x = 0;

    std::vector<Eigen::VectorXd*> _intermediate_x;
    std::vector<Eigen::VectorXd*> _intermediate_u;

    const VectorVertex* _x1 = nullptr;
    const VectorVertex* _u1 = nullptr;
    const VectorVertex* _x2 = nullptr;
    const VectorVertex* _u2 = nullptr;
    const ScalarVertex* _dt = nullptr;

 public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

class ConstControlCombinedCompressedCollocationEdge : public MixedEdge<VectorVertex, VectorVertex, ScalarVertex, VectorVertex>
{
 public:
    using Ptr  = std::shared_ptr<ConstControlCombinedCompressedCollocationEdge>;
    using UPtr = std::unique_ptr<ConstControlCombinedCompressedCollocationEdge>;

    explicit ConstControlCombinedCompressedCollocationEdge(SystemDynamicsInterface::Ptr dynamics, StageCost::Ptr stage_cost,
                                                           StageEqualityConstraint::Ptr stage_eq, StageInequalityConstraint::Ptr stage_ineq,
                                                           bool eval_intermediate_constr, int k, VectorVertex& x1, VectorVertex& u1, ScalarVertex& dt,
                                                           VectorVertex& x2)
        : MixedEdge<VectorVertex, VectorVertex, ScalarVertex, VectorVertex>(x1, u1, dt, x2),
          _dynamics(dynamics),
          _stage_cost(stage_cost),
          _stage_eq(stage_eq),
          _stage_ineq(stage_ineq),
          _k(k),
          _eval_intermediate_constr(eval_intermediate_constr)
    {
        int x_dim = _dynamics->getStateDimension();
        _dynamics_quadrature.resize(x_dim);

        _dim_obj      = (_stage_cost && _stage_cost->getIntegralStateControlTermDimension(_k) > 0) ? 1 : 0;  // we have only a single obejtive value
        _dim_int_eq   = _stage_eq ? _stage_eq->getIntegralStateControlTermDimension(k) : 0;
        _dim_eq       = x_dim + _dim_int_eq;
        _dim_int_ineq = _stage_ineq ? _stage_ineq->getIntegralStateControlTermDimension(k) : 0;
        _dim_ineq     = _dim_int_ineq;
        _dim_dyn      = _dynamics->getStateDimension();

        _x1 = static_cast<const VectorVertex*>(_vertices[0]);
        _u1 = static_cast<const VectorVertex*>(_vertices[1]);
        _dt = static_cast<const ScalarVertex*>(_vertices[2]);
        _x2 = static_cast<const VectorVertex*>(_vertices[3]);
    }

    virtual ~ConstControlCombinedCompressedCollocationEdge()
    {
        for (int i = 0; i < (int)_intermediate_x.size(); ++i) delete _intermediate_x[i];
    }

    // implements interface method
    int getObjectiveDimension() const override { return _dim_obj; }
    // implements interface method
    int getEqualityDimension() const override { return _dim_eq; }
    int getInequalityDimension() const override { return _dim_ineq; }

    // implement in child class:
    bool isObjectiveLinear() const override { return false; }
    bool isEqualityLinear() const override { return false; }
    bool isInequalityLinear() const override { return false; }

    // bool providesJacobian() const override { return false; }

    bool isObjectiveLeastSquaresForm() const override { return false; }

    void precompute() override
    {
        assert(_collocation);
        assert(_dynamics);

        assert(_x1->getDimension() == _dynamics->getStateDimension());
        assert(_u1->getDimension() == _dynamics->getInputDimension());
        assert(_x2->getDimension() == _dynamics->getStateDimension());

        _collocation->quadrature(_x1->values(), _u1->values(), _x2->values(), _dt->value(), *_dynamics, _dynamics_quadrature, _intermediate_x);
    }

    void computeObjectiveValues(Eigen::Ref<Eigen::VectorXd> obj_values) override
    {
        // evaluate quadrature formula with previously computed intermediate states
        auto integrand = [this](const Eigen::VectorXd& x, Eigen::Ref<Eigen::VectorXd> result) {
            _stage_cost->computeIntegralStateControlTerm(_k, x, _u1->values(), result);
        };
        _collocation->quadrature(_x1->values(), _x2->values(), _dt->value(), _intermediate_x, integrand, obj_values);
    }

    void computeEqualityValues(Eigen::Ref<Eigen::VectorXd> eq_values) override
    {
        eq_values.head(_dynamics_quadrature.size()).noalias() = _x2->values() - _x1->values() - _dynamics_quadrature;
        if (_dim_int_eq > 0)
        {
            // evaluate quadrature formula with previously computed intermediate states
            auto integrand = [this](const Eigen::VectorXd& x, Eigen::Ref<Eigen::VectorXd> result) {
                _stage_eq->computeIntegralStateControlTerm(_k, x, _u1->values(), result);
            };
            _collocation->quadrature(_x1->values(), _x2->values(), _dt->value(), _intermediate_x, integrand,
                                     eq_values.segment(_dim_dyn, _dim_int_eq));
        }

        // in case we consider intermediate constraints
        if (_dim_nonint_eq > 0)
        {
            int cur_idx = _dim_dyn + _dim_int_eq;

            for (int i = 0; i < (int)_intermediate_x.size(); ++i)
            {
                _stage_eq->computeConcatenatedNonIntegralStateTerms(_k, *_intermediate_x[i], _u1->values(), _dt->value(),
                                                                    eq_values.segment(cur_idx, _dim_nonint_eq));
                cur_idx += _dim_nonint_eq;
            }
            assert(cur_idx == getEqualityDimension());
        }
    }
    void computeInequalityValues(Eigen::Ref<Eigen::VectorXd> ineq_values) override
    {
        if (_dim_int_ineq > 0)
        {
            // evaluate quadrature formula with previously computed intermediate states
            auto integrand = [this](const Eigen::VectorXd& x, Eigen::Ref<Eigen::VectorXd> result) {
                _stage_ineq->computeIntegralStateControlTerm(_k, x, _u1->values(), result);
            };
            _collocation->quadrature(_x1->values(), _x2->values(), _dt->value(), _intermediate_x, integrand, ineq_values);
        }

        // in case we consider intermediate constraints, the following dimensions are > 0
        int cur_idx = _dim_int_eq;
        if (_dim_nonint_ineq > 0)
        {
            for (int i = 0; i < (int)_intermediate_x.size(); ++i)
            {
                _stage_ineq->computeConcatenatedNonIntegralStateTerms(_k, *_intermediate_x[i], _u1->values(), _dt->value(),
                                                                      ineq_values.segment(cur_idx, _dim_nonint_ineq));
                cur_idx += _dim_nonint_eq;
            }
        }

        // lower bounds
        if (_dim_lb_x > 0)
        {
            for (int i = 0; i < _intermediate_x.size(); ++i)
            {
                for (int j = 0; j < _intermediate_x[i]->size(); ++j)
                {
                    if (_x1->hasFiniteLowerBound(j))
                    {
                        ineq_values(cur_idx++) = _x1->lowerBound()[j] - (*_intermediate_x[i])[j];  // x_min - x <= 0
                    }
                }
            }
        }

        // upper bounds
        if (_dim_ub_x > 0)
        {
            for (int i = 0; i < _intermediate_x.size(); ++i)
            {
                for (int j = 0; j < _intermediate_x[i]->size(); ++j)
                {
                    if (_x1->hasFiniteUpperBound(j))
                    {
                        ineq_values(cur_idx++) = (*_intermediate_x[i])[j] - _x1->upperBound()[j];  // x - x_max <= 0
                    }
                }
            }
        }
        assert(cur_idx == getInequalityDimension());
    }

    void setCollocationMethod(QuadratureCollocationInterface::Ptr quadrature)
    {
        assert(!quadrature->isIntermediateControlSubjectToOptim() && quadrature->isSupportingCompressedStatesMode());
        _collocation = quadrature;
        _intermediate_x.resize(_collocation->getNumIntermediatePoints(), new Eigen::VectorXd(_dim_dyn));

        activateIntermediateConstraints();
    }

 protected:
    void activateIntermediateConstraints()
    {
        // get number inequalities in case _eval_intermediate_constr is turned on
        if (_eval_intermediate_constr)
        {
            int num_interm_states = _collocation->getNumIntermediatePoints();
            _dim_nonint_eq        = _stage_eq ? _stage_eq->getConcatenatedNonIntegralStateTermDimension(_k) : 0;
            _dim_nonint_ineq      = _stage_ineq ? _stage_ineq->getConcatenatedNonIntegralStateTermDimension(_k) : 0;

            // consider number of finite bounds for all vertices
            // we apply the bounds of vertex x1
            _dim_lb_x = _x1->getNumberFiniteLowerBounds(false);
            _dim_ub_x = _x1->getNumberFiniteUpperBounds(false);

            // update overall dimensions
            _dim_eq   = _dim_int_eq + _dim_dyn + num_interm_states * _dim_nonint_eq;
            _dim_ineq = _dim_int_ineq + num_interm_states * (_dim_nonint_ineq + _dim_lb_x + _dim_ub_x);
        }
    }

 private:
    SystemDynamicsInterface::Ptr _dynamics;
    QuadratureCollocationInterface::Ptr _collocation;

    StageCost::Ptr _stage_cost;
    StageEqualityConstraint::Ptr _stage_eq;
    StageInequalityConstraint::Ptr _stage_ineq;

    int _dim_obj      = 1;
    int _dim_eq       = 0;
    int _dim_int_eq   = 0;
    int _dim_ineq     = 0;
    int _dim_int_ineq = 0;
    int _dim_dyn      = 0;

    int _dim_nonint_eq   = 0;
    int _dim_nonint_ineq = 0;

    int _dim_lb_x = 0;
    int _dim_ub_x = 0;

    int _k = 0;

    bool _eval_intermediate_constr = false;

    Eigen::VectorXd _dynamics_quadrature;
    std::vector<Eigen::VectorXd*> _intermediate_x;

    const VectorVertex* _x1 = nullptr;
    const VectorVertex* _u1 = nullptr;
    const VectorVertex* _x2 = nullptr;
    const ScalarVertex* _dt = nullptr;

 public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

class CombinedCompressedCollocationEdge : public MixedEdge<VectorVertex, VectorVertex, ScalarVertex, VectorVertex, VectorVertex>
{
 public:
    using Ptr  = std::shared_ptr<CombinedCompressedCollocationEdge>;
    using UPtr = std::unique_ptr<CombinedCompressedCollocationEdge>;

    explicit CombinedCompressedCollocationEdge(SystemDynamicsInterface::Ptr dynamics, StageCost::Ptr stage_cost,
                                               StageEqualityConstraint::Ptr stage_eq, StageInequalityConstraint::Ptr stage_ineq,
                                               bool eval_intermediate_constr, int k, VectorVertex& x1, VectorVertex& u1, ScalarVertex& dt,
                                               VectorVertex& u2, VectorVertex& x2)
        : MixedEdge<VectorVertex, VectorVertex, ScalarVertex, VectorVertex, VectorVertex>(x1, u1, dt, u2, x2),
          _dynamics(dynamics),
          _stage_cost(stage_cost),
          _stage_eq(stage_eq),
          _stage_ineq(stage_ineq),
          _k(k),
          _eval_intermediate_constr(eval_intermediate_constr)
    {
        int x_dim = _dynamics->getStateDimension();
        _dynamics_quadrature.resize(x_dim);

        _dim_obj      = (_stage_cost && _stage_cost->getIntegralStateControlTermDimension(_k) > 0) ? 1 : 0;  // we have only a single obejtive value
        _dim_int_eq   = _stage_eq ? _stage_eq->getIntegralStateControlTermDimension(k) : 0;
        _dim_eq       = x_dim + _dim_int_eq;
        _dim_int_ineq = _stage_ineq ? _stage_ineq->getIntegralStateControlTermDimension(k) : 0;
        _dim_ineq     = _dim_int_ineq;
        _dim_dyn      = _dynamics->getStateDimension();

        _x1 = static_cast<const VectorVertex*>(_vertices[0]);
        _u1 = static_cast<const VectorVertex*>(_vertices[1]);
        _dt = static_cast<const ScalarVertex*>(_vertices[2]);
        _u2 = static_cast<const VectorVertex*>(_vertices[3]);
        _x2 = static_cast<const VectorVertex*>(_vertices[4]);
    }

    virtual ~CombinedCompressedCollocationEdge()
    {
        for (int i = 0; i < (int)_intermediate_x.size(); ++i) delete _intermediate_x[i];
        for (int i = 0; i < (int)_intermediate_u.size(); ++i) delete _intermediate_u[i];
    }

    // implements interface method
    int getObjectiveDimension() const override { return _dim_obj; }
    // implements interface method
    int getEqualityDimension() const override { return _dim_eq; }
    int getInequalityDimension() const override { return _dim_ineq; }

    // implement in child class:
    bool isObjectiveLinear() const override { return false; }
    bool isEqualityLinear() const override { return false; }
    bool isInequalityLinear() const override { return false; }

    // bool providesJacobian() const override { return false; }

    bool isObjectiveLeastSquaresForm() const override { return false; }

    void precompute() override
    {
        assert(_collocation);
        assert(_dynamics);

        assert(_x1->getDimension() == _dynamics->getStateDimension());
        assert(_u1->getDimension() == _dynamics->getInputDimension());
        assert(_u2->getDimension() == _dynamics->getInputDimension());
        assert(_x2->getDimension() == _dynamics->getStateDimension());

        _collocation->quadrature(_x1->values(), _u1->values(), _x2->values(), _u2->values(), _dt->value(), *_dynamics, _dynamics_quadrature,
                                 _intermediate_x, _intermediate_u);
    }

    void computeObjectiveValues(Eigen::Ref<Eigen::VectorXd> obj_values) override
    {
        // evaluate quadrature formula with previously computed intermediate states
        auto integrand = [this](const Eigen::VectorXd& x, const Eigen::VectorXd& u, Eigen::Ref<Eigen::VectorXd> result) {
            _stage_cost->computeIntegralStateControlTerm(_k, x, u, result);
        };
        _collocation->quadrature(_x1->values(), _u1->values(), _x2->values(), _u2->values(), _dt->value(), _intermediate_x, _intermediate_u,
                                 integrand, obj_values);
    }

    void computeEqualityValues(Eigen::Ref<Eigen::VectorXd> eq_values) override
    {
        eq_values.head(_dynamics_quadrature.size()).noalias() = _x2->values() - _x1->values() - _dynamics_quadrature;
        if (_dim_int_eq > 0)
        {
            // evaluate quadrature formula with previously computed intermediate states
            auto integrand = [this](const Eigen::VectorXd& x, const Eigen::VectorXd& u, Eigen::Ref<Eigen::VectorXd> result) {
                _stage_eq->computeIntegralStateControlTerm(_k, x, u, result);
            };
            _collocation->quadrature(_x1->values(), _u1->values(), _x2->values(), _u2->values(), _dt->value(), _intermediate_x, _intermediate_u,
                                     integrand, eq_values.segment(_dim_dyn, _dim_int_eq));
        }

        // in case we consider intermediate constraints
        if (_dim_nonint_eq > 0)
        {
            int cur_idx = _dim_dyn + _dim_int_eq;
            for (int i = 0; i < (int)_intermediate_x.size(); ++i)
            {
                _stage_eq->computeConcatenatedNonIntegralStateTerms(_k, *_intermediate_x[i], *_intermediate_u[i], _dt->value(),
                                                                    eq_values.segment(cur_idx, _dim_nonint_eq));
                cur_idx += _dim_nonint_eq;
            }
            assert(cur_idx == getEqualityDimension());
        }
    }
    void computeInequalityValues(Eigen::Ref<Eigen::VectorXd> ineq_values) override
    {
        if (_dim_int_ineq > 0)
        {
            // evaluate quadrature formula with previously computed intermediate states
            auto integrand = [this](const Eigen::VectorXd& x, const Eigen::VectorXd& u, Eigen::Ref<Eigen::VectorXd> result) {
                _stage_ineq->computeIntegralStateControlTerm(_k, x, u, result);
            };
            _collocation->quadrature(_x1->values(), _u1->values(), _x2->values(), _u2->values(), _dt->value(), _intermediate_x, _intermediate_u,
                                     integrand, ineq_values);
        }

        // in case we consider intermediate constraints, the following dimensions are > 0
        int cur_idx = _dim_int_eq;
        if (_dim_nonint_ineq > 0)
        {
            for (int i = 0; i < (int)_intermediate_x.size(); ++i)
            {
                _stage_ineq->computeConcatenatedNonIntegralStateTerms(_k, *_intermediate_x[i], *_intermediate_u[i], _dt->value(),
                                                                      ineq_values.segment(cur_idx, _dim_nonint_ineq));
                cur_idx += _dim_nonint_eq;
            }
        }

        // lower bounds
        if (_dim_lb_x > 0)
        {
            for (int i = 0; i < (int)_intermediate_x.size(); ++i)
            {
                for (int j = 0; j < _intermediate_x[i]->size(); ++j)
                {
                    if (_x1->hasFiniteLowerBound(j))
                    {
                        ineq_values(cur_idx++) = _x1->lowerBound()[j] - (*_intermediate_x[i])[j];  // x_min - x <= 0
                    }
                }
            }
        }

        // upper bounds
        if (_dim_ub_x > 0)
        {
            for (int i = 0; i < (int)_intermediate_x.size(); ++i)
            {
                for (int j = 0; j < _intermediate_x[i]->size(); ++j)
                {
                    if (_x1->hasFiniteUpperBound(j))
                    {
                        ineq_values(cur_idx++) = (*_intermediate_x[i])[j] - _x1->upperBound()[j];  // x - x_max <= 0
                    }
                }
            }
        }
        assert(cur_idx == getInequalityDimension());
    }

    void setCollocationMethod(QuadratureCollocationInterface::Ptr quadrature)
    {
        assert(!quadrature->isIntermediateControlSubjectToOptim() && quadrature->isSupportingCompressedStatesMode());

        _collocation = quadrature;
        _intermediate_x.resize(_collocation->getNumIntermediatePoints(), new Eigen::VectorXd(_dim_dyn));
        _intermediate_u.resize(_collocation->getNumIntermediatePoints(), new Eigen::VectorXd(_dynamics->getInputDimension()));

        activateIntermediateConstraints();
    }

 protected:
    void activateIntermediateConstraints()
    {
        // get number inequalities in case _eval_intermediate_constr is turned on
        // TODO(roesmann) currently we only evaluate intermediate states, since u is usually just a linear function
        if (_eval_intermediate_constr)
        {
            int num_interm_states = _collocation->getNumIntermediatePoints();
            _dim_nonint_eq        = _stage_eq ? _stage_eq->getConcatenatedNonIntegralStateTermDimension(_k) : 0;
            _dim_nonint_ineq      = _stage_ineq ? _stage_ineq->getConcatenatedNonIntegralStateTermDimension(_k) : 0;

            // consider number of finite bounds for all vertices
            // we apply the bounds of vertex x1
            _dim_lb_x = _x1->getNumberFiniteLowerBounds(false);
            _dim_ub_x = _x1->getNumberFiniteUpperBounds(false);

            // update overall dimensions
            _dim_eq   = _dim_int_eq + _dim_dyn + num_interm_states * _dim_nonint_eq;
            _dim_ineq = _dim_int_ineq + num_interm_states * (_dim_nonint_ineq + _dim_lb_x + _dim_ub_x);
        }
    }

 private:
    SystemDynamicsInterface::Ptr _dynamics;
    QuadratureCollocationInterface::Ptr _collocation;

    StageCost::Ptr _stage_cost;
    StageEqualityConstraint::Ptr _stage_eq;
    StageInequalityConstraint::Ptr _stage_ineq;

    int _dim_obj      = 1;
    int _dim_eq       = 0;
    int _dim_int_eq   = 0;
    int _dim_ineq     = 0;
    int _dim_int_ineq = 0;
    int _dim_dyn      = 0;

    int _dim_nonint_eq   = 0;
    int _dim_nonint_ineq = 0;

    int _dim_lb_x = 0;
    int _dim_ub_x = 0;

    int _k = 0;

    bool _eval_intermediate_constr = false;

    Eigen::VectorXd _dynamics_quadrature;
    std::vector<Eigen::VectorXd*> _intermediate_x;
    std::vector<Eigen::VectorXd*> _intermediate_u;

    const VectorVertex* _x1 = nullptr;
    const VectorVertex* _u1 = nullptr;
    const VectorVertex* _u2 = nullptr;
    const VectorVertex* _x2 = nullptr;
    const ScalarVertex* _dt = nullptr;

 public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

class CombinedCompressedCollocationMultipleControlsEdge : public MixedEdge<>
{
 public:
    using Ptr  = std::shared_ptr<CombinedCompressedCollocationMultipleControlsEdge>;
    using UPtr = std::unique_ptr<CombinedCompressedCollocationMultipleControlsEdge>;

    explicit CombinedCompressedCollocationMultipleControlsEdge(SystemDynamicsInterface::Ptr dynamics, QuadratureCollocationInterface::Ptr quadrature,
                                                               StageCost::Ptr stage_cost, StageEqualityConstraint::Ptr stage_eq,
                                                               StageInequalityConstraint::Ptr stage_ineq, bool eval_intermediate_constr, int k)
        : MixedEdge<>(5 + quadrature->getNumIntermediatePoints()),
          _dynamics(dynamics),
          _collocation(quadrature),
          _stage_cost(stage_cost),
          _stage_eq(stage_eq),
          _stage_ineq(stage_ineq),
          _k(k),
          _eval_intermediate_constr(eval_intermediate_constr),
          _num_intermediate_points(quadrature->getNumIntermediatePoints())
    {
        assert(quadrature->isIntermediateControlSubjectToOptim());
    }

    void finalize()
    {
        int x_dim = _dynamics->getStateDimension();
        _dynamics_quadrature.resize(x_dim);

        _dim_obj      = (_stage_cost && _stage_cost->getIntegralStateControlTermDimension(_k) > 0) ? 1 : 0;  // we have only a single obejtive value
        _dim_int_eq   = _stage_eq ? _stage_eq->getIntegralStateControlTermDimension(_k) : 0;
        _dim_eq       = x_dim + _dim_int_eq;
        _dim_int_ineq = _stage_ineq ? _stage_ineq->getIntegralStateControlTermDimension(_k) : 0;
        _dim_ineq     = _dim_int_ineq;
        _dim_dyn      = _dynamics->getStateDimension();

        _x1 = static_cast<const VectorVertex*>(_vertices[0]);
        _u1 = static_cast<const VectorVertex*>(_vertices[1]);
        _dt = static_cast<const ScalarVertex*>(_vertices[2]);
        _u2 = static_cast<const VectorVertex*>(_vertices[3]);
        _x2 = static_cast<const VectorVertex*>(_vertices[4]);

        int idx = 5;
        _intermediate_u.clear();
        for (int i = 0; i < _num_intermediate_points; ++i)
        {
            VectorVertex* interm_u = static_cast<VectorVertex*>(_vertices[idx + i]);
            _intermediate_u.push_back(&interm_u->values());
        }

        clearInternalBuffer();

        for (int i = 0; i < _num_intermediate_points; ++i)
        {
            _intermediate_x.push_back(new Eigen::VectorXd(_dim_dyn));
        }

        activateIntermediateConstraints();
    }

    virtual ~CombinedCompressedCollocationMultipleControlsEdge() { clearInternalBuffer(); }

    // implements interface method
    int getObjectiveDimension() const override { return _dim_obj; }
    // implements interface method
    int getEqualityDimension() const override { return _dim_eq; }
    int getInequalityDimension() const override { return _dim_ineq; }

    // implement in child class:
    bool isObjectiveLinear() const override { return false; }
    bool isEqualityLinear() const override { return false; }
    bool isInequalityLinear() const override { return false; }

    // bool providesJacobian() const override { return false; }

    bool isObjectiveLeastSquaresForm() const override { return false; }

    void precompute() override
    {
        assert(_collocation);
        assert(_dynamics);

        assert(_x1->getDimension() == _dynamics->getStateDimension());
        assert(_u1->getDimension() == _dynamics->getInputDimension());
        assert(_u2->getDimension() == _dynamics->getInputDimension());
        assert(_x2->getDimension() == _dynamics->getStateDimension());

        // perform quadrature with given control mid points
        _collocation->quadrature(_x1->values(), _u1->values(), _x2->values(), _u2->values(), _dt->value(), *_dynamics, _intermediate_u,
                                 _dynamics_quadrature, _intermediate_x);
    }

    void computeObjectiveValues(Eigen::Ref<Eigen::VectorXd> obj_values) override
    {
        // evaluate quadrature formula with previously computed intermediate states
        auto integrand = [this](const Eigen::VectorXd& x, const Eigen::VectorXd& u, Eigen::Ref<Eigen::VectorXd> result) {
            _stage_cost->computeIntegralStateControlTerm(_k, x, u, result);
        };
        _collocation->quadrature(_x1->values(), _u1->values(), _x2->values(), _u2->values(), _dt->value(), _intermediate_x, _intermediate_u,
                                 integrand, obj_values);
    }

    void computeEqualityValues(Eigen::Ref<Eigen::VectorXd> eq_values) override
    {
        eq_values.head(_dynamics_quadrature.size()).noalias() = _x2->values() - _x1->values() - _dynamics_quadrature;
        if (_dim_int_eq > 0)
        {
            // evaluate quadrature formula with previously computed intermediate states
            auto integrand = [this](const Eigen::VectorXd& x, const Eigen::VectorXd& u, Eigen::Ref<Eigen::VectorXd> result) {
                _stage_eq->computeIntegralStateControlTerm(_k, x, u, result);
            };
            _collocation->quadrature(_x1->values(), _u1->values(), _x2->values(), _u2->values(), _dt->value(), _intermediate_x, _intermediate_u,
                                     integrand, eq_values.segment(_dim_dyn, _dim_int_eq));
        }

        // in case we consider intermediate constraints
        if (_dim_nonint_eq > 0)
        {
            int cur_idx = _dim_dyn + _dim_int_eq;
            for (int i = 0; i < (int)_intermediate_x.size(); ++i)
            {
                _stage_eq->computeConcatenatedNonIntegralStateTerms(_k, *_intermediate_x[i], *_intermediate_u[i], _dt->value(),
                                                                    eq_values.segment(cur_idx, _dim_nonint_eq));
                cur_idx += _dim_nonint_eq;
            }
            assert(cur_idx == getEqualityDimension());
        }
    }
    void computeInequalityValues(Eigen::Ref<Eigen::VectorXd> ineq_values) override
    {
        if (_dim_int_ineq > 0)
        {
            // evaluate quadrature formula with previously computed intermediate states
            auto integrand = [this](const Eigen::VectorXd& x, const Eigen::VectorXd& u, Eigen::Ref<Eigen::VectorXd> result) {
                _stage_ineq->computeIntegralStateControlTerm(_k, x, u, result);
            };
            _collocation->quadrature(_x1->values(), _u1->values(), _x2->values(), _u2->values(), _dt->value(), _intermediate_x, _intermediate_u,
                                     integrand, ineq_values);
        }

        // in case we consider intermediate constraints, the following dimensions are > 0
        int cur_idx = _dim_int_eq;
        if (_dim_nonint_ineq > 0)
        {
            for (int i = 0; i < (int)_intermediate_x.size(); ++i)
            {
                _stage_ineq->computeConcatenatedNonIntegralStateTerms(_k, *_intermediate_x[i], *_intermediate_u[i], _dt->value(),
                                                                      ineq_values.segment(cur_idx, _dim_nonint_ineq));
                cur_idx += _dim_nonint_eq;
            }
        }

        // lower bounds
        if (_dim_lb_x > 0)
        {
            for (int i = 0; i < (int)_intermediate_x.size(); ++i)
            {
                for (int j = 0; j < _intermediate_x[i]->size(); ++j)
                {
                    if (_x1->hasFiniteLowerBound(j))
                    {
                        ineq_values(cur_idx++) = _x1->lowerBound()[j] - (*_intermediate_x[i])[j];  // x_min - x <= 0
                    }
                }
            }
        }

        // upper bounds
        if (_dim_ub_x > 0)
        {
            for (int i = 0; i < (int)_intermediate_x.size(); ++i)
            {
                for (int j = 0; j < _intermediate_x[i]->size(); ++j)
                {
                    if (_x1->hasFiniteUpperBound(j))
                    {
                        ineq_values(cur_idx++) = (*_intermediate_x[i])[j] - _x1->upperBound()[j];  // x - x_max <= 0
                    }
                }
            }
        }
        assert(cur_idx == getInequalityDimension());
    }

 protected:
    void activateIntermediateConstraints()
    {
        // get number inequalities in case _eval_intermediate_constr is turned on
        // TODO(roesmann) currently we only evaluate intermediate states, since u is usually just a linear function
        if (_eval_intermediate_constr)
        {
            int num_interm_states = _collocation->getNumIntermediatePoints();
            _dim_nonint_eq        = _stage_eq ? _stage_eq->getConcatenatedNonIntegralStateTermDimension(_k) : 0;
            _dim_nonint_ineq      = _stage_ineq ? _stage_ineq->getConcatenatedNonIntegralStateTermDimension(_k) : 0;

            // consider number of finite bounds for all vertices
            // we apply the bounds of vertex x1
            _dim_lb_x = _x1->getNumberFiniteLowerBounds(false);
            _dim_ub_x = _x1->getNumberFiniteUpperBounds(false);

            // update overall dimensions
            _dim_eq   = _dim_int_eq + _dim_dyn + num_interm_states * _dim_nonint_eq;
            _dim_ineq = _dim_int_ineq + num_interm_states * (_dim_nonint_ineq + _dim_lb_x + _dim_ub_x);
        }
    }

    void clearInternalBuffer()
    {
        for (int i = 0; i < _intermediate_x.size(); ++i) delete _intermediate_x[i];
        _intermediate_x.clear();
    }

 private:
    SystemDynamicsInterface::Ptr _dynamics;
    QuadratureCollocationInterface::Ptr _collocation;

    StageCost::Ptr _stage_cost;
    StageEqualityConstraint::Ptr _stage_eq;
    StageInequalityConstraint::Ptr _stage_ineq;

    int _dim_obj      = 1;
    int _dim_eq       = 0;
    int _dim_int_eq   = 0;
    int _dim_ineq     = 0;
    int _dim_int_ineq = 0;
    int _dim_dyn      = 0;

    int _dim_nonint_eq   = 0;
    int _dim_nonint_ineq = 0;

    int _dim_lb_x = 0;
    int _dim_ub_x = 0;

    int _k = 0;

    bool _eval_intermediate_constr = false;

    int _num_intermediate_points;

    Eigen::VectorXd _dynamics_quadrature;
    std::vector<Eigen::VectorXd*> _intermediate_x;
    std::vector<Eigen::VectorXd*> _intermediate_u;

    const VectorVertex* _x1 = nullptr;
    const VectorVertex* _u1 = nullptr;
    const VectorVertex* _u2 = nullptr;
    const VectorVertex* _x2 = nullptr;
    const ScalarVertex* _dt = nullptr;

 public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

class UncompressedCollocationEdge : public Edge<>
{
 public:
    using Ptr  = std::shared_ptr<UncompressedCollocationEdge>;
    using UPtr = std::unique_ptr<UncompressedCollocationEdge>;

    explicit UncompressedCollocationEdge(SystemDynamicsInterface::Ptr dynamics, QuadratureCollocationInterface::Ptr quadrature)
        : Edge<>(5 + quadrature->getNumIntermediatePoints() * (1 + (int)quadrature->isIntermediateControlSubjectToOptim())),
          _dynamics(dynamics),
          _collocation(quadrature),
          _num_intermediate_points(quadrature->getNumIntermediatePoints())
    {
    }

    void finalize()
    {
        _dim_x = _dynamics->getStateDimension();

        _x1 = static_cast<const VectorVertex*>(_vertices[0]);
        _u1 = static_cast<const VectorVertex*>(_vertices[1]);
        _dt = static_cast<const ScalarVertex*>(_vertices[2]);
        _u2 = static_cast<const VectorVertex*>(_vertices[3]);
        _x2 = static_cast<const VectorVertex*>(_vertices[4]);

        int idx = 5;
        _intermediate_x.clear();
        for (int i = 0; i < _num_intermediate_points; ++i)
        {
            VectorVertex* interm_x = static_cast<VectorVertex*>(_vertices[idx + i]);
            _intermediate_x.push_back(&interm_x->values());
        }

        idx += _num_intermediate_points;

        _dimension = _dim_x + _num_intermediate_points * _dim_x;

        clearInternalBuffer();

        if (_collocation->isIntermediateControlSubjectToOptim())
        {
            _intermediate_u_internal_buf = false;
            for (int i = 0; i < _num_intermediate_points; ++i)
            {
                VectorVertex* interm_u = static_cast<VectorVertex*>(_vertices[idx + i]);
                _intermediate_u.push_back(&interm_u->values());
            }
        }
        else if (_num_intermediate_points > 0)
        {
            _intermediate_u_internal_buf = true;
            for (int i = 0; i < _num_intermediate_points; ++i)
            {
                _intermediate_u.push_back(new Eigen::VectorXd(_dynamics->getInputDimension()));
            }
        }
    }

    // implements interface method
    int getDimension() const override { return _dimension; }
    // implement in child class:
    bool isLinear() const override { return false; }

    // bool providesJacobian() const override { return false; }

    bool isLeastSquaresForm() const override { return false; }

    void computeValues(Eigen::Ref<Eigen::VectorXd> values) override
    {
        assert(_collocation);
        assert(_dynamics);
        assert(_x1->getDimension() == _dynamics->getStateDimension());
        assert(_u1->getDimension() == _dynamics->getInputDimension());
        assert(_x2->getDimension() == _dynamics->getStateDimension());
        assert(_u2->getDimension() == _dynamics->getInputDimension());

        if (_intermediate_u_internal_buf)
        {
            _collocation->computeIntermediateControls(_u1->values(), _u2->values(), _intermediate_u);
        }

        // preform collocation with given mid points:
        _collocation->quadrature(_x1->values(), _u1->values(), _x2->values(), _u2->values(), _dt->value(), *_dynamics, _intermediate_x,
                                 _intermediate_u, values.head(_dim_x), values.tail(_dim_x * _num_intermediate_points));

        values.head(_dim_x).noalias() -= (_x2->values() - _x1->values());
    }

 protected:
    void clearInternalBuffer()
    {
        if (_intermediate_u_internal_buf)
        {
            for (int i = 0; i < _intermediate_u.size(); ++i) delete _intermediate_u[i];
        }
        _intermediate_u.clear();
    }

 private:
    SystemDynamicsInterface::Ptr _dynamics;
    QuadratureCollocationInterface::Ptr _collocation;

    Eigen::VectorXd _dynamics_quadrature;
    Eigen::VectorXd _midpoint_error;

    int _dim_x     = 0;
    int _dimension = 0;

    int _num_intermediate_points = 0;
    std::vector<Eigen::VectorXd*> _intermediate_x;
    std::vector<Eigen::VectorXd*> _intermediate_u;
    bool _intermediate_u_internal_buf = false;

    const VectorVertex* _x1 = nullptr;
    const VectorVertex* _u1 = nullptr;
    const VectorVertex* _x2 = nullptr;
    const VectorVertex* _u2 = nullptr;
    const ScalarVertex* _dt = nullptr;

 public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

// vertices
// x1
// u1
// dt
// x2
// u2
// intermediate_x1
// ..
// intermediate_u1
// ..
class CombinedUncompressedCollocationEdge : public MixedEdge<>
{
 public:
    using Ptr  = std::shared_ptr<CombinedUncompressedCollocationEdge>;
    using UPtr = std::unique_ptr<CombinedUncompressedCollocationEdge>;

    explicit CombinedUncompressedCollocationEdge(SystemDynamicsInterface::Ptr dynamics, QuadratureCollocationInterface::Ptr quadrature,
                                                 StageCost::Ptr stage_cost, StageEqualityConstraint::Ptr stage_eq,
                                                 StageInequalityConstraint::Ptr stage_ineq, bool eval_intermediate_constr, int k)
        : MixedEdge<>(5 + quadrature->getNumIntermediatePoints() * (1 + (int)quadrature->isIntermediateControlSubjectToOptim())),
          _dynamics(dynamics),
          _collocation(quadrature),
          _stage_cost(stage_cost),
          _stage_eq(stage_eq),
          _stage_ineq(stage_ineq),
          _k(k),
          _eval_intermediate_constr(eval_intermediate_constr),
          _num_intermediate_points(quadrature->getNumIntermediatePoints())
    {
    }

    virtual ~CombinedUncompressedCollocationEdge() { clearInternalBuffer(); }

    void finalize()
    {
        int x_dim = _dynamics->getStateDimension();

        _dim_obj    = (_stage_cost && _stage_cost->getIntegralStateControlTermDimension(_k) > 0) ? 1 : 0;  // we have only a single obejtive value
        _dim_int_eq = _stage_eq ? _stage_eq->getIntegralStateControlTermDimension(_k) : 0;
        // Note, we have additional equations for intermediate points in uncompressed mode:
        _dim_eq       = x_dim + _dim_int_eq + x_dim * _num_intermediate_points;
        _dim_int_ineq = _stage_ineq ? _stage_ineq->getIntegralStateControlTermDimension(_k) : 0;
        _dim_ineq     = _dim_int_ineq;
        _dim_dyn      = _dynamics->getStateDimension();

        _dynamics_quadrature.resize(x_dim);
        _midpoint_error.resize(x_dim * _num_intermediate_points);

        _x1 = static_cast<const VectorVertex*>(_vertices[0]);
        _u1 = static_cast<const VectorVertex*>(_vertices[1]);
        _dt = static_cast<const ScalarVertex*>(_vertices[2]);
        _u2 = static_cast<const VectorVertex*>(_vertices[3]);
        _x2 = static_cast<const VectorVertex*>(_vertices[4]);

        int idx = 5;
        _intermediate_x.clear();
        for (int i = 0; i < _num_intermediate_points; ++i)
        {
            VectorVertex* interm_x = static_cast<VectorVertex*>(_vertices[idx + i]);
            _intermediate_x.push_back(&interm_x->values());
        }

        idx += _num_intermediate_points;

        clearInternalBuffer();

        if (_collocation->isIntermediateControlSubjectToOptim())
        {
            _intermediate_u_internal_buf = false;
            for (int i = 0; i < _num_intermediate_points; ++i)
            {
                VectorVertex* interm_u = static_cast<VectorVertex*>(_vertices[idx + i]);
                _intermediate_u.push_back(&interm_u->values());
            }
        }
        else if (_num_intermediate_points > 0)
        {
            _intermediate_u_internal_buf = true;
            for (int i = 0; i < _num_intermediate_points; ++i)
            {
                _intermediate_u.push_back(new Eigen::VectorXd(_dynamics->getInputDimension()));
            }
        }

        activateIntermediateConstraints();
    }

    // implements interface method
    int getObjectiveDimension() const override { return _dim_obj; }
    // implements interface method
    int getEqualityDimension() const override { return _dim_eq; }
    int getInequalityDimension() const override { return _dim_ineq; }

    // implement in child class:
    bool isObjectiveLinear() const override { return false; }
    bool isEqualityLinear() const override { return false; }
    bool isInequalityLinear() const override { return false; }

    // bool providesJacobian() const override { return false; }

    bool isObjectiveLeastSquaresForm() const override { return false; }

    void precompute() override
    {
        assert(_collocation);
        assert(_dynamics);

        assert(_x1->getDimension() == _dynamics->getStateDimension());
        assert(_u1->getDimension() == _dynamics->getInputDimension());
        assert(_u2->getDimension() == _dynamics->getInputDimension());
        assert(_x2->getDimension() == _dynamics->getStateDimension());

        if (_intermediate_u_internal_buf)
        {
            _collocation->computeIntermediateControls(_u1->values(), _u2->values(), _intermediate_u);
        }

        // preform collocation with given mid points:
        _collocation->quadrature(_x1->values(), _u1->values(), _x2->values(), _u2->values(), _dt->value(), *_dynamics, _intermediate_x,
                                 _intermediate_u, _dynamics_quadrature, _midpoint_error);
    }

    void computeObjectiveValues(Eigen::Ref<Eigen::VectorXd> obj_values) override
    {
        // evaluate quadrature formula with previously computed intermediate states
        auto integrand = [this](const Eigen::VectorXd& x, const Eigen::VectorXd& u, Eigen::Ref<Eigen::VectorXd> result) {
            _stage_cost->computeIntegralStateControlTerm(_k, x, u, result);
        };
        _collocation->quadrature(_x1->values(), _u1->values(), _x2->values(), _u2->values(), _dt->value(), _intermediate_x, _intermediate_u,
                                 integrand, obj_values);
    }

    void computeEqualityValues(Eigen::Ref<Eigen::VectorXd> eq_values) override
    {
        eq_values.head(_dim_dyn).noalias() = _x2->values() - _x1->values() - _dynamics_quadrature;
        eq_values.segment(_dim_dyn, _midpoint_error.size()).noalias() = _midpoint_error;

        int idx = _midpoint_error.size() + _dim_dyn;

        if (_dim_int_eq > 0)
        {
            // evaluate quadrature formula with previously computed intermediate states
            auto integrand = [this](const Eigen::VectorXd& x, const Eigen::VectorXd& u, Eigen::Ref<Eigen::VectorXd> result) {
                _stage_eq->computeIntegralStateControlTerm(_k, x, u, result);
            };
            _collocation->quadrature(_x1->values(), _u1->values(), _x2->values(), _u2->values(), _dt->value(), _intermediate_x, _intermediate_u,
                                     integrand, eq_values.segment(idx, _dim_int_eq));
            idx += _dim_int_eq;
        }

        // in case we consider intermediate constraints
        if (_dim_nonint_eq > 0)
        {
            for (int i = 0; i < (int)_intermediate_x.size(); ++i)
            {
                if (i < _intermediate_u.size())
                {
                    _stage_eq->computeConcatenatedNonIntegralStateTerms(_k, *_intermediate_x[i], _u1->values(), _dt->value(),
                                                                        eq_values.segment(idx, _dim_nonint_eq));
                }
                idx += _dim_nonint_eq;
            }
            assert(idx == getEqualityDimension());
        }
    }
    void computeInequalityValues(Eigen::Ref<Eigen::VectorXd> ineq_values) override
    {
        if (_dim_int_ineq > 0)
        {
            // evaluate quadrature formula with previously computed intermediate states
            auto integrand = [this](const Eigen::VectorXd& x, const Eigen::VectorXd& u, Eigen::Ref<Eigen::VectorXd> result) {
                _stage_ineq->computeIntegralStateControlTerm(_k, x, u, result);
            };
            _collocation->quadrature(_x1->values(), _u1->values(), _x2->values(), _u2->values(), _dt->value(), _intermediate_x, _intermediate_u,
                                     integrand, ineq_values);
        }

        // in case we consider intermediate constraints, the following dimensions are > 0
        int cur_idx = _dim_int_eq;
        if (_dim_nonint_ineq > 0)
        {
            for (int i = 0; i < _collocation->getNumIntermediatePoints(); ++i)
            {
                _stage_ineq->computeConcatenatedNonIntegralStateTerms(_k, _x1->values(), _u1->values(), _dt->value(),
                                                                      ineq_values.segment(cur_idx, _dim_nonint_ineq));
                cur_idx += _dim_nonint_eq;
            }
        }
        assert(cur_idx == getInequalityDimension());
    }

 protected:
    void activateIntermediateConstraints()
    {
        // get number inequalities in case _eval_intermediate_constr is turned on
        // TODO(roesmann) currently we only evaluate intermediate states, since u is usually just a linear function
        if (_eval_intermediate_constr)
        {
            int num_interm_states = _collocation->getNumIntermediatePoints();
            _dim_nonint_eq        = _stage_eq ? _stage_eq->getConcatenatedNonIntegralStateTermDimension(_k) : 0;
            _dim_nonint_ineq      = _stage_ineq ? _stage_ineq->getConcatenatedNonIntegralStateTermDimension(_k) : 0;

            // update overall dimensions
            _dim_eq   = _dim_int_eq + _dim_dyn + _dim_dyn * _num_intermediate_points + num_interm_states * _dim_nonint_eq;
            _dim_ineq = _dim_int_ineq + num_interm_states * _dim_nonint_ineq;
        }
    }

    void clearInternalBuffer()
    {
        if (_intermediate_u_internal_buf)
        {
            for (int i = 0; i < _intermediate_u.size(); ++i) delete _intermediate_u[i];
        }
        _intermediate_u.clear();
    }

 private:
    SystemDynamicsInterface::Ptr _dynamics;
    QuadratureCollocationInterface::Ptr _collocation;

    StageCost::Ptr _stage_cost;
    StageEqualityConstraint::Ptr _stage_eq;
    StageInequalityConstraint::Ptr _stage_ineq;

    int _dim_obj      = 1;
    int _dim_eq       = 0;
    int _dim_int_eq   = 0;
    int _dim_ineq     = 0;
    int _dim_int_ineq = 0;
    int _dim_dyn      = 0;

    int _dim_nonint_eq   = 0;
    int _dim_nonint_ineq = 0;

    int _k = 0;

    bool _eval_intermediate_constr = false;

    Eigen::VectorXd _dynamics_quadrature;
    Eigen::VectorXd _midpoint_error;

    int _num_intermediate_points = 0;
    std::vector<Eigen::VectorXd*> _intermediate_x;
    std::vector<Eigen::VectorXd*> _intermediate_u;
    bool _intermediate_u_internal_buf = false;

    const VectorVertex* _x1 = nullptr;
    const VectorVertex* _u1 = nullptr;
    const VectorVertex* _u2 = nullptr;
    const VectorVertex* _x2 = nullptr;
    const ScalarVertex* _dt = nullptr;

 public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

}  // namespace corbo

#endif  // SRC_OPTIMAL_CONTROL_INCLUDE_CORBO_OPTIMAL_CONTROL_STRUCTURED_OCP_EDGES_COLLOCATION_EDGES_H_
