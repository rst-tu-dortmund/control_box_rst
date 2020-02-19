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

#ifndef SRC_OPTIMAL_CONTROL_INCLUDE_CORBO_OPTIMAL_CONTROL_FUNCTIONS_HYBRID_COST_H_
#define SRC_OPTIMAL_CONTROL_INCLUDE_CORBO_OPTIMAL_CONTROL_FUNCTIONS_HYBRID_COST_H_

#include <corbo-optimal-control/functions/minimum_time.h>
#include <corbo-optimal-control/functions/quadratic_cost.h>

#include <algorithm>
#include <cmath>
#include <memory>

namespace corbo {

class MinTimeQuadraticGainScheduled : public StageCost
{
 public:
    MinTimeQuadraticGainScheduled() = default;

    StageCost::Ptr getInstance() const override { return std::make_shared<MinTimeQuadraticGainScheduled>(); }

    bool hasNonIntegralTerms(int k) const override
    {
        return (_gain_quadratic > 1e-2 && _quad_cost.hasNonIntegralTerms(k)) || (_gain_to > 1e-2 && _min_time.hasNonIntegralTerms(k));
    }
    bool hasIntegralTerms(int k) const override { return _gain_quadratic > 1e-2 && _quad_cost.hasIntegralTerms(k); }

    int getNonIntegralDtTermDimension(int k) const override { return _gain_to > 1e-2 ? _min_time.getNonIntegralDtTermDimension(k) : 0; }
    int getNonIntegralStateTermDimension(int k) const override { return _gain_quadratic > 1e-2 ? _quad_cost.getNonIntegralStateTermDimension(k) : 0; }
    int getNonIntegralControlTermDimension(int k) const override
    {
        return _gain_quadratic > 1e-2 ? _quad_cost.getNonIntegralControlTermDimension(k) : 0;
    }
    int getIntegralStateControlTermDimension(int k) const override
    {
        return _gain_quadratic > 1e-2 ? _quad_cost.getIntegralStateControlTermDimension(k) : 0;
    }

    bool isLsqFormNonIntegralStateTerm(int k) const override { return _quad_cost.isLsqFormNonIntegralStateTerm(k); }
    bool isLsqFormNonIntegralControlTerm(int k) const override { return _quad_cost.isLsqFormNonIntegralControlTerm(k); }

    bool update(int n, double t, ReferenceTrajectoryInterface& xref, ReferenceTrajectoryInterface& uref, ReferenceTrajectoryInterface* sref,
                bool single_dt, const Eigen::VectorXd& x0, StagePreprocessor::Ptr stage_preprocessor, const std::vector<double>& dts,
                const DiscretizationGridInterface* grid) override
    {
        double prev_gain_quad = _gain_quadratic;
        double prev_gain_to   = _gain_to;

        if (_gain_quadratic > 1e-2)
        {
            // revert scaling to reset default weights
            // TODO(roesmann) maybe better backup the default matrices
            _quad_cost.scaleCurrentWeightQ(1.0 / _gain_quadratic);
            _quad_cost.scaleCurrentWeightR(1.0 / _gain_quadratic);
        }

        bool changed = false;

        if (_min_time.update(n, t, xref, uref, sref, single_dt, x0, stage_preprocessor, dts, grid)) changed = true;
        if (_quad_cost.update(n, t, xref, uref, sref, single_dt, x0, stage_preprocessor, dts, grid)) changed = true;

        // compute distance to goal
        Eigen::VectorXd xd = xref.getReferenceCached(n) - x0;
        double dist        = xd.transpose() * _quad_cost.getWeightQ() * xd;  // TODO(roesmann): use qf?

        // calculate the tan hyperbolica to approximate a sigmoid function
        // shift the tanh by "3" in order to align it approx. at the origin, such that f(0)=0;
        double aux = 0.5 * std::tanh(_gamma * dist - 3);

        _gain_to        = aux + 0.5;   // normalize to [0,1]
        _gain_quadratic = -aux + 0.5;  // normalize to [0,1]

        if (_gain_quadratic > 1e-2)
        {
            // consider as active
            _quad_cost.scaleCurrentWeightQ(_gain_quadratic);
            _quad_cost.scaleCurrentWeightR(_gain_quadratic);
        }

        // update time-optimal weight
        if (single_dt)
            _min_time.setCustomWeight((double)(n - 1) * _gain_to);
        else
            _min_time.setCustomWeight(_gain_to);  // default weight is 1

        if (changed) return true;

        if (_gain_to > 1e-2 && prev_gain_to <= 1e-2)
        {
            return true;
        }
        if (_gain_to <= 1e-2 && prev_gain_to > 1e-2)
        {
            return true;
        }
        if (_gain_quadratic > 1e-2 && prev_gain_quad <= 1e-2)
        {
            return true;
        }
        if (_gain_quadratic <= 1e-2 && prev_gain_quad > 1e-2)
        {
            return true;
        }
        return false;
    }

    void computeNonIntegralDtTerm(int k, double dt, Eigen::Ref<Eigen::VectorXd> cost) const override
    {
        _min_time.computeNonIntegralDtTerm(k, dt, cost);
    }

    void computeNonIntegralStateTerm(int k, const Eigen::Ref<const Eigen::VectorXd>& x_k, Eigen::Ref<Eigen::VectorXd> cost) const override
    {
        _quad_cost.computeNonIntegralStateTerm(k, x_k, cost);
    }

    void computeNonIntegralControlTerm(int k, const Eigen::Ref<const Eigen::VectorXd>& u_k, Eigen::Ref<Eigen::VectorXd> cost) const override
    {
        _quad_cost.computeNonIntegralControlTerm(k, u_k, cost);
    }

    void computeIntegralStateControlTerm(int k, const Eigen::Ref<const Eigen::VectorXd>& x_k, const Eigen::Ref<const Eigen::VectorXd>& u_k,
                                         Eigen::Ref<Eigen::VectorXd> cost) const override
    {
        _quad_cost.computeIntegralStateControlTerm(k, x_k, u_k, cost);
    }

#ifdef MESSAGE_SUPPORT
    bool fromMessage(const messages::StageCost& message, std::stringstream* issues) override
    {
        bool success = true;

        const messages::MinTimeQuadraticGainScheduled& msg = message.min_time_quad_form_gain_scheduled();

        success = success && _min_time.fromMessage(msg.minimum_time(), issues);
        success = success && _quad_cost.fromMessage(msg.quadratic_form(), issues);

        _gamma = msg.scale();

        return success;
    }

    void toMessage(messages::StageCost& message) const override
    {
        messages::MinTimeQuadraticGainScheduled* msg = message.mutable_min_time_quad_form_gain_scheduled();

        _min_time.toMessage(*msg->mutable_minimum_time());
        _quad_cost.toMessage(*msg->mutable_quadratic_form());

        msg->set_scale(_gamma);
    }
#endif

 protected:
    MinimumTime _min_time;
    QuadraticFormCost _quad_cost;

    double _gamma = 0.1;

    double _gain_to        = 1;
    double _gain_quadratic = 0;
};
FACTORY_REGISTER_STAGE_COST(MinTimeQuadraticGainScheduled)

class MinTimeQuadratic : public StageCost
{
 public:
    MinTimeQuadratic() = default;

    StageCost::Ptr getInstance() const override { return std::make_shared<MinTimeQuadratic>(); }

    bool hasNonIntegralTerms(int k) const override
    {
        return (k >= _quad_k_min && _quad_cost.hasNonIntegralTerms(k)) || _min_time.hasNonIntegralTerms(k);
    }
    bool hasIntegralTerms(int k) const override { return k >= _quad_k_min && _quad_cost.hasIntegralTerms(k); }

    int getNonIntegralDtTermDimension(int k) const override { return _min_time.getNonIntegralDtTermDimension(k); }
    int getNonIntegralStateTermDimension(int k) const override { return k >= _quad_k_min ? _quad_cost.getNonIntegralStateTermDimension(k) : 0; }
    int getNonIntegralControlTermDimension(int k) const override { return k >= _quad_k_min ? _quad_cost.getNonIntegralControlTermDimension(k) : 0; }
    int getIntegralStateControlTermDimension(int k) const override
    {
        return k >= _quad_k_min ? _quad_cost.getIntegralStateControlTermDimension(k) : 0;
    }

    bool isLsqFormNonIntegralStateTerm(int k) const override { return _quad_cost.isLsqFormNonIntegralStateTerm(k); }
    bool isLsqFormNonIntegralControlTerm(int k) const override { return _quad_cost.isLsqFormNonIntegralControlTerm(k); }

    bool update(int n, double t, ReferenceTrajectoryInterface& xref, ReferenceTrajectoryInterface& uref, ReferenceTrajectoryInterface* sref,
                bool single_dt, const Eigen::VectorXd& x0, StagePreprocessor::Ptr stage_preprocessor, const std::vector<double>& dts,
                const DiscretizationGridInterface* grid) override
    {
        int quad_k_min = 0;
        if (_only_last_n > 0)
        {
            quad_k_min = std::max(n - _only_last_n, 0);  // TODO(roesmann) or plus 1
        }

        bool changed = false;
        if (quad_k_min != _quad_k_min)
        {
            changed     = true;
            _quad_k_min = quad_k_min;
        }

        if (_min_time.update(n, t, xref, uref, sref, single_dt, x0, stage_preprocessor, dts, grid)) changed = true;
        if (_quad_cost.update(n, t, xref, uref, sref, single_dt, x0, stage_preprocessor, dts, grid)) changed = true;

        return changed;
    }

    void computeNonIntegralDtTerm(int k, double dt, Eigen::Ref<Eigen::VectorXd> cost) const override
    {
        _min_time.computeNonIntegralDtTerm(k, dt, cost);
    }

    void computeNonIntegralStateTerm(int k, const Eigen::Ref<const Eigen::VectorXd>& x_k, Eigen::Ref<Eigen::VectorXd> cost) const override
    {
        assert(k >= _quad_k_min);
        _quad_cost.computeNonIntegralStateTerm(k, x_k, cost);
    }

    void computeNonIntegralControlTerm(int k, const Eigen::Ref<const Eigen::VectorXd>& u_k, Eigen::Ref<Eigen::VectorXd> cost) const override
    {
        assert(k >= _quad_k_min);
        _quad_cost.computeNonIntegralControlTerm(k, u_k, cost);
    }

    void computeIntegralStateControlTerm(int k, const Eigen::Ref<const Eigen::VectorXd>& x_k, const Eigen::Ref<const Eigen::VectorXd>& u_k,
                                         Eigen::Ref<Eigen::VectorXd> cost) const override
    {
        assert(k >= _quad_k_min);
        _quad_cost.computeIntegralStateControlTerm(k, x_k, u_k, cost);
    }

#ifdef MESSAGE_SUPPORT
    bool fromMessage(const messages::StageCost& message, std::stringstream* issues) override
    {
        bool success = true;

        const messages::MinTimeQuadratic& msg = message.min_time_quad_form();

        success = success && _min_time.fromMessage(msg.minimum_time(), issues);
        success = success && _quad_cost.fromMessage(msg.quadratic_form(), issues);

        _only_last_n = msg.only_last_n();

        return success;
    }

    void toMessage(messages::StageCost& message) const override
    {
        messages::MinTimeQuadratic* msg = message.mutable_min_time_quad_form();

        _min_time.toMessage(*msg->mutable_minimum_time());
        _quad_cost.toMessage(*msg->mutable_quadratic_form());

        msg->set_only_last_n(_only_last_n);
    }
#endif

 protected:
    MinimumTime _min_time;
    QuadraticFormCost _quad_cost;

    int _only_last_n = 0;

    int _quad_k_min = 0;
};
FACTORY_REGISTER_STAGE_COST(MinTimeQuadratic)

}  // namespace corbo

#endif  // SRC_OPTIMAL_CONTROL_INCLUDE_CORBO_OPTIMAL_CONTROL_FUNCTIONS_HYBRID_COST_H_
