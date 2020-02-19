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

#ifndef SRC_OPTIMAL_CONTROL_INCLUDE_CORBO_OPTIMAL_CONTROL_FUNCTIONS_MINIMUM_TIME_H_
#define SRC_OPTIMAL_CONTROL_INCLUDE_CORBO_OPTIMAL_CONTROL_FUNCTIONS_MINIMUM_TIME_H_

#include <corbo-optimal-control/functions/stage_functions.h>

#include <cmath>
#include <memory>

namespace corbo {

class MinimumTime : public StageCost
{
 public:
    using Ptr = std::shared_ptr<MinimumTime>;

    MinimumTime() = default;

    MinimumTime(bool lsq_form) : _lsq_form(lsq_form) {}

    StageCost::Ptr getInstance() const override { return std::make_shared<MinimumTime>(); }

    bool hasNonIntegralTerms(int k) const override { return true; }
    bool hasIntegralTerms(int k) const override { return false; }

    int getNonIntegralDtTermDimension(int k) const override { return (k == 0 || !_single_dt) ? 1 : 0; }
    bool isLsqFormNonIntegralDtTerm(int k) const override { return _lsq_form; }

    bool update(int n, double /*t*/, ReferenceTrajectoryInterface& /*xref*/, ReferenceTrajectoryInterface& /*uref*/,
                ReferenceTrajectoryInterface* /*sref*/, bool single_dt, const Eigen::VectorXd& x0, StagePreprocessor::Ptr stage_preprocessor,
                const std::vector<double>& dts, const DiscretizationGridInterface* /*grid*/) override
    {
        if (!_custom_weight)
        {
            _single_dt = single_dt;
            if (single_dt)
                _weight = _lsq_form ? std::sqrt(n - 1) : (n - 1);
            else
            {
                _weight = _lsq_form ? std::sqrt(n - 1) : 1.0;
                //   _weight = 1.0;
            }
        }
        return false;
    }

    void computeNonIntegralDtTerm(int k, double dt, Eigen::Ref<Eigen::VectorXd> cost) const override
    {
        if (!_single_dt || k == 0)
            cost[0] = _weight * dt;
        else
        {
            PRINT_DEBUG_NAMED("this method should not be called in single_dt mode and k>0");
        }
    }

    void setCustomWeight(double weight)
    {
        _custom_weight = true;
        _weight        = weight;
    }
    double getCustomWeight() const { return _weight; }

    void setLsqForm(bool lsq_form) { _lsq_form = lsq_form; }

#ifdef MESSAGE_SUPPORT
    virtual bool fromMessage(const messages::MinimumTime& message, std::stringstream* issues)
    {
        _lsq_form = message.lsq_form();
        return true;
    }
    virtual void toMessage(messages::MinimumTime& message) const { message.set_lsq_form(_lsq_form); }

    bool fromMessage(const messages::StageCost& message, std::stringstream* issues) override { return fromMessage(message.minimum_time(), issues); }
    void toMessage(messages::StageCost& message) const override { toMessage(*message.mutable_minimum_time()); }
#endif

 protected:
    bool _lsq_form      = false;
    double _weight      = 1.0;
    bool _single_dt     = false;
    bool _custom_weight = false;
};
FACTORY_REGISTER_STAGE_COST(MinimumTime)

class MinimumTimeRegularized : public StageCost
{
 public:
    using Ptr = std::shared_ptr<MinimumTimeRegularized>;

    MinimumTimeRegularized() = default;

    StageCost::Ptr getInstance() const override { return std::make_shared<MinimumTimeRegularized>(); }

    bool hasNonIntegralTerms(int k) const override { return true; }
    bool hasIntegralTerms(int k) const override { return false; }

    int getNonIntegralDtTermDimension(int k) const override { return 1; }
    bool isLsqFormNonIntegralDtTerm(int k) const override { return false; }

    bool update(int n, double /*t*/, ReferenceTrajectoryInterface& /*xref*/, ReferenceTrajectoryInterface& /*uref*/,
                ReferenceTrajectoryInterface* /*sref*/, bool /*single_dt*/, const Eigen::VectorXd& /*x0*/,
                StagePreprocessor::Ptr /*stage_preprocessor*/, const std::vector<double>& /*dts*/,
                const DiscretizationGridInterface* /*grid*/) override
    {
        return false;
    }

    void computeNonIntegralDtTerm(int k, double dt, Eigen::Ref<Eigen::VectorXd> cost) const override { cost[0] = dt + _reg_factor * dt * dt; }

    void setRegularizationFactor(double reg_factor) { _reg_factor = reg_factor; }

#ifdef MESSAGE_SUPPORT
    virtual bool fromMessage(const messages::MinimumTimeRegularized& message, std::stringstream* issues)
    {
        _reg_factor = message.reg_factor();
        return true;
    }
    virtual void toMessage(messages::MinimumTimeRegularized& message) const { message.set_reg_factor(_reg_factor); }

    bool fromMessage(const messages::StageCost& message, std::stringstream* issues) override
    {
        return fromMessage(message.minimum_time_regularized(), issues);
    }
    void toMessage(messages::StageCost& message) const override { toMessage(*message.mutable_minimum_time_regularized()); }
#endif

 protected:
    double _reg_factor = 0.0;
};
FACTORY_REGISTER_STAGE_COST(MinimumTimeRegularized)

}  // namespace corbo

#endif  // SRC_OPTIMAL_CONTROL_INCLUDE_CORBO_OPTIMAL_CONTROL_FUNCTIONS_MINIMUM_TIME_H_
