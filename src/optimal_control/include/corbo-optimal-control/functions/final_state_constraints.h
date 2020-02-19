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

#ifndef SRC_OPTIMAL_CONTROL_INCLUDE_CORBO_OPTIMAL_CONTROL_FUNCTIONS_FINAL_STATE_CONSTRAINTS_H_
#define SRC_OPTIMAL_CONTROL_INCLUDE_CORBO_OPTIMAL_CONTROL_FUNCTIONS_FINAL_STATE_CONSTRAINTS_H_

#include <corbo-optimal-control/functions/final_state_cost.h>

#include <corbo-core/reference_trajectory.h>
#include <corbo-optimal-control/functions/stage_functions.h>
#include <corbo-systems/system_dynamics_interface.h>

#include <memory>

namespace corbo {

class TerminalBall : public FinalStageConstraint
{
 public:
    using Ptr      = std::shared_ptr<TerminalBall>;
    using ConstPtr = std::shared_ptr<const TerminalBall>;

    TerminalBall() = default;

    TerminalBall(const Eigen::Ref<const Eigen::MatrixXd>& S, double gamma)
    {
        setWeightS(S);
        setGamma(gamma);
    }

    FinalStageConstraint::Ptr getInstance() const override { return std::make_shared<TerminalBall>(); }

    bool isEqualityConstraint() const override { return false; }

    int getNonIntegralStateTermDimension(int k) const override { return 1; }

    bool setWeightS(const Eigen::Ref<const Eigen::MatrixXd>& S);
    bool setWeightS(const Eigen::DiagonalMatrix<double, -1>& S);
    const Eigen::MatrixXd& getWeightS() const { return _S; }

    void setGamma(double gamma) { _gamma = gamma; }
    double getGamma() { return _gamma; }

    void computeNonIntegralStateTerm(int k, const Eigen::Ref<const Eigen::VectorXd>& x_k, Eigen::Ref<Eigen::VectorXd> cost) const override;

    bool update(int n, double t, ReferenceTrajectoryInterface& xref, ReferenceTrajectoryInterface& uref, ReferenceTrajectoryInterface* sref,
                bool single_dt, const Eigen::VectorXd& x0, FinalStageCost::ConstPtr final_stage_cost, StagePreprocessor::Ptr stage_preprocessor,
                const std::vector<double>& dts, const DiscretizationGridInterface* /*grid*/) override
    {
        _x_ref = &xref;

        _zero_x_ref = _x_ref->isZero();

        return false;
    }

    bool checkParameters(int state_dim, int control_dim, FinalStageCost::ConstPtr final_stage_cost, std::stringstream* issues) const override;

#ifdef MESSAGE_SUPPORT
    bool fromMessage(const messages::FinalStageConstraints& message, std::stringstream* issues) override;
    void toMessage(messages::FinalStageConstraints& message) const override;
#endif

 protected:
    Eigen::MatrixXd _S;
    Eigen::DiagonalMatrix<double, -1> _S_diag;
    double _gamma = 0.0;

    const ReferenceTrajectoryInterface* _x_ref = nullptr;
    bool _zero_x_ref                           = false;

    bool _diagonal_mode               = false;
    bool _diagonal_mode_intentionally = false;
};
FACTORY_REGISTER_FINAL_STAGE_CONSTRAINT(TerminalBall)

class TerminalBallInheritFromCost : public TerminalBall
{
 public:
    using Ptr      = std::shared_ptr<TerminalBallInheritFromCost>;
    using ConstPtr = std::shared_ptr<const TerminalBallInheritFromCost>;

    TerminalBallInheritFromCost() = default;

    FinalStageConstraint::Ptr getInstance() const override { return std::make_shared<TerminalBallInheritFromCost>(); }

    int getNonIntegralStateTermDimension(int k) const override { return 1; }

    void setWeightS(const Eigen::Ref<const Eigen::MatrixXd>& S) = delete;
    bool setWeightS(const Eigen::DiagonalMatrix<double, -1>& S) = delete;

    bool update(int n, double t, ReferenceTrajectoryInterface& xref, ReferenceTrajectoryInterface& uref, ReferenceTrajectoryInterface* sref,
                bool single_dt, const Eigen::VectorXd& x0, FinalStageCost::ConstPtr final_stage_cost, StagePreprocessor::Ptr stage_preprocessor,
                const std::vector<double>& dts, const DiscretizationGridInterface* /*grid*/) override;

    void computeNonIntegralStateTerm(int k, const Eigen::Ref<const Eigen::VectorXd>& x_k, Eigen::Ref<Eigen::VectorXd> cost) const override;

    bool checkParameters(int state_dim, int control_dim, FinalStageCost::ConstPtr final_stage_cost, std::stringstream* issues) const override;

#ifdef MESSAGE_SUPPORT
    bool fromMessage(const messages::FinalStageConstraints& message, std::stringstream* issues) override;
    void toMessage(messages::FinalStageConstraints& message) const override;
#endif

    double _gamma = 0.0;
};
FACTORY_REGISTER_FINAL_STAGE_CONSTRAINT(TerminalBallInheritFromCost)

class TerminalEqualityConstraint : public FinalStageConstraint
{
 public:
    using Ptr      = std::shared_ptr<TerminalEqualityConstraint>;
    using ConstPtr = std::shared_ptr<const TerminalEqualityConstraint>;

    TerminalEqualityConstraint() = default;

    TerminalEqualityConstraint(const Eigen::Ref<const Eigen::VectorXd>& xref) { setXRef(xref); }

    FinalStageConstraint::Ptr getInstance() const override { return std::make_shared<TerminalEqualityConstraint>(); }

    bool isEqualityConstraint() const override { return true; }

    int getNonIntegralStateTermDimension(int k) const override { return _xref.size(); }

    void setXRef(const Eigen::Ref<const Eigen::VectorXd>& xref) { _xref = xref; }
    const Eigen::VectorXd& getXRef() const { return _xref; }

    void computeNonIntegralStateTerm(int k, const Eigen::Ref<const Eigen::VectorXd>& x_k, Eigen::Ref<Eigen::VectorXd> cost) const override
    {
        assert(x_k.size() == _xref.size());
        assert(cost.size() == _xref.size());
        cost = x_k - _xref;
    }

    bool update(int n, double t, ReferenceTrajectoryInterface& xref, ReferenceTrajectoryInterface& uref, ReferenceTrajectoryInterface* sref,
                bool single_dt, const Eigen::VectorXd& x0, FinalStageCost::ConstPtr final_stage_cost, StagePreprocessor::Ptr stage_preprocessor,
                const std::vector<double>& dts, const DiscretizationGridInterface* /*grid*/) override
    {
        return false;
    }

    bool checkParameters(int state_dim, int control_dim, FinalStageCost::ConstPtr final_stage_cost, std::stringstream* issues) const override
    {
        if (state_dim != _xref.size())
        {
            if (issues)
                *issues << "TerminalEqualityConstraint: Dimension of xref (" << _xref.size() << ") does not coincide with state dimension ("
                        << state_dim << ")." << std::endl;
        }
        return true;
    }

#ifdef MESSAGE_SUPPORT
    bool fromMessage(const messages::TerminalEqualityConstraint& message, std::stringstream* issues)
    {
        _xref = Eigen::Map<const Eigen::VectorXd>(message.x_ref().data(), message.x_ref_size());
        return true;
    }
    void toMessage(messages::TerminalEqualityConstraint& message) const
    {
        message.mutable_x_ref()->Resize(_xref.size(), 0);
        Eigen::Map<Eigen::VectorXd>(message.mutable_x_ref()->mutable_data(), _xref.size()) = _xref;
    }

    bool fromMessage(const messages::FinalStageConstraints& message, std::stringstream* issues) override
    {
        return fromMessage(message.terminal_equality_constraint(), issues);
    }
    void toMessage(messages::FinalStageConstraints& message) const override { toMessage(*message.mutable_terminal_equality_constraint()); }
#endif

 protected:
    Eigen::VectorXd _xref;
};
FACTORY_REGISTER_FINAL_STAGE_CONSTRAINT(TerminalEqualityConstraint)

class TerminalPartialEqualityConstraint : public TerminalEqualityConstraint
{
 public:
    using Ptr      = std::shared_ptr<TerminalPartialEqualityConstraint>;
    using ConstPtr = std::shared_ptr<const TerminalPartialEqualityConstraint>;

    TerminalPartialEqualityConstraint() = default;
    explicit TerminalPartialEqualityConstraint(const Eigen::Ref<const Eigen::Matrix<bool, -1, 1>>& active_components)
    {
        setActiveComponents(active_components);
    }

    FinalStageConstraint::Ptr getInstance() const override { return std::make_shared<TerminalPartialEqualityConstraint>(); }

    bool isEqualityConstraint() const override { return true; }

    int getNonIntegralStateTermDimension(int k) const override { return _num_active; }

    void setActiveComponents(const Eigen::Ref<const Eigen::Matrix<bool, -1, 1>>& active_components)
    {
        _active_components = active_components;
        _num_active        = _active_components.count();
    }

    void setXRef(const Eigen::Ref<const Eigen::VectorXd>& xref, const Eigen::Ref<const Eigen::Matrix<bool, -1, 1>>& xfixed)
    {
        _xref = xref;
        setActiveComponents(xfixed);
    }
    const Eigen::VectorXd& getXRef() const { return _xref; }
    const Eigen::Matrix<bool, -1, 1>& getXFixed() const { return _active_components; }

    void computeNonIntegralStateTerm(int k, const Eigen::Ref<const Eigen::VectorXd>& x_k, Eigen::Ref<Eigen::VectorXd> cost) const override
    {
        assert(x_k.size() == _xref.size());
        assert(cost.size() == _num_active);
        assert(_active_components.size() == _xref.size());
        assert(_num_active == _active_components.count());
        int idx = 0;
        for (int i = 0; i < _active_components.size(); ++i)
        {
            if (_active_components[i])
            {
                cost[idx] = x_k[i] - _xref[i];
                ++idx;
            }
        }
        assert(idx == _num_active);
    }

    bool update(int n, double t, ReferenceTrajectoryInterface& xref, ReferenceTrajectoryInterface& uref, ReferenceTrajectoryInterface* sref,
                bool single_dt, const Eigen::VectorXd& x0, FinalStageCost::ConstPtr final_stage_cost, StagePreprocessor::Ptr stage_preprocessor,
                const std::vector<double>& dts, const DiscretizationGridInterface* /*grid*/) override
    {
        return false;
    }

    bool checkParameters(int state_dim, int control_dim, FinalStageCost::ConstPtr final_stage_cost, std::stringstream* issues) const override
    {
        if (state_dim != _xref.size())
        {
            if (issues)
                *issues << "TerminalEqualityConstraint: Dimension of xref (" << _xref.size() << ") does not coincide with state dimension ("
                        << state_dim << ")." << std::endl;
        }
        if (state_dim != _active_components.size())
        {
            if (issues)
                *issues << "TerminalEqualityConstraint: Dimension of active_components (" << _active_components.size()
                        << ") does not coincide with state dimension (" << state_dim << ")." << std::endl;
        }
        return true;
    }

#ifdef MESSAGE_SUPPORT
    bool fromMessage(const messages::TerminalPartialEqualityConstraint& message, std::stringstream* issues)
    {
        _xref = Eigen::Map<const Eigen::VectorXd>(message.x_ref().data(), message.x_ref_size());

        setActiveComponents(Eigen::Map<const Eigen::Matrix<bool, -1, 1>>(message.active_components().data(), message.active_components_size()));
        return true;
    }
    void toMessage(messages::TerminalPartialEqualityConstraint& message) const
    {
        message.mutable_x_ref()->Resize(_xref.size(), 0);
        Eigen::Map<Eigen::VectorXd>(message.mutable_x_ref()->mutable_data(), _xref.size()) = _xref;

        // xf fixed states
        if (_active_components.size() > 0)
        {
            message.mutable_active_components()->Resize(_active_components.size(), false);
            Eigen::Map<Eigen::Matrix<bool, -1, 1>>(message.mutable_active_components()->mutable_data(), _active_components.size()) =
                _active_components;
        }
    }

    bool fromMessage(const messages::FinalStageConstraints& message, std::stringstream* issues) override
    {
        return fromMessage(message.terminal_partial_eq_constr(), issues);
    }
    void toMessage(messages::FinalStageConstraints& message) const override { toMessage(*message.mutable_terminal_partial_eq_constr()); }
#endif

 protected:
    Eigen::Matrix<bool, -1, 1> _active_components;
    int _num_active = 0;
};
FACTORY_REGISTER_FINAL_STAGE_CONSTRAINT(TerminalPartialEqualityConstraint)

}  // namespace corbo

#endif  // SRC_OPTIMAL_CONTROL_INCLUDE_CORBO_OPTIMAL_CONTROL_FUNCTIONS_FINAL_STATE_CONSTRAINTS_H_
