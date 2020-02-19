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

#ifndef SRC_CONTROLLERS_INCLUDE_CORBO_CONTROLLERS_DUAL_MODE_CONTROLLER_H_
#define SRC_CONTROLLERS_INCLUDE_CORBO_CONTROLLERS_DUAL_MODE_CONTROLLER_H_

#include <corbo-controllers/controller_interface.h>
#include <corbo-controllers/predictive_controller.h>

#include <memory>

namespace corbo {

/**
 * @brief Dual mode controller
 *
 * @ingroup controllers
 *
 * Use predictive controller until some conditions are satisfied (terminal region, min dt)
 * and then switch to an LQR.
 *
 * @see ControllerInterface LqrController PredictiveController
 *
 * @author Christoph Rösmann (christoph.roesmann@tu-dortmund.de)
 */
class DualModeController : public ControllerInterface
{
 public:
    using Ptr = std::shared_ptr<DualModeController>;

    DualModeController() = default;

    // implements interface method
    int getControlInputDimension() const override { return _pred_controller.getControlInputDimension(); }
    // implements interface method
    int getStateDimension() const override { return _pred_controller.getStateDimension(); }
    // implements interface method
    bool hasPiecewiseConstantControls() const override
    {
        return _local_ctrl_active ? _local_controller->hasPiecewiseConstantControls() : _pred_controller.hasPiecewiseConstantControls();
    }
    // implements interface method
    bool providesFutureControls() const override
    {
        return _local_ctrl_active ? _local_controller->providesFutureControls() : _pred_controller.providesFutureControls();
    }  // TODO(roesmann): how should we handle this with the LQR?
    // implements interface method
    bool providesFutureStates() const override
    {
        return _local_ctrl_active ? _local_controller->providesFutureStates() : _pred_controller.providesFutureStates();
    }  // TODO(roesmann): how should we handle this with the LQR?

    // implements interface method
    ControllerInterface::Ptr getInstance() const override { return std::make_shared<DualModeController>(); }
    static Ptr getInstanceStatic() { return std::make_shared<DualModeController>(); }

    // implements interface method
    bool initialize(const StateVector& x, ReferenceTrajectoryInterface& expected_xref, ReferenceTrajectoryInterface& expected_uref,
                    const Duration& expected_dt, const Time& t, ReferenceTrajectoryInterface* sref = nullptr) override;

    // implements interface method
    bool step(const StateVector& x, ReferenceTrajectoryInterface& xref, ReferenceTrajectoryInterface& uref, const Duration& dt, const Time& t,
              TimeSeries::Ptr u_sequence, TimeSeries::Ptr x_sequence, SignalTargetInterface* signal_target = nullptr,
              ReferenceTrajectoryInterface* sref = nullptr, ReferenceTrajectoryInterface* xinit = nullptr,
              ReferenceTrajectoryInterface* uinit = nullptr, const std::string& ns = "") override;

    // implements interface method
    double getControlDuration() const override
    {
        return _local_ctrl_active ? _local_controller->getControlDuration() : _pred_controller.getControlDuration();
    }

    bool setWeightS(const Eigen::Ref<const Eigen::MatrixXd>& S);

    void setLocalController(ControllerInterface::Ptr local_controller) { _local_controller = local_controller; }

    // implements interface method
    bool supportsAsynchronousControl() const override { return true; }

    // implements interface method
    void getAvailableSignals(SignalTargetInterface& signal_target, const std::string& ns = "") const override;

    // implements interface method
    void reset() override;

    void sendSignals(double t, SignalTargetInterface& signal_target, const std::string& ns = "") const override;

    ControllerStatistics::Ptr getStatistics() const override
    {
        return std::make_shared<ControllerStatistics>(_statistics);
    }  // TODO(roesmann): make_shared?

#ifdef MESSAGE_SUPPORT
    void toMessage(corbo::messages::DualModeController& message) const;
    void fromMessage(const corbo::messages::DualModeController& message, std::stringstream* issues = nullptr);

    // implements interface method
    void toMessage(corbo::messages::Controller& message) const override { toMessage(*message.mutable_dual_mode_controller()); }
    // implements interface method
    void fromMessage(const corbo::messages::Controller& message, std::stringstream* issues = nullptr) override
    {
        fromMessage(message.dual_mode_controller(), issues);
    }
#endif

 protected:
    bool isInsideInTerminalBall(const Eigen::Ref<const Eigen::VectorXd>& x0, const Eigen::Ref<const Eigen::VectorXd>& xf) const;

 private:
    ControllerStatistics _statistics;

    PredictiveController _pred_controller;
    ControllerInterface::Ptr _local_controller;

    bool _switch_dt            = false;
    bool _switch_terminal_ball = false;

    bool _local_ctrl_active = false;

    Eigen::MatrixXd _S;
    double _gamma = 0.0;

    double _min_dt = 0;

    bool _first_run = true;

    bool _initialized = false;
};

FACTORY_REGISTER_CONTROLLER(DualModeController)

}  // namespace corbo

#endif  // SRC_CONTROLLERS_INCLUDE_CORBO_CONTROLLERS_DUAL_MODE_CONTROLLER_H_
