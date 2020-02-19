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

#ifndef SRC_CONTROLLERS_INCLUDE_CORBO_CONTROLLERS_PREDICTIVE_CONTROLLER_H_
#define SRC_CONTROLLERS_INCLUDE_CORBO_CONTROLLERS_PREDICTIVE_CONTROLLER_H_

#include <corbo-controllers/controller_interface.h>
#include <corbo-optimal-control/optimal_control_problem_interface.h>

#include <memory>

namespace corbo {

/**
 * @brief Predictive controller
 *
 * @ingroup controllers
 *
 * Implementation of a predictive controller that accepts a generic
 * optimal control problem (OCP) that is solved within each step command.
 *
 * The OCP is invoked repeatedly up to a user-specified number of iterations:
 * refer to numOcpIterations().
 *
 * @see ControllerInterface LqrController PidController
 *
 * @author Christoph Rösmann (christoph.roesmann@tu-dortmund.de)
 */
class PredictiveController : public ControllerInterface
{
 public:
    using Ptr = std::shared_ptr<PredictiveController>;

    PredictiveController();

    // implements interface method
    int getControlInputDimension() const override { return _ocp ? _ocp->getControlInputDimension() : 0; }
    // implements interface method
    int getStateDimension() const override { return _ocp ? _ocp->getStateDimension() : 0; }
    // implements interface method
    bool hasPiecewiseConstantControls() const override { return _ocp ? _ocp->isConstantControlAction() : false; }
    // implements interface method
    bool providesFutureControls() const override { return _ocp ? _ocp->providesFutureControls() : false; }
    // implements interface method
    bool providesFutureStates() const override { return _ocp ? _ocp->providesFutureStates() : false; }

    // implements interface method
    ControllerInterface::Ptr getInstance() const override { return std::make_shared<PredictiveController>(); }
    static Ptr getInstanceStatic() { return std::make_shared<PredictiveController>(); }

    void setOptimalControlProblem(OptimalControlProblemInterface::Ptr ocp) { _ocp = ocp; }
    OptimalControlProblemInterface::Ptr getOptimalControlProblem() { return _ocp; }

    // implements interface method
    bool initialize(const StateVector& x, ReferenceTrajectoryInterface& expected_xref, ReferenceTrajectoryInterface& expected_uref,
                    const Duration& expected_dt, const Time& t, ReferenceTrajectoryInterface* sref = nullptr) override;

    // implements interface method
    bool step(const StateVector& x, ReferenceTrajectoryInterface& xref, ReferenceTrajectoryInterface& uref, const Duration& dt, const Time& t,
              TimeSeries::Ptr u_sequence, TimeSeries::Ptr x_sequence, SignalTargetInterface* signal_target = nullptr,
              ReferenceTrajectoryInterface* sref = nullptr, ReferenceTrajectoryInterface* xinit = nullptr,
              ReferenceTrajectoryInterface* uinit = nullptr, const std::string& ns = "") override;

    // implements interface method
    double getControlDuration() const override { return _ocp ? _ocp->getFirstDt() : 0.0; }

    // implements interface method
    bool supportsAsynchronousControl() const override { return true; }

    // implements interface method
    void getAvailableSignals(SignalTargetInterface& signal_target, const std::string& ns = "") const override;

    // implements interface method
    void reset() override;

    const int& getNumOcpIterations() const { return _num_ocp_iterations; }
    void setNumOcpIterations(int ocp_iter) { _num_ocp_iterations = ocp_iter; }

    const bool& isPublishPrediction() const { return _publish_prediction; }
    void setPublishPrediction(bool publish) { _publish_prediction = publish; }

    void setOutputControlSequenceLenght(bool activate) { _output_control_sequence = activate; }
    void setOutputStateSequenceLenght(bool activate) { _output_state_sequence = activate; }

    void setAutoUpdatePreviousControl(bool enable) { _auto_update_prev_control = enable; }
    bool getAutoUpdatePreviousControl() const { return _auto_update_prev_control; }

    void sendSignals(double t, SignalTargetInterface& signal_target, const std::string& ns = "") const override;

    ControllerStatistics::Ptr getStatistics() const override
    {
        return std::make_shared<ControllerStatistics>(_statistics);
    }  // TODO(roesmann): make_shared?

#ifdef MESSAGE_SUPPORT
    void toMessage(corbo::messages::PredictiveController& message) const;
    void fromMessage(const corbo::messages::PredictiveController& message, std::stringstream* issues = nullptr);

    // implements interface method
    void toMessage(corbo::messages::Controller& message) const override { toMessage(*message.mutable_predictive_controller()); }
    // implements interface method
    void fromMessage(const corbo::messages::Controller& message, std::stringstream* issues = nullptr) override
    {
        fromMessage(message.predictive_controller(), issues);
    }
#endif

 protected:
    OptimalControlProblemInterface::Ptr _ocp;
    TimeSeries::Ptr _x_ts;
    TimeSeries::Ptr _u_ts;

    ControllerStatistics _statistics;

    bool _auto_update_prev_control = true;

    int _num_ocp_iterations  = 1;
    bool _publish_prediction = true;
    // double _publish_prediction_dt = -1;
    bool _output_control_sequence = false;
    bool _output_state_sequence   = false;

    bool _initialized = false;
};

FACTORY_REGISTER_CONTROLLER(PredictiveController)

}  // namespace corbo

#endif  // SRC_CONTROLLERS_INCLUDE_CORBO_CONTROLLERS_PREDICTIVE_CONTROLLER_H_
