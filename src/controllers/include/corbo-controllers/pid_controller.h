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

#ifndef SRC_CONTROLLERS_INCLUDE_CORBO_CONTROLLERS_PID_CONTROLLER_H_
#define SRC_CONTROLLERS_INCLUDE_CORBO_CONTROLLERS_PID_CONTROLLER_H_

#include <corbo-controllers/controller_interface.h>
#include <memory>

namespace corbo {

/**
 * @brief PID controller
 *
 * @ingroup controllers
 *
 * Implementation of a discrete-time PID controller.
 *
 * @see ControllerInterface LqrController PredictiveController
 *
 * @author Christoph Rösmann (christoph.roesmann@tu-dortmund.de)
 *
 * @todo Anti-windup implementation
 */
class PidController : public ControllerInterface
{
 public:
    PidController();

    // implements interface method
    int getControlInputDimension() const override { return _num_parallel_pid; }
    // implements interface method
    int getStateDimension() const override { return _num_parallel_pid; }
    // implements interface method
    bool hasPiecewiseConstantControls() const override { return false; }
    // implements interface method
    bool providesFutureControls() const override { return false; }
    // implements interface method
    bool providesFutureStates() const override { return false; }

    // implements interface method
    Ptr getInstance() const override { return std::make_shared<PidController>(); }
    //! Return a newly created shared instance of the implemented class (static method)
    static Ptr getInstanceStatic() { return std::make_shared<PidController>(); }

    const double& getGainP() const { return _p_gain; }
    //! Set proportional gain
    void setGainP(double p_gain) { _p_gain = p_gain; }
    //! Set integral gain
    const double& getGainI() const { return _i_gain; }
    void setGainI(double i_gain) { _i_gain = i_gain; }
    //! Set differential gain
    const double& getGainD() const { return _d_gain; }
    void setGainD(double d_gain) { _d_gain = d_gain; }

    void setNumParallelPid(int num_parallel_pid)
    {
        _num_parallel_pid = num_parallel_pid;
        _p_error.resize(_num_parallel_pid, 0.0);
        _i_error.resize(_num_parallel_pid, 0.0);
        _d_error.resize(_num_parallel_pid, 0.0);
    }

    // implements interface method
    bool step(const StateVector& x, ReferenceTrajectoryInterface& xref, ReferenceTrajectoryInterface& uref, const Duration& dt, const Time& t,
              TimeSeries::Ptr u_sequence, TimeSeries::Ptr x_sequence, SignalTargetInterface* signal_target = nullptr,
              ReferenceTrajectoryInterface* sref = nullptr, ReferenceTrajectoryInterface* xinit = nullptr,
              ReferenceTrajectoryInterface* uinit = nullptr, const std::string& ns = "") override;

    // implements interface method
    void getAvailableSignals(SignalTargetInterface& signal_target, const std::string& ns = "") const override;

    ControllerStatistics::Ptr getStatistics() const override { return {}; }

#ifdef MESSAGE_SUPPORT
    // implements interface method
    void toMessage(corbo::messages::Controller& message) const override;
    // implements interface method
    void fromMessage(const corbo::messages::Controller& message, std::stringstream* issues = nullptr) override;
#endif

    // implements interface method
    void reset() override;

    //! Specify whether the control error should be published via signal
    void publishError(bool active) { _publish_error = active; }

 private:
    // parameters
    double _p_gain = 0.2;
    double _i_gain = 0.0;
    double _d_gain = 0.0;

    // internal states
    std::vector<double> _p_error;
    std::vector<double> _i_error;
    std::vector<double> _d_error;
    // double _control_error_last = 0;

    int _num_parallel_pid = 1;

    bool _publish_error = true;
};

FACTORY_REGISTER_CONTROLLER(PidController)

}  // namespace corbo

#endif  // SRC_CONTROLLERS_INCLUDE_CORBO_CONTROLLERS_PID_CONTROLLER_H_
