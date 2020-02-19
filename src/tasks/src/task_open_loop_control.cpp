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

#include <corbo-communication/utilities.h>
#include <corbo-core/console.h>
#include <corbo-core/time.h>
#include <corbo-core/time_series.h>
#include <corbo-tasks/task_open_loop_control.h>
#include <string>

namespace corbo {

OpenLoopControlTask::OpenLoopControlTask()
{
    // default reference
    _xreference = std::make_shared<StaticReference>(Eigen::Matrix<double, 1, 1>::Ones());
    _ureference = std::make_shared<ZeroReference>(1);
}

void OpenLoopControlTask::getAvailableSignals(const Environment& environment, SignalTargetInterface& signal_target, const std::string& ns) const
{
    if (environment.getController())
    {
        environment.getController()->getAvailableSignals(signal_target);

        if (environment.getController()->getControlInputDimension() != property::INHERITED)
            signal_target.registerMeasurement(ns + "control_input", environment.getController()->getControlInputDimension(), {}, true);
    }

    if (environment.getPlant())
    {
        environment.getPlant()->getAvailableSignals(signal_target);
        if (environment.getPlant()->getOutputDimension() != property::INHERITED)
            signal_target.registerMeasurement(ns + "plant_output", environment.getPlant()->getOutputDimension());
    }

    if (environment.getObserver())
    {
        environment.getObserver()->getAvailableSignals(signal_target);
        if (environment.getObserver()->getStateDimension() != property::INHERITED)
            signal_target.registerMeasurement(ns + "observed_states", environment.getObserver()->getStateDimension());
    }
}

void OpenLoopControlTask::performTask(Environment& environment, SignalTargetInterface* signal_target, std::string* msg, const std::string& ns)
{
    // environment->reset();

    // initialize variables
    Time t(0);
    Duration dt(_dt);

    ControllerInterface* controller = environment.getController().get();
    PlantInterface* plant           = environment.getPlant().get();
    ObserverInterface* observer     = environment.getObserver().get();

    using ControlVector = Eigen::VectorXd;
    using StateVector   = Eigen::VectorXd;
    using OutputVector  = Eigen::VectorXd;

    TimeSeries::Ptr u_sequence = std::make_shared<TimeSeries>();
    TimeSeries::Ptr x_sequence = std::make_shared<TimeSeries>();

    ControlVector u(controller->getControlInputDimension());
    OutputVector y(plant->getOutputDimension());
    StateVector x(controller->getStateDimension());

    if (!verify(environment, msg)) return;

    // initialize modules
    if (!controller->initialize(x, *_xreference, *_ureference, dt, t)) PRINT_FATAL("Controller initialization failed.");
    if (!plant->initialize()) PRINT_FATAL("Plant initialization failed.");

    // request open-loop control action
    if (!plant->output(y, t, signal_target, ns)) PRINT_ERROR("OpenLoopControlTask::performTask(): error while retreiving plant output.");
    if (!observer->observe(y, x, Duration(0), t, signal_target, ns)) PRINT_ERROR("OpenLoopControlTask::performTask(): observer error.");
    Time time_pre_step = Time::now();

    // if (!controller->step(x, *_xreference, *_ureference, Duration(0), t, u, signal_target))
    bool controller_success =
        controller->step(x, *_xreference, *_ureference, Duration(0), t, u_sequence, x_sequence, signal_target, nullptr, nullptr, nullptr, ns);
    if (!controller_success)
    {
        u.setZero();
        PRINT_ERROR("OpenLoopControlTask::performTask(): controller error.");
    }

    PRINT_INFO("Controller CPU time: " << (Time::now() - time_pre_step).toSec() * 1000 << " ms.");

    TimeSeries::Interpolation u_interp_strategy =
        controller->hasPiecewiseConstantControls() ? TimeSeries::Interpolation::ZeroOrderHold : TimeSeries::Interpolation::Linear;

    // send controller signals
    if (signal_target) controller->sendSignals(t.toSec(), *signal_target, ns);

    if (!u_sequence || u_sequence->getTimeDimension() < 1)
    {
        PRINT_ERROR("OpenLoopControlTask::performTask(): control does not provide any open-loop control sequence. Canceling.");
        return;
    }

    // set final time
    Time tf(u_sequence->getTime().back() + 1e-12);

    Rate rate(_realtime_sync ? dt : Duration(1e-5));  // wait a small amount also if realtime sync is disabled

    int t_idx = 0;  // just for the case (_dt <= 0)

    // perform actual open-loop task
    while (t <= tf && ok())
    {
        // plant output
        if (!plant->output(y, t, signal_target, ns)) PRINT_ERROR("OpenLoopControlTask::performTask(): error while retreiving plant output.");
        if (signal_target) signal_target->sendMeasurement(ns + "plant_output", t.toSec(), y);

        // observer
        if (!observer->observe(y, x, Duration(0), t, signal_target)) PRINT_ERROR("OpenLoopControlTask::performTask(): observer error.");
        if (signal_target) signal_target->sendMeasurement(ns + "observed_states", t.toSec(), x);

        // get current u
        if (controller_success && !u_sequence->getValuesInterpolate(t.toSec(), u, u_interp_strategy, TimeSeries::Extrapolation::ZeroOrderHold))
            PRINT_ERROR("OpenLoopControlTask::performTask(): control sequence interpolation error.");

        if (signal_target) signal_target->sendMeasurement(ns + "control_input", t.toSec(), u);

        // if dt<=0 -> inherit from open loop sequence
        if (_dt <= 0 && t_idx < u_sequence->getTimeDimension() - 1)
        {
            dt.fromSec(u_sequence->getTimeRef()[t_idx + 1] - u_sequence->getTimeRef()[t_idx]);
            ++t_idx;
        }

        // control plant
        plant->control(u, dt, t, signal_target, ns);
        // plant->control(u_sequence, x_sequence, dt, t, signal_target);
        // this does not work in open-loop, since u_seq is not pruned w.r.t. current
        // TODO(roesmann): add prune method to TimeSeries object

        if (!rate.sleep() && _realtime_sync)
        {
            PRINT_WARNING("OpenLoopControlTask(): rate exceeded (" << rate.lastCycleTime().toSec() << "s/" << dt << "s).");
        }
        t += dt;
    }

    plant->stop();

    signal_target->sendMeasurement(ns + "ctrl_succeess", 0.0, {(double)controller_success});
}

bool OpenLoopControlTask::verify(const Environment& environment, std::string* msg) const
{
    bool ret_val = true;

    if (msg) msg->clear();

    if (environment.getController() && !environment.getController()->providesFutureControls())
    {
        ret_val = false;
        if (msg) *msg += "The provided controller does not support open loop control tasks.\n";
        return ret_val;
    }

    // check if all objects are set
    if (!_xreference)
    {
        ret_val = false;
        if (msg) *msg += "State reference trajectory not specified for OpenLoopControlTask\n";
    }

    // check if all objects are set
    if (!_ureference)
    {
        ret_val = false;
        if (msg) *msg += "Control reference trajectory not specified for OpenLoopControlTask\n";
    }

    // verify environment
    std::string environment_msg;
    ret_val = ret_val && environment.verify(&environment_msg);
    if (msg) *msg += environment_msg;
    if (ret_val == false) return false;  // we need all objects in environment allocated!

    // check reference dimension
    if (environment.getController()->getStateDimension() != _xreference->getDimension())
    {
        ret_val = false;
        if (msg)
            *msg += "State reference trajectory dimension (" + std::to_string(_xreference->getDimension()) +
                    ") does not match controller state dimension (" + std::to_string(environment.getController()->getStateDimension()) + ").\n";
    }
    if (environment.getController()->getControlInputDimension() != _ureference->getDimension())
    {
        ret_val = false;
        if (msg)
            *msg += "Control reference trajectory dimension (" + std::to_string(_ureference->getDimension()) +
                    ") does not match controller control input dimension (" +
                    std::to_string(environment.getController()->getControlInputDimension()) + ").\n";
    }

    if (environment.getPlant()->requiresFutureControls() && !environment.getController()->providesFutureControls())
    {
        ret_val = false;
        if (msg) *msg += "Controller does not support control sequences, that are required by the plant";
    }

    if (environment.getPlant()->requiresFutureStates() && !environment.getController()->providesFutureStates())
    {
        ret_val = false;
        if (msg) *msg += "Controller does not support state sequences, that are required by the plant";
    }

    return ret_val;
}

#ifdef MESSAGE_SUPPORT
void OpenLoopControlTask::toMessage(corbo::messages::OpenLoopControlTask& message) const
{
    message.set_dt(_dt);
    message.set_realtime_sync(_realtime_sync);
    if (_xreference) _xreference->toMessage(*message.mutable_xreference());
    if (_ureference) _ureference->toMessage(*message.mutable_ureference());
}
void OpenLoopControlTask::fromMessage(const corbo::messages::OpenLoopControlTask& message, std::stringstream* issues)
{
    _dt            = message.dt();
    _realtime_sync = message.realtime_sync();

    // xreference
    _xreference.reset();
    if (message.has_xreference())
    {
        // construct object
        std::string type;
        if (util::get_oneof_field_type(message.xreference(), "reference", type, false))
        {
            ReferenceTrajectoryInterface::Ptr xreference = ReferenceTrajectoryFactory::instance().create(type);
            // import parameters
            if (xreference)
            {
                xreference->fromMessage(message.xreference(), issues);
                setStateReference(xreference);
            }
        }
    }

    // ureference
    _ureference.reset();
    if (message.has_ureference())
    {
        // construct object
        std::string type;
        if (util::get_oneof_field_type(message.ureference(), "reference", type, false))
        {
            ReferenceTrajectoryInterface::Ptr ureference = ReferenceTrajectoryFactory::instance().create(type);
            // import parameters
            if (ureference)
            {
                ureference->fromMessage(message.ureference(), issues);
                setControlReference(ureference);
            }
        }
    }
}
#endif

}  // namespace corbo
