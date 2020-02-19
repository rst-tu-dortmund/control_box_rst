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
#include <corbo-tasks/task_closed_loop_control.h>
#include <string>

namespace corbo {

ClosedLoopControlTask::ClosedLoopControlTask()
{
    // default reference
    _xreference = std::make_shared<StaticReference>(Eigen::Matrix<double, 1, 1>::Ones());
    _ureference = std::make_shared<ZeroReference>(1);
}

void ClosedLoopControlTask::getAvailableSignals(const Environment& environment, SignalTargetInterface& signal_target, const std::string& ns) const
{
    if (environment.getController())
    {
        environment.getController()->getAvailableSignals(signal_target);

        if (environment.getController()->getControlInputDimension() != property::INHERITED)
        {
            signal_target.registerMeasurement(ns + "control_input", environment.getController()->getControlInputDimension(), {}, true);
        }
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

    if (_compensate_cpu_time || _compensate_dead_time)
    {
        signal_target.registerMeasurement(ns + "comp_states", environment.getObserver()->getStateDimension());
    }

    if (_xreference)
    {
        signal_target.registerMeasurement(ns + "reference/x", _xreference->getDimension());
    }
    if (_ureference)
    {
        signal_target.registerMeasurement(ns + "reference/u", _ureference->getDimension());
    }
}

void ClosedLoopControlTask::performTask(Environment& environment, SignalTargetInterface* signal_target, std::string* msg, const std::string& ns)
{
    // environment->reset();

    if (_dt <= 0 && !environment.getController()->supportsAsynchronousControl())
    {
        PRINT_ERROR("ClosedLoopControlTask::performTask(): dt <= 0 selected, but current controller does not support asynchronous control.");
        return;
    }

    // initialize variables
    Time t(0);
    Time tf(_sim_time);
    Duration dt(_dt);

    ControllerInterface* controller = environment.getController().get();
    PlantInterface* plant           = environment.getPlant().get();
    ObserverInterface* observer     = environment.getObserver().get();

    using ControlVector = Eigen::VectorXd;
    using StateVector   = Eigen::VectorXd;
    using OutputVector  = Eigen::VectorXd;

    TimeSeries::Ptr u_sequence = std::make_shared<TimeSeries>();
    TimeSeries::Ptr x_sequence = std::make_shared<TimeSeries>();
    u_sequence->add(0, ControlVector::Zero(plant->getOutputDimension()));

    OutputVector y(plant->getOutputDimension());
    StateVector x(controller->getStateDimension());

    Eigen::VectorXd xref(_xreference->getDimension());
    Eigen::VectorXd uref(_ureference->getDimension());

    if (!verify(environment, msg)) return;

    // initialize modules
    if (!controller->initialize(x, *_xreference, *_ureference, dt, t))
    {
        PRINT_ERROR("Controller initialization failed.");
        return;
    }
    if (!plant->initialize())
    {
        PRINT_ERROR("Plant initialization failed.");
        return;
    }

    if ((_compensate_cpu_time || _compensate_dead_time) && !_compensator.initialize())
    {
        PRINT_ERROR("Compensator initialization failed.");
        return;
    }

    Rate rate(_realtime_sync ? dt : Duration(1e-6));  // wait a small amount also if realtime sync is disabled

    Time t_measure_x;
    Time compensator_meas_start;
    Duration cpu_time(0);
    Duration cpu_time_cache(0);

    double comp_dt  = 0;
    double deadtime = 0;

    if (_compensate_dead_time) deadtime = _compensator.getDeadTime();

    std::vector<std::pair<double, Eigen::VectorXd>> useq_predict;
    if (_time_value_buffer.isEmpty()) _time_value_buffer.setInitialValue(ControlVector::Zero(plant->getOutputDimension()));

    // perform actual closed-loop task
    Time t_wall_start = Time::now();

    while (t <= tf && ok())
    {
        Duration last_dt = (t == Time(0)) ? Duration(0) : dt;

        if (_use_wall_time) t = (Time::now() - t_wall_start).toTime();

        // publish reference signals
        if (signal_target)
        {
            _xreference->getReference(t, xref);
            _ureference->getReference(t, uref);
            signal_target->sendMeasurement(ns + "reference/x", t.toSec(), xref);
            signal_target->sendMeasurement(ns + "reference/u", t.toSec(), uref);
        }

        // plant output
        if (_use_wall_time) t = (Time::now() - t_wall_start).toTime();
        if (!plant->output(y, t, signal_target, ns)) PRINT_ERROR("ClosedLoopControl::performTask(): error while retreiving plant output.");
        if (signal_target) signal_target->sendMeasurement(ns + "plant_output", t.toSec(), y);

        t_measure_x            = t;
        compensator_meas_start = Time::now();

        // observer
        if (_use_wall_time) t = (Time::now() - t_wall_start).toTime();
        if (!observer->observe(y, x, last_dt, t, signal_target, ns)) PRINT_ERROR("ClosedLoopControl::performTask(): observer error.");
        if (signal_target) signal_target->sendMeasurement(ns + "observed_states", t.toSec(), x);

        // compensator (we might want to compensate for long controller CPU times)

        if (_compensate_cpu_time)
        {
            comp_dt = (_computation_delay < 0) ? cpu_time.toSec() : _computation_delay;
        }

        if (_compensate_cpu_time || _compensate_dead_time)
        {
            _time_value_buffer.getValues(t_measure_x.toSec() - deadtime, comp_dt + deadtime, useq_predict);
            _compensator.predict(x, useq_predict, comp_dt + deadtime, x);

            // TODO(roesmann) seems not to work in gui
            // if (signal_target) signal_target->sendMeasurement("comp_states", t_measure_x.toSec() + comp_dt, x);
        }

        // controller
        if (_use_wall_time) t = (Time::now() - t_wall_start).toTime();
        if (!controller->step(x, *_xreference, *_ureference, last_dt, t, u_sequence, x_sequence, signal_target, nullptr, nullptr, nullptr, ns))
        {
            PRINT_ERROR("ClosedLoopControl::performTask(): controller error.");
            u_sequence->clear();
            u_sequence->add(t.toSec(), Eigen::VectorXd::Zero(controller->getControlInputDimension()));
        }

        // if dt<=0 -> inherit from controller (asynchronous control mode)
        if (_dt <= 0)
        {
            if (controller->getControlDuration() >= _min_dt)
            {
                if (controller->getControlDuration() <= _max_dt)
                {
                    dt.fromSec(controller->getControlDuration());
                    rate.updateCycleTime(dt);
                }
                else
                {
                    PRINT_WARNING_ONCE(
                        "ClosedLoopControl::performTask(): asychnronous control mode: controller returned dt>_max_dt. Setting dt=dt_max. This "
                        "warning is printed once.");
                    dt.fromSec(_max_dt);
                }
            }
            else
            {
                PRINT_WARNING_ONCE(
                    "ClosedLoopControl::performTask(): asychnronous control mode: controller returned dt<_min_dt. Setting dt=dt_min. This warning is "
                    "printed once.");
                dt.fromSec(_min_dt);
            }
        }

        // control plant
        if (_use_wall_time) t = (Time::now() - t_wall_start).toTime();
        plant->control(u_sequence, x_sequence, dt, t, signal_target, ns);

        if (_compensate_cpu_time || _compensate_dead_time)
        {
            _time_value_buffer.appendValues(t.toSec(), u_sequence->getValuesMap(0));
        }

        // save cpu time required for the controller step (used for compensation in the next step)
        cpu_time       = Time::now() - compensator_meas_start;
        cpu_time_cache = cpu_time;

        if (_use_wall_time && _compensate_cpu_time && _computation_delay < 0)
        {
            if (_computation_delay_filter) cpu_time.fromSec(_computation_delay_filter->filter(t.toSec(), cpu_time.toSec()));
        }

        if (signal_target)
        {
            if (_use_wall_time) t.fromSec(t_measure_x.toSec() + deadtime + cpu_time_cache.toSec());  // TODO(kraemer) correct?
            controller->sendSignals(t.toSec(), *signal_target, ns);
            if (u_sequence && !u_sequence->isEmpty()) signal_target->sendMeasurement(ns + "control_input", t.toSec(), u_sequence->getValuesMap(0));
        }

        if (!rate.sleep() && _realtime_sync)
        {
            PRINT_WARNING("ClosedLoopControlTask(): rate exceeded (" << rate.lastCycleTime().toSec() << "s/" << dt << "s).");
        }
        if (!_use_wall_time) t += dt;
    }

    plant->stop();
}

bool ClosedLoopControlTask::verify(const Environment& environment, std::string* msg) const
{
    bool ret_val = true;

    if (msg) msg->clear();

    // check if all objects are set
    if (!_xreference)
    {
        ret_val = false;
        if (msg) *msg += "State reference trajectory not specified for ClosedLoopControlTask\n";
    }

    // check if all objects are set
    if (!_ureference)
    {
        ret_val = false;
        if (msg) *msg += "Control reference trajectory not specified for ClosedLoopControlTask\n";
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

void ClosedLoopControlTask::reset()
{
    _time_value_buffer.reset();
    if (_computation_delay_filter) _computation_delay_filter->reset();
}

#ifdef MESSAGE_SUPPORT
void ClosedLoopControlTask::toMessage(corbo::messages::ClosedLoopControlTask& message) const
{
    message.set_sim_time(_sim_time);
    message.set_dt(_dt);
    message.set_realtime_sync(_realtime_sync);
    message.set_use_wall_time(_use_wall_time);
    message.set_min_dt(_min_dt);
    message.set_max_dt(_max_dt);
    if (_xreference) _xreference->toMessage(*message.mutable_xreference());
    if (_ureference) _ureference->toMessage(*message.mutable_ureference());

    message.set_compensate_dead_time(_compensate_dead_time);
    message.set_compensate_cpu_time(_compensate_cpu_time);
    message.set_computation_delay(_computation_delay);

    _compensator.toMessage(*message.mutable_compensator());
    if (_computation_delay_filter) _computation_delay_filter->toMessage(*message.mutable_computation_delay_filter());
}
void ClosedLoopControlTask::fromMessage(const corbo::messages::ClosedLoopControlTask& message, std::stringstream* issues)
{
    _sim_time      = message.sim_time();
    _dt            = message.dt();
    _realtime_sync = message.realtime_sync();
    _use_wall_time = message.use_wall_time();
    _min_dt        = message.min_dt();
    _max_dt        = message.max_dt();

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

    // compensator
    _computation_delay    = message.computation_delay();
    _compensate_dead_time = message.compensate_dead_time();
    _compensate_cpu_time  = message.compensate_cpu_time();

    _compensator.fromMessage(message.compensator(), issues);

    _computation_delay_filter.reset();

    if (message.has_computation_delay_filter() && !message.computation_delay_filter().has_no_filter())
    {
        // construct object
        std::string type;
        if (util::get_oneof_field_type(message.computation_delay_filter(), "filter", type, false))
        {
            FilterInterface::Ptr filter = FilterFactory::instance().create(type);
            // import parameters
            if (filter)
            {
                filter->fromMessage(message.computation_delay_filter(), issues);
                _computation_delay_filter = filter;
            }
        }
    }
}
#endif

}  // namespace corbo
