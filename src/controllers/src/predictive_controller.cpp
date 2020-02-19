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

#include <corbo-controllers/predictive_controller.h>

#include <corbo-communication/utilities.h>
#include <corbo-core/console.h>

namespace corbo {

PredictiveController::PredictiveController() {}

bool PredictiveController::initialize(const StateVector& x, ReferenceTrajectoryInterface& expected_xref, ReferenceTrajectoryInterface& expected_uref,
                                      const Duration& expected_dt, const Time& t, ReferenceTrajectoryInterface* sref)
{
    if (!_ocp || !_ocp->initialize())
    {
        PRINT_ERROR("PredictiveController::initialize(): no ocp specified or ocp initialization failed.");
        return false;
    }
    _initialized = true;
    return true;
}

bool PredictiveController::step(const ControllerInterface::StateVector& x, ReferenceTrajectoryInterface& xref, ReferenceTrajectoryInterface& uref,
                                const Duration& dt, const Time& t, TimeSeries::Ptr u_sequence, TimeSeries::Ptr x_sequence,
                                SignalTargetInterface* signal_target, ReferenceTrajectoryInterface* sref, ReferenceTrajectoryInterface* xinit,
                                ReferenceTrajectoryInterface* uinit, const std::string& ns)
{
    if (!_initialized)
    {
        if (!initialize(x, xref, uref, dt, t, sref)) return false;
    }

    ControlVector u;

    if (!_ocp) return false;

    if (_auto_update_prev_control) _ocp->setPreviousControlInputDt(dt.toSec());

    bool success = false;

    Time t_pre = Time::now();

    for (int i = 0; i < _num_ocp_iterations; ++i) success = _ocp->compute(x, xref, uref, sref, t, i == 0, signal_target, xinit, uinit, ns);

    _statistics.step_time = Time::now() - t_pre;

    success = success && _ocp->getFirstControlInput(u);

    if (_auto_update_prev_control) _ocp->setPreviousControlInput(u, dt.toSec());  // cache control input for next step call

    _ocp->getTimeSeries(x_sequence, u_sequence);
    _x_ts = x_sequence;
    _u_ts = u_sequence;

    return success;
}

void PredictiveController::getAvailableSignals(SignalTargetInterface& signal_target, const std::string& ns) const
{
    if (_publish_prediction && _ocp)
    {
        signal_target.registerTimeSeries(ns + "prediction/x", _ocp->getStateDimension(), false);
        signal_target.registerTimeSeries(ns + "prediction/u", _ocp->getControlInputDimension(), _ocp->isConstantControlAction());
        signal_target.registerMeasurement(ns + "prediction/n", 1, {}, false);
        signal_target.registerMeasurement(ns + "prediction/first_dt", 1, {}, false);
        signal_target.registerMeasurement(ns + "prediction/objective", 1, {}, false);
        signal_target.registerMeasurement(ns + "prediction/cpu_time", 1, {}, true);
    }
}

void PredictiveController::reset()
{
    if (_ocp) _ocp->reset();
}

void PredictiveController::sendSignals(double t, SignalTargetInterface& signal_target, const std::string& ns) const
{
    if (_publish_prediction)
    {
        if (_x_ts && _x_ts->getTimeDimension() > 0)
        {
            _x_ts->setTimeFromStart(t);
            signal_target.sendTimeSeries(ns + "prediction/x", _x_ts);
        }

        if (_u_ts && _u_ts->getTimeDimension() > 0)
        {
            _u_ts->setTimeFromStart(t);
            signal_target.sendTimeSeries(ns + "prediction/u", _u_ts);
        }

        signal_target.sendMeasurement(ns + "prediction/n", t, {(double)_ocp->getN()});
        signal_target.sendMeasurement(ns + "prediction/first_dt", t, {_ocp->getFirstDt()});
        signal_target.sendMeasurement(ns + "prediction/objective", t, {_ocp->getCurrentObjectiveValue()});
        signal_target.sendMeasurement(ns + "prediction/cpu_time", t, {_statistics.step_time.toSec()});
    }
}

#ifdef MESSAGE_SUPPORT
void PredictiveController::toMessage(corbo::messages::PredictiveController& message) const
{
    if (_ocp) _ocp->toMessage(*message.mutable_optimal_control_problem());
    message.set_num_ocp_iterations(_num_ocp_iterations);
    message.set_auto_update_prev_control(_auto_update_prev_control);
    message.set_publish_prediction(_publish_prediction);
}
void PredictiveController::fromMessage(const corbo::messages::PredictiveController& message, std::stringstream* issues)
{
    // ocp
    if (!message.has_optimal_control_problem())
    {
        if (issues) *issues << "PredictiveController: no optimal control problem specified.\n";
        return;
    }

    _num_ocp_iterations       = message.num_ocp_iterations();
    _auto_update_prev_control = message.auto_update_prev_control();
    _publish_prediction       = message.publish_prediction();

    // construct object
    std::string type;
    util::get_oneof_field_type(message.optimal_control_problem(), "optimal_control_problem", type, false);
    OptimalControlProblemInterface::Ptr ocp = OptimalControlProgramFactory::instance().create(type);
    // import parameters
    if (ocp)
    {
        ocp->fromMessage(message.optimal_control_problem(), issues);
        setOptimalControlProblem(ocp);
    }
    else
    {
        if (issues) *issues << "PredictiveController: unknown optimal control problem specified.\n";
        return;
    }
}
#endif

}  // namespace corbo
