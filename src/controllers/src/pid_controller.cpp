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

#include <corbo-controllers/pid_controller.h>
#include <corbo-core/console.h>

#include <algorithm>

namespace corbo {

PidController::PidController() { setNumParallelPid(_num_parallel_pid); }

bool PidController::step(const ControllerInterface::StateVector& x, ReferenceTrajectoryInterface& xref, ReferenceTrajectoryInterface& uref,
                         const Duration& dt, const Time& t, TimeSeries::Ptr u_sequence, TimeSeries::Ptr x_sequence,
                         SignalTargetInterface* signal_target, ReferenceTrajectoryInterface* sref, ReferenceTrajectoryInterface* xinit,
                         ReferenceTrajectoryInterface* uinit, const std::string& ns)
{
    assert(x.rows() == xref.getDimension() && "Dimension mismatch in controller: current state x and reference");

    ReferenceTrajectoryInterface::OutputVector yref;
    xref.getReference(t, yref);

    if (yref.rows() != _num_parallel_pid)
    {
        PRINT_WARNING("PidController number of parallel PID controllers does not match dimension of the state vector.");
        return false;
    }
    if (uref.getDimension() != _num_parallel_pid)
    {
        PRINT_WARNING("PidController: number of parallel PID controllers does not match number of required control inputs");
        return false;
    }

    assert(_num_parallel_pid == _p_error.size());
    assert(_num_parallel_pid == _i_error.size());
    assert(_num_parallel_pid == _d_error.size());

    ControlVector u;

    for (int i = 0; i < _num_parallel_pid; ++i)
    {
        double p_error_last = _p_error[i];

        _p_error[i] = yref[i] - x[i];
        _d_error[i] = dt.toSec() > 0 ? (_p_error[i] - p_error_last) / dt.toSec() : 0.0;
        _i_error[i] += dt.toSec() * _p_error[i];

        u.resize(_num_parallel_pid);
        u[i] = _p_gain * _p_error[i] + _i_gain * _i_error[i] + _d_gain * _d_error[i];

        if (signal_target && _publish_error)
        {
            signal_target->sendMeasurement(ns + "controller/error/p", t.toSec(), _p_error);
            signal_target->sendMeasurement(ns + "controller/error/i", t.toSec(), _i_error);
            signal_target->sendMeasurement(ns + "controller/error/d", t.toSec(), _d_error);
        }
    }

    u_sequence->clear();
    x_sequence->clear();
    u_sequence->add(0.0, u);
    x_sequence->add(0.0, x);

    return true;
}

void PidController::getAvailableSignals(SignalTargetInterface& signal_target, const std::string& ns) const
{
    if (_publish_error)
    {
        signal_target.registerMeasurement(ns + "controller/error/p", _num_parallel_pid);
        signal_target.registerMeasurement(ns + "controller/error/i", _num_parallel_pid);
        signal_target.registerMeasurement(ns + "controller/error/d", _num_parallel_pid);
    }
}

void PidController::reset()
{
    std::fill(_p_error.begin(), _p_error.end(), 0.0);
    std::fill(_d_error.begin(), _d_error.end(), 0.0);
    std::fill(_i_error.begin(), _i_error.end(), 0.0);
}

#ifdef MESSAGE_SUPPORT
void PidController::toMessage(corbo::messages::Controller& message) const
{
    message.mutable_pid_controller()->set_p_gain(_p_gain);
    message.mutable_pid_controller()->set_i_gain(_i_gain);
    message.mutable_pid_controller()->set_d_gain(_d_gain);

    message.mutable_pid_controller()->set_num_parallel_pid(_num_parallel_pid);

    message.mutable_pid_controller()->set_publish_error(_publish_error);
}

void PidController::fromMessage(const corbo::messages::Controller& message, std::stringstream* issues)
{
    _p_gain = message.pid_controller().p_gain();
    _i_gain = message.pid_controller().i_gain();
    _d_gain = message.pid_controller().d_gain();

    setNumParallelPid(message.pid_controller().num_parallel_pid());

    _publish_error = message.pid_controller().publish_error();
}
#endif

}  // namespace corbo
