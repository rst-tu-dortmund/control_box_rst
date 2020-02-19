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

#include <corbo-plants/simulated_plant_threaded.h>

#include <corbo-core/global.h>
#include <corbo-core/time.h>

#include <memory>
#include <thread>

namespace corbo {

SimulatedPlantThreaded::SimulatedPlantThreaded() : SimulatedPlant() {}

SimulatedPlantThreaded::SimulatedPlantThreaded(SystemDynamicsInterface::Ptr dynamics, SystemOutputInterface::Ptr output)
    : SimulatedPlant(dynamics, output)
{
}

SimulatedPlantThreaded::~SimulatedPlantThreaded()
{
    _stop_thread = true;
    if (_worker_thread.joinable())
    {
        _worker_thread.join();
    }
}

bool SimulatedPlantThreaded::initialize()
{
    if (!_dynamics || !_output || !_integrator) return false;

    // initialize control vector with zeros
    _control.setZero(_dynamics->getInputDimension());

    _stop_thread = false;

    _worker_thread = std::thread(&SimulatedPlantThreaded::simulate, this);
    return true;
}

void SimulatedPlantThreaded::stop()
{
    _stop_thread = true;
    _worker_thread.join();
}

void SimulatedPlantThreaded::reset()
{
    _stop_thread = true;
    _worker_thread.join();

    SimulatedPlant::reset();
}

void SimulatedPlantThreaded::simulate()
{
    StateVector x(_current_state.rows());
    StateVector next_state(x.rows());
    ControlVector u;

    Duration dt(1.0 / _sim_rate);
    Time t(0);

    std::vector<std::pair<double, ControlVector>> delayed_u_seq;

    if (_time_value_buffer.isEmpty()) _time_value_buffer.setInitialValue(ControlVector::Zero(getInputDimension()));

    Rate sim_rate(_sim_rate);
    while (!_stop_thread && ok())
    {
        {
            std::lock_guard<std::mutex> locker(_control_mutex);
            u = _control;  // Get most recent control vector
        }
        {
            std::lock_guard<std::mutex> locker(_current_state_mutex);
            x = _current_state;  // Get current state vector
        }

        // disturb input if desired
        ControlVector u_disturbed(u.rows());
        if (_input_disturbance)
            _input_disturbance->disturb(t, u, u_disturbed);
        else
            u_disturbed = u;

        // we might have deadtimes
        delayed_u_seq.clear();

        _time_value_buffer.appendValues(t.toSec(), u_disturbed);
        _time_value_buffer.getValues(t.toSec() - _dynamics->getDeadTime(), dt.toSec(), delayed_u_seq);

        for (int i = 0; i < delayed_u_seq.size(); ++i)
        {
            double cur_dt = delayed_u_seq[i].first;

            if (_dynamics->isContinuousTime())
            {
                // _integrator->integrate(_current_state, u, dt.toSec(), *_dynamics, next_state);
                _integrator->solveIVP(x, delayed_u_seq[i].second, cur_dt, *_dynamics, next_state);
            }
            else
            {
                // TODO(roesmann): we need to check how discrete states behave with the deadtime simulator!!!
                PRINT_WARNING_COND(_dynamics->getDeadTime() > 0, "Discrete-time systems with deadtime not yet tested/fully implemented...");
                _dynamics->dynamics(x, delayed_u_seq[i].second, next_state);
            }
            x = next_state;  // TODO(roesmann): we can avoid this if the integrator would not be affacted by alias
            t += Duration(cur_dt);
        }

        // write values back to shared cache
        {
            std::lock_guard<std::mutex> locker(_current_state_mutex);
            _current_state = x;  // Get current state vector
        }

        // TODO(roesmann): should we use the lastCycleTime() to update our simulation dt?
        // This might introduce new effects, so let's define at least the following warning:
        if (!sim_rate.sleep())
        {
            PRINT_WARNING_NAMED("Rate exceeded (" << sim_rate.lastCycleTime().toSec() * 1000.0 << "ms/" << dt.toSec() * 1000.0
                                                  << "ms). Simulation results are probably useless. You might reduce sim_rate.");
        }
    }
}

bool SimulatedPlantThreaded::control(const TimeSeries::ConstPtr& u_sequence, const TimeSeries::ConstPtr& x_sequence, const Duration& dt,
                                     const Time& t, SignalTargetInterface* /*signal_target*/, const std::string& /*ns*/)
{
    if (u_sequence->getTimeDimension() < 1)
    {
        PRINT_ERROR_NAMED("u_sequence is empty.");
        return false;
    }

    {
        std::lock_guard<std::mutex> locker(_control_mutex);
        _control = u_sequence->getValuesMap(0);  // Get most recent control vector
    }

    return true;
}

bool SimulatedPlantThreaded::output(OutputVector& output, const Time& t, SignalTargetInterface* signal_target, const std::string& ns)
{
    if (!_output) return false;

    {
        std::lock_guard<std::mutex> locker(_current_state_mutex);
        _output->output(_current_state, output);

        if (signal_target) signal_target->sendMeasurement(ns + "plant/state", t.toSec(), _current_state);
    }

    // disturb output if desired
    if (_output_disturbance) _output_disturbance->disturb(t, output, output);

    return true;
}

bool SimulatedPlantThreaded::setState(const Eigen::Ref<const Eigen::VectorXd>& state)
{
    std::lock_guard<std::mutex> l(_current_state_mutex);
    bool success = setInitialState(state);
    return success;
}

#ifdef MESSAGE_SUPPORT
void SimulatedPlantThreaded::toMessage(corbo::messages::SimulatedPlantThreaded& message) const
{
    SimulatedPlant::toMessage(*message.mutable_simulated_plant());

    message.set_sim_rate(_sim_rate);
}

void SimulatedPlantThreaded::fromMessage(const corbo::messages::SimulatedPlantThreaded& message, std::stringstream* issues)
{
    SimulatedPlant::fromMessage(message.simulated_plant(), issues);

    _sim_rate = message.sim_rate();
    if (_sim_rate <= 0)
    {
        if (issues) *issues << "SimulatedPlantThreaded:: sim_rate must be positive" << std::endl;
        _sim_rate = 10;
    }
}
#endif

}  // namespace corbo
