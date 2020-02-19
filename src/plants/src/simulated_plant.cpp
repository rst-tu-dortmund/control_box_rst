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
#include <corbo-numerics/explicit_integrators.h>
#include <corbo-plants/simulated_plant.h>
#include <memory>

namespace corbo {

SimulatedPlant::SimulatedPlant()
{
    // default integrator
    setIntegrator(std::make_shared<IntegratorExplicitEuler>());
}

SimulatedPlant::SimulatedPlant(SystemDynamicsInterface::Ptr dynamics, SystemOutputInterface::Ptr output)
{
    setSystemDynamics(dynamics);  // this also sets the initial state
    setOutputFunction(output);

    // default integrator
    setIntegrator(std::make_shared<IntegratorExplicitEuler>());
}

int SimulatedPlant::getOutputDimension() const
{
    if (!_output) return 0;
    if (_output->getOutputDimension() == property::INHERITED)
    {
        if (_dynamics) return _dynamics->getStateDimension();
        return 0;
    }

    return _output->getOutputDimension();
}

/*
bool SimulatedPlant::control(const ControlVector& u, const Duration& dt, const Time& t, SignalTargetInterface* signal_target)
{
    if (!_dynamics || !_output || !_integrator) return false;

    StateVector next_state(_current_state.rows());

    // disturb input if desired
    ControlVector u_disturbed(u.rows());
    if (_input_disturbance)
        _input_disturbance->disturb(t, u, u_disturbed);
    else
        u_disturbed = u;

    if (_dynamics->isContinuousTime())
    {
        // _integrator->integrate(_current_state, u, dt.toSec(), *_dynamics, next_state);
        _integrator->solveIVPWithInnerDt(_current_state, u_disturbed, dt.toSec(), *_dynamics, next_state);
    }
    else
    {
        _dynamics->dynamics(_current_state, u_disturbed, next_state);
    }
    _current_state = next_state;

    if (signal_target) signal_target->sendMeasurement("plant/state", t.toSec(), _current_state);

    return true;
}
*/

bool SimulatedPlant::control(const TimeSeries::ConstPtr& u_sequence, const TimeSeries::ConstPtr& x_sequence, const Duration& dt, const Time& t,
                             SignalTargetInterface* signal_target, const std::string& ns)
{
    if (!_dynamics || !_output || !_integrator) return false;

    StateVector next_state(_current_state.rows());
    ControlVector u(getInputDimension());

    if (u_sequence->getTimeDimension() < 1)
    {
        PRINT_ERROR_NAMED("u_sequence is empty.");
        return false;
    }

    u = u_sequence->getValuesMap(0);  // Get most recent control vector

    // disturb input if desired
    ControlVector u_disturbed(u.rows());
    if (_input_disturbance)
        _input_disturbance->disturb(t, u, u_disturbed);
    else
        u_disturbed = u;

    // we might have deadtimes so use buffer
    std::vector<std::pair<double, ControlVector>> delayed_u_seq;

    if (_time_value_buffer.isEmpty()) _time_value_buffer.setInitialValue(ControlVector::Zero(getInputDimension()));

    _time_value_buffer.appendValues(t.toSec(), u_disturbed);
    _time_value_buffer.getValues(t.toSec() - _dynamics->getDeadTime(), dt.toSec(), delayed_u_seq);

    double cur_t = t.toSec();
    for (int i = 0; i < delayed_u_seq.size(); ++i)
    {
        double cur_dt = delayed_u_seq[i].first;

        if (_dynamics->isContinuousTime())
        {
            // _integrator->integrate(_current_state, u, dt.toSec(), *_dynamics, next_state);
            _integrator->solveIVP(_current_state, delayed_u_seq[i].second, cur_dt, *_dynamics, next_state);
        }
        else
        {
            // TODO(roesmann): we need to check how discrete states behave with the deadtime simulator!!!
            PRINT_WARNING_COND(_dynamics->getDeadTime() > 0, "Discrete-time systems with deadtime not yet tested/fully implemented...");
            _dynamics->dynamics(_current_state, delayed_u_seq[i].second, next_state);
        }
        _current_state = next_state;

        // disturb state if desired
        if (_state_disturbance) _state_disturbance->disturb(t, _current_state, _current_state);

        cur_t += cur_dt;

        if (signal_target) signal_target->sendMeasurement(ns + "plant/state", cur_t, _current_state);  // TODO(roesmann): check time stamp!
    }
    return true;
}

bool SimulatedPlant::output(OutputVector& output, const Time& t, SignalTargetInterface* /*signal_target*/, const std::string& /*ns*/)
{
    if (!_output) return false;

    _output->output(_current_state, output);

    // disturb output if desired
    if (_output_disturbance) _output_disturbance->disturb(t, output, output);

    return true;
}

void SimulatedPlant::setSystemDynamics(SystemDynamicsInterface::Ptr dynamics)
{
    _dynamics = dynamics;
    // Set initial state and control (important: we set the dimensions here as well!)
    setInitialState(StateVector::Zero(dynamics->getStateDimension()));
    _current_control = ControlVector::Zero(dynamics->getInputDimension());
}

void SimulatedPlant::setOutputFunction(SystemOutputInterface::Ptr output) { _output = output; }

void SimulatedPlant::setIntegrator(NumericalIntegratorExplicitInterface::Ptr integrator) { _integrator = integrator; }

void SimulatedPlant::setInputDisturbance(DisturbanceInterface::Ptr disturbance) { _input_disturbance = disturbance; }

void SimulatedPlant::setOutputDisturbance(DisturbanceInterface::Ptr disturbance) { _output_disturbance = disturbance; }

void SimulatedPlant::setStateDisturbance(DisturbanceInterface::Ptr disturbance) { _state_disturbance = disturbance; }

bool SimulatedPlant::setInitialState(const StateVector& x_init)
{
    if (x_init.rows() != _dynamics->getStateDimension())
    {
        // PRINT_ERROR("SimulatedPlant::setInitialState(): dimension mismatch between x_init and _dynamics->getStateDimension()");
        return false;
    }
    _current_state = x_init;
    _initial_state = _current_state;
    return true;
}

void SimulatedPlant::getAvailableSignals(SignalTargetInterface& signal_target, const std::string& ns) const
{
    if (_dynamics) signal_target.registerMeasurement(ns + "plant/state", _dynamics->getStateDimension());
}

void SimulatedPlant::reset()
{
    _current_state = _initial_state;
    _time_value_buffer.reset();

    if (_dynamics) _dynamics->reset();
    if (_output) _output->reset();
    if (_input_disturbance) _input_disturbance->reset();
    if (_state_disturbance) _state_disturbance->reset();
    if (_output_disturbance) _output_disturbance->reset();

    _current_control.setZero();
}

#ifdef MESSAGE_SUPPORT
void SimulatedPlant::toMessage(messages::SimulatedPlant& message) const
{
    if (_dynamics) _dynamics->toMessage(*message.mutable_system_dynamics());
    if (_output) _output->toMessage(*message.mutable_output_function());
    if (_integrator) _integrator->toMessage(*message.mutable_integrator());

    if (_input_disturbance) _input_disturbance->toMessage(*message.mutable_input_disturbance());
    if (_output_disturbance) _output_disturbance->toMessage(*message.mutable_output_disturbance());
    if (_state_disturbance) _state_disturbance->toMessage(*message.mutable_state_disturbance());

    message.mutable_x0()->Resize(_current_state.rows(), 0);
    Eigen::Map<StateVector>(message.mutable_x0()->mutable_data(), _current_state.rows()) = _current_state;
}

void SimulatedPlant::fromMessage(const corbo::messages::SimulatedPlant& message, std::stringstream* issues)
{
    _dynamics.reset();
    _current_state.resize(0);
    _current_control.resize(0);
    _output.reset();
    _integrator.reset();

    // system dynamics
    if (message.has_system_dynamics())
    {
        // construct object
        std::string type;
        util::get_oneof_field_type_expand_isolated(message.system_dynamics(), "system_dynamics", type, false, 1);
        SystemDynamicsInterface::Ptr dynamics = SystemDynamicsFactory::instance().create(type);
        // import parameters
        if (dynamics)
        {
            StateVector x0;
            dynamics->fromMessage(message.system_dynamics(), issues);
            setSystemDynamics(dynamics);
        }
        else if (issues)
        {
            *issues << "SimulatedPlant: No system dynamics specified.\n";
            return;
        }
    }

    // output function
    if (message.has_output_function())
    {
        // construct object
        std::string type;
        util::get_oneof_field_type_expand_isolated(message.output_function(), "output_function", type, false, 1);
        SystemOutputInterface::Ptr output = OutputFunctionFactory::instance().create(type);
        // import parameters
        if (output)
        {
            output->fromMessage(message.output_function(), issues);
            setOutputFunction(output);
        }
    }

    // integrator
    if (message.has_integrator())
    {
        // construct object
        std::string type;
        util::get_oneof_field_type(message.integrator(), "explicit_integrator", type, false);
        NumericalIntegratorExplicitInterface::Ptr integrator = create_from_factory<NumericalIntegratorExplicitInterface>(type);
        // import parameters
        if (integrator)
        {
            integrator->fromMessage(message.integrator(), issues);
            setIntegrator(integrator);
        }
    }

    // start state
    if (message.x0_size() > 0 && _dynamics)
    {
        int dim = message.x0_size();
        Eigen::VectorXd x0(dim);
        for (int i = 0; i < dim; ++i) x0[i] = message.x0(i);
        if (!setInitialState(x0))
        {
            *issues << "Size of state x0 (" << dim << ") does not comply with system dynamics dimension (" << _dynamics->getStateDimension()
                    << ").\n";
        }
    }
    else if (issues)
        *issues << "SimulatedPlant: dimension of x0 must be larger then zero.\n";

    // input disturbance
    if (message.has_input_disturbance() && !message.input_disturbance().has_no_disturbance())
    {
        // construct object
        std::string type;
        util::get_oneof_field_type(message.input_disturbance(), "disturbance", type, false);
        DisturbanceInterface::Ptr disturbance = create_from_factory<DisturbanceInterface>(type);
        // import parameters
        if (disturbance)
        {
            disturbance->fromMessage(message.input_disturbance(), issues);
            if (!disturbance->checkParameters(getInputDimension(), issues)) return;
            setInputDisturbance(disturbance);
        }
    }

    // output disturbance
    if (message.has_output_disturbance() && !message.output_disturbance().has_no_disturbance())
    {
        // construct object
        std::string type;
        util::get_oneof_field_type(message.output_disturbance(), "disturbance", type, false);
        DisturbanceInterface::Ptr disturbance = create_from_factory<DisturbanceInterface>(type);
        // import parameters
        if (disturbance)
        {
            disturbance->fromMessage(message.output_disturbance(), issues);
            if (!disturbance->checkParameters(getOutputDimension(), issues)) return;
            setOutputDisturbance(disturbance);
        }
    }

    // state disturbance
    if (message.has_state_disturbance() && !message.state_disturbance().has_no_disturbance())
    {
        // construct object
        std::string type;
        util::get_oneof_field_type(message.state_disturbance(), "disturbance", type, false);
        DisturbanceInterface::Ptr disturbance = create_from_factory<DisturbanceInterface>(type);
        // import parameters
        if (disturbance)
        {
            disturbance->fromMessage(message.state_disturbance(), issues);
            if (!disturbance->checkParameters(_dynamics->getStateDimension(), issues)) return;
            setStateDisturbance(disturbance);
        }
        else if (issues)
        {
            *issues << "SimulatedPlant: Cannot set state disturbance model.\n";
        }
    }
}
#endif

}  // namespace corbo
