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

#include <corbo-systems/one_step_predictor.h>

#include <corbo-communication/utilities.h>

namespace corbo {

double OneStepPredictor::getDeadTime() { return _dynamics ? _dynamics->getDeadTime() : 0.0; }

bool OneStepPredictor::initialize()
{
    if (!_dynamics || !_integrator) return false;

    _initialized = true;

    return true;
}

void OneStepPredictor::predict(const Eigen::Ref<const StateVector>& x0, std::vector<std::pair<double, ControlVector>> u_seq, double dt,
                               Eigen::Ref<StateVector> x1)
{
    if (dt <= 0)
    {
        x1 = x0;
        return;
    }

    if (!_initialized) initialize();

    if (x1.rows() != x0.rows())
    {
        PRINT_ERROR_NAMED("size of x1 does not match size of x0. Aborting.");
        return;
    }

    Eigen::VectorXd x = x0;

    for (int i = 0; i < u_seq.size(); ++i)
    {
        double cur_dt = u_seq[i].first;

        if (_dynamics->isContinuousTime())
        {
            // _integrator->integrate(_current_state, u, dt.toSec(), *_dynamics, next_state);
            _integrator->solveIVP(x, u_seq[i].second, cur_dt, *_dynamics, x1);
        }
        else
        {
            // TODO(roesmann): we need to check how discrete states behave with the deadtime simulator!!!
            PRINT_WARNING_COND(_dynamics->getDeadTime() > 0, "Discrete-time systems with deadtime not yet tested/fully implemented...");
            _dynamics->dynamics(x, u_seq[i].second, x1);
        }
        x = x1;  // TODO(roesmann): place for improving efficiency (this is not required in the last step...)
    }
}

void OneStepPredictor::setSystemDynamics(SystemDynamicsInterface::Ptr dynamics) { _dynamics = dynamics; }

#ifdef MESSAGE_SUPPORT
void OneStepPredictor::toMessage(messages::OneStepPredictor& message) const
{
    if (_dynamics) _dynamics->toMessage(*message.mutable_system_dynamics());
    if (_integrator) _integrator->toMessage(*message.mutable_integrator());
}

void OneStepPredictor::fromMessage(const messages::OneStepPredictor& message, std::stringstream* issues)
{
    _dynamics.reset();
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
            dynamics->fromMessage(message.system_dynamics(), issues);
            setSystemDynamics(dynamics);
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
}
#endif

}  // namespace corbo
