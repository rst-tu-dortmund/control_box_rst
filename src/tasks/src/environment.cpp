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

#include <corbo-tasks/environment.h>

#include <corbo-core/console.h>

#include <string>

namespace corbo {

Environment::Environment(ControllerInterface::Ptr controller, ObserverInterface::Ptr observer, PlantInterface::Ptr plant)
{
    setController(controller);
    setObserver(observer);
    setPlant(plant);
}

void Environment::setController(ControllerInterface::Ptr controller) { _controller = controller; }

void Environment::setObserver(ObserverInterface::Ptr observer) { _observer = observer; }

void Environment::setPlant(PlantInterface::Ptr plant) { _plant = plant; }

bool Environment::verify(std::string* msg) const
{
    if (msg) msg->clear();
    bool ret_val = true;
    // first check pointer
    if (!hasController())
    {
        ret_val = false;
        if (msg) *msg += "Controller not specified.\n";
    }

    if (!hasObserver())
    {
        ret_val = false;
        if (msg) *msg += "Observer not specified.\n";
    }

    if (!hasPlant())
    {
        ret_val = false;
        if (msg) *msg += "Plant not specified.\n";
    }
    if (!ret_val) return false;

    // now check dimensions
    int controller_state_dim   = getController()->getStateDimension();
    int controller_control_dim = getController()->getControlInputDimension();

    int observer_output_dim = getObserver()->getOutputDimension();
    int observer_state_dim  = getObserver()->getStateDimension();

    int plant_input_dim  = getPlant()->getInputDimension();
    int plant_output_dim = getPlant()->getOutputDimension();

    // control input vs. plant input
    if (controller_control_dim != plant_input_dim)
    {
        ret_val = false;
        if (msg)
        {
            *msg += "Contol input dimension (" + std::to_string(controller_control_dim) + ") does not match plant intput dimension (" +
                    std::to_string(plant_input_dim) + ").\n";
        }
    }

    // plant output vs observer input (output_vector)
    if (plant_output_dim != observer_output_dim && observer_output_dim != property::INHERITED)
    {
        ret_val = false;
        if (msg)
        {
            *msg += "Plant output dimension (" + std::to_string(plant_output_dim) + ") does not match the input dimension of the observer (" +
                    std::to_string(observer_output_dim) + ").\n";
        }
    }

    // observer state vs controller state
    if (observer_state_dim != controller_state_dim && observer_state_dim != property::INHERITED)
    {
        ret_val = false;
        if (msg)
        {
            *msg += "Observer state dimension (" + std::to_string(observer_state_dim) + ") does not match the state dimension of the controller (" +
                    std::to_string(controller_state_dim) + ").\n";
        }
    }

    // check if inheritance is appropriate

    if (observer_output_dim == property::INHERITED && observer_state_dim != property::INHERITED)
    {
        if (observer_state_dim != controller_state_dim)
        {
            ret_val = false;
            if (msg)
            {
                *msg += "Observer input dimension is inherited from the plant output dimension (" + std::to_string(plant_output_dim) +
                        ") and does not match the state dimension of the observer (" + std::to_string(observer_state_dim) + ").\n";
            }
        }
    }

    if (observer_output_dim == property::INHERITED && observer_state_dim == property::INHERITED)
    {
        if (plant_output_dim != controller_state_dim)
        {
            ret_val = false;
            if (msg)
            {
                *msg += "Observer state dimension is inherited from the plant output dimension (" + std::to_string(plant_output_dim) +
                        ") and does not match the state dimension of the controller (" + std::to_string(controller_state_dim) + ").\n";
            }
        }
    }

    if (observer_output_dim != property::INHERITED && observer_state_dim == property::INHERITED)
    {
        if (observer_state_dim != controller_state_dim)
        {
            ret_val = false;
            if (msg)
            {
                *msg += "Observer state dimension is inherited from the observer input dimension (" + std::to_string(observer_output_dim) +
                        ") and does not match the state dimension of the controller (" + std::to_string(controller_state_dim) + ").\n";
            }
        }
    }
    return ret_val;
}

void Environment::reset()
{
    if (_controller) _controller->reset();
    if (_observer) _observer->reset();
    if (_plant) _plant->reset();
}
}  // namespace corbo
