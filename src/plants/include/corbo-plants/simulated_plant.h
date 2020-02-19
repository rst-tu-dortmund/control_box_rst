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

#ifndef SRC_PLANTS_INCLUDE_CORBO_PLANTS_SIMULATED_PLANT_H_
#define SRC_PLANTS_INCLUDE_CORBO_PLANTS_SIMULATED_PLANT_H_

#include <corbo-plants/plant_interface.h>

#include <corbo-numerics/integrator_interface.h>
#include <corbo-plants/disturbance_interface.h>
#include <corbo-systems/output_function_interface.h>
#include <corbo-systems/system_dynamics_interface.h>
#include <corbo-systems/time_value_buffer.h>

#include <memory>

namespace corbo {

/**
 * @brief Adapter class for plant simulation
 *
 * @ingroup plants
 *
 * This class provides a PlantInterface implementation for the
 * simulation of plants defined in terms of a state space model,
 * in particular a SystemDynamicsInterface and a SystemOutputInterface.
 *
 * In the continuous-time case, a numerical integrator solves
 * the initial value problem w.r.t. given control inputs.
 *
 * @see SimulatedPlantThreaded PlantInterface
 *
 * @author Christoph Rösmann (christoph.roesmann@tu-dortmund.de)
 */
class SimulatedPlant : public PlantInterface
{
 public:
    using StateVector = Eigen::VectorXd;

    //! Default constructor
    SimulatedPlant();
    //! Construct simulated plant with system dynamics and output function
    SimulatedPlant(SystemDynamicsInterface::Ptr dynamics, SystemOutputInterface::Ptr output);

    // implements interface method
    Ptr getInstance() const override { return std::make_shared<SimulatedPlant>(); }
    // implements interface method
    int getInputDimension() const override { return _dynamics ? _dynamics->getInputDimension() : 0; }
    // implements interface method
    int getOutputDimension() const override;
    // implements interface method
    bool requiresFutureControls() const override { return false; }
    // implements interface method
    bool requiresFutureStates() const override { return false; }

    // implements interface method
    // bool control(const ControlVector& u, const Duration& dt, const Time& t, SignalTargetInterface* signal_target = nullptr) override;

    // implements interface method
    bool control(const TimeSeries::ConstPtr& u_sequence, const TimeSeries::ConstPtr& x_sequence, const Duration& dt, const Time& t,
                 SignalTargetInterface* signal_target = nullptr, const std::string& ns = "") override;

    // implements interface method
    bool output(OutputVector& output, const Time& t, SignalTargetInterface* signal_target = nullptr, const std::string& ns = "") override;

    // implements interface method
    void reset() override;

    // implements interface method
    bool setState(const Eigen::Ref<const Eigen::VectorXd>& state) override { return setInitialState(state); }

    //! Set the system dynamics of the simulated plant
    void setSystemDynamics(SystemDynamicsInterface::Ptr dynamics);
    //! Set the output function of the simulated plant
    void setOutputFunction(SystemOutputInterface::Ptr output);
    //! Set a numerical integrator for continuous-time dynamics
    void setIntegrator(NumericalIntegratorExplicitInterface::Ptr integrator);
    //! Set plant input disturbance model
    void setInputDisturbance(DisturbanceInterface::Ptr disturbance);
    //! Set plant output disturbance model
    void setOutputDisturbance(DisturbanceInterface::Ptr disturbance);
    //! Set plant state disturbance model
    void setStateDisturbance(DisturbanceInterface::Ptr disturbance);

    //! Specify an initial state x0 [SystemDynamicsInterface::getStateDimension() x 1]
    bool setInitialState(const StateVector& x_init);

    // implements interface method
    void getAvailableSignals(SignalTargetInterface& signal_target, const std::string& ns = "") const override;

#ifdef MESSAGE_SUPPORT
    //! Export to message
    void toMessage(messages::SimulatedPlant& message) const;
    //! Import from message
    void fromMessage(const messages::SimulatedPlant& message, std::stringstream* issues = nullptr);

    // implements interface method
    void toMessage(messages::Plant& message) const override { toMessage(*message.mutable_simulated_plant()); }
    // implements interface method
    void fromMessage(const messages::Plant& message, std::stringstream* issues = nullptr) override { fromMessage(message.simulated_plant(), issues); }
#endif

 protected:
    SystemDynamicsInterface::Ptr _dynamics;
    SystemOutputInterface::Ptr _output;

    NumericalIntegratorExplicitInterface::Ptr _integrator;

    DisturbanceInterface::Ptr _input_disturbance;
    DisturbanceInterface::Ptr _output_disturbance;
    DisturbanceInterface::Ptr _state_disturbance;

    TimeValueBuffer _time_value_buffer;

    StateVector _initial_state;  // just a cache to allow resetting
    StateVector _current_state;
    ControlVector _current_control;
};

FACTORY_REGISTER_PLANT(SimulatedPlant)

}  // namespace corbo

#endif  // SRC_PLANTS_INCLUDE_CORBO_PLANTS_SIMULATED_PLANT_H_
