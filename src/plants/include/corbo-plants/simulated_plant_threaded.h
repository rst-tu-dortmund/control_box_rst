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

#ifndef SRC_PLANTS_INCLUDE_CORBO_PLANTS_SIMULATED_PLANT_THREADED_H_
#define SRC_PLANTS_INCLUDE_CORBO_PLANTS_SIMULATED_PLANT_THREADED_H_

#include <corbo-plants/simulated_plant.h>

#include <memory>
#include <mutex>
#include <thread>

namespace corbo {

/**
 * @brief Adapter class for plant simulation (threaded version)
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
 * @see SimulatedPlant PlantInterface
 *
 * @author Christoph Rösmann (christoph.roesmann@tu-dortmund.de)
 */
class SimulatedPlantThreaded : public SimulatedPlant
{
 public:
    using StateVector = Eigen::VectorXd;

    //! Default constructor
    SimulatedPlantThreaded();
    //! Construct simulated plant with system dynamics and output function
    SimulatedPlantThreaded(SystemDynamicsInterface::Ptr dynamics, SystemOutputInterface::Ptr output);

    virtual ~SimulatedPlantThreaded();

    // implements interface method
    Ptr getInstance() const override { return std::make_shared<SimulatedPlantThreaded>(); }

    // implements interface method
    bool initialize() override;

    // implements interface method
    void stop() override;

    // implements interface method
    void reset() override;

    // implements interface method
    bool control(const TimeSeries::ConstPtr& u_sequence, const TimeSeries::ConstPtr& x_sequence, const Duration& dt, const Time& t,
                 SignalTargetInterface* signal_target = nullptr, const std::string& ns = "") override;

    // implements interface method
    bool output(OutputVector& output, const Time& t, SignalTargetInterface* signal_target = nullptr, const std::string& ns = "") override;

    // implements interface method
    bool setState(const Eigen::Ref<const Eigen::VectorXd>& state) override;

    //! Specify rate for the simulation thread
    void setSimulationRate(double sim_rate) { _sim_rate = sim_rate; }

#ifdef MESSAGE_SUPPORT
    //! Export to message
    void toMessage(messages::SimulatedPlantThreaded& message) const;
    //! Import from message
    void fromMessage(const messages::SimulatedPlantThreaded& message, std::stringstream* issues = nullptr);

    // implements interface method
    void toMessage(messages::Plant& message) const override { toMessage(*message.mutable_simulated_plant_threaded()); }
    // implements interface method
    void fromMessage(const messages::Plant& message, std::stringstream* issues = nullptr) override
    {
        fromMessage(message.simulated_plant_threaded(), issues);
    }
#endif

 protected:
    void simulate();

 private:
    double _sim_rate = 10;

    std::thread _worker_thread;
    bool _stop_thread = false;

    std::mutex _control_mutex;
    ControlVector _control;

    std::mutex _current_state_mutex;
};

FACTORY_REGISTER_PLANT(SimulatedPlantThreaded)

}  // namespace corbo

#endif  // SRC_PLANTS_INCLUDE_CORBO_PLANTS_SIMULATED_PLANT_THREADED_H_
