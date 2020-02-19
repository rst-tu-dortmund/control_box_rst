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

#ifndef SRC_PLANTS_INCLUDE_CORBO_PLANTS_PLANT_INTERFACE_H_
#define SRC_PLANTS_INCLUDE_CORBO_PLANTS_PLANT_INTERFACE_H_

#include <corbo-core/factory.h>
#include <corbo-core/signal_target_interface.h>
#include <corbo-core/types.h>

#ifdef MESSAGE_SUPPORT
#include <corbo-communication/messages/plants/plant.pb.h>
#endif

#include <memory>

namespace corbo {

/**
 * @brief Interface class for plants
 *
 * @ingroup plants
 *
 * This class specifies methods that are required to be implemented by specific
 * plants in order to allow their general utilization in a variety of control tasks.
 *
 * @remark This interface is provided with factory support (PlantFactoryPlantFactory).
 *
 * @see SimulatedPlant
 *
 * @author Christoph Rösmann (christoph.roesmann@tu-dortmund.de)
 */
class PlantInterface
{
 public:
    using Ptr           = std::shared_ptr<PlantInterface>;
    using ControlVector = Eigen::VectorXd;
    using StateVector   = Eigen::VectorXd;
    using OutputVector  = Eigen::VectorXd;

    //! Virtual destructor
    virtual ~PlantInterface() {}

    //! Return a newly created shared instance of the implemented class
    virtual Ptr getInstance() const = 0;

    //! Get access to the associated factory
    static Factory<PlantInterface>& getFactory() { return Factory<PlantInterface>::instance(); }

    //! Return the plant input dimension (u)
    virtual int getInputDimension() const = 0;
    //! Return the plant output dimension (y)
    virtual int getOutputDimension() const = 0;

    virtual bool requiresFutureControls() const = 0;

    virtual bool requiresFutureStates() const = 0;

    //! Initialize plant
    virtual bool initialize() { return true; }

    //! Stop plant (you might probably use this to set the plant into a safe setpoint)
    virtual void stop() {}

    virtual void reset() {}

    /**
     * @brief Send commands to plant
     *
     * @param[in]  u   Control command that is commanded to the plant [getInputDimension() x 1]
     * @param[in]  dt  Specify the intended duration for the signal (usually only relevant for simulation)
     * @param[in]  t   Current time stamp (can be sim-time or system-time, but compatible to state and control references)
     * @param[in,out] signal_target  Target for occuring signals [optional]
     *
     * @returns true if transmission of commands was successful, otherwise false.
     */
    virtual bool control(const ControlVector& u, const Duration& dt, const Time& t, SignalTargetInterface* signal_target = nullptr,
                         const std::string& ns = "");

    /**
     * @brief Send commands to plant
     *
     * @param[in]  u_sequence   Sequence of controls that are commanded to the plant
     * @param[in]  x_sequence   Sequence of states that are commanded to the plant
     * @param[in]  dt  Specify the intended duration for the signal (usually only relevant for simulation)
     * @param[in]  t   Current time stamp (can be sim-time or system-time, but compatible to state and control references)
     * @param[in,out] signal_target  Target for occuring signals [optional]
     *
     * @returns true if transmission of commands was successful, otherwise false.
     */
    virtual bool control(const TimeSeries::ConstPtr& u_sequence, const TimeSeries::ConstPtr& x_sequence, const Duration& dt, const Time& t,
                         SignalTargetInterface* signal_target = nullptr, const std::string& ns = "") = 0;

    /**
     * @brief Retrieve current plant output (measurements)
     *
     * @param[out]     output         Plant output will be written to this vector [getOutputDimension() x 1]
     * @param[in]  t   Current time stamp (can be sim-time or system-time, but compatible to state and control references)
     * @param[in,out]  signal_target  Target for occuring signals [optional]
     *
     * @returns true if receiving of measurements was successful, otherwise false.
     */
    virtual bool output(OutputVector& output, const Time& t, SignalTargetInterface* signal_target = nullptr, const std::string& ns = "") = 0;

    /**
     * @brief Retrieve available signals from the plant
     *
     * Register a-priori known signals at the signal target.
     * Registration is optional.
     * Note, during control() or output() execution further signals might occur without
     * registration (in case the they are not known in advance or the implementation lacks a proper registration).
     * @param[in,out] signal_target   Target for occuring signals [optional]
     */
    virtual void getAvailableSignals(SignalTargetInterface& signal_target, const std::string& ns = "") const {}

    /**
     * @brief Set/move plant to a desired state (if possible)
     * @param[in] Desired state vector
     * @return true if successful, false otherwise
     */
    virtual bool setState(const Eigen::Ref<const Eigen::VectorXd>& state) { return false; }

#ifdef MESSAGE_SUPPORT
    //! Export plant settings to message
    virtual void toMessage(messages::Plant& message) const {}
    //! Import plant settings from message
    virtual void fromMessage(const messages::Plant& message, std::stringstream* issues = nullptr) {}
#endif
};

using PlantFactory = Factory<PlantInterface>;
#define FACTORY_REGISTER_PLANT(type) FACTORY_REGISTER_OBJECT(type, PlantInterface)

}  // namespace corbo

#endif  // SRC_PLANTS_INCLUDE_CORBO_PLANTS_PLANT_INTERFACE_H_
