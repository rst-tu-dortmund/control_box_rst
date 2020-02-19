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

#ifndef SRC_CONTROLLERS_INCLUDE_CORBO_CONTROLLERS_CONTROLLER_INTERFACE_H_
#define SRC_CONTROLLERS_INCLUDE_CORBO_CONTROLLERS_CONTROLLER_INTERFACE_H_

#include <corbo-controllers/statistics.h>
#include <corbo-core/factory.h>
#include <corbo-core/reference_trajectory.h>
#include <corbo-core/signal_target_interface.h>
#include <corbo-core/time.h>
#include <corbo-core/types.h>
#include <corbo-systems/output_function_interface.h>

#ifdef MESSAGE_SUPPORT
#include <corbo-communication/messages/controllers/controllers.pb.h>
#endif

#include <memory>

namespace corbo {

/**
 * @brief Interface class for controllers
 *
 * @ingroup controllers
 *
 * This class specifies methods that are required to be implemented by specific
 * controllers in order to allow their general utilization in a variety of control tasks.
 *
 * @remark This interface is provided with factory support (ControllerFactory).
 *
 * @see PidController LqrController PredictiveController
 *
 * @author Christoph Rösmann (christoph.roesmann@tu-dortmund.de)
 */
class ControllerInterface
{
 public:
    using Ptr           = std::shared_ptr<ControllerInterface>;
    using UPtr          = std::unique_ptr<ControllerInterface>;
    using StateVector   = Eigen::VectorXd;
    using ControlVector = Eigen::VectorXd;

    //! Virtual destructor
    virtual ~ControllerInterface() {}

    //! Return a newly created shared instance of the implemented class
    virtual Ptr getInstance() const = 0;

    //! Get access to the associated factory
    static Factory<ControllerInterface>& getFactory() { return Factory<ControllerInterface>::instance(); }

    //! Return the control input dimension
    virtual int getControlInputDimension() const = 0;
    /**
     * @brief Return the dimension of the required plant state/output
     *
     * Depending on the controller type, the state dimension can also be just the plant output or
     * the complete measured or observed state vector.
     * @returns the controllers state dimension
     */
    virtual int getStateDimension() const = 0;

    virtual bool providesFutureControls() const = 0;

    virtual bool providesFutureStates() const = 0;

    /**
     * @brief Return true if the controller returns piecewise constant control pieces
     *
     * @return true if the controls are piecewise constant
     */
    virtual bool hasPiecewiseConstantControls() const = 0;

    /**
     * @brief Initialize the controller
     *
     * Initialization should be optional but it should facilitate memory allocation and trajectory initialization
     * respectively warm starting.
     * @param[in] x                 Current plant state [getStateDimension() x 1]
     * @param[in] expected_xref     State reference (writable in order to allow to the reference object to precompute caches)
     * @param[in] expected_uref     Control reference (writable in order to allow to the reference object to precompute caches)
     * @param[in] expected_dt       Expected sampling interval length (controller rate)
     * @param[in] t                 Current time stamp (can be sim-time or system-time, but compatible to state and control references)
     * @return true if initialization was successful, false otherwise.
     */
    virtual bool initialize(const StateVector& x, ReferenceTrajectoryInterface& expected_xref, ReferenceTrajectoryInterface& expected_uref,
                            const Duration& expected_dt, const Time& t, ReferenceTrajectoryInterface* expected_sref = nullptr)
    {
        return true;
    }

    /**
     * @brief Perform actual controller step / control law computation
     *
     * @param[in]     x              Current plant state [getStateDimension() x 1]
     * @param[in]     xref           State reference (writable in order to allow to the reference object to precompute caches)
     * @param[in]     uref           Control reference (writable in order to allow to the reference object to precompute caches)
     * @param[in]     dt             Last sampling interval length (does not (!) mean that the subsequent one will be identical)
     * @param[in]     t              Current time stamp (can be sim-time or system-time, but compatible to state and control references)
     * @param[out]    u              Computed control input for the plant [getControlInputDimension() x 1]
     * @param[in,out] signal_target  Target for occuring signals [optional]
     * @return true if step was successful and computation of \c u was successful, false otherwise.
     */
    virtual bool step(const StateVector& x, ReferenceTrajectoryInterface& xref, ReferenceTrajectoryInterface& uref, const Duration& dt, const Time& t,
                      ControlVector& u, SignalTargetInterface* signal_target = nullptr, const std::string& ns = "");

    virtual bool step(const StateVector& x, ReferenceTrajectoryInterface& xref, ReferenceTrajectoryInterface& uref, const Duration& dt, const Time& t,
                      TimeSeries::Ptr u_sequence, TimeSeries::Ptr x_sequence, SignalTargetInterface* signal_target = nullptr,
                      ReferenceTrajectoryInterface* sref = nullptr, ReferenceTrajectoryInterface* xinit = nullptr,
                      ReferenceTrajectoryInterface* uinit = nullptr, const std::string& ns = "") = 0;

    //! Return the duration for which the control u obtained from step() is valid (useful for asynchronous control)
    virtual double getControlDuration() const { return 0.0; }
    //! Specify whether the controllers step function is independent of dt and getControlDuration() returns a valid value
    virtual bool supportsAsynchronousControl() const { return false; }

    //! Reset internal controller state and caches
    virtual void reset() = 0;

    /**
     * @brief Retrieve available signals from the controller
     *
     * Register a-priori known signals at the signal target.
     * Registration is optional.
     * Note, during step() execution further signals might occur without
     * registration (in case the they are not known in advance or the implementation lacks a proper registration).
     * @param[in,out] signal_target   Target for occuring signals [optional]
     */
    virtual void getAvailableSignals(SignalTargetInterface& signal_target, const std::string& ns = "") const {}

    virtual void sendSignals(double t, SignalTargetInterface& signal_target, const std::string& ns = "") const {}

    virtual ControllerStatistics::Ptr getStatistics() const { return {}; }

#ifdef MESSAGE_SUPPORT
    //! Export controller to message
    virtual void toMessage(corbo::messages::Controller& message) const {}
    //! Import controller from message
    virtual void fromMessage(const corbo::messages::Controller& message, std::stringstream* issues = nullptr) {}
#endif
};

using ControllerFactory = Factory<ControllerInterface>;
#define FACTORY_REGISTER_CONTROLLER(type) FACTORY_REGISTER_OBJECT(type, ControllerInterface)

}  // namespace corbo

#endif  // SRC_CONTROLLERS_INCLUDE_CORBO_CONTROLLERS_CONTROLLER_INTERFACE_H_
