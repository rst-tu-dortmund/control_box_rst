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

#ifndef SRC_OBSERVERS_INCLUDE_CORBO_OBSERVERS_OBSERVER_INTERFACE_H_
#define SRC_OBSERVERS_INCLUDE_CORBO_OBSERVERS_OBSERVER_INTERFACE_H_

#include <corbo-core/factory.h>
#include <corbo-core/signal_target_interface.h>
#include <corbo-core/time.h>
#include <corbo-core/types.h>

#ifdef MESSAGE_SUPPORT
#include <corbo-communication/messages/observers/observers.pb.h>
#endif

#include <memory>

namespace corbo {

/**
 * @brief Interface class for observers
 *
 * @ingroup observers
 *
 * This class specifies methods that are required to be implemented by specific
 * observers in order to allow their general utilization in a variety of control tasks.
 *
 * An observer takes the measured plant output as input and
 * computes an estimate of a complete system state (w.r.t. a system model).
 *
 * @remark This interface is provided with factory support (ObserverFactory).
 *
 * @see NoObserver
 *
 * @author Christoph Rösmann (christoph.roesmann@tu-dortmund.de)
 */
class ObserverInterface
{
 public:
    using Ptr          = std::shared_ptr<ObserverInterface>;
    using UPtr         = std::unique_ptr<ObserverInterface>;
    using StateVector  = Eigen::VectorXd;
    using OutputVector = Eigen::VectorXd;

    //! Virtual destructor
    virtual ~ObserverInterface() {}

    //! Return a newly created shared instance of the implemented class
    virtual Ptr getInstance() const = 0;

    //! Get access to the associated factory
    static Factory<ObserverInterface>& getFactory() { return Factory<ObserverInterface>::instance(); }

    /**
     * @brief Return the dimension of the supported plant output
     *
     * This is not the output of the observer (actually it is the input)
     * but the output of the plant.
     *
     * @returns the supported plants output dimension
     */
    virtual int getOutputDimension() const = 0;

    /**
     * @brief Return the dimension of the observed state
     * @returns the state dimension
     */
    virtual int getStateDimension() const = 0;

    /**
     * @brief Perform actual observer step / state estimation
     *
     * @param[in]     y              Current plant output [getOutputDimension() x 1]
     * @param[out]    x              Estimated state vector [getStateDimension() x 1]
     * @param[in]     dt             Last sampling interval length (does not (!) mean that the subsequent one will be identical)
     * @param[in]     t              Current time stamp (can be sim-time or system-time, but compatible to state and control references)
     * @param[in,out] signal_target  Target for occuring signals [optional]
     * @return true if step was successful and computation of \c u was successful, false otherwise.
     */
    virtual bool observe(const OutputVector& y, StateVector& x, const Duration& dt, const Time& t,
                         SignalTargetInterface* signal_target = nullptr, const std::string& ns = "") = 0;

    //! Reset internal observer state and caches
    virtual void reset() {}

    /**
     * @brief Retrieve available signals from the observer
     *
     * Register a-priori known signals at the signal target.
     * Registration is optional.
     * Note, during observe() execution further signals might occur without
     * registration (in case the they are not known in advance or the implementation lacks a proper registration).
     * @param[in,out] signal_target   Target for occuring signals [optional]
     */
    virtual void getAvailableSignals(SignalTargetInterface& signal_target, const std::string& ns = "") const {}

#ifdef MESSAGE_SUPPORT
    //! Export observer to message
    virtual void toMessage(corbo::messages::Observer& message) const {}
    //! Import observer from message
    virtual void fromMessage(const corbo::messages::Observer& message, std::stringstream* issues = nullptr) {}
#endif
};

using ObserverFactory = Factory<ObserverInterface>;
#define FACTORY_REGISTER_OBSERVER(type) FACTORY_REGISTER_OBJECT(type, ObserverInterface)

/**
 * @brief Dummy observer that feed forwards the plant output
 *
 * @ingroup observers
 *
 * This dummy observer taks the plant output \f$ y \f$ directly directly as
 * state \f$ x \f$, in particular $\f x = y \f$.
 *
 * @see ObserverInterface
 *
 * @author Christoph Rösmann (christoph.roesmann@tu-dortmund.de)
 */
class NoObserver : public ObserverInterface
{
 public:
    // implements interface method
    Ptr getInstance() const override { return std::make_shared<NoObserver>(); }

    // implements interface method
    int getOutputDimension() const override { return property::INHERITED; }
    // implements interface method
    int getStateDimension() const override { return property::INHERITED; }
    // implements interface method
    bool observe(const OutputVector& y, StateVector& x, const Duration& /*dt*/, const Time& /*t*/,
                 SignalTargetInterface* signal_target = nullptr, const std::string& ns = "") override
    {
        x = y;
        return true;
    }

#ifdef MESSAGE_SUPPORT
    // implements interface method
    void toMessage(corbo::messages::Observer& message) const override { message.mutable_no_observer(); }
#endif
};

FACTORY_REGISTER_OBSERVER(NoObserver)

}  // namespace corbo

#endif  // SRC_OBSERVERS_INCLUDE_CORBO_OBSERVERS_OBSERVER_INTERFACE_H_
