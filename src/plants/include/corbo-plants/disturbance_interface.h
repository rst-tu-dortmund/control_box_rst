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

#ifndef SRC_PLANTS_INCLUDE_CORBO_PLANTS_DISTURBANCE_INTERFACE_H_
#define SRC_PLANTS_INCLUDE_CORBO_PLANTS_DISTURBANCE_INTERFACE_H_

#include <corbo-core/factory.h>
#include <corbo-core/time.h>
#include <corbo-core/types.h>

#ifdef MESSAGE_SUPPORT
#include <corbo-communication/messages/plants/disturbance.pb.h>
#endif

#include <memory>

namespace corbo {

class DisturbanceInterface
{
 public:
    using Ptr = std::shared_ptr<DisturbanceInterface>;

    //! Virtual destructor
    virtual ~DisturbanceInterface() {}

    //! Return a newly created shared instance of the implemented class
    virtual Ptr getInstance() const = 0;

    /**
     * @brief Modify values according to the underlying disturbance model
     *
     * @remarks This method allows alias between input and output parameters
     *
     * @param[in] t                  Current time stamp
     * @param[in] values             Value vector to be disturbed
     * @param[out] disturbed_values  Disturbed value vector (alias-safe) [must be preallocated]
     */
    virtual void disturb(const Time& t, const Eigen::Ref<const Eigen::VectorXd>& values, Eigen::Ref<Eigen::VectorXd> disturbed_values) = 0;

    /**
     * @brief Check the underlying parameter configuration for validity
     *
     * This method might be useful if the class has been configured via a message from another class
     * and the internal dimensions must fulfil some requirements.
     *
     * @param[in]  values_dim  Expected dimension of the value vector to be disturbed
     * @param[out] issues      Issue related messages are forwarded to this stream (optional)
     */
    virtual bool checkParameters(int values_dim, std::stringstream* issues) const { return true; }

    //! reset internal state
    virtual void reset() = 0;

    //! Get access to the associated factory
    static Factory<DisturbanceInterface>& getFactory() { return Factory<DisturbanceInterface>::instance(); }

#ifdef MESSAGE_SUPPORT
    //! Export plant settings to message
    virtual void toMessage(messages::Disturbance& message) const {}
    //! Import plant settings from message
    virtual void fromMessage(const messages::Disturbance& message, std::stringstream* issues = nullptr) {}
#endif
};

using DisturbanceFactory = Factory<DisturbanceInterface>;
#define FACTORY_REGISTER_DISTURBANCE(type) FACTORY_REGISTER_OBJECT(type, DisturbanceInterface)

}  // namespace corbo

#endif  // SRC_PLANTS_INCLUDE_CORBO_PLANTS_DISTURBANCE_INTERFACE_H_
