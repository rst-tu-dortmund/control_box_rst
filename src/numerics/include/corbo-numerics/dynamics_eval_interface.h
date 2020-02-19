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

#ifndef SRC_NUMERICS_INCLUDE_CORBO_NUMERICS_DYNAMICS_EVAL_INTERFACE_H_
#define SRC_NUMERICS_INCLUDE_CORBO_NUMERICS_DYNAMICS_EVAL_INTERFACE_H_

#include <corbo-core/factory.h>
#include <corbo-core/range.h>
#include <corbo-core/types.h>
#include <corbo-systems/system_dynamics_interface.h>

#ifdef MESSAGE_SUPPORT
#include <corbo-communication/messages/numerics/dynamics_eval.pb.h>
#endif

#include <functional>
#include <memory>

namespace corbo {

/**
 * @brief Interface class for solving and evaluating dynamics
 *
 * @ingroup numerics
 * *
 * @remark This interface is provided with factory support (DynamicsEvalInterface).
 *
 * @see ForwardDifferences CentralDifferences NumericalIntegratorExplicitInterface
 *      NumericalIntegratorExplicitInterface
 *
 * @author Christoph Rösmann (christoph.roesmann@tu-dortmund.de)
 */
class DynamicsEvalInterface
{
 public:
    using Ptr  = std::shared_ptr<DynamicsEvalInterface>;
    using UPtr = std::unique_ptr<DynamicsEvalInterface>;

    using StateVector = Eigen::VectorXd;
    using InputVector = Eigen::VectorXd;

    //! Virtual destructor
    virtual ~DynamicsEvalInterface() {}

    //! Return a newly allocated instances of the inherited class.
    virtual Ptr getInstance() const = 0;

    //! Get access to the accociated factory
    static Factory<DynamicsEvalInterface>& getFactory() { return Factory<DynamicsEvalInterface>::instance(); }

    /**
     * @brief Compute error between two consecutive (discrete) states
     *
     * \f[
     *      e =  f(x,u,t) - \dot{x}
     * \f].
     *
     * @remarks We require no alias between \c error and all other parameters!
     *
     * @param[in]  x1      Initial state vector [SystemDynamicsInterface::getStateDimension() x 1]
     * @param[in]  u1      Constant control input vector [SystemDynamicsInterface::getInputDimension() x 1]
     * @param[in]  x2      Final state vector [SystemDynamicsInterface::getStateDimension() x 1]
     * @param[in]  dt      Time interval length
     * @param[in]  system  System dynamics object
     * @param[out] error   Resulting error [SystemDynamicsInterface::getStateDimension() x 1] (must be preallocated)
     */
    virtual void computeEqualityConstraint(const StateVector& x1, const InputVector& u1, const StateVector& x2, double dt,
                                           const SystemDynamicsInterface& system, Eigen::Ref<Eigen::VectorXd> error) = 0;

    virtual bool interpolate(const Eigen::Ref<const Eigen::VectorXd>& x1, const Eigen::Ref<const Eigen::VectorXd>& u1,
                             const Eigen::Ref<const Eigen::VectorXd>& x2, const Eigen::Ref<const Eigen::VectorXd>& u2, double dt,
                             const SystemDynamicsInterface& system, const Range& range, std::vector<Eigen::VectorXd>& states,
                             std::vector<Eigen::VectorXd>& controls) = 0;

    virtual bool interpolate(const std::vector<const Eigen::VectorXd*>& x, const std::vector<const Eigen::VectorXd*>& u, double dt,
                             const SystemDynamicsInterface& system, const Range& range, std::vector<Eigen::VectorXd>& states,
                             std::vector<Eigen::VectorXd>& controls)
    {
        if (x.size() > 1 && u.size() > 1)
            return interpolate(*x[0], *u[0], *x.back(), *u.back(), dt, system, range, states, controls);
        else if (x.size() > 1 && u.size() == 1)
            return interpolate(*x[0], *u[0], *x.back(), *u[0], dt, system, range, states, controls);
        PRINT_WARNING_NAMED("There is no appropriate overload for this method...");
        return false;
    }

#ifdef MESSAGE_SUPPORT
    //! Export selected class and parameters to a given message
    virtual void toMessage(corbo::messages::DynamicsEval& message) const {}
    //! Import selected class and parameters from a given message (optionally pass any issues to the caller).
    virtual void fromMessage(const corbo::messages::DynamicsEval& message, std::stringstream* issues = nullptr) {}
#endif
};

using DynamicsEvalFactory = Factory<DynamicsEvalInterface>;
#define FACTORY_REGISTER_DYNAMICS_EVAL(type) FACTORY_REGISTER_OBJECT_ID(type, DynamicsEvalInterface, 0)

}  // namespace corbo

#endif  // SRC_NUMERICS_INCLUDE_CORBO_NUMERICS_DYNAMICS_EVAL_INTERFACE_H_
