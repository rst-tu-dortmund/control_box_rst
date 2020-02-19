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

#ifndef SRC_SYSTEMS_INCLUDE_CORBO_SYSTEMS_SYSTEM_DYNAMICS_INTERFACE_H_
#define SRC_SYSTEMS_INCLUDE_CORBO_SYSTEMS_SYSTEM_DYNAMICS_INTERFACE_H_

#include <corbo-core/factory.h>
#include <corbo-core/time.h>
#include <corbo-core/types.h>

#ifdef MESSAGE_SUPPORT
#include <corbo-communication/messages/systems/system_dynamics.pb.h>
#endif

#include <memory>

namespace corbo {

class FiniteDifferencesInterface;  // forward declaration

/**
 * @brief Interface class for system dynamic equations
 *
 * @ingroup systems
 *
 * This class specifies methods that are required to be implemented by specific
 * subclasses in order to allow their general utilization in a variety of control tasks.
 *
 * System dynamics can be either defined in continuous-time, e.g.
 * \f[
 *      \dot{x} = f(x, u)
 * \f]
 * or in discrete-time, e.g.
 * \f[
 *      x_{k+1} = f(x_k, u_k)
 * \f].
 * Subclasses must overload isContinuousTime() appropriately.
 *
 * @remark This interface is provided with factory support (SystemDynamicsFactory).
 *
 * @see SystemOutputInterface PlantInterface
 *
 * @author Christoph Rösmann (christoph.roesmann@tu-dortmund.de)
 */
class SystemDynamicsInterface
{
 public:
    using Ptr           = std::shared_ptr<SystemDynamicsInterface>;
    using StateVector   = Eigen::VectorXd;
    using ControlVector = Eigen::VectorXd;

    //! Default constructor
    SystemDynamicsInterface();

    //! Default destructor
    virtual ~SystemDynamicsInterface() = default;

    //! Return a newly created shared instance of the implemented class
    virtual Ptr getInstance() const = 0;

    /**
     * @brief Check if the system dynamics are defined in continuous-time
     *
     * Continous-time equations are defined as \f$ \dot{x} = f(x,u) \f$
     * and discrete-time equations as \f$ x_{k+1} = f(x_k,u_k) \f$.
     * @returns true if the system is specified in continuous time, false otherwise.
     */
    virtual bool isContinuousTime() const = 0;

    /**
     * @brief Check if the system dynamics are linear
     *
     * Linear system can be written in the form
     * \f$ \dot{x} = Ax + Bu \f$ (continuous-time) or
     * \f$ x_{k+1} = A x_k + B u_k  \f$ (discrete-time).
     * Consequently, getLinearA() and getLinearB() are independet of x0 and u0.
     * @returns true if the system behavior is linear, false otherwise.
     */
    virtual bool isLinear() const = 0;

    //! Return the plant input dimension (u)
    virtual int getInputDimension() const = 0;
    //! Return state dimension (x)
    virtual int getStateDimension() const = 0;
    //! Return deadtime which might be taken into account by controllers/simulators if supported
    virtual double getDeadTime() const { return _deadtime; }

    /**
     * @brief Evaluate the system dynamics equation
     *
     * Implement this method to specify the actual system dynamics equation
     * \f$ \dot{x} = f(x,u) \f$ (continuous-time) or
     * \f$ x_{k+1} = f(x_k, u_k) \f$ (discrete-time).
     *
     * @param[in]  x   State vector x [getStateDimension() x 1]
     * @param[in]  u   Control input vector u [getInputDimension() x 1]
     * @param[out] f   Resulting function value \f$ \dot{x} \f$ respectively \f$ x_{k+1} \f$
     *                 (f must be preallocated with dimension [getStateDimension() x 1])
     */
    virtual void dynamics(const Eigen::Ref<const StateVector>& x, const Eigen::Ref<const ControlVector>& u, Eigen::Ref<StateVector> f) const = 0;

    /**
     * @brief Return linear system matrix A (linearized system dynamics)
     *
     * Linearizes the system dynamics equation \f$ \dot{x} = f(x,u) \f$ (cont. case)
     * specified with dynamics() using finite differences (refer to setLinearizationMethod()).
     * The linearizes system dynamics at x=x0 and u=u0 are defined as:
     * \f[
     *     \dot{x} = A x + B u
     * \f]
     * And for the discrete-time case:
     * \f[
     *     x_{k+1} = A x_k + B u_k
     * \f].
     * This method returns matrix \f$ A \f$. For matrix \f$ B \f$ refer to getLinearB().
     *
     * @remarks This method might be overriden in order to specify analytically derived matrices.
     *
     * @param[in]   x0    State at which the system dynamics should be linearized [getStateDimension() x 1]
     * @param[in]   u0    Control input at which the dynamics should be linearized [getInputDimension() x 1]
     * @param[out]  A     Linear system matrix A [getStateDimension() x getStateDimension()] (must be preallocated)
     */
    virtual void getLinearA(const StateVector& x0, const ControlVector& u0, Eigen::MatrixXd& A) const;

    /**
     * @brief Return linear control input matrix B (linearized system dynamics)
     *
     * Refer to the description of method getLinearA().
     *
     * @remarks This method might be overriden in order to specify analytically derived matrices.
     *
     * @param[in]   x0    State at which the system dynamics should be linearized [getStateDimension() x 1]
     * @param[in]   u0    Control input at which the dynamics should be linearized [getInputDimension() x 1]
     * @param[out]  B     Linear input matrix B [getStateDimension() x getInputDimension()] (must be preallocated)
     */
    virtual void getLinearB(const StateVector& x0, const ControlVector& u0, Eigen::MatrixXd& B) const;

    /**
     * @brief Set linearization method in case getLinearA() or getLinearB() are not overriden
     * @param[in] lin_method  shared instance of a FiniteDifferencesInterface implementation
     */
    void setLinearizationMethod(std::shared_ptr<FiniteDifferencesInterface> lin_method);

    virtual void reset() {}

#ifdef MESSAGE_SUPPORT
    //! Export system dynamics to message
    virtual void toMessage(messages::SystemDynamics& message) const;
    //! Import system dynamics from message
    virtual void fromMessage(const messages::SystemDynamics& message, std::stringstream* issues = nullptr);
#endif

 private:
    std::shared_ptr<FiniteDifferencesInterface> _linearization_method;
    double _deadtime = 0.0;
};

using SystemDynamicsFactory = Factory<SystemDynamicsInterface>;
#define FACTORY_REGISTER_SYSTEM_DYNAMICS(type) FACTORY_REGISTER_OBJECT(type, SystemDynamicsInterface)

}  // namespace corbo

#endif  // SRC_SYSTEMS_INCLUDE_CORBO_SYSTEMS_SYSTEM_DYNAMICS_INTERFACE_H_
