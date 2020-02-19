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

#ifndef SRC_NUMERICS_INCLUDE_CORBO_NUMERICS_INTEGRATOR_INTERFACE_H_
#define SRC_NUMERICS_INCLUDE_CORBO_NUMERICS_INTEGRATOR_INTERFACE_H_

#include <corbo-core/factory.h>
#include <corbo-core/time.h>
#include <corbo-core/types.h>
#include <corbo-numerics/dynamics_eval_interface.h>
#include <corbo-systems/system_dynamics_interface.h>

#ifdef MESSAGE_SUPPORT
#include <corbo-communication/messages/numerics/explicit_integrators.pb.h>
#endif

#include <cmath>
#include <functional>
#include <memory>

namespace corbo {

///**
// * @brief Interface for numerical integrators (explicit and implicit)
// *
// * @ingroup numerics
// *
// * Explicit numerical integrators in the context of ordinary or partial
// * differential equations determine the integral value only
// * based on the previous or current time instance, e.g.:
// * \f[
// *      x(t_0 + \Delta t) = \int_{t_0}^{t_0+\Delta T} f(x(t_0), t) dt
// * \f].
// *
// * Implicit methods require the solution of a (nonlinear) equation, since
// * they also depend on the solution itself:
// * \f[
// *      x(t_0 + \Delta t) = \int_{t_0}^{t_0+\Delta T} f(x(t_0+ \Delta T), x(t_0), t) dt
// * \f].
// *
// * This generic interface mainly provides access to the underlying quadrature rules.
// * These might be used as part of an optimization or might be solved via root-finding directly.
// * For a direct use of integrators you might refer to the subclass NumericalIntegratorExplicitInterface.
// *
// * @remark This interface is provided with factory support.
// *
// * @see NumericalIntegratorExplicitInterface FiniteDifferencesInterface
// *
// * @author Christoph Rösmann (christoph.roesmann@tu-dortmund.de)
// *
// * @todo This interface is not yet completed (messages, subclasses, ...)
// */
// class NumericalIntegratorInterface : public DynamicsEvalInterface
// {
// public:
//    using Ptr  = std::shared_ptr<NumericalIntegratorInterface>;
//    using UPtr = std::unique_ptr<NumericalIntegratorInterface>;
//
//    using StateVector = Eigen::VectorXd;
//    using InputVector = Eigen::VectorXd;
//    using UnaryFunction = const std::function<void(const Eigen::VectorXd&, Eigen::Ref<Eigen::VectorXd>)>;
//
//    //! Virtual destructor
//    virtual ~NumericalIntegratorInterface() {}
//    //! Return a newly created shared instance of the implemented class
//    virtual DynamicsEvalInterface::Ptr getInstance() const override = 0;
//
//    /**
//     * @brief Allocate memory for a given state dimension
//     *
//     * @remarks This method is optional and not guaranteed to be called, but if so this might speed up computation.
//     * @param[in] state_dim  Expected dimension of the state vector
//     */
//     virtual void initialize(int state_dim) {}
//
//    /**
//     * @brief Compute equality constraint containing the implicit quadrature/interpolation formula (system dynamics specialization)
//     *
//     * \f[
//     *      e =  \int_{t=0}^{t=\Delta T} f(x(t), u_1) dt - x_2 = 0
//     * \f].
//     *
//     * @param[in]  x1      Initial state vector [SystemDynamicsInterface::getStateDimension() x 1]
//     * @param[in]  u1      Constant control input vector [SystemDynamicsInterface::getInputDimension() x 1]
//     * @param[in]  x2      Final state vector [SystemDynamicsInterface::getStateDimension() x 1]
//     * @param[in]  dt      Time interval length
//     * @param[in]  system  System dynamics object
//     * @param[out] error   Resulting error [SystemDynamicsInterface::getStateDimension() x 1] (must be preallocated)
//     */
//    void computeEqualityConstraint(const StateVector& x1, const InputVector& u1, const StateVector& x2, double dt, const SystemDynamicsInterface&
//    system,
//                                           Eigen::Ref<Eigen::VectorXd> error) override = 0;
//
//    /**
//     * @brief Compute equality constraint containing the implicit quadrature/interpolation formula (unary function specialization)
//     *
//     * \f[
//     *      e =  \int_{t=0}^{t=\Delta T} f(x(t)) dt - x_2 = 0
//     * \f].
//     *
//     * @param[in]  x1      Initial function argument
//     * @param[in]  x2      Final function argument
//     * @param[in]  dt      Time interval length
//     * @param[in]  fun     f(x) to be integrated in the interval [0, dt]
//     * @param[out] error   Resulting error [SystemDynamicsInterface::getStateDimension() x 1] (must be preallocated)
//     */
//    virtual void computeEqualityConstraint(const Eigen::Ref<const Eigen::VectorXd>& x1, const Eigen::Ref<const Eigen::VectorXd>& x2, double dt,
//                                           const UnaryFunction& fun, Eigen::Ref<Eigen::VectorXd> error) = 0;
//
// #ifdef MESSAGE_SUPPORT
//    virtual void toMessage(messages::Integrator& message) const
//    {}
//    virtual void fromMessage(const messages::Integrator& message, std::stringstream* issues = nullptr) {}
// #endif
// };
//
// #define FACTORY_REGISTER_INTEGRATOR(type) FACTORY_REGISTER_DYNAMICS_EVAL(type)

/**
 * @brief Interface for explicit numerical integrators
 *
 * @ingroup numerics
 *
 * Explicit numerical integrators in the context of ordinary or partial
 * differential equations determine the integral value only
 * based on the previous or current time instance, e.g.:
 * \f[
 *      x(t_0 + \Delta t) = \int_{t_0}^{t_0+\Delta T} f(x(t_0), t) dt
 * \f].
 *
 * @remark This interface is provided with factory support (DynamicsEval).
 *
 * @see NumericalIntegratorInterface FiniteDifferencesInterface
 *
 * @author Christoph Rösmann (christoph.roesmann@tu-dortmund.de)
 */
class NumericalIntegratorExplicitInterface : public DynamicsEvalInterface
{
 public:
    using Ptr  = std::shared_ptr<NumericalIntegratorExplicitInterface>;
    using UPtr = std::unique_ptr<NumericalIntegratorExplicitInterface>;

    using UnaryFunction = const std::function<void(const Eigen::VectorXd&, Eigen::Ref<Eigen::VectorXd>)>;
    // using DynamicsEvalInterface::StateVector;
    // using DynamicsEvalInterface::InputVector;

    //! Virtual destructor
    virtual ~NumericalIntegratorExplicitInterface() {}
    //! Return a newly created shared instance of the implemented class
    DynamicsEvalInterface::Ptr getInstance() const override = 0;

    //! Return the convergence order
    virtual int getConvergenceOrder() const = 0;

    /**
     * @brief Allocate memory for a given state dimension
     *
     * @remarks This method is optional and not guaranteed to be called, but if so this might speed up computation.
     * @param[in] state_dim  Expected dimension of the state vector
     */
    virtual void initialize(int state_dim) {}

    /**
     * @brief Solution of the initial value problem
     *
     * \f[
     *      x(t=\Delta T) = \int_{t=0}^{t=\Delta T} f(x(t)) dt
     * \f]
     * with $ x(t=0) = x_1 $.
     *
     * @param[in]  x1      Initial state vector [SystemDynamicsInterface::getStateDimension() x 1]
     * @param[in]  u1      Constant control input vector [SystemDynamicsInterface::getInputDimension() x 1]
     * @param[in]  dt      Time interval length
     * @param[in]  system  System dynamics object
     * @param[out] x2      Resulting state vector [SystemDynamicsInterface::getStateDimension() x 1] (must be preallocated)
     */
    virtual void solveIVP(const Eigen::VectorXd& x1, double dt, const UnaryFunction& fun, Eigen::Ref<Eigen::VectorXd> x2) = 0;

    /**
     * @brief Solution of the initial value problem
     *
     * \f[
     *      x(t=\Delta T) = \int_{t=0}^{t=\Delta T} f(x(t), u_1) dt
     * \f]
     * with $ x(t=0) = x_1 $.
     *
     * @param[in]  x1      Initial state vector [SystemDynamicsInterface::getStateDimension() x 1]
     * @param[in]  u1      Constant control input vector [SystemDynamicsInterface::getInputDimension() x 1]
     * @param[in]  dt      Time interval length
     * @param[in]  system  System dynamics object
     * @param[out] x2      Resulting state vector [SystemDynamicsInterface::getStateDimension() x 1] (must be preallocated)
     */
    virtual void solveIVP(const StateVector& x1, const InputVector& u1, double dt, const SystemDynamicsInterface& system,
                          Eigen::Ref<Eigen::VectorXd> x2) = 0;
    
    void computeEqualityConstraint(const StateVector& x1, const InputVector& u1, const StateVector& x2, double dt,
                                   const SystemDynamicsInterface& system, Eigen::Ref<Eigen::VectorXd> error) override
    {
        solveIVP(x1, u1, dt, system, error);
        error -= x2;
    }

    void computeEqualityConstraint(const Eigen::Ref<const Eigen::VectorXd>& x1, const Eigen::Ref<const Eigen::VectorXd>& x2, double dt,
                                   const UnaryFunction& fun, Eigen::Ref<Eigen::VectorXd> error)
    {
        solveIVP(x1, dt, fun, error);
        error -= x2;
    }
    
    bool interpolate(const Eigen::Ref<const Eigen::VectorXd>& x1, const Eigen::Ref<const Eigen::VectorXd>& u1,
                     const Eigen::Ref<const Eigen::VectorXd>& /*x2*/, const Eigen::Ref<const Eigen::VectorXd>& /*u2*/, double dt,
                     const SystemDynamicsInterface& system, const Range& range, std::vector<Eigen::VectorXd>& states,
                     std::vector<Eigen::VectorXd>& controls) override
    {
        if (range.getStart() < 0 || range.getEnd() > dt)
        {
            PRINT_ERROR_NAMED("range must be in the interval [0, dt]");
            return false;
        }
        if (dt < 1e-15)
        {
            states.push_back(x1);
            controls.push_back(u1);
            return true;
        }

        Eigen::VectorXd x_cur = x1;
        Eigen::VectorXd x_next(x1.size());

        int n = range.getNumInRange();
        for (int i = 0; i < n; ++i)
        {
            states.push_back(x_cur);
            controls.push_back(u1);  // constant controls per default

            if (i < n - 1)
            {
                solveIVP(x_cur, u1, range.getStep(), system, x_next);
            }
            x_cur = x_next;
        }

        if (range.includeEnd())
        {
            solveIVP(x_cur, u1, range.getRemainder(), system, x_next);
        }
        return true;
    }

#ifdef MESSAGE_SUPPORT
    //! Export integrator settings to message
    virtual void toMessage(messages::ExplicitIntegrator& message) const {}
    //! Import integrator settings from message
    virtual void fromMessage(const messages::ExplicitIntegrator& message, std::stringstream* issues = nullptr) {}

    //! Export selected class and parameters to a Integrator message
    // void toMessage(corbo::messages::Integrator& message) const override { toMessage(*message.mutable_explicit_integrator()); }
    //! Import selected class and parameters from a Integrator message (optionally pass any issues to the caller).
    // void fromMessage(const corbo::messages::Integrator& message, std::stringstream* issues = nullptr) override
    //{
    //    if (message.has_explicit_integrator()) fromMessage(message.explicit_integrator());
    //}

    //! Export selected class and parameters to a DynamicsEval message
    void toMessage(corbo::messages::DynamicsEval& message) const override { toMessage(*message.mutable_integrator()); }
    //! Import selected class and parameters from a DynamicsEval message (optionally pass any issues to the caller).
    void fromMessage(const corbo::messages::DynamicsEval& message, std::stringstream* issues = nullptr) override
    {
        if (message.has_integrator()) fromMessage(message.integrator());
    }
#endif
};

// #define FACTORY_REGISTER_EXPLICIT_INTEGRATOR(type) FACTORY_REGISTER_INTEGRATOR(type)
#define FACTORY_REGISTER_EXPLICIT_INTEGRATOR(type) FACTORY_REGISTER_DYNAMICS_EVAL(type)

// class NumericalIntegratorImplicitInterface : public NumericalIntegratorInterface
// {
// public:
//    using Ptr  = std::shared_ptr<NumericalIntegratorImplicitInterface>;
//    using UPtr = std::unique_ptr<NumericalIntegratorImplicitInterface>;
//
//    //! Virtual destructor
//    virtual ~NumericalIntegratorImplicitInterface() {}
//    //! Return a newly created shared instance of the implemented class
//    virtual DynamicsEvalInterface::Ptr getInstance() const override = 0;
//
//    /**
//     * @brief Allocate memory for a given state dimension
//     *
//     * @remarks This method is optional and not guaranteed to be called, but if so this might speed up computation.
//     * @param[in] state_dim  Expected dimension of the state vector
//     */
//    //virtual void initialize(int state_dim) {}
//
//
//    virtual void quadrature(const StateVector& x1, const InputVector& u1, const StateVector& x2, double dt, const SystemDynamicsInterface& system,
//                           Eigen::Ref<Eigen::VectorXd> integral_value) = 0;
//
//    // integral of x(t) over the interval [0, dt]
//    virtual void quadrature(const Eigen::Ref<const Eigen::VectorXd>& x1, const Eigen::Ref<const Eigen::VectorXd>& x2, double dt, const
//    UnaryFunction& fun,
//                            Eigen::Ref<Eigen::VectorXd> integral_value) = 0;
//
//    void computeEqualityConstraint(const StateVector& x1, const InputVector& u1, const StateVector& x2, double dt, const SystemDynamicsInterface&
//    system,
//                                   Eigen::Ref<Eigen::VectorXd> error) override
//    {
//        quadrature(x1, u1, x2, dt, system, error);
//        error -= x2;
//    }
//
//    void computeEqualityConstraint(const Eigen::Ref<const Eigen::VectorXd>& x1, const Eigen::Ref<const Eigen::VectorXd>& x2, double dt,
//                                   const UnaryFunction& fun, Eigen::Ref<Eigen::VectorXd> error) override
//    {
//        quadrature(x1, x2, dt, fun, error);
//        error -= x2;
//    }
//
// #ifdef MESSAGE_SUPPORT
//    //! Export integrator settings to message
//    virtual void toMessage(messages::ImplicitIntegrator& message) const {}
//    //! Import integrator settings from message
//    virtual void fromMessage(const messages::ImplicitIntegrator& message, std::stringstream* issues = nullptr) {}
//
//    //! Export selected class and parameters to a Integrator message
//    void toMessage(corbo::messages::Integrator& message) const override { toMessage(*message.mutable_implicit_integrator()); }
//    //! Import selected class and parameters from a Integrator message (optionally pass any issues to the caller).
//    void fromMessage(const corbo::messages::Integrator& message, std::stringstream* issues = nullptr) override
//    {
//        if (message.has_explicit_integrator()) fromMessage(message.implicit_integrator());
//    }
//
//    //! Export selected class and parameters to a DynamicsEval message
//    //void toMessage(corbo::messages::DynamicsEval& message) const override { toMessage(*message.mutable_integrator()); }
//    //! Import selected class and parameters from a DynamicsEval message (optionally pass any issues to the caller).
//    //void fromMessage(const corbo::messages::DynamicsEval& message, std::stringstream* issues = nullptr) override
//    //{
//    //    if (message.has_integrator()) fromMessage(message.integrator());
//    //}
// #endif
// };
//
// #define FACTORY_REGISTER_IMPLICIT_INTEGRATOR(type) FACTORY_REGISTER_INTEGRATOR(type)

}  // namespace corbo

#endif  // SRC_NUMERICS_INCLUDE_CORBO_NUMERICS_INTEGRATOR_INTERFACE_H_
