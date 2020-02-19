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

#ifndef SRC_NUMERICS_INCLUDE_CORBO_NUMERICS_QUADRATURE_INTERFACE_H_
#define SRC_NUMERICS_INCLUDE_CORBO_NUMERICS_QUADRATURE_INTERFACE_H_

#include <corbo-core/factory.h>
#include <corbo-core/time.h>
#include <corbo-core/time_series.h>
#include <corbo-core/types.h>
#include <corbo-numerics/collocation_interface.h>

#ifdef MESSAGE_SUPPORT
#include <corbo-communication/messages/numerics/quadrature.pb.h>
#endif

#include <functional>
#include <memory>

namespace corbo {

class QuadratureCollocationInterface : public CollocationInterface
{
 public:
    using Ptr  = std::shared_ptr<QuadratureCollocationInterface>;
    using UPtr = std::unique_ptr<QuadratureCollocationInterface>;

    using UnaryFunction  = const std::function<void(const Eigen::VectorXd&, Eigen::Ref<Eigen::VectorXd>)>;
    using BinaryFunction = const std::function<void(const Eigen::VectorXd&, const Eigen::VectorXd&, Eigen::Ref<Eigen::VectorXd>)>;

    //! Virtual destructor
    virtual ~QuadratureCollocationInterface() {}
    //! Return a newly created shared instance of the implemented class
    DynamicsEvalInterface::Ptr getInstance() const override = 0;

    // useful for collocation
    virtual bool isSupportingCompressedStatesMode() const    = 0;
    virtual bool isIntermediateControlSubjectToOptim() const = 0;
    virtual int getNumIntermediatePoints() const             = 0;

    // TODO(roesmann): PROPOSAL:
    virtual int getNumIntermediateControls() const = 0;
    virtual int getNumIntermediateStates() const   = 0;

    /**
     * @brief Allocate memory for a given state dimension
     *
     * @remarks This method is optional and not guaranteed to be called, but if so this might speed up computation.
     * @param[in] state_dim  Expected dimension of the state vector
     */
    // virtual void initialize(int state_dim) {}

    virtual void quadrature(const std::vector<const Eigen::VectorXd*>& x1_points, const std::vector<const Eigen::VectorXd*>& u1_points,
                            const Eigen::VectorXd& u2, const Eigen::VectorXd& x2, double dt, const SystemDynamicsInterface& system,
                            Eigen::Ref<Eigen::VectorXd> integral_value) = 0;
    
    // integral of f(x(t), u(t)) over the interval [0, dt]
    virtual void quadrature(const std::vector<const Eigen::VectorXd*>& x1_points, const std::vector<const Eigen::VectorXd*>& u1_points,
                            const Eigen::VectorXd& u2, const Eigen::VectorXd& x2, double dt,
                            const BinaryFunction& fun, Eigen::Ref<Eigen::VectorXd> integral_value) = 0;

    virtual bool interpolate(const std::vector<const Eigen::VectorXd*>& x1_points, const std::vector<const Eigen::VectorXd*>& u1_points,
                             const Eigen::VectorXd& u2, const Eigen::VectorXd& x2, double t1, double dt, double dt_interp,
                             const SystemDynamicsInterface& system, bool skip_u2_x2, TimeSeries& ts_x, TimeSeries& ts_u) = 0;

    virtual void quadrature(const StateVector& x1, const InputVector& u1, const StateVector& x2, double dt, const SystemDynamicsInterface& system,
                            Eigen::Ref<Eigen::VectorXd> integral_value) = 0;

    virtual void quadrature(const StateVector& x1, const InputVector& u1, const StateVector& x2, const InputVector& u2, double dt,
                            const SystemDynamicsInterface& system, Eigen::Ref<Eigen::VectorXd> integral_value) = 0;

    virtual void quadrature(const StateVector& x1, const InputVector& u1, const StateVector& x2, double dt, const SystemDynamicsInterface& system,
                            Eigen::Ref<Eigen::VectorXd> integral_value, std::vector<Eigen::VectorXd*>& intermediate_x) = 0;

    virtual void quadrature(const StateVector& x1, const InputVector& u1, const StateVector& x2, const InputVector& u2, double dt,
                            const SystemDynamicsInterface& system, const std::vector<Eigen::VectorXd*>& intermediate_x,
                            const std::vector<Eigen::VectorXd*>& intermediate_u, Eigen::Ref<Eigen::VectorXd> integral_value,
                            Eigen::Ref<Eigen::VectorXd> interm_x_error) = 0;

    virtual void quadrature(const StateVector& x1, const InputVector& u1, const StateVector& x2, const InputVector& u2, double dt,
                            const SystemDynamicsInterface& system, const std::vector<Eigen::VectorXd*>& intermediate_u,
                            Eigen::Ref<Eigen::VectorXd> integral_value, std::vector<Eigen::VectorXd*>& intermediate_x) = 0;

    virtual void quadrature(const StateVector& x1, const InputVector& u1, const StateVector& x2, const InputVector& u2, double dt,
                            const SystemDynamicsInterface& system, Eigen::Ref<Eigen::VectorXd> integral_value,
                            std::vector<Eigen::VectorXd*>& intermediate_x, std::vector<Eigen::VectorXd*>& intermediate_u) = 0;

    // integral of f(x(t)) over the interval [0, dt]
    virtual void quadrature(const Eigen::Ref<const Eigen::VectorXd>& x1, const Eigen::Ref<const Eigen::VectorXd>& x2, double dt,
                            const UnaryFunction& fun, Eigen::Ref<Eigen::VectorXd> integral_value) = 0;

    // integral of f(x(t)) over the interval [0, dt] with precomputed intermediate points
    virtual void quadrature(const Eigen::Ref<const Eigen::VectorXd>& x1, const Eigen::Ref<const Eigen::VectorXd>& x2, double dt,
                            const std::vector<Eigen::VectorXd*>& intermediate_x, const UnaryFunction& fun,
                            Eigen::Ref<Eigen::VectorXd> integral_value) = 0;

    // integral of f(x(t)) over the interval [0, dt] and return intermediate points (preallocate the intermediate_points container!!!)
    virtual void quadrature(const Eigen::Ref<const Eigen::VectorXd>& x1, const Eigen::Ref<const Eigen::VectorXd>& x2, double dt,
                            const UnaryFunction& fun, Eigen::Ref<Eigen::VectorXd> integral_value, std::vector<Eigen::VectorXd*>& intermediate_x) = 0;

    // integral of f(x(t), u(t)) over the interval [0, dt]
    virtual void quadrature(const Eigen::Ref<const Eigen::VectorXd>& x1, const Eigen::Ref<const Eigen::VectorXd>& u1,
                            const Eigen::Ref<const Eigen::VectorXd>& x2, const Eigen::Ref<const Eigen::VectorXd>& u2, double dt,
                            const BinaryFunction& fun, Eigen::Ref<Eigen::VectorXd> integral_value) = 0;

    // integral of f(x(t), u(t)) over the interval [0, dt] with precomputed intermediate points
    virtual void quadrature(const Eigen::Ref<const Eigen::VectorXd>& x1, const Eigen::Ref<const Eigen::VectorXd>& u1,
                            const Eigen::Ref<const Eigen::VectorXd>& x2, const Eigen::Ref<const Eigen::VectorXd>& u2, double dt,
                            const std::vector<Eigen::VectorXd*>& intermediate_x, const std::vector<Eigen::VectorXd*>& intermediate_u,
                            const BinaryFunction& fun, Eigen::Ref<Eigen::VectorXd> integral_value) = 0;

    // integral of f(x(t), u(t)) over the interval [0, dt] and return intermediate points (preallocate the intermediate_points container!!!)
    virtual void quadrature(const Eigen::Ref<const Eigen::VectorXd>& x1, const Eigen::Ref<const Eigen::VectorXd>& u1,
                            const Eigen::Ref<const Eigen::VectorXd>& x2, const Eigen::Ref<const Eigen::VectorXd>& u2, double dt,
                            const BinaryFunction& fun, Eigen::Ref<Eigen::VectorXd> integral_value, std::vector<Eigen::VectorXd*>& intermediate_x,
                            std::vector<Eigen::VectorXd*>& intermediate_u) = 0;

    virtual bool computeIntermediateControls(const Eigen::Ref<const Eigen::VectorXd>& u1, const Eigen::Ref<const Eigen::VectorXd>& u2,
                                             std::vector<Eigen::VectorXd*>& intermediate_u) = 0;

    void computeEqualityConstraint(const StateVector& x1, const InputVector& u1, const StateVector& x2, double dt,
                                   const SystemDynamicsInterface& system, Eigen::Ref<Eigen::VectorXd> error) override
    {
        quadrature(x1, u1, x2, dt, system, error);
        error += x1 - x2;
    }

    void computeEqualityConstraint(const StateVector& x1, const InputVector& u1, const StateVector& x2, const InputVector& u2, double dt,
                                   const SystemDynamicsInterface& system, Eigen::Ref<Eigen::VectorXd> error)
    {
        quadrature(x1, u1, x2, u2, dt, system, error);
        error += x1 - x2;
    }

    void computeEqualityConstraint(const Eigen::Ref<const Eigen::VectorXd>& x1, const Eigen::Ref<const Eigen::VectorXd>& x2, double dt,
                                   const UnaryFunction& fun, Eigen::Ref<Eigen::VectorXd> error)
    {
        quadrature(x1, x2, dt, fun, error);
        error += x1 - x2;
    }

    bool interpolate(const Eigen::Ref<const Eigen::VectorXd>& x1, const Eigen::Ref<const Eigen::VectorXd>& u1,
                     const Eigen::Ref<const Eigen::VectorXd>& x2, const Eigen::Ref<const Eigen::VectorXd>& u2, double dt,
                     const SystemDynamicsInterface& system, const Range& range, std::vector<Eigen::VectorXd>& states,
                     std::vector<Eigen::VectorXd>& controls) override = 0;

#ifdef MESSAGE_SUPPORT
    //! Export integrator settings to message
    virtual void toMessage(messages::Quadrature& message) const {}
    //! Import integrator settings from message
    virtual void fromMessage(const messages::Quadrature& message, std::stringstream* issues = nullptr) {}

    //! Export selected class and parameters to a Integrator message
    void toMessage(corbo::messages::Collocation& message) const override { toMessage(*message.mutable_quadrature()); }
    //! Import selected class and parameters from a Integrator message (optionally pass any issues to the caller).
    void fromMessage(const corbo::messages::Collocation& message, std::stringstream* issues = nullptr) override
    {
        if (message.has_quadrature()) fromMessage(message.quadrature());
    }

//! Export selected class and parameters to a DynamicsEval message
// void toMessage(corbo::messages::DynamicsEval& message) const override { toMessage(*message.mutable_integrator()); }
//! Import selected class and parameters from a DynamicsEval message (optionally pass any issues to the caller).
// void fromMessage(const corbo::messages::DynamicsEval& message, std::stringstream* issues = nullptr) override
//{
//    if (message.has_integrator()) fromMessage(message.integrator());
//}
#endif
};

#define FACTORY_REGISTER_QUADRATURE_COLLOCATION(type) FACTORY_REGISTER_COLLOCATION(type)

}  // namespace corbo

#endif  // SRC_NUMERICS_INCLUDE_CORBO_NUMERICS_QUADRATURE_INTERFACE_H_
