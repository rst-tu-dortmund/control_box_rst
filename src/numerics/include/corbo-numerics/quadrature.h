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

#ifndef SRC_NUMERICS_INCLUDE_CORBO_NUMERICS_QUADRATURE_H_
#define SRC_NUMERICS_INCLUDE_CORBO_NUMERICS_QUADRATURE_H_

#include <corbo-numerics/quadrature_interface.h>

#include <cmath>
#include <memory>

namespace corbo {

/**
 * @brief Rectangle/midpoint rule (approximates function as constant between two points)
 *
 * @ingroup numerics
 *
 * \f[
 *      \int_0^{\Delta T} f(t) = \Delta T f(0.5 \Delta T)
 * \f].
 *
 * @see NumericalIntegratorInterface NumericalIntegratorImplicitInterface
 *
 * @author Christoph Rösmann (christoph.roesmann@tu-dortmund.de)
 */
class QuadratureRectangleRule : public QuadratureCollocationInterface
{
 public:
    // Implements interface method
    DynamicsEvalInterface::Ptr getInstance() const override { return std::make_shared<QuadratureRectangleRule>(); }

    // Implements interface method
    bool isSupportingCompressedStatesMode() const override { return true; }
    bool isIntermediateControlSubjectToOptim() const override { return false; }
    int getNumIntermediatePoints() const override { return 0; }

    int getNumIntermediateControls() const override { return 0; }
    int getNumIntermediateStates() const override { return 0; }

    // Implements interface method
    void quadrature(const std::vector<const Eigen::VectorXd*>& x1_points, const std::vector<const Eigen::VectorXd*>& u1_points,
                    const Eigen::VectorXd& u2, const Eigen::VectorXd& x2, double dt, const SystemDynamicsInterface& system,
                    Eigen::Ref<Eigen::VectorXd> integral_value) override
    {
        assert(x1_points.size() == 1);
        assert(u1_points.size() == 1);

        system.dynamics(0.5 * (*x1_points.front() + x2), 0.5 * (*u1_points.front() + u2), integral_value);
        integral_value *= dt;
    }
    
    void quadrature(const std::vector<const Eigen::VectorXd*>& x1_points, const std::vector<const Eigen::VectorXd*>& u1_points,
                            const Eigen::VectorXd& u2, const Eigen::VectorXd& x2, double dt,
                            const BinaryFunction& fun, Eigen::Ref<Eigen::VectorXd> integral_value) override
    {
        assert(x1_points.size() == 1);
        assert(u1_points.size() == 1);
        
        fun(0.5 * (*x1_points.front() + x2), 0.5 * (*u1_points.front() + u2), integral_value);
        integral_value *= dt;
    }

    bool interpolate(const std::vector<const Eigen::VectorXd*>& x1_points, const std::vector<const Eigen::VectorXd*>& u1_points,
                     const Eigen::VectorXd& u2, const Eigen::VectorXd& x2, double t1, double dt, double dt_interp,
                     const SystemDynamicsInterface& system, bool skip_u2_x2, TimeSeries& ts_x, TimeSeries& ts_u) override
    {
        assert(x1_points.size() == 1);
        assert(u1_points.size() == 1);

        Eigen::VectorXd xc = 0.5 * (*x1_points.front() + x2);
        Eigen::VectorXd uc = 0.5 * (*u1_points.front() + u2);

        if (dt_interp <= 0 || dt_interp >= dt)
        {
            ts_x.add(t1, xc);
            ts_u.add(t1, uc);
        }
        else
        {
            double t = t1;
            int n    = std::floor(dt / dt_interp);
            for (int i = 0; i < n; ++i)
            {
                ts_x.add(t, xc);
                ts_u.add(t, uc);
                t += dt_interp;
            }
        }

        if (!skip_u2_x2)
        {
            ts_x.add(t1 + dt, xc);
            ts_u.add(t1 + dt, uc);
        }
        return true;
    }

    // Implements interface method
    void quadrature(const StateVector& x1, const InputVector& u1, const StateVector& x2, double dt, const SystemDynamicsInterface& system,
                    Eigen::Ref<Eigen::VectorXd> integral_value) override
    {
        assert(integral_value.size() == x1.size());
        assert(dt > 0 && "dt must be greater then zero!");

        system.dynamics(0.5 * (x1 + x2), u1, integral_value);
        integral_value *= dt;
    }

    // Implements interface method
    void quadrature(const StateVector& x1, const InputVector& u1, const StateVector& x2, double dt, const SystemDynamicsInterface& system,
                    Eigen::Ref<Eigen::VectorXd> integral_value, std::vector<Eigen::VectorXd*>& intermediate_x) override
    {
        assert(integral_value.size() == x1.size());
        assert(dt > 0 && "dt must be greater then zero!");

        system.dynamics(0.5 * (x1 + x2), u1, integral_value);
        integral_value *= dt;
    }

    // Implements interface method
    void quadrature(const StateVector& x1, const InputVector& u1, const StateVector& x2, const InputVector& u2, double dt,
                    const SystemDynamicsInterface& system, Eigen::Ref<Eigen::VectorXd> integral_value) override
    {
        assert(integral_value.size() == x1.size());
        assert(dt > 0 && "dt must be greater then zero!");

        system.dynamics(0.5 * (x1 + x2), 0.5 * (u1 + u2), integral_value);
        integral_value *= dt;
    }

    void quadrature(const StateVector& x1, const InputVector& u1, const StateVector& x2, const InputVector& u2, double dt,
                    const SystemDynamicsInterface& system, const std::vector<Eigen::VectorXd*>& intermediate_x,
                    const std::vector<Eigen::VectorXd*>& intermediate_u, Eigen::Ref<Eigen::VectorXd> integral_value,
                    Eigen::Ref<Eigen::VectorXd> interm_x_error) override
    {
        system.dynamics(0.5 * (x1 + x2), 0.5 * (u1 + u2), integral_value);
        integral_value *= dt;
    }

    void quadrature(const StateVector& x1, const InputVector& u1, const StateVector& x2, const InputVector& u2, double dt,
                    const SystemDynamicsInterface& system, const std::vector<Eigen::VectorXd*>& intermediate_u,
                    Eigen::Ref<Eigen::VectorXd> integral_value, std::vector<Eigen::VectorXd*>& intermediate_x) override
    {
        assert(integral_value.size() == x1.size());
        assert(dt > 0 && "dt must be greater then zero!");

        system.dynamics(0.5 * (x1 + x2), 0.5 * (u1 + u2), integral_value);
        integral_value *= dt;
    }

    // Implements interface method
    void quadrature(const StateVector& x1, const InputVector& u1, const StateVector& x2, const InputVector& u2, double dt,
                    const SystemDynamicsInterface& system, Eigen::Ref<Eigen::VectorXd> integral_value, std::vector<Eigen::VectorXd*>& intermediate_x,
                    std::vector<Eigen::VectorXd*>& intermediate_u) override
    {
        assert(integral_value.size() == x1.size());
        assert(dt > 0 && "dt must be greater then zero!");

        system.dynamics(0.5 * (x1 + x2), 0.5 * (u1 + u2), integral_value);
        integral_value *= dt;
    }

    // Implements interface method
    void quadrature(const Eigen::Ref<const Eigen::VectorXd>& x1, const Eigen::Ref<const Eigen::VectorXd>& x2, double dt, const UnaryFunction& fun,
                    Eigen::Ref<Eigen::VectorXd> integral_value) override
    {
        assert(dt > 0 && "dt must be greater then zero!");

        fun(0.5 * (x1 + x2), integral_value);
        integral_value *= dt;
    }

    // Implements interface method
    void quadrature(const Eigen::Ref<const Eigen::VectorXd>& x1, const Eigen::Ref<const Eigen::VectorXd>& x2, double dt,
                    const std::vector<Eigen::VectorXd*>& intermediate_x, const UnaryFunction& fun,
                    Eigen::Ref<Eigen::VectorXd> integral_value) override
    {
        fun(0.5 * (x1 + x2), integral_value);
        integral_value *= dt;
    }

    // Implements interface method
    void quadrature(const Eigen::Ref<const Eigen::VectorXd>& x1, const Eigen::Ref<const Eigen::VectorXd>& x2, double dt, const UnaryFunction& fun,
                    Eigen::Ref<Eigen::VectorXd> integral_value, std::vector<Eigen::VectorXd*>& intermediate_x) override
    {
        assert(intermediate_x.size() == 1);
        assert(intermediate_x[0]->size() == x1.size());
        intermediate_x[0]->noalias() = 0.5 * (x1 + x2);
        fun(*intermediate_x[0], integral_value);
        integral_value *= dt;
    }

    // Implements interface method
    void quadrature(const Eigen::Ref<const Eigen::VectorXd>& x1, const Eigen::Ref<const Eigen::VectorXd>& u1,
                    const Eigen::Ref<const Eigen::VectorXd>& x2, const Eigen::Ref<const Eigen::VectorXd>& u2, double dt, const BinaryFunction& fun,
                    Eigen::Ref<Eigen::VectorXd> integral_value) override
    {
        assert(dt > 0 && "dt must be greater then zero!");

        fun(0.5 * (x1 + x2), 0.5 * (u1 + u2), integral_value);
        integral_value *= dt;
    }

    // Implements interface method
    void quadrature(const Eigen::Ref<const Eigen::VectorXd>& x1, const Eigen::Ref<const Eigen::VectorXd>& u1,
                    const Eigen::Ref<const Eigen::VectorXd>& x2, const Eigen::Ref<const Eigen::VectorXd>& u2, double dt,
                    const std::vector<Eigen::VectorXd*>& intermediate_x, const std::vector<Eigen::VectorXd*>& intermediate_u,
                    const BinaryFunction& fun, Eigen::Ref<Eigen::VectorXd> integral_value) override
    {
        fun(0.5 * (x1 + x2), 0.5 * (u1 + u2), integral_value);
        integral_value *= dt;
    }

    // Implements interface method
    void quadrature(const Eigen::Ref<const Eigen::VectorXd>& x1, const Eigen::Ref<const Eigen::VectorXd>& u1,
                    const Eigen::Ref<const Eigen::VectorXd>& x2, const Eigen::Ref<const Eigen::VectorXd>& u2, double dt, const BinaryFunction& fun,
                    Eigen::Ref<Eigen::VectorXd> integral_value, std::vector<Eigen::VectorXd*>& intermediate_x,
                    std::vector<Eigen::VectorXd*>& intermediate_u) override
    {
        fun(0.5 * (x1 + x2), 0.5 * (u1 + u2), integral_value);
        integral_value *= dt;
    }

    bool interpolate(const Eigen::Ref<const Eigen::VectorXd>& x1, const Eigen::Ref<const Eigen::VectorXd>& u1,
                     const Eigen::Ref<const Eigen::VectorXd>& x2, const Eigen::Ref<const Eigen::VectorXd>& u2, double dt,
                     const SystemDynamicsInterface& /*system*/, const Range& range, std::vector<Eigen::VectorXd>& states,
                     std::vector<Eigen::VectorXd>& controls) override
    {
        if (range.getStart() < 0 || range.getEnd() > dt)
        {
            PRINT_ERROR_NAMED("range must be in the interval [0, dt]");
            return false;
        }

        Eigen::VectorXd xcenter = 0.5 * (x1 + x2);
        Eigen::VectorXd ucenter = 0.5 * (u1 + u2);

        if (dt < 1e-15)
        {
            states.push_back(x1);
            controls.push_back(ucenter);
            return true;
        }

        int n = range.getNumInRange();

        double t = range.getStart();
        for (int i = 0; i < n; ++i, t += range.getStep())
        {
            double tau = (t - range.getStart()) / dt;
            states.push_back(x1 + tau * (x2 - x1));
            controls.push_back(ucenter);
        }
        if (range.includeEnd())
        {
            states.push_back(x2);
            controls.push_back(ucenter);
        }
        return true;
    }

    bool computeIntermediateControls(const Eigen::Ref<const Eigen::VectorXd>& u1, const Eigen::Ref<const Eigen::VectorXd>& u2,
                                     std::vector<Eigen::VectorXd*>& intermediate_u) override
    {
        return true;
    }

#ifdef MESSAGE_SUPPORT
    // Implements interface method
    void toMessage(messages::Quadrature& message) const override
    {
        QuadratureCollocationInterface::toMessage(message);
        message.mutable_rectangle_rule();
    }
#endif
};

FACTORY_REGISTER_QUADRATURE_COLLOCATION(QuadratureRectangleRule)

/**
 * @brief Trapezoidal rule (approximates function as linear function between two points)
 *
 * @ingroup numerics
 *
 * \f[
 *      \int_0^{\Delta T} f(t) = \Delta T 0.5 * (f(0) + f(\Delta T))
 * \f].
 *
 * @see NumericalIntegratorInterface NumericalIntegratorImplicitInterface
 *
 * @author Christoph Rösmann (christoph.roesmann@tu-dortmund.de)
 */
class QuadratureTrapezoidalRule : public QuadratureCollocationInterface
{
 public:
    // Implements interface method
    DynamicsEvalInterface::Ptr getInstance() const override { return std::make_shared<QuadratureTrapezoidalRule>(); }

    // Implements interface method
    bool isSupportingCompressedStatesMode() const override { return true; }
    bool isIntermediateControlSubjectToOptim() const override { return false; }
    int getNumIntermediatePoints() const override { return 0; }

    int getNumIntermediateControls() const override { return 0; }
    int getNumIntermediateStates() const override { return 0; }

    // Implements interface method
    void quadrature(const std::vector<const Eigen::VectorXd*>& x1_points, const std::vector<const Eigen::VectorXd*>& u1_points,
                    const Eigen::VectorXd& u2, const Eigen::VectorXd& x2, double dt, const SystemDynamicsInterface& system,
                    Eigen::Ref<Eigen::VectorXd> integral_value) override
    {
        assert(x1_points.size() == 1);
        assert(u1_points.size() == 1);

        Eigen::VectorXd f1(x2.size());
        system.dynamics(*x1_points.front(), *u1_points.front(), f1);
        system.dynamics(x2, u2, integral_value);
        integral_value = dt * 0.5 * (f1 + integral_value);
    }

    void quadrature(const std::vector<const Eigen::VectorXd*>& x1_points, const std::vector<const Eigen::VectorXd*>& u1_points,
                    const Eigen::VectorXd& u2, const Eigen::VectorXd& x2, double dt,
                    const BinaryFunction& fun, Eigen::Ref<Eigen::VectorXd> integral_value) override
    {
        assert(x1_points.size() == 1);
        assert(u1_points.size() == 1);
        
        Eigen::VectorXd f1(x2.size());
        fun(*x1_points.front(), *u1_points.front(), f1);
        fun(x2, u2, integral_value);
        integral_value = dt * 0.5 * (f1 + integral_value);
    }
    
    bool interpolate(const std::vector<const Eigen::VectorXd*>& x1_points, const std::vector<const Eigen::VectorXd*>& u1_points,
                     const Eigen::VectorXd& u2, const Eigen::VectorXd& x2, double t1, double dt, double dt_interp,
                     const SystemDynamicsInterface& system, bool skip_u2_x2, TimeSeries& ts_x, TimeSeries& ts_u) override
    {
        assert(x1_points.size() == 1);
        assert(u1_points.size() == 1);

        if (dt_interp <= 0 || dt_interp >= dt)
        {
            ts_x.add(t1, *x1_points[0]);
            ts_u.add(t1, *u1_points[0]);
        }
        else if (dt > 0)
        {
            Eigen::VectorXd xvel = (x2 - *x1_points[0]) / dt;
            Eigen::VectorXd uvel = (u2 - *u1_points[0]) / dt;
            int n                = std::floor(dt / dt_interp);
            int t = t1;
            for (int i = 0; i < n; ++i)
            {
                ts_x.add(t, *x1_points[0] + (double)i * dt_interp * xvel);
                ts_u.add(t, *u1_points[0] + (double)i * dt_interp * uvel);
                t += dt_interp;
            }
        }

        if (!skip_u2_x2)
        {
            ts_x.add(t1 + dt, x2);
            ts_u.add(t1 + dt, u2);
        }
        return true;
    }

    // Implements interface method
    void quadrature(const StateVector& x1, const InputVector& u1, const StateVector& x2, double dt, const SystemDynamicsInterface& system,
                    Eigen::Ref<Eigen::VectorXd> integral_value) override
    {
        assert(integral_value.size() == x1.size());
        assert(dt > 0 && "dt must be greater then zero!");

        Eigen::VectorXd f1(x1.size());
        system.dynamics(x1, u1, f1);
        system.dynamics(x2, u1, integral_value);
        integral_value = dt * 0.5 * (f1 + integral_value);
    }

    // Implements interface method
    void quadrature(const StateVector& x1, const InputVector& u1, const StateVector& x2, double dt, const SystemDynamicsInterface& system,
                    Eigen::Ref<Eigen::VectorXd> integral_value, std::vector<Eigen::VectorXd*>& /*intermediate_x*/) override
    {
        assert(integral_value.size() == x1.size());
        assert(dt > 0 && "dt must be greater then zero!");

        Eigen::VectorXd f1(x1.size());
        system.dynamics(x1, u1, f1);
        system.dynamics(x2, u1, integral_value);
        integral_value = dt * 0.5 * (f1 + integral_value);
    }

    // Implements interface method
    void quadrature(const StateVector& x1, const InputVector& u1, const StateVector& x2, const InputVector& u2, double dt,
                    const SystemDynamicsInterface& system, Eigen::Ref<Eigen::VectorXd> integral_value) override
    {
        assert(integral_value.size() == x1.size());
        assert(dt > 0 && "dt must be greater then zero!");

        Eigen::VectorXd f1(x1.size());
        system.dynamics(x1, u1, f1);
        system.dynamics(x2, u2, integral_value);
        integral_value = dt * 0.5 * (f1 + integral_value);
    }

    void quadrature(const StateVector& x1, const InputVector& u1, const StateVector& x2, const InputVector& u2, double dt,
                    const SystemDynamicsInterface& system, const std::vector<Eigen::VectorXd*>& /*intermediate_x*/,
                    const std::vector<Eigen::VectorXd*>& /*intermediate_u*/, Eigen::Ref<Eigen::VectorXd> integral_value,
                    Eigen::Ref<Eigen::VectorXd> /*interm_x_error*/) override
    {
        assert(dt > 0 && "dt must be greater then zero!");

        Eigen::VectorXd f1(integral_value.size());
        system.dynamics(x1, u1, f1);
        system.dynamics(x2, u2, integral_value);
        integral_value = dt * 0.5 * (f1 + integral_value);
    }

    void quadrature(const StateVector& x1, const InputVector& u1, const StateVector& x2, const InputVector& u2, double dt,
                    const SystemDynamicsInterface& system, const std::vector<Eigen::VectorXd*>& /*intermediate_u*/,
                    Eigen::Ref<Eigen::VectorXd> integral_value, std::vector<Eigen::VectorXd*>& /*intermediate_x*/) override
    {
        assert(integral_value.size() == x1.size());
        assert(dt > 0 && "dt must be greater then zero!");

        Eigen::VectorXd f1(x1.size());
        system.dynamics(x1, u1, f1);
        system.dynamics(x2, u2, integral_value);
        integral_value = dt * 0.5 * (f1 + integral_value);
    }

    // Implements interface method
    void quadrature(const StateVector& x1, const InputVector& u1, const StateVector& x2, const InputVector& u2, double dt,
                    const SystemDynamicsInterface& system, Eigen::Ref<Eigen::VectorXd> integral_value,
                    std::vector<Eigen::VectorXd*>& /*intermediate_x*/, std::vector<Eigen::VectorXd*>& /*intermediate_u*/) override
    {
        assert(integral_value.size() == x1.size());
        assert(dt > 0 && "dt must be greater then zero!");

        Eigen::VectorXd f1(x1.size());
        system.dynamics(x1, u1, f1);
        system.dynamics(x2, u2, integral_value);
        integral_value = dt * 0.5 * (f1 + integral_value);
    }
    // Implements interface method
    void quadrature(const Eigen::Ref<const Eigen::VectorXd>& x1, const Eigen::Ref<const Eigen::VectorXd>& x2, double dt, const UnaryFunction& fun,
                    Eigen::Ref<Eigen::VectorXd> integral_value) override
    {
        assert(dt > 0 && "dt must be greater then zero!");

        Eigen::VectorXd f1(integral_value.size());
        fun(x1, f1);
        fun(x2, integral_value);
        integral_value = dt * 0.5 * (f1 + integral_value);
    }

    // Implements interface method
    void quadrature(const Eigen::Ref<const Eigen::VectorXd>& x1, const Eigen::Ref<const Eigen::VectorXd>& x2, double dt,
                    const std::vector<Eigen::VectorXd*>& /*intermediate_x*/, const UnaryFunction& fun,
                    Eigen::Ref<Eigen::VectorXd> integral_value) override
    {
        assert(dt > 0 && "dt must be greater then zero!");

        Eigen::VectorXd f1(integral_value.size());
        fun(x1, f1);
        fun(x2, integral_value);
        integral_value = dt * 0.5 * (f1 + integral_value);
    }

    // Implements interface method
    void quadrature(const Eigen::Ref<const Eigen::VectorXd>& x1, const Eigen::Ref<const Eigen::VectorXd>& x2, double dt, const UnaryFunction& fun,
                    Eigen::Ref<Eigen::VectorXd> integral_value, std::vector<Eigen::VectorXd*>& /*intermediate_x*/) override
    {
        assert(dt > 0 && "dt must be greater then zero!");

        Eigen::VectorXd f1(integral_value.size());
        fun(x1, f1);
        fun(x2, integral_value);
        integral_value = dt * 0.5 * (f1 + integral_value);
    }

    // Implements interface method
    void quadrature(const Eigen::Ref<const Eigen::VectorXd>& x1, const Eigen::Ref<const Eigen::VectorXd>& u1,
                    const Eigen::Ref<const Eigen::VectorXd>& x2, const Eigen::Ref<const Eigen::VectorXd>& u2, double dt, const BinaryFunction& fun,
                    Eigen::Ref<Eigen::VectorXd> integral_value) override
    {
        assert(dt > 0 && "dt must be greater then zero!");

        Eigen::VectorXd f1(integral_value.size());
        fun(x1, u1, f1);
        fun(x2, u2, integral_value);
        integral_value = dt * 0.5 * (f1 + integral_value);
    }

    // Implements interface method
    void quadrature(const Eigen::Ref<const Eigen::VectorXd>& x1, const Eigen::Ref<const Eigen::VectorXd>& u1,
                    const Eigen::Ref<const Eigen::VectorXd>& x2, const Eigen::Ref<const Eigen::VectorXd>& u2, double dt,
                    const std::vector<Eigen::VectorXd*>& /*intermediate_x*/, const std::vector<Eigen::VectorXd*>& /*intermediate_u*/,
                    const BinaryFunction& fun, Eigen::Ref<Eigen::VectorXd> integral_value) override
    {
        assert(dt > 0 && "dt must be greater then zero!");

        Eigen::VectorXd f1(integral_value.size());
        fun(x1, u1, f1);
        fun(x2, u2, integral_value);
        integral_value = dt * 0.5 * (f1 + integral_value);
    }

    // Implements interface method
    void quadrature(const Eigen::Ref<const Eigen::VectorXd>& x1, const Eigen::Ref<const Eigen::VectorXd>& u1,
                    const Eigen::Ref<const Eigen::VectorXd>& x2, const Eigen::Ref<const Eigen::VectorXd>& u2, double dt, const BinaryFunction& fun,
                    Eigen::Ref<Eigen::VectorXd> integral_value, std::vector<Eigen::VectorXd*>& /*intermediate_x*/,
                    std::vector<Eigen::VectorXd*>& /*intermediate_u*/) override
    {
        assert(dt > 0 && "dt must be greater then zero!");

        Eigen::VectorXd f1(integral_value.size());
        fun(x1, u1, f1);
        fun(x2, u2, integral_value);
        integral_value = dt * 0.5 * (f1 + integral_value);
    }

    bool interpolate(const Eigen::Ref<const Eigen::VectorXd>& x1, const Eigen::Ref<const Eigen::VectorXd>& u1,
                     const Eigen::Ref<const Eigen::VectorXd>& x2, const Eigen::Ref<const Eigen::VectorXd>& u2, double dt,
                     const SystemDynamicsInterface& system, const Range& range, std::vector<Eigen::VectorXd>& states,
                     std::vector<Eigen::VectorXd>& controls) override
    {
        if (range.getStart() < 0 || range.getEnd() > dt)
        {
            PRINT_ERROR_NAMED("range must be in the interval [0, dt]");
            return false;
        }

        states.push_back(x1);
        controls.push_back(u1);

        if (dt < 1e-15) return true;

        Eigen::VectorXd f1(system.getStateDimension());
        Eigen::VectorXd f2(system.getStateDimension());
        system.dynamics(x1, u1, f1);
        system.dynamics(x2, u2, f2);

        int n    = range.getNumInRange();
        double t = range.getStart() + range.getStep();  // we already added x1 and u1
        for (int i = 1; i < n; ++i, t += range.getStep())
        {
            double tau = t - range.getStart();
            states.push_back(x1 + tau * f1 + tau * tau / (2.0 * dt) * (f2 - f1));  // quadratic state trajectoy
            controls.push_back(u1 + tau / dt * (u2 - u1));                         // linear controls
        }
        if (range.includeEnd())
        {
            states.push_back(x2);
            controls.push_back(u2);
        }
        return true;
    }

    bool computeIntermediateControls(const Eigen::Ref<const Eigen::VectorXd>& u1, const Eigen::Ref<const Eigen::VectorXd>& u2,
                                     std::vector<Eigen::VectorXd*>& intermediate_u) override
    {
        return true;
    }

#ifdef MESSAGE_SUPPORT
    // Implements interface method
    void toMessage(messages::Quadrature& message) const override
    {
        QuadratureCollocationInterface::toMessage(message);
        message.mutable_trapezoidal_rule();
    }
#endif
};

FACTORY_REGISTER_QUADRATURE_COLLOCATION(QuadratureTrapezoidalRule)

/**
 * @brief Hermite-Simpson rule (approximates function as quadratic polynomial between two points)
 *
 * @ingroup numerics
 *
 * \f[
 *      \int_0^{\Delta T} f(t) = \frac{\Delta T}{6} (f(0) + 4 f(0.5 \Delta T) + f(\Delta T)
 * \f].
 *
 * Hermite-Simpson collocation approximates the system dynamics by
 * quadratic polynomials between two consecutive points.
 * This version approximations the control trajectory by linear segments.
 * Refer to the QuadratureHermiteSimpson for a quadratic control interpolation.
 * Refer to
 * - https://mec560sbu.github.io/2016/09/30/direct_collocation/
 * - http://epubs.siam.org/doi/pdf/10.1137/16M1062569
 *
 * @see NumericalIntegratorInterface NumericalIntegratorImplicitInterface
 *
 * @author Christoph Rösmann (christoph.roesmann@tu-dortmund.de)
 */
class QuadratureHermiteSimpsonLinearControl : public QuadratureCollocationInterface
{
 public:
    // Implements interface method
    DynamicsEvalInterface::Ptr getInstance() const override { return std::make_shared<QuadratureHermiteSimpsonLinearControl>(); }

    // Implements interface method
    bool isSupportingCompressedStatesMode() const override { return true; }
    bool isIntermediateControlSubjectToOptim() const override { return false; }
    int getNumIntermediatePoints() const override { return 1; }

    int getNumIntermediateControls() const override { return 0; }
    int getNumIntermediateStates() const override { return 0; }

    // Implements interface method
    void quadrature(const std::vector<const Eigen::VectorXd*>& x1_points, const std::vector<const Eigen::VectorXd*>& u1_points,
                    const Eigen::VectorXd& u2, const Eigen::VectorXd& x2, double dt, const SystemDynamicsInterface& system,
                    Eigen::Ref<Eigen::VectorXd> integral_value) override
    {
        assert(x1_points.size() == 1);
        assert(u1_points.size() == 1);

        Eigen::VectorXd f1(x2.size());
        system.dynamics(*x1_points.front(), *u1_points.front(), f1);

        Eigen::VectorXd f2(x2.size());
        system.dynamics(x2, u2, f2);

        Eigen::VectorXd x_center = 0.5 * (*x1_points.front() + x2) + dt / 8.0 * (f1 - f2);

        system.dynamics(x_center, 0.5 * (*u1_points.front() + u2), integral_value);

        integral_value = dt / 6.0 * (f1 + 4.0 * integral_value + f2);
    }
    
    void quadrature(const std::vector<const Eigen::VectorXd*>& x1_points, const std::vector<const Eigen::VectorXd*>& u1_points,
                    const Eigen::VectorXd& u2, const Eigen::VectorXd& x2, double dt,
                    const BinaryFunction& fun, Eigen::Ref<Eigen::VectorXd> integral_value) override
    {
        assert(x1_points.size() == 1);
        assert(u1_points.size() == 1);
        
        Eigen::VectorXd f1(x2.size());
        fun(*x1_points.front(), *u1_points.front(), f1);
        
        Eigen::VectorXd f2(x2.size());
        fun(x2, u2, f2);
        
        Eigen::VectorXd x_center = 0.5 * (*x1_points.front() + x2) + dt / 8.0 * (f1 - f2);
        
        fun(x_center, 0.5 * (*u1_points.front() + u2), integral_value);
        
        integral_value = dt / 6.0 * (f1 + 4.0 * integral_value + f2);
    }

    bool interpolate(const std::vector<const Eigen::VectorXd*>& x1_points, const std::vector<const Eigen::VectorXd*>& u1_points,
                     const Eigen::VectorXd& u2, const Eigen::VectorXd& x2, double t1, double dt, double dt_interp,
                     const SystemDynamicsInterface& system, bool skip_u2_x2, TimeSeries& ts_x, TimeSeries& ts_u) override
    {
        assert(x1_points.size() == 1);
        assert(u1_points.size() == 1);

        if (dt_interp <= 0 || dt_interp >= dt)
        {
            ts_x.add(t1, *x1_points[0]);
            ts_u.add(t1, *u1_points[0]);
        }
        else if (dt > 0)
        {
            Eigen::VectorXd f1(x2.size());
            system.dynamics(*x1_points.front(), *u1_points.front(), f1);

            Eigen::VectorXd f2(x2.size());
            system.dynamics(x2, u2, f2);

            Eigen::VectorXd fc(x2.size());
            Eigen::VectorXd xc = 0.5 * (*x1_points.front() + x2) + dt / 8.0 * (f1 - f2);
            system.dynamics(xc, 0.5 * (*u1_points.front() + u2), fc);

            Eigen::VectorXd uvel = (u2 - *u1_points[0]) / dt;
            int n                = std::floor(dt / dt_interp);
            double t             = t1;
            for (int i = 0; i < n; ++i)
            {
                double tau = (t - t1);

                // cubic interpolation for states (since quadratic interplation for xdot)
                Eigen::VectorXd x_interp = *x1_points.front() + f1 * tau + 0.5 / dt * tau * tau * (-3.0 * f1 + 4.0 * fc - f2) +
                                           1.0 / 3.0 / (dt * dt) * std::pow(tau, 3) * (2.0 * f1 - 4.0 * fc + 2.0 * f2);

                ts_x.add(t1 + tau, x_interp);
                ts_u.add(t1 + tau, *u1_points[0] + (double)i * dt_interp * uvel);

                t += dt_interp;
            }
        }

        if (!skip_u2_x2)
        {
            ts_x.add(t1 + dt, x2);
            ts_u.add(t1 + dt, u2);
        }
        return true;
    }

    void quadrature(const StateVector& x1, const InputVector& u1, const StateVector& x2, double dt, const SystemDynamicsInterface& system,
                    Eigen::Ref<Eigen::VectorXd> integral_value) override
    {
        assert(integral_value.size() == x1.size());
        assert(dt > 0 && "dt must be greater then zero!");

        Eigen::VectorXd f1(x1.size());
        system.dynamics(x1, u1, f1);

        Eigen::VectorXd f2(x1.size());
        system.dynamics(x2, u1, f2);

        Eigen::VectorXd x_center = 0.5 * (x1 + x2) + dt / 8.0 * (f1 - f2);

        system.dynamics(x_center, u1, integral_value);

        integral_value = dt / 6.0 * (f1 + 4.0 * integral_value + f2);
    }

    void quadrature(const StateVector& x1, const InputVector& u1, const StateVector& x2, double dt, const SystemDynamicsInterface& system,
                    Eigen::Ref<Eigen::VectorXd> integral_value, std::vector<Eigen::VectorXd*>& intermediate_x) override
    {
        assert(integral_value.size() == x1.size());
        assert(dt > 0 && "dt must be greater then zero!");
        assert(intermediate_x.size() == 1);
        assert(intermediate_x[0]->size() == x1.size());

        Eigen::VectorXd f1(x1.size());
        system.dynamics(x1, u1, f1);

        Eigen::VectorXd f2(x1.size());
        system.dynamics(x2, u1, f2);

        intermediate_x[0]->noalias() = 0.5 * (x1 + x2) + dt / 8.0 * (f1 - f2);

        system.dynamics(*intermediate_x[0], u1, integral_value);

        integral_value = dt / 6.0 * (f1 + 4.0 * integral_value + f2);
    }

    // Implements interface method
    void quadrature(const StateVector& x1, const InputVector& u1, const StateVector& x2, const InputVector& u2, double dt,
                    const SystemDynamicsInterface& system, Eigen::Ref<Eigen::VectorXd> integral_value) override
    {
        assert(integral_value.size() == x1.size());
        assert(dt > 0 && "dt must be greater then zero!");

        Eigen::VectorXd f1(x1.size());
        system.dynamics(x1, u1, f1);

        Eigen::VectorXd f2(x1.size());
        system.dynamics(x2, u2, f2);

        Eigen::VectorXd x_center = 0.5 * (x1 + x2) + dt / 8.0 * (f1 - f2);

        system.dynamics(x_center, 0.5 * (u1 + u2), integral_value);

        integral_value = dt / 6.0 * (f1 + 4.0 * integral_value + f2);
    }

    void quadrature(const StateVector& x1, const InputVector& u1, const StateVector& x2, const InputVector& u2, double dt,
                    const SystemDynamicsInterface& system, const std::vector<Eigen::VectorXd*>& intermediate_x,
                    const std::vector<Eigen::VectorXd*>& intermediate_u, Eigen::Ref<Eigen::VectorXd> integral_value,
                    Eigen::Ref<Eigen::VectorXd> interm_x_error) override
    {
        assert(intermediate_x.size() == 1);
        assert(intermediate_x[0]->size() == x1.size());
        assert(intermediate_u.size() == 1);
        assert(intermediate_u[0]->size() == u1.size());
        assert(dt > 0 && "dt must be greater then zero!");

        Eigen::VectorXd f1(integral_value.size());
        system.dynamics(x1, u1, f1);

        Eigen::VectorXd f2(integral_value.size());
        system.dynamics(x2, u2, f2);

        system.dynamics(*intermediate_x[0], *intermediate_u[0], integral_value);

        integral_value = dt / 6.0 * (f1 + 4.0 * integral_value + f2);

        assert(interm_x_error.size() == x1.size());
        interm_x_error = *intermediate_x[0] - (0.5 * (x1 + x2) + dt / 8.0 * (f1 - f2));
    }

    void quadrature(const StateVector& x1, const InputVector& u1, const StateVector& x2, const InputVector& u2, double dt,
                    const SystemDynamicsInterface& system, const std::vector<Eigen::VectorXd*>& intermediate_u,
                    Eigen::Ref<Eigen::VectorXd> integral_value, std::vector<Eigen::VectorXd*>& intermediate_x) override
    {
        assert(integral_value.size() == x1.size());
        assert(dt > 0 && "dt must be greater then zero!");
        assert(intermediate_x.size() == 1);
        assert(intermediate_x[0]->size() == x1.size());
        assert(intermediate_u.size() == 1);
        assert(intermediate_u[0]->size() == u1.size());

        Eigen::VectorXd f1(x1.size());
        system.dynamics(x1, u1, f1);

        Eigen::VectorXd f2(x1.size());
        system.dynamics(x2, u2, f2);

        intermediate_x[0]->noalias() = 0.5 * (x1 + x2) + dt / 8.0 * (f1 - f2);

        system.dynamics(*intermediate_x[0], *intermediate_u[0], integral_value);

        integral_value = dt / 6.0 * (f1 + 4.0 * integral_value + f2);
    }

    void quadrature(const StateVector& x1, const InputVector& u1, const StateVector& x2, const InputVector& u2, double dt,
                    const SystemDynamicsInterface& system, Eigen::Ref<Eigen::VectorXd> integral_value, std::vector<Eigen::VectorXd*>& intermediate_x,
                    std::vector<Eigen::VectorXd*>& intermediate_u) override
    {
        assert(integral_value.size() == x1.size());
        assert(dt > 0 && "dt must be greater then zero!");
        assert(intermediate_x.size() == 1);
        assert(intermediate_x[0]->size() == x1.size());
        assert(intermediate_u.size() == 1);
        assert(intermediate_u[0]->size() == u1.size());

        Eigen::VectorXd f1(x1.size());
        system.dynamics(x1, u1, f1);

        Eigen::VectorXd f2(x1.size());
        system.dynamics(x2, u2, f2);

        intermediate_x[0]->noalias() = 0.5 * (x1 + x2) + dt / 8.0 * (f1 - f2);
        intermediate_u[0]->noalias() = 0.5 * (u1 + u2);

        system.dynamics(*intermediate_x[0], *intermediate_u[0], integral_value);

        integral_value = dt / 6.0 * (f1 + 4.0 * integral_value + f2);
    }

    // Implements interface method
    void quadrature(const Eigen::Ref<const Eigen::VectorXd>& x1, const Eigen::Ref<const Eigen::VectorXd>& x2, double dt, const UnaryFunction& fun,
                    Eigen::Ref<Eigen::VectorXd> integral_value) override
    {
        assert(dt > 0 && "dt must be greater then zero!");

        Eigen::VectorXd f1(integral_value.size());
        fun(x1, f1);

        Eigen::VectorXd f2(integral_value.size());
        fun(x2, f2);

        Eigen::VectorXd x_center = 0.5 * (x1 + x2) + dt / 8.0 * (f1 - f2);

        fun(x_center, integral_value);

        integral_value = dt / 6.0 * (f1 + 4.0 * integral_value + f2);
    }

    // Implements interface method
    void quadrature(const Eigen::Ref<const Eigen::VectorXd>& x1, const Eigen::Ref<const Eigen::VectorXd>& x2, double dt,
                    const std::vector<Eigen::VectorXd*>& intermediate_x, const UnaryFunction& fun,
                    Eigen::Ref<Eigen::VectorXd> integral_value) override
    {
        assert(intermediate_x.size() == 1);
        assert(intermediate_x[0]->size() == x1.size());
        assert(dt > 0 && "dt must be greater then zero!");

        Eigen::VectorXd f1(integral_value.size());
        fun(x1, f1);

        Eigen::VectorXd f2(integral_value.size());
        fun(x2, f2);

        fun(*intermediate_x[0], integral_value);

        integral_value = dt / 6.0 * (f1 + 4.0 * integral_value + f2);
    }

    // Implements interface method
    void quadrature(const Eigen::Ref<const Eigen::VectorXd>& x1, const Eigen::Ref<const Eigen::VectorXd>& x2, double dt, const UnaryFunction& fun,
                    Eigen::Ref<Eigen::VectorXd> integral_value, std::vector<Eigen::VectorXd*>& intermediate_x) override
    {
        assert(intermediate_x.size() == 1);
        assert(intermediate_x[0]->size() == x1.size());
        assert(dt > 0 && "dt must be greater then zero!");

        Eigen::VectorXd f1(integral_value.size());
        fun(x1, f1);

        Eigen::VectorXd f2(integral_value.size());
        fun(x2, f2);

        *intermediate_x[0] = 0.5 * (x1 + x2) + dt / 8.0 * (f1 - f2);

        fun(*intermediate_x[0], integral_value);

        integral_value = dt / 6.0 * (f1 + 4.0 * integral_value + f2);
    }

    // Implements interface method
    void quadrature(const Eigen::Ref<const Eigen::VectorXd>& x1, const Eigen::Ref<const Eigen::VectorXd>& u1,
                    const Eigen::Ref<const Eigen::VectorXd>& x2, const Eigen::Ref<const Eigen::VectorXd>& u2, double dt, const BinaryFunction& fun,
                    Eigen::Ref<Eigen::VectorXd> integral_value) override
    {
        assert(dt > 0 && "dt must be greater then zero!");

        Eigen::VectorXd f1(integral_value.size());
        fun(x1, u1, f1);

        Eigen::VectorXd f2(integral_value.size());
        fun(x2, u2, f2);

        Eigen::VectorXd x_center = 0.5 * (x1 + x2) + dt / 8.0 * (f1 - f2);

        fun(x_center, 0.5 * (u1 + u2), integral_value);

        integral_value = dt / 6.0 * (f1 + 4.0 * integral_value + f2);
    }

    // Implements interface method
    void quadrature(const Eigen::Ref<const Eigen::VectorXd>& x1, const Eigen::Ref<const Eigen::VectorXd>& u1,
                    const Eigen::Ref<const Eigen::VectorXd>& x2, const Eigen::Ref<const Eigen::VectorXd>& u2, double dt,
                    const std::vector<Eigen::VectorXd*>& intermediate_x, const std::vector<Eigen::VectorXd*>& intermediate_u,
                    const BinaryFunction& fun, Eigen::Ref<Eigen::VectorXd> integral_value) override
    {
        assert(intermediate_x.size() == 1);
        assert(intermediate_x[0]->size() == x1.size());
        assert(intermediate_u.size() == 1);
        assert(intermediate_u[0]->size() == u1.size());
        assert(dt > 0 && "dt must be greater then zero!");

        Eigen::VectorXd f1(integral_value.size());
        fun(x1, u1, f1);

        Eigen::VectorXd f2(integral_value.size());
        fun(x2, u2, f2);

        fun(*intermediate_x[0], *intermediate_u[0], integral_value);

        integral_value = dt / 6.0 * (f1 + 4.0 * integral_value + f2);
    }

    // Implements interface method
    void quadrature(const Eigen::Ref<const Eigen::VectorXd>& x1, const Eigen::Ref<const Eigen::VectorXd>& u1,
                    const Eigen::Ref<const Eigen::VectorXd>& x2, const Eigen::Ref<const Eigen::VectorXd>& u2, double dt, const BinaryFunction& fun,
                    Eigen::Ref<Eigen::VectorXd> integral_value, std::vector<Eigen::VectorXd*>& intermediate_x,
                    std::vector<Eigen::VectorXd*>& intermediate_u) override
    {
        assert(intermediate_x.size() == 1);
        assert(intermediate_x[0]->size() == x1.size());
        assert(intermediate_u.size() == 1);
        assert(intermediate_u[0]->size() == u1.size());
        assert(dt > 0 && "dt must be greater then zero!");

        Eigen::VectorXd f1(integral_value.size());
        fun(x1, u1, f1);

        Eigen::VectorXd f2(integral_value.size());
        fun(x2, u2, f2);

        *intermediate_x[0] = 0.5 * (x1 + x2) + dt / 8.0 * (f1 - f2);
        *intermediate_u[0] = 0.5 * (u1 + u2);

        fun(*intermediate_x[0], *intermediate_u[0], integral_value);

        integral_value = dt / 6.0 * (f1 + 4.0 * integral_value + f2);
    }

    bool interpolate(const Eigen::Ref<const Eigen::VectorXd>& x1, const Eigen::Ref<const Eigen::VectorXd>& u1,
                     const Eigen::Ref<const Eigen::VectorXd>& x2, const Eigen::Ref<const Eigen::VectorXd>& u2, double dt,
                     const SystemDynamicsInterface& system, const Range& range, std::vector<Eigen::VectorXd>& states,
                     std::vector<Eigen::VectorXd>& controls) override
    {
        if (range.getStart() < 0 || range.getEnd() > dt)
        {
            PRINT_ERROR_NAMED("range must be in the interval [0, dt]");
            return false;
        }

        states.push_back(x1);
        controls.push_back(u1);

        if (dt < 1e-15) return true;

        Eigen::VectorXd f1(system.getStateDimension());
        system.dynamics(x1, u1, f1);

        Eigen::VectorXd f2(system.getStateDimension());
        system.dynamics(x2, u2, f2);

        Eigen::VectorXd xc = 0.5 * (x1 + x2) + dt / 8.0 * (f1 - f2);
        Eigen::VectorXd uc = 0.5 * (u1 + u2);

        Eigen::VectorXd fc(system.getStateDimension());
        system.dynamics(xc, uc, fc);

        int n    = range.getNumInRange();
        double t = range.getStart() + range.getStep();  // x1 and u1 are already added
        for (int i = 1; i < n; ++i, t += range.getStep())
        {
            double tau = (t - range.getStart());

            // quadratic interpolation
            //            Eigen::VectorXd u_interp =
            //                2.0 / dt_square * ((tau - dt / 2.0) * (tau - dt) * u1 - 2.0 * tau * (tau - dt) * u_center + tau * (tau - dt / 2.0) *
            //                u2);

            // linear interpolation for controls
            controls.push_back(u1 + tau / dt * (u2 - u1));

            // cubic interpolation for states (since quadratic interplation for xdot)
            Eigen::VectorXd x_interp = x1 + f1 * tau + 0.5 / dt * tau * tau * (-3.0 * f1 + 4.0 * fc - f2) +
                                       1.0 / 3.0 / (dt * dt) * std::pow(tau, 3) * (2.0 * f1 - 4.0 * fc + 2.0 * f2);
            states.push_back(x_interp);
        }

        if (range.includeEnd())
        {
            states.push_back(x2);
            controls.push_back(u2);
        }
        return true;
    }

    bool computeIntermediateControls(const Eigen::Ref<const Eigen::VectorXd>& u1, const Eigen::Ref<const Eigen::VectorXd>& u2,
                                     std::vector<Eigen::VectorXd*>& intermediate_u) override
    {
        assert(intermediate_u.size() == 1);
        assert(intermediate_u[0]->size() == u1.size());
        *intermediate_u[0] = 0.5 * (u1 + u2);
        return true;
    }

#ifdef MESSAGE_SUPPORT
    // Implements interface method
    void toMessage(messages::Quadrature& message) const override
    {
        QuadratureCollocationInterface::toMessage(message);
        message.mutable_hermite_simpson_linear_u();
    }
#endif
};

FACTORY_REGISTER_QUADRATURE_COLLOCATION(QuadratureHermiteSimpsonLinearControl)

/**
 * @brief Hermite-Simpson rule (approximates function as quadratic polynomial between two points)
 *
 * @ingroup numerics
 *
 * \f[
 *      \int_0^{\Delta T} f(t) = \frac{\Delta T}{6} (f(0) + 4 f(0.5 \Delta T) + f(\Delta T)
 * \f].
 *
 * Hermite-Simpson collocation approximates the system dynamics and control trajectory by
 * quadratic polynomials between two consecutive points.
 * Hence, the resulting state trajectory is composed of cubic hermite splines.
 * Refer to
 * - https://mec560sbu.github.io/2016/09/30/direct_collocation/
 * - http://epubs.siam.org/doi/pdf/10.1137/16M1062569
 *
 * @see NumericalIntegratorInterface NumericalIntegratorImplicitInterface
 *
 * @author Christoph Rösmann (christoph.roesmann@tu-dortmund.de)
 */
class QuadratureHermiteSimpson : public QuadratureCollocationInterface
{
 public:
    // Implements interface method
    DynamicsEvalInterface::Ptr getInstance() const override { return std::make_shared<QuadratureHermiteSimpson>(); }

    // Implements interface method
    bool isSupportingCompressedStatesMode() const override { return true; }
    bool isIntermediateControlSubjectToOptim() const override { return true; }
    int getNumIntermediatePoints() const override { return 1; }

    int getNumIntermediateControls() const override { return 1; }
    int getNumIntermediateStates() const override { return 0; }

    // Implements interface method
    void quadrature(const std::vector<const Eigen::VectorXd*>& x1_points, const std::vector<const Eigen::VectorXd*>& u1_points,
                    const Eigen::VectorXd& u2, const Eigen::VectorXd& x2, double dt, const SystemDynamicsInterface& system,
                    Eigen::Ref<Eigen::VectorXd> integral_value) override
    {
        assert(x1_points.size() == 1);
        assert(u1_points.size() == 2);

        Eigen::VectorXd f1(integral_value.size());
        system.dynamics(*x1_points.front(), *u1_points[0], f1);

        Eigen::VectorXd f2(integral_value.size());
        system.dynamics(x2, u2, f2);

        Eigen::VectorXd x_mid = 0.5 * (*x1_points.front() + x2) + dt / 8.0 * (f1 - f2);
        system.dynamics(x_mid, *u1_points[1], integral_value);

        integral_value = dt / 6.0 * (f1 + 4.0 * integral_value + f2);
    }

    void quadrature(const std::vector<const Eigen::VectorXd*>& x1_points, const std::vector<const Eigen::VectorXd*>& u1_points,
                    const Eigen::VectorXd& u2, const Eigen::VectorXd& x2, double dt,
                    const BinaryFunction& fun, Eigen::Ref<Eigen::VectorXd> integral_value) override
    {
        assert(x1_points.size() == 1);
        assert(u1_points.size() == 2);
        
        Eigen::VectorXd f1(integral_value.size());
        fun(*x1_points.front(), *u1_points[0], f1);
        
        Eigen::VectorXd f2(integral_value.size());
        fun(x2, u2, f2);
        
        Eigen::VectorXd x_mid = 0.5 * (*x1_points.front() + x2) + dt / 8.0 * (f1 - f2);
        fun(x_mid, *u1_points[1], integral_value);
        
        integral_value = dt / 6.0 * (f1 + 4.0 * integral_value + f2);
    }
    
    bool interpolate(const std::vector<const Eigen::VectorXd*>& x1_points, const std::vector<const Eigen::VectorXd*>& u1_points,
                     const Eigen::VectorXd& u2, const Eigen::VectorXd& x2, double t1, double dt, double dt_interp,
                     const SystemDynamicsInterface& system, bool skip_u2_x2, TimeSeries& ts_x, TimeSeries& ts_u) override
    {
        assert(x1_points.size() == 1);
        assert(u1_points.size() == 2);

        if (dt_interp <= 0 || dt_interp >= dt)
        {
            ts_x.add(t1, *x1_points[0]);
            ts_u.add(t1, *u1_points[0]);
        }
        else if (dt > 0)
        {
            Eigen::VectorXd f1(x2.size());
            system.dynamics(*x1_points.front(), *u1_points.front(), f1);

            Eigen::VectorXd f2(x2.size());
            system.dynamics(x2, u2, f2);

            Eigen::VectorXd fc(x2.size());
            Eigen::VectorXd xc = 0.5 * (*x1_points.front() + x2) + dt / 8.0 * (f1 - f2);
            system.dynamics(xc, *u1_points[1], fc);

            int n    = std::floor(dt / dt_interp);
            double t = t1;
            for (int i = 0; i < n; ++i)
            {
                double tau = (t - t1);

                // cubic interpolation for states (since quadratic interplation for xdot)
                Eigen::VectorXd x_interp = *x1_points.front() + f1 * tau + 0.5 / dt * tau * tau * (-3.0 * f1 + 4.0 * fc - f2) +
                                           1.0 / 3.0 / (dt * dt) * std::pow(tau, 3) * (2.0 * f1 - 4.0 * fc + 2.0 * f2);

                Eigen::VectorXd u_interp = *u1_points[0] + tau / dt * (-3.0 * *u1_points[0] + 4.0 * *u1_points[1] - u2) +
                                           2.0 * tau * tau / (dt * dt) * (*u1_points[0] - 2.0 * *u1_points[1] + u2);

                ts_x.add(t1 + tau, x_interp);
                ts_u.add(t1 + tau, u_interp);

                t += dt_interp;
            }
        }

        if (!skip_u2_x2)
        {
            ts_x.add(t1 + dt, x2);
            ts_u.add(t1 + dt, u2);
        }
        return true;
    }

    void quadrature(const StateVector& x1, const InputVector& u1, const StateVector& x2, double dt, const SystemDynamicsInterface& system,
                    Eigen::Ref<Eigen::VectorXd> integral_value) override
    {
        assert(integral_value.size() == x1.size());
        assert(dt > 0 && "dt must be greater then zero!");

        Eigen::VectorXd f1(x1.size());
        system.dynamics(x1, u1, f1);

        Eigen::VectorXd f2(x1.size());
        system.dynamics(x2, u1, f2);

        Eigen::VectorXd x_center = 0.5 * (x1 + x2) + dt / 8.0 * (f1 - f2);

        system.dynamics(x_center, u1, integral_value);

        integral_value = dt / 6.0 * (f1 + 4.0 * integral_value + f2);
    }

    void quadrature(const StateVector& x1, const InputVector& u1, const StateVector& x2, double dt, const SystemDynamicsInterface& system,
                    Eigen::Ref<Eigen::VectorXd> integral_value, std::vector<Eigen::VectorXd*>& intermediate_x) override
    {
        assert(integral_value.size() == x1.size());
        assert(dt > 0 && "dt must be greater then zero!");
        assert(intermediate_x.size() == 1);
        assert(intermediate_x[0]->size() == x1.size());

        Eigen::VectorXd f1(x1.size());
        system.dynamics(x1, u1, f1);

        Eigen::VectorXd f2(x1.size());
        system.dynamics(x2, u1, f2);

        intermediate_x[0]->noalias() = 0.5 * (x1 + x2) + dt / 8.0 * (f1 - f2);

        system.dynamics(*intermediate_x[0], u1, integral_value);

        integral_value = dt / 6.0 * (f1 + 4.0 * integral_value + f2);
    }

    // Implements interface method
    void quadrature(const StateVector& x1, const InputVector& u1, const StateVector& x2, const InputVector& u2, double dt,
                    const SystemDynamicsInterface& system, Eigen::Ref<Eigen::VectorXd> integral_value) override
    {
        PRINT_WARNING_NAMED(
            "This function overload is not supported by HermiteSimpson-Collocation since we require the midpoint between u1 and u2 as extra "
            "argument.");
    }

    void quadrature(const StateVector& x1, const InputVector& u1, const StateVector& x2, const InputVector& u2, double dt,
                    const SystemDynamicsInterface& system, const std::vector<Eigen::VectorXd*>& intermediate_x,
                    const std::vector<Eigen::VectorXd*>& intermediate_u, Eigen::Ref<Eigen::VectorXd> integral_value,
                    Eigen::Ref<Eigen::VectorXd> interm_x_error) override
    {
        assert(intermediate_x.size() == 1);
        assert(intermediate_x[0]->size() == x1.size());
        assert(intermediate_u.size() == 1);
        assert(intermediate_u[0]->size() == u1.size());
        assert(dt > 0 && "dt must be greater then zero!");

        Eigen::VectorXd f1(integral_value.size());
        system.dynamics(x1, u1, f1);

        Eigen::VectorXd f2(integral_value.size());
        system.dynamics(x2, u2, f2);

        system.dynamics(*intermediate_x[0], *intermediate_u[0], integral_value);

        integral_value = dt / 6.0 * (f1 + 4.0 * integral_value + f2);

        assert(interm_x_error.size() == x1.size());
        interm_x_error = *intermediate_x[0] - (0.5 * (x1 + x2) + dt / 8.0 * (f1 - f2));
    }

    void quadrature(const StateVector& x1, const InputVector& u1, const StateVector& x2, const InputVector& u2, double dt,
                    const SystemDynamicsInterface& system, const std::vector<Eigen::VectorXd*>& intermediate_u,
                    Eigen::Ref<Eigen::VectorXd> integral_value, std::vector<Eigen::VectorXd*>& intermediate_x) override
    {
        assert(intermediate_x.size() == 1);
        assert(intermediate_x[0]->size() == x1.size());
        assert(intermediate_u.size() == 1);
        assert(intermediate_u[0]->size() == u1.size());
        assert(dt > 0 && "dt must be greater then zero!");

        Eigen::VectorXd f1(integral_value.size());
        system.dynamics(x1, u1, f1);

        Eigen::VectorXd f2(integral_value.size());
        system.dynamics(x2, u2, f2);

        intermediate_x[0]->noalias() = 0.5 * (x1 + x2) + dt / 8.0 * (f1 - f2);
        system.dynamics(*intermediate_x[0], *intermediate_u[0], integral_value);

        integral_value = dt / 6.0 * (f1 + 4.0 * integral_value + f2);
    }

    void quadrature(const StateVector& x1, const InputVector& u1, const StateVector& x2, const InputVector& u2, double dt,
                    const SystemDynamicsInterface& system, Eigen::Ref<Eigen::VectorXd> integral_value, std::vector<Eigen::VectorXd*>& intermediate_x,
                    std::vector<Eigen::VectorXd*>& intermediate_u) override
    {
        PRINT_WARNING_NAMED(
            "This function overload is not supported by HermiteSimpson-Collocation since we require the midpoint between u1 and u2 as extra "
            "argument.");
    }

    // Implements interface method
    void quadrature(const Eigen::Ref<const Eigen::VectorXd>& x1, const Eigen::Ref<const Eigen::VectorXd>& x2, double dt, const UnaryFunction& fun,
                    Eigen::Ref<Eigen::VectorXd> integral_value) override
    {
        assert(dt > 0 && "dt must be greater then zero!");

        Eigen::VectorXd f1(integral_value.size());
        fun(x1, f1);

        Eigen::VectorXd f2(integral_value.size());
        fun(x2, f2);

        Eigen::VectorXd x_center = 0.5 * (x1 + x2) + dt / 8.0 * (f1 - f2);

        fun(x_center, integral_value);

        integral_value = dt / 6.0 * (f1 + 4.0 * integral_value + f2);
    }

    // Implements interface method
    void quadrature(const Eigen::Ref<const Eigen::VectorXd>& x1, const Eigen::Ref<const Eigen::VectorXd>& x2, double dt,
                    const std::vector<Eigen::VectorXd*>& intermediate_x, const UnaryFunction& fun,
                    Eigen::Ref<Eigen::VectorXd> integral_value) override
    {
        assert(intermediate_x.size() == 1);
        assert(intermediate_x[0]->size() == x1.size());
        assert(dt > 0 && "dt must be greater then zero!");

        Eigen::VectorXd f1(integral_value.size());
        fun(x1, f1);

        Eigen::VectorXd f2(integral_value.size());
        fun(x2, f2);

        fun(*intermediate_x[0], integral_value);

        integral_value = dt / 6.0 * (f1 + 4.0 * integral_value + f2);
    }

    // Implements interface method
    void quadrature(const Eigen::Ref<const Eigen::VectorXd>& x1, const Eigen::Ref<const Eigen::VectorXd>& x2, double dt, const UnaryFunction& fun,
                    Eigen::Ref<Eigen::VectorXd> integral_value, std::vector<Eigen::VectorXd*>& intermediate_x) override
    {
        assert(intermediate_x.size() == 1);
        assert(intermediate_x[0]->size() == x1.size());
        assert(dt > 0 && "dt must be greater then zero!");

        Eigen::VectorXd f1(integral_value.size());
        fun(x1, f1);

        Eigen::VectorXd f2(integral_value.size());
        fun(x2, f2);

        *intermediate_x[0] = 0.5 * (x1 + x2) + dt / 8.0 * (f1 - f2);

        fun(*intermediate_x[0], integral_value);

        integral_value = dt / 6.0 * (f1 + 4.0 * integral_value + f2);
    }

    // Implements interface method
    void quadrature(const Eigen::Ref<const Eigen::VectorXd>& x1, const Eigen::Ref<const Eigen::VectorXd>& u1,
                    const Eigen::Ref<const Eigen::VectorXd>& x2, const Eigen::Ref<const Eigen::VectorXd>& u2, double dt, const BinaryFunction& fun,
                    Eigen::Ref<Eigen::VectorXd> integral_value) override
    {
        PRINT_WARNING_NAMED(
            "This function overload is not supported by HermiteSimpson-Collocation since we require the midpoint between u1 and u2 as extra "
            "argument.");
    }

    // Implements interface method
    void quadrature(const Eigen::Ref<const Eigen::VectorXd>& x1, const Eigen::Ref<const Eigen::VectorXd>& u1,
                    const Eigen::Ref<const Eigen::VectorXd>& x2, const Eigen::Ref<const Eigen::VectorXd>& u2, double dt,
                    const std::vector<Eigen::VectorXd*>& intermediate_x, const std::vector<Eigen::VectorXd*>& intermediate_u,
                    const BinaryFunction& fun, Eigen::Ref<Eigen::VectorXd> integral_value) override
    {
        assert(intermediate_x.size() == 1);
        assert(intermediate_x[0]->size() == x1.size());
        assert(intermediate_u.size() == 1);
        assert(intermediate_u[0]->size() == u1.size());
        assert(dt > 0 && "dt must be greater then zero!");

        Eigen::VectorXd f1(integral_value.size());
        fun(x1, u1, f1);

        Eigen::VectorXd f2(integral_value.size());
        fun(x2, u2, f2);

        fun(*intermediate_x[0], *intermediate_u[0], integral_value);

        integral_value = dt / 6.0 * (f1 + 4.0 * integral_value + f2);
    }

    // Implements interface method
    void quadrature(const Eigen::Ref<const Eigen::VectorXd>& x1, const Eigen::Ref<const Eigen::VectorXd>& u1,
                    const Eigen::Ref<const Eigen::VectorXd>& x2, const Eigen::Ref<const Eigen::VectorXd>& u2, double dt, const BinaryFunction& fun,
                    Eigen::Ref<Eigen::VectorXd> integral_value, std::vector<Eigen::VectorXd*>& intermediate_x,
                    std::vector<Eigen::VectorXd*>& intermediate_u) override
    {
        PRINT_WARNING_NAMED(
            "This function overload is not supported by HermiteSimpson-Collocation since we require the midpoint between u1 and u2 as extra "
            "argument.");
    }

    bool interpolate(const Eigen::Ref<const Eigen::VectorXd>& x1, const Eigen::Ref<const Eigen::VectorXd>& u1,
                     const Eigen::Ref<const Eigen::VectorXd>& x2, const Eigen::Ref<const Eigen::VectorXd>& u2, double dt,
                     const SystemDynamicsInterface& system, const Range& range, std::vector<Eigen::VectorXd>& states,
                     std::vector<Eigen::VectorXd>& controls) override
    {
        return false;
    }

    bool interpolate(const std::vector<const Eigen::VectorXd*>& x, const std::vector<const Eigen::VectorXd*>& u, double dt,
                     const SystemDynamicsInterface& system, const Range& range, std::vector<Eigen::VectorXd>& states,
                     std::vector<Eigen::VectorXd>& controls) override
    {
        if ((x.size() != 2 && x.size() != 3) || u.size() != 3)
        {
            PRINT_DEBUG_NAMED("x-dim must be 2 in compressed mode and 3 in uncompressed mode. U-dim must be 3");
            return false;
        }

        if (range.getStart() < 0 || range.getEnd() > dt)
        {
            PRINT_ERROR_NAMED("range must be in the interval [0, dt]");
            return false;
        }

        states.push_back(*x[0]);
        controls.push_back(*u[0]);

        if (dt < 1e-15) return true;

        Eigen::VectorXd f1(system.getStateDimension());
        system.dynamics(*x[0], *u[0], f1);

        Eigen::VectorXd f2(system.getStateDimension());
        system.dynamics(*x.back(), *u[2], f2);

        Eigen::VectorXd fc(system.getStateDimension());
        if (x.size() == 2)
        {
            Eigen::VectorXd xc = 0.5 * (*x[0] + *x[1]) + dt / 8.0 * (f1 - f2);
            system.dynamics(xc, *u[1], fc);
        }
        else
        {
            assert(x.size() == 3);
            system.dynamics(*x[1], *u[1], fc);
        }

        int n    = range.getNumInRange();
        double t = range.getStart() + range.getStep();  // x1 and u1 are already added
        for (int i = 1; i < n; ++i, t += range.getStep())
        {
            double tau = (t - range.getStart());

            // quadratic interpolation
            // Eigen::VectorXd u_interp =
            //    2.0 / (dt * dt) * ((tau - dt / 2.0) * (tau - dt) * *u[0] - 2.0 * tau * (tau - dt) * *u[1] + tau * (tau - dt / 2.0) * *u[2]);
            Eigen::VectorXd u_interp = *u[0] + tau / dt * (-3 * *u[0] + 4 * *u[1] - *u[2]) + 2 * tau * tau / (dt * dt) * (*u[0] - 2 * *u[1] + *u[2]);

            // linear interpolation for controls
            controls.push_back(u_interp);

            // cubic interpolation for states (since quadratic interplation for xdot)
            Eigen::VectorXd x_interp = *x[0] + f1 * tau + 0.5 / dt * tau * tau * (-3.0 * f1 + 4.0 * fc - f2) +
                                       1.0 / 3.0 / (dt * dt) * std::pow(tau, 3) * (2.0 * f1 - 4.0 * fc + 2.0 * f2);
            states.push_back(x_interp);
        }

        if (range.includeEnd())
        {
            states.push_back(*x.back());
            controls.push_back(*u[2]);
        }
        return true;
    }

    bool computeIntermediateControls(const Eigen::Ref<const Eigen::VectorXd>& u1, const Eigen::Ref<const Eigen::VectorXd>& u2,
                                     std::vector<Eigen::VectorXd*>& intermediate_u) override
    {
        return false;
    }

#ifdef MESSAGE_SUPPORT
    // Implements interface method
    void toMessage(messages::Quadrature& message) const override
    {
        QuadratureCollocationInterface::toMessage(message);
        message.mutable_hermite_simpson();
    }
#endif
};

FACTORY_REGISTER_QUADRATURE_COLLOCATION(QuadratureHermiteSimpson)

}  // namespace corbo

#endif  // SRC_NUMERICS_INCLUDE_CORBO_NUMERICS_QUADRATURE_H_
