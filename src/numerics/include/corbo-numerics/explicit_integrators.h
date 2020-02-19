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

#ifndef SRC_NUMERICS_INCLUDE_CORBO_NUMERICS_EXPLICIT_INTEGRATORS_H_
#define SRC_NUMERICS_INCLUDE_CORBO_NUMERICS_EXPLICIT_INTEGRATORS_H_

#include <corbo-communication/utilities.h>
#include <corbo-numerics/integrator_interface.h>
#include <memory>

namespace corbo {

/**
 * @brief Forward euler (explicit euler) integration
 *
 * @ingroup numerics
 *
 * \f[
 *      x_2 = x_1 + \Delta T f(x,u)
 * \f].
 *
 * @see NumericalIntegratorInterface NumericalIntegratorExplicitInterface
 *
 * @author Christoph Rösmann (christoph.roesmann@tu-dortmund.de)
 */
class IntegratorExplicitEuler : public NumericalIntegratorExplicitInterface
{
 public:
    // Implements interface method
    DynamicsEvalInterface::Ptr getInstance() const override { return std::make_shared<IntegratorExplicitEuler>(); }

    // Implements interface method
    int getConvergenceOrder() const override { return 1; }

    // Implements interface method
    void solveIVP(const Eigen::VectorXd& x1, double dt, const std::function<void(const Eigen::VectorXd&, Eigen::Ref<Eigen::VectorXd>)>& fun,
                  Eigen::Ref<Eigen::VectorXd> x2) override
    {
        fun(x1, x2);
        x2 *= dt;
        x2 += x1;
    }

    // Implements interface method
    void solveIVP(const StateVector& x1, const InputVector& u1, double dt, const SystemDynamicsInterface& system,
                  Eigen::Ref<Eigen::VectorXd> x2) override
    {
        system.dynamics(x1, u1, x2);
        x2 *= dt;
        x2 += x1;
    }

#ifdef MESSAGE_SUPPORT
    // Implements interface method
    void toMessage(messages::ExplicitIntegrator& message) const override
    {
        NumericalIntegratorExplicitInterface::toMessage(message);
        message.mutable_explicit_euler();
    }
#endif
};

FACTORY_REGISTER_EXPLICIT_INTEGRATOR(IntegratorExplicitEuler)

/**
 * @brief 2th Order Runge-Kutta Integrator (explicit)
 *
 * @ingroup numerics
 *
 * Refer to http://www.mathematik.uni-stuttgart.de/studium/infomat/Numerische-Mathematik-II-SS11/Matlab/ode_4.pdf
 *
 * @see NumericalIntegratorInterface NumericalIntegratorExplicitInterface
 *
 * @author Alexander Westermann (alexander.westermann@tu-dortmund.de)
 */
class IntegratorExplicitRungeKutta2 : public NumericalIntegratorExplicitInterface
{
 public:
    // Implements interface method
    DynamicsEvalInterface::Ptr getInstance() const override { return std::make_shared<IntegratorExplicitRungeKutta2>(); }

    // Implements interface method
    int getConvergenceOrder() const override { return 2; }

    void initialize(int state_dim) override
    {
        _k1.resize(state_dim);
        _k2.resize(state_dim);
    }

    // Implements interface method
    void solveIVP(const Eigen::VectorXd& x1, double dt, const std::function<void(const Eigen::VectorXd&, Eigen::Ref<Eigen::VectorXd>)>& fun,
                  Eigen::Ref<Eigen::VectorXd> x2) override
    {
        if (x1.size() != _k1.size()) initialize(x1.size());

        fun(x1, _k1);  // k1
        _k1 *= dt;
        fun(x1 + _k1, _k2);  // k2
        _k2 *= dt;

        x2 = x1 + (_k1 + _k2) / 2.0;
    }

    // Implements interface method
    void solveIVP(const StateVector& x1, const InputVector& u1, double dt, const SystemDynamicsInterface& system,
                  Eigen::Ref<Eigen::VectorXd> x2) override
    {
        if (x1.size() != _k1.size()) initialize(x1.size());

        system.dynamics(x1, u1, _k1);  // k1
        _k1 *= dt;
        system.dynamics(x1 + _k1, u1, _k2);  // k2
        _k2 *= dt;

        x2 = x1 + (_k1 + _k2) / 2.0;
    }

#ifdef MESSAGE_SUPPORT
    // Implements interface method
    void toMessage(messages::ExplicitIntegrator& message) const override
    {
        NumericalIntegratorExplicitInterface::toMessage(message);
        message.mutable_runge_kutta_2();
    }
#endif

 private:
    Eigen::VectorXd _k1;
    Eigen::VectorXd _k2;
};

FACTORY_REGISTER_EXPLICIT_INTEGRATOR(IntegratorExplicitRungeKutta2)

/**
 * @brief 3th Order Runge-Kutta Integrator (explicit)
 *
 * @ingroup numerics
 *
 * Refer to https://www.math.uni-hamburg.de/home/hofmann/lehrveranstaltungen/sommer05/prosem/Vortrag12.pdf
 *
 * @see NumericalIntegratorInterface NumericalIntegratorExplicitInterface
 *
 * @author Alexander Westermann (alexander.westermann@tu-dortmund.de)
 */
class IntegratorExplicitRungeKutta3 : public NumericalIntegratorExplicitInterface
{
 public:
    // Implements interface method
    DynamicsEvalInterface::Ptr getInstance() const override { return std::make_shared<IntegratorExplicitRungeKutta3>(); }

    // Implements interface method
    int getConvergenceOrder() const override { return 3; }

    void initialize(int state_dim) override
    {
        _k1.resize(state_dim);
        _k2.resize(state_dim);
        _k3.resize(state_dim);
    }

    // Implements interface method
    void solveIVP(const Eigen::VectorXd& x1, double dt, const std::function<void(const Eigen::VectorXd&, Eigen::Ref<Eigen::VectorXd>)>& fun,
                  Eigen::Ref<Eigen::VectorXd> x2) override
    {
        if (x1.size() != _k1.size()) initialize(x1.size());

        fun(x1, _k1);  // k1
        _k1 *= dt;
        fun(x1 + (_k1 / 2.0), _k2);  // k2
        _k2 *= dt;
        fun(x1 - _k1 + 2.0 * _k2, _k3);  // k3
        _k3 *= dt;

        x2 = x1 + (_k1 + 4.0 * _k2 + _k3) / 6.0;
    }

    // Implements interface method
    void solveIVP(const StateVector& x1, const InputVector& u1, double dt, const SystemDynamicsInterface& system,
                  Eigen::Ref<Eigen::VectorXd> x2) override
    {
        if (x1.size() != _k1.size()) initialize(x1.size());

        system.dynamics(x1, u1, _k1);  // k1
        _k1 *= dt;
        system.dynamics(x1 + (_k1 / 2.0), u1, _k2);  // k2
        _k2 *= dt;
        system.dynamics(x1 - _k1 + 2.0 * _k2, u1, _k3);  // k3
        _k3 *= dt;

        x2 = x1 + (_k1 + 4.0 * _k2 + _k3) / 6.0;
    }

#ifdef MESSAGE_SUPPORT
    // Implements interface method
    void toMessage(messages::ExplicitIntegrator& message) const override
    {
        NumericalIntegratorExplicitInterface::toMessage(message);
        message.mutable_runge_kutta_3();
    }
#endif

 private:
    Eigen::VectorXd _k1;
    Eigen::VectorXd _k2;
    Eigen::VectorXd _k3;
};

FACTORY_REGISTER_EXPLICIT_INTEGRATOR(IntegratorExplicitRungeKutta3)

/**
 * @brief 4th Order Runge-Kutta Integrator (explicit)
 *
 * @ingroup numerics
 *
 * Also called as classical Runge-Kutta method.
 * Refer to http://en.wikipedia.org/wiki/Runge–Kutta_methods
 *
 * @see NumericalIntegratorInterface NumericalIntegratorExplicitInterface
 *
 * @author Christoph Rösmann (christoph.roesmann@tu-dortmund.de)
 */
class IntegratorExplicitRungeKutta4 : public NumericalIntegratorExplicitInterface
{
 public:
    // Implements interface method
    DynamicsEvalInterface::Ptr getInstance() const override { return std::make_shared<IntegratorExplicitRungeKutta4>(); }

    // Implements interface method
    int getConvergenceOrder() const override { return 4; }

    void initialize(int state_dim) override
    {
        _k1.resize(state_dim);
        _k2.resize(state_dim);
        _k3.resize(state_dim);
        _k4.resize(state_dim);
    }

    // Implements interface method
    void solveIVP(const Eigen::VectorXd& x1, double dt, const std::function<void(const Eigen::VectorXd&, Eigen::Ref<Eigen::VectorXd>)>& fun,
                  Eigen::Ref<Eigen::VectorXd> x2) override
    {
        if (x1.size() != _k1.size()) initialize(x1.size());

        fun(x1, _k1);  // k1
        _k1 *= dt;
        fun(x1 + _k1 / 2.0, _k2);  // k2
        _k2 *= dt;
        fun(x1 + _k2 / 2.0, _k3);  // k3
        _k3 *= dt;
        fun(x1 + _k3, _k4);  // k4
        _k4 *= dt;

        x2 = x1 + (_k1 + 2.0 * _k2 + 2.0 * _k3 + _k4) / 6.0;  // TODO(roesmann) Is it possible to express k4 with x2 -> less temporaries
    }

    // Implements interface method
    void solveIVP(const StateVector& x1, const InputVector& u1, double dt, const SystemDynamicsInterface& system,
                  Eigen::Ref<Eigen::VectorXd> x2) override
    {
        if (x1.size() != _k1.size()) initialize(x1.size());

        system.dynamics(x1, u1, _k1);  // k1
        _k1 *= dt;
        system.dynamics(x1 + _k1 / 2.0, u1, _k2);  // k2
        _k2 *= dt;
        system.dynamics(x1 + _k2 / 2.0, u1, _k3);  // k3
        _k3 *= dt;
        system.dynamics(x1 + _k3, u1, _k4);  // k4
        _k4 *= dt;

        x2 = x1 + (_k1 + 2.0 * _k2 + 2.0 * _k3 + _k4) / 6.0;  // TODO(roesmann) Is it possible to express k4 with x2 -> less temporaries
    }

#ifdef MESSAGE_SUPPORT
    // Implements interface method
    void toMessage(messages::ExplicitIntegrator& message) const override
    {
        NumericalIntegratorExplicitInterface::toMessage(message);
        message.mutable_runge_kutta_4();
    }
#endif

 private:
    Eigen::VectorXd _k1;
    Eigen::VectorXd _k2;
    Eigen::VectorXd _k3;
    Eigen::VectorXd _k4;
};

FACTORY_REGISTER_EXPLICIT_INTEGRATOR(IntegratorExplicitRungeKutta4)

/**
 * @brief 5th Order Runge-Kutta Integrator (explicit, slightly modified)
 *
 * @ingroup numerics
 *
 * For more information please refer to
 * http://www.jstor.org/stable/pdfplus/2027775.pdf?acceptTC=true&jpdConfirm=true (Equation 2)
 *
 * @see NumericalIntegratorInterface NumericalIntegratorExplicitInterface IntegratorExplicitRungeKutta4
 *
 * @author Christoph Rösmann (christoph.roesmann@tu-dortmund.de)
 */
class IntegratorExplicitRungeKutta5 : public NumericalIntegratorExplicitInterface
{
 public:
    // Implements interface method
    DynamicsEvalInterface::Ptr getInstance() const override { return std::make_shared<IntegratorExplicitRungeKutta5>(); }

    // Implements interface method
    int getConvergenceOrder() const override { return 5; }

    void initialize(int state_dim) override
    {
        _k1.resize(state_dim);
        _k2.resize(state_dim);
        _k3.resize(state_dim);
        _k4.resize(state_dim);
        _k5.resize(state_dim);
        _k6.resize(state_dim);
    }

    // Implements interface method
    void solveIVP(const Eigen::VectorXd& x1, double dt, const std::function<void(const Eigen::VectorXd&, Eigen::Ref<Eigen::VectorXd>)>& fun,
                  Eigen::Ref<Eigen::VectorXd> x2) override
    {
        if (x1.size() != _k1.size()) initialize(x1.size());

        fun(x1, _k1);  // k1
        _k1 *= dt;
        fun(x1 + 4.0 * _k1 / 11.0, _k2);  // k2
        _k2 *= dt;
        fun(x1 + (9.0 * _k1 + 11.0 * _k2) / 50.0, _k3);  // k3
        _k3 *= dt;
        fun(x1 + (-11.0 * _k2 + 15.0 * _k3) / 4.0, _k4);  // k4
        _k4 *= dt;
        fun(x1 + ((81.0 + 9.0 * std::sqrt(6.0)) * _k1 + (255.0 - 55.0 * std::sqrt(6.0)) * _k3 + (24.0 - 14.0 * std::sqrt(6.0)) * _k4) / 600.0,
            _k5);  // k5
        _k5 *= dt;
        fun(x1 + ((81.0 - 9.0 * std::sqrt(6.0)) * _k1 + (255.0 + 55.0 * std::sqrt(6.0)) * _k3 + (24.0 + 14.0 * std::sqrt(6.0)) * _k4) / 600.0,
            _k6);  // k6
        _k6 *= dt;

        x2 = x1 + (4.0 * _k1 + (16.0 + std::sqrt(6.0)) * _k5 + (16.0 - std::sqrt(6.0)) * _k6) / 36.0;
    }

    // Implements interface method
    void solveIVP(const StateVector& x1, const InputVector& u1, double dt, const SystemDynamicsInterface& system,
                  Eigen::Ref<Eigen::VectorXd> x2) override
    {
        if (x1.size() != _k1.size()) initialize(x1.size());

        system.dynamics(x1, u1, _k1);  // k1
        _k1 *= dt;
        system.dynamics(x1 + 4.0 * _k1 / 11.0, u1, _k2);  // k2
        _k2 *= dt;
        system.dynamics(x1 + (9.0 * _k1 + 11.0 * _k2) / 50.0, u1, _k3);  // k3
        _k3 *= dt;
        system.dynamics(x1 + (-11.0 * _k2 + 15.0 * _k3) / 4.0, u1, _k4);  // k4
        _k4 *= dt;
        system.dynamics(
            x1 + ((81.0 + 9.0 * std::sqrt(6.0)) * _k1 + (255.0 - 55.0 * std::sqrt(6.0)) * _k3 + (24.0 - 14.0 * std::sqrt(6.0)) * _k4) / 600.0, u1,
            _k5);  // k5
        _k5 *= dt;
        system.dynamics(
            x1 + ((81.0 - 9.0 * std::sqrt(6.0)) * _k1 + (255.0 + 55.0 * std::sqrt(6.0)) * _k3 + (24.0 + 14.0 * std::sqrt(6.0)) * _k4) / 600.0, u1,
            _k6);  // k6
        _k6 *= dt;

        x2 = x1 + (4.0 * _k1 + (16.0 + std::sqrt(6.0)) * _k5 + (16.0 - std::sqrt(6.0)) * _k6) / 36.0;
    }

#ifdef MESSAGE_SUPPORT
    // Implements interface method
    void toMessage(messages::ExplicitIntegrator& message) const override
    {
        NumericalIntegratorExplicitInterface::toMessage(message);
        message.mutable_runge_kutta_5();
    }
#endif

 private:
    Eigen::VectorXd _k1;
    Eigen::VectorXd _k2;
    Eigen::VectorXd _k3;
    Eigen::VectorXd _k4;
    Eigen::VectorXd _k5;
    Eigen::VectorXd _k6;
};

FACTORY_REGISTER_EXPLICIT_INTEGRATOR(IntegratorExplicitRungeKutta5)

/**
 * @brief 6th Order Runge-Kutta Integrator (explicit)
 *
 * @ingroup numerics
 *
 * For more information please refer to
 * https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19680027281.pdf

 *
 * @see NumericalIntegratorInterface NumericalIntegratorExplicitInterface IntegratorExplicitRungeKutta6
 *
 * @author Alexander Westermann (alexander.westermann@tu-dortmund.de)
 */
class IntegratorExplicitRungeKutta6 : public NumericalIntegratorExplicitInterface
{
 public:
    // Implements interface method
    DynamicsEvalInterface::Ptr getInstance() const override { return std::make_shared<IntegratorExplicitRungeKutta6>(); }

    // Implements interface method
    int getConvergenceOrder() const override { return 6; }

    void initialize(int state_dim) override
    {
        _k1.resize(state_dim);
        _k2.resize(state_dim);
        _k3.resize(state_dim);
        _k4.resize(state_dim);
        _k5.resize(state_dim);
        _k6.resize(state_dim);
        _k7.resize(state_dim);
        _k8.resize(state_dim);
    }

    // Implements interface method
    void solveIVP(const Eigen::VectorXd& x1, double dt, const std::function<void(const Eigen::VectorXd&, Eigen::Ref<Eigen::VectorXd>)>& fun,
                  Eigen::Ref<Eigen::VectorXd> x2) override
    {
        if (x1.size() != _k1.size()) initialize(x1.size());

        fun(x1, _k1);  // k1
        _k1 *= dt;
        fun(x1 + 2.0 * _k1 / 33.0, _k2);  // k2
        _k2 *= dt;
        fun(x1 + 4.0 * _k2 / 33.0, _k3);  // k3
        _k3 *= dt;
        fun(x1 + (_k1 + 3.0 * _k3) / 22.0, _k4);  // k4
        _k4 *= dt;
        fun(x1 + (43.0 * _k1 - 165.0 * _k3 + 144.0 * _k4) / 64.0, _k5);  // k5
        _k5 *= dt;
        fun(x1 + (-4053483.0 * _k1 + 16334703.0 * _k3 - 12787632.0 * _k4 + 1057536.0 * _k5) / 826686.0, _k6);  // k6
        _k6 *= dt;
        fun(x1 + (169364139.0 * _k1 - 663893307.0 * _k3 + 558275718.0 * _k4 - 29964480.0 * _k5 + 35395542.0 * _k6) / 80707214.0, _k7);  // k7
        _k7 *= dt;
        fun(x1 + (-733.0 * _k1 + 3102.0 * _k3) / 176.0 - (335763.0 * _k4 / 23296.0) + (216.0 * _k5 / 77.0) - (4617.0 * _k6 / 2816.0) +
                (7203.0 * _k7 / 9152.0),
            _k8);  // k8
        _k8 *= dt;

        x2 = x1 + (336336.0 * _k1 + 1771561.0 * _k4 + 1916928.0 * _k5 + 597051.0 * _k6 + 1411788.0 * _k7 + 256256.0 * _k8) / 6289920.0;
    }

    // Implements interface method
    void solveIVP(const StateVector& x1, const InputVector& u1, double dt, const SystemDynamicsInterface& system,
                  Eigen::Ref<Eigen::VectorXd> x2) override
    {
        if (x1.size() != _k1.size()) initialize(x1.size());

        system.dynamics(x1, u1, _k1);  // k1
        _k1 *= dt;
        system.dynamics(x1 + 2.0 * _k1 / 33.0, u1, _k2);  // k2
        _k2 *= dt;
        system.dynamics(x1 + 4.0 * _k2 / 33.0, u1, _k3);  // k3
        _k3 *= dt;
        system.dynamics(x1 + (_k1 + 3.0 * _k3) / 22.0, u1, _k4);  // k4
        _k4 *= dt;
        system.dynamics(x1 + (43.0 * _k1 - 165.0 * _k3 + 144.0 * _k4) / 64.0, u1, _k5);  // k5
        _k5 *= dt;
        system.dynamics(x1 + (-4053483.0 * _k1 + 16334703.0 * _k3 - 12787632.0 * _k4 + 1057536.0 * _k5) / 826686.0, u1, _k6);  // k6
        _k6 *= dt;
        system.dynamics(x1 + (169364139.0 * _k1 - 663893307.0 * _k3 + 558275718.0 * _k4 - 29964480.0 * _k5 + 35395542.0 * _k6) / 80707214.0, u1,
                        _k7);  // k7
        _k7 *= dt;
        system.dynamics(x1 + (-733.0 * _k1 + 3102.0 * _k3) / 176.0 - (335763.0 * _k4 / 23296.0) + (216.0 * _k5 / 77.0) - (4617.0 * _k6 / 2816.0) +
                            (7203.0 * _k7 / 9152.0),
                        u1, _k8);  // k8
        _k8 *= dt;

        x2 = x1 + (336336.0 * _k1 + 1771561.0 * _k4 + 1916928.0 * _k5 + 597051.0 * _k6 + 1411788.0 * _k7 + 256256.0 * _k8) / 6289920.0;
    }

#ifdef MESSAGE_SUPPORT
    // Implements interface method
    void toMessage(messages::ExplicitIntegrator& message) const override
    {
        NumericalIntegratorExplicitInterface::toMessage(message);
        message.mutable_runge_kutta_6();
    }
#endif

 private:
    Eigen::VectorXd _k1;
    Eigen::VectorXd _k2;
    Eigen::VectorXd _k3;
    Eigen::VectorXd _k4;
    Eigen::VectorXd _k5;
    Eigen::VectorXd _k6;
    Eigen::VectorXd _k7;
    Eigen::VectorXd _k8;
};

FACTORY_REGISTER_EXPLICIT_INTEGRATOR(IntegratorExplicitRungeKutta6)

/**
 * @brief 7th Order Runge-Kutta Integrator (explicit)
 *
 * @ingroup numerics
 *
 * For more information please refer to
 * https://link.springer.com/content/pdf/10.1007%2FBF02234758.pdf
 *
 * @see NumericalIntegratorInterface NumericalIntegratorExplicitInterface IntegratorExplicitRungeKutta7
 *
 * @author Alexander Westermann (alexander.westermann@tu-dortmund.de)
 */
class IntegratorExplicitRungeKutta7 : public NumericalIntegratorExplicitInterface
{
 public:
    // Implements interface method
    DynamicsEvalInterface::Ptr getInstance() const override { return std::make_shared<IntegratorExplicitRungeKutta7>(); }

    // Implements interface method
    int getConvergenceOrder() const override { return 7; }

    void initialize(int state_dim) override
    {
        _k1.resize(state_dim);
        _k2.resize(state_dim);
        _k3.resize(state_dim);
        _k4.resize(state_dim);
        _k5.resize(state_dim);
        _k6.resize(state_dim);
        _k7.resize(state_dim);
        _k8.resize(state_dim);
        _k9.resize(state_dim);
        _k10.resize(state_dim);
        _k11.resize(state_dim);
    }

    // Implements interface method
    void solveIVP(const Eigen::VectorXd& x1, double dt, const std::function<void(const Eigen::VectorXd&, Eigen::Ref<Eigen::VectorXd>)>& fun,
                  Eigen::Ref<Eigen::VectorXd> x2) override
    {
        if (x1.size() != _k1.size()) initialize(x1.size());

        fun(x1, _k1);  // k1
        _k1 *= dt;
        fun(x1 + 2.0 * _k1 / 27.0, _k2);  // k2
        _k2 *= dt;
        fun(x1 + (_k1 + 3.0 * _k2) / 36.0, _k3);  // k3
        _k3 *= dt;
        fun(x1 + (_k1 + 3.0 * _k3) / 24.0, _k4);  // k4
        _k4 *= dt;
        fun(x1 + (80.0 * _k1 - 300.0 * _k3 + 300.0 * _k4) / 192.0, _k5);  // k5
        _k5 *= dt;
        fun(x1 + (_k1 + 5.0 * _k4 + 4.0 * _k5) / 20.0, _k6);  // k6
        _k6 *= dt;
        fun(x1 + (-25.0 * _k1 + 125.0 * _k4 - 260.0 * _k5 + 250.0 * _k6) / 108.0, _k7);  // k7
        _k7 *= dt;
        fun(x1 + (93.0 * _k1 + 244.0 * _k5 - 200.0 * _k6 + 13.0 * _k7) / 900.0, _k8);  // k8
        _k8 *= dt;
        fun(x1 + (1080.0 * _k1 - 4770.0 * _k4 + 8448.0 * _k5 - 6420.0 * _k6 + 402.0 * _k7 + 1620.0 * _k8) / 540.0, _k9);  // k9
        _k9 *= dt;
        fun(x1 + (-12285.0 * _k1 + 3105.0 * _k4 - 105408.0 * _k5 + 83970.0 * _k6 - 4617.0 * _k7 + 41310.0 * _k8 - 1215.0 * _k9) / 14580.0,
            _k10);  // k10
        _k10 *= dt;
        fun(x1 + (2383.0 * _k1 - 8525.0 * _k4 + 17984.0 * _k5 - 15050.0 * _k6 + 2133.0 * _k7 + 2250.0 * _k8 + 1125.0 * _k9 + 1800.0 * _k10) / 4100.0,
            _k11);  // k11
        _k11 *= dt;

        x2 = x1 + (41.0 * _k1 + 272.0 * _k6 + 216.0 * _k7 + 216.0 * _k8 + 27.0 * _k9 + 27.0 * _k10 + 41.0 * _k11) / 840.0;
    }

    // Implements interface method
    void solveIVP(const StateVector& x1, const InputVector& u1, double dt, const SystemDynamicsInterface& system,
                  Eigen::Ref<Eigen::VectorXd> x2) override
    {
        if (x1.size() != _k1.size()) initialize(x1.size());

        system.dynamics(x1, u1, _k1);  // k1
        _k1 *= dt;
        system.dynamics(x1 + 2.0 * _k1 / 27.0, u1, _k2);  // k2
        _k2 *= dt;
        system.dynamics(x1 + (_k1 + 3.0 * _k2) / 36.0, u1, _k3);  // k3
        _k3 *= dt;
        system.dynamics(x1 + (_k1 + 3.0 * _k3) / 24.0, u1, _k4);  // k4
        _k4 *= dt;
        system.dynamics(x1 + (80.0 * _k1 - 300.0 * _k3 + 300.0 * _k4) / 192.0, u1, _k5);  // k5
        _k5 *= dt;
        system.dynamics(x1 + (_k1 + 5.0 * _k4 + 4.0 * _k5) / 20.0, u1, _k6);  // k6
        _k6 *= dt;
        system.dynamics(x1 + (-25.0 * _k1 + 125.0 * _k4 - 260.0 * _k5 + 250.0 * _k6) / 108.0, u1, _k7);  // k7
        _k7 *= dt;
        system.dynamics(x1 + (93.0 * _k1 + 244.0 * _k5 - 200.0 * _k6 + 13.0 * _k7) / 900.0, u1, _k8);  // k8
        _k8 *= dt;
        system.dynamics(x1 + (12.0 * _k1 - 53.0 * _k4) / 6.0 + (1408.0 * _k5 - 1070.0 * _k6 + 67.0 * _k7 + 270.0 * _k8) / 90.0, u1, _k9);  // k9
        _k9 *= dt;
        system.dynamics(x1 + (-12285.0 * _k1 + 3105.0 * _k4 - 105408.0 * _k5 + 83970.0 * _k6 - 4617.0 * _k7 + 41310.0 * _k8 - 1215.0 * _k9) / 14580.0,
                        u1, _k10);  // k10
        _k10 *= dt;
        system.dynamics(
            x1 + (2383.0 * _k1 - 8525.0 * _k4 + 17984.0 * _k5 - 15050.0 * _k6 + 2133.0 * _k7 + 2250.0 * _k8 + 1125.0 * _k9 + 1800.0 * _k10) / 4100.0,
            u1, _k11);  // k11
        _k11 *= dt;

        x2 = x1 + (41.0 * _k1 + 272.0 * _k6 + 216.0 * _k7 + 216.0 * _k8 + 27.0 * _k9 + 27.0 * _k10 + 41.0 * _k11) / 840.0;
    }

#ifdef MESSAGE_SUPPORT
    // Implements interface method
    void toMessage(messages::ExplicitIntegrator& message) const override
    {
        NumericalIntegratorExplicitInterface::toMessage(message);
        message.mutable_runge_kutta_7();
    }
#endif

 private:
    Eigen::VectorXd _k1;
    Eigen::VectorXd _k2;
    Eigen::VectorXd _k3;
    Eigen::VectorXd _k4;
    Eigen::VectorXd _k5;
    Eigen::VectorXd _k6;
    Eigen::VectorXd _k7;
    Eigen::VectorXd _k8;
    Eigen::VectorXd _k9;
    Eigen::VectorXd _k10;
    Eigen::VectorXd _k11;
};

FACTORY_REGISTER_EXPLICIT_INTEGRATOR(IntegratorExplicitRungeKutta7)

/**
 * @brief Adaptive-Step-Size-Control
 *
 * @ingroup numerics
 *
 * For more information please refer to
 * Lars Grüne, Jürgen Pannek: "Nonlinear Model Predictive Control" - chapter 9.3: "Adaptive Step Size Control -
 * Algorithm 9.7
 *
 * @see NumericalIntegratorInterface NumericalIntegratorExplicitInterface IntegratorAdaptiveStepSize
 *
 * @author Alexander Westermann (alexander.westermann@tu-dortmund.de)
 */
class IntegratorAdaptiveStepSize : public NumericalIntegratorExplicitInterface
{
 public:
    // Implements interface method
    DynamicsEvalInterface::Ptr getInstance() const override { return std::make_shared<IntegratorAdaptiveStepSize>(); }

    // Implements interface method
    int getConvergenceOrder() const override { return _p2; }

    void initialize(int state_dim) override
    {
        _integrator1->initialize(state_dim);
        _integrator2->initialize(state_dim);
    }

    // Implements interface method
    void solveIVP(const Eigen::VectorXd& x1, double dt, const std::function<void(const Eigen::VectorXd&, Eigen::Ref<Eigen::VectorXd>)>& fun,
                  Eigen::Ref<Eigen::VectorXd> x2) override
    {
        // start-values:
        _tau   = 0.0;
        _tau2  = 0.0;
        _h_new = dt;

        // initialization
        _help2 = x1;
        x2     = x1;
        _help  = _help2;

        // start-loop:
        while (std::abs(_tau - dt) > std::numeric_limits<double>::epsilon())
        {
            _h = _h_new;
            if ((_tau + _h) > dt)
            {
                _h = dt - _tau;
            }

            while (true)
            {
                _tau2 = _tau + _h;

                // solve the integration-problem
                _integrator1->solveIVP(_help2, _h, fun, _help);
                _integrator2->solveIVP(_help2, _h, fun, x2);

                // evaluating the error in the in the discrete  l²-norm
                _epsilon_hat = 0.0;
                _help        = _help - x2;
                for (int k = 0; k < _help.size(); k++)
                {
                    _epsilon_hat = _epsilon_hat + pow(_help(k), 2.0);
                }
                _epsilon_hat = std::pow(_epsilon_hat, 0.5);

                // evaluation of the new step-size
                _h_new = pow(0.9 * (_tol / _epsilon_hat), (1.0 / (double(_p1) + 1.0))) * _h;
                if (_epsilon_hat > _tol)
                {
                    _h = _h_new;
                }
                else
                {
                    break;
                }
            }

            // update parameters
            _tau   = _tau2;
            _help2 = x2;
        }
    }

    // Implements interface method
    void solveIVP(const StateVector& x1, const InputVector& u1, double dt, const SystemDynamicsInterface& system,
                  Eigen::Ref<Eigen::VectorXd> x2) override
    {
        // start-values:
        _tau   = 0.0;
        _tau2  = 0.0;
        _h_new = dt;

        // initialization
        _help2 = x1;
        x2     = x1;
        _help  = _help2;

        // start-loop:
        while (std::abs(_tau - dt) > std::numeric_limits<double>::epsilon())
        {
            _h = _h_new;
            if ((_tau + _h) > dt)
            {
                _h = dt - _tau;
            }

            while (true)
            {
                _tau2 = _tau + _h;

                // solve the integration-problem
                _integrator1->solveIVP(_help2, u1, _h, system, _help);
                _integrator2->solveIVP(_help2, u1, _h, system, x2);

                // evaluating the error in the in the discrete  l²-norm
                _epsilon_hat = 0.0;
                _help        = _help - x2;
                for (int k = 0; k < _help.size(); k++)
                {
                    _epsilon_hat = _epsilon_hat + pow(_help(k), 2.0);
                }
                _epsilon_hat = pow(_epsilon_hat, 0.5);

                // evaluation of the new step-size
                _h_new = std::pow(0.9 * (_tol / _epsilon_hat), (1.0 / (double(_p1) + 1.0))) * _h;
                if (_epsilon_hat > _tol)
                {
                    _h = _h_new;
                }
                else
                {
                    break;
                }
            }

            // update parameters
            _tau   = _tau2;
            _help2 = x2;
        }
    }

#ifdef MESSAGE_SUPPORT
    void toMessage(messages::IntegratorAdaptiveStepSize& message) const
    {
        message.set_tol(_tol);
        if (_integrator1) _integrator1->toMessage(*message.mutable_integrator1());
        if (_integrator2) _integrator2->toMessage(*message.mutable_integrator2());
    }
    void fromMessage(const messages::IntegratorAdaptiveStepSize& message, std::stringstream* issues = nullptr)
    {
        _tol = message.tol();
        // create integrator1
        std::string type;
        util::get_oneof_field_type(message.integrator1(), "explicit_integrator", type, false);
        NumericalIntegratorExplicitInterface::Ptr int1 = create_from_factory<NumericalIntegratorExplicitInterface>(type);
        // import parameters
        if (int1)
        {
            int1->fromMessage(message.integrator1(), issues);
            _integrator1 = int1;
            _p1          = int1->getConvergenceOrder();
        }
        else
        {
            if (issues) *issues << "IntegratorAdaptiveStepSize: unknown integrator1 specified.\n";
            return;
        }
        // create integrator2
        type = "";
        util::get_oneof_field_type(message.integrator2(), "explicit_integrator", type, false);
        NumericalIntegratorExplicitInterface::Ptr int2 = create_from_factory<NumericalIntegratorExplicitInterface>(type);
        // import parameters
        if (int2)
        {
            int2->fromMessage(message.integrator2(), issues);
            _integrator2 = int2;
            _p2          = int2->getConvergenceOrder();
        }
        else
        {
            if (issues) *issues << "IntegratorAdaptiveStepSize: unknown integrator2 specified.\n";
            return;
        }
        // error-message if _p1 >= _p2
        if (_p1 >= _p2 && issues)
            *issues << "The chosen methods are not applicable! The order of the first method has to be lower than the order of the second method!\n";
    }

    void toMessage(messages::ExplicitIntegrator& message) const override
    {
        NumericalIntegratorExplicitInterface::toMessage(message);
        toMessage(*message.mutable_adaptive_step_size());
    }
    void fromMessage(const messages::DynamicsEval& message, std::stringstream* issues = nullptr) override
    {
        NumericalIntegratorExplicitInterface::fromMessage(message, issues);
        fromMessage(message.integrator().adaptive_step_size(), issues);
    }
#endif
 protected:
    double _tol;
    NumericalIntegratorExplicitInterface::Ptr _integrator1 = std::make_shared<IntegratorExplicitRungeKutta2>();
    NumericalIntegratorExplicitInterface::Ptr _integrator2 = std::make_shared<IntegratorExplicitRungeKutta7>();

    int _p1 = 0;
    int _p2 = 0;

 private:
    double _epsilon_hat;
    double _tau;
    double _tau2;
    double _h;
    double _h_new;
    Eigen::VectorXd _help;
    Eigen::VectorXd _help2;
};

FACTORY_REGISTER_EXPLICIT_INTEGRATOR(IntegratorAdaptiveStepSize)

/**
 * @brief Multi-stage integration fixed step
 *
 * @ingroup numerics
 *
 * @see NumericalIntegratorInterface NumericalIntegratorExplicitInterface
 *
 * @author Christoph Rösmann (christoph.roesmann@tu-dortmund.de)
 */
class IntegratorMultiStageFixedStep : public NumericalIntegratorExplicitInterface
{
 public:
    // Implements interface method
    DynamicsEvalInterface::Ptr getInstance() const override { return std::make_shared<IntegratorMultiStageFixedStep>(); }

    // Implements interface method
    int getConvergenceOrder() const override { return _integrator->getConvergenceOrder(); }

    // Implements interface method
    void solveIVP(const StateVector& x1, const InputVector& u1, double dt, const SystemDynamicsInterface& system,
                  Eigen::Ref<Eigen::VectorXd> x2) override
    {
        assert(_integrator);
        if (dt <= _inner_integration_dt || _inner_integration_dt <= 0)
        {
            _integrator->solveIVP(x1, u1, dt, system, x2);
            return;
        }

        int n               = dt / _inner_integration_dt;
        double dt_remainder = std::fmod(dt, _inner_integration_dt);

        _x_sub.resize(x1.size());
        _integrator->solveIVP(x1, u1, _inner_integration_dt, system, _x_sub);
        for (int i = 1; i < n; ++i)
        {
            _integrator->solveIVP(_x_sub, u1, _inner_integration_dt, system, x2);
            _x_sub = x2;
        }
        if (dt_remainder > 0) _integrator->solveIVP(_x_sub, u1, dt_remainder, system, x2);
    }

    // Implements interface method
    void solveIVP(const Eigen::VectorXd& x1, double dt, const UnaryFunction& fun, Eigen::Ref<Eigen::VectorXd> x2) override
    {
        assert(_integrator);
        if (dt <= _inner_integration_dt || _inner_integration_dt <= 0)
        {
            _integrator->solveIVP(x1, dt, fun, x2);
            return;
        }

        int n               = dt / _inner_integration_dt;
        double dt_remainder = std::fmod(dt, _inner_integration_dt);

        _x_sub.resize(x1.size());
        _integrator->solveIVP(x1, _inner_integration_dt, fun, _x_sub);
        for (int i = 1; i < n; ++i)
        {
            _integrator->solveIVP(_x_sub, _inner_integration_dt, fun, x2);
            _x_sub = x2;
        }
        if (dt_remainder > 1e-8)
            _integrator->solveIVP(_x_sub, dt_remainder, fun, x2);  // TODO(roesmann): maybe use std::numerical_limits::epsilon or round_error
    }

    void setSubIntegrator(NumericalIntegratorExplicitInterface::Ptr integrator) { _integrator = integrator; }

#ifdef MESSAGE_SUPPORT
    void fromMessage(const messages::IntegratorMultiStageFixedStep& message, std::stringstream* issues = nullptr)
    {
        _inner_integration_dt = message.inner_integration_dt();

        // create integrator
        std::string type;
        util::get_oneof_field_type(message.integrator(), "explicit_integrator", type, false);
        NumericalIntegratorExplicitInterface::Ptr integrator = create_from_factory<NumericalIntegratorExplicitInterface>(type);
        // import parameters
        if (integrator)
        {
            integrator->fromMessage(message.integrator(), issues);
            _integrator = integrator;
        }
        else
        {
            if (issues) *issues << "IntegratorMultiStageFixedStep: unknown integrator specified.\n";
            return;
        }
    }

    void toMessage(messages::IntegratorMultiStageFixedStep& message) const
    {
        message.set_inner_integration_dt(_inner_integration_dt);

        if (_integrator) _integrator->toMessage(*message.mutable_integrator());
    }

    // Implements interface method
    void fromMessage(const messages::ExplicitIntegrator& message, std::stringstream* issues = nullptr) override
    {
        NumericalIntegratorExplicitInterface::fromMessage(message);
        fromMessage(message.multi_stage_fixed_step());
    }

    // Implements interface method
    void toMessage(messages::ExplicitIntegrator& message) const override
    {
        NumericalIntegratorExplicitInterface::toMessage(message);
        toMessage(*message.mutable_multi_stage_fixed_step());
    }
#endif

 private:
    NumericalIntegratorExplicitInterface::Ptr _integrator;
    double _inner_integration_dt = 0.1;

    Eigen::VectorXd _x_sub;
};

FACTORY_REGISTER_EXPLICIT_INTEGRATOR(IntegratorMultiStageFixedStep)

/**
 * @brief Multi-stage integration scaled
 *
 * @ingroup numerics
 *
 * @see NumericalIntegratorInterface NumericalIntegratorExplicitInterface
 *
 * @author Christoph Rösmann (christoph.roesmann@tu-dortmund.de)
 */
class IntegratorMultiStageScaled : public NumericalIntegratorExplicitInterface
{
 public:
    // Implements interface method
    DynamicsEvalInterface::Ptr getInstance() const override { return std::make_shared<IntegratorMultiStageScaled>(); }

    // Implements interface method
    int getConvergenceOrder() const override { return _integrator->getConvergenceOrder(); }

    // Implements interface method
    void solveIVP(const StateVector& x1, const InputVector& u1, double dt, const SystemDynamicsInterface& system,
                  Eigen::Ref<Eigen::VectorXd> x2) override
    {
        assert(_n > 0);

        double inner_dt = dt / (double)_n;

        _x_sub.resize(x1.size());
        _integrator->solveIVP(x1, u1, inner_dt, system, _x_sub);
        for (int i = 1; i < _n; ++i)
        {
            _integrator->solveIVP(_x_sub, u1, inner_dt, system, x2);
            _x_sub = x2;
        }
    }

    // Implements interface method
    void solveIVP(const Eigen::VectorXd& x1, double dt, const UnaryFunction& fun, Eigen::Ref<Eigen::VectorXd> x2) override
    {
        assert(_n > 0);

        double inner_dt = dt / (double)_n;

        _x_sub.resize(x1.size());
        _integrator->solveIVP(x1, inner_dt, fun, _x_sub);
        for (int i = 1; i < _n; ++i)
        {
            _integrator->solveIVP(_x_sub, inner_dt, fun, x2);
            _x_sub = x2;
        }
    }

    void setSubIntegrator(NumericalIntegratorExplicitInterface::Ptr integrator) { _integrator = integrator; }
    void setN(int n) { _n = n; }

#ifdef MESSAGE_SUPPORT
    void fromMessage(const messages::IntegratorMultiStageScaled& message, std::stringstream* issues = nullptr)
    {
        _n = message.n();

        // create integrator
        std::string type;
        util::get_oneof_field_type(message.integrator(), "explicit_integrator", type, false);
        NumericalIntegratorExplicitInterface::Ptr integrator = create_from_factory<NumericalIntegratorExplicitInterface>(type);
        // import parameters
        if (integrator)
        {
            integrator->fromMessage(message.integrator(), issues);
            _integrator = integrator;
        }
        else
        {
            if (issues) *issues << "IntegratorMultiStageScaled: unknown integrator specified.\n";
            return;
        }
    }

    void toMessage(messages::IntegratorMultiStageScaled& message) const
    {
        message.set_n(_n);

        if (_integrator) _integrator->toMessage(*message.mutable_integrator());
    }

    // Implements interface method
    void fromMessage(const messages::ExplicitIntegrator& message, std::stringstream* issues = nullptr) override
    {
        NumericalIntegratorExplicitInterface::fromMessage(message);
        fromMessage(message.multi_stage_scaled());
    }

    // Implements interface method
    void toMessage(messages::ExplicitIntegrator& message) const override
    {
        NumericalIntegratorExplicitInterface::toMessage(message);
        toMessage(*message.mutable_multi_stage_scaled());
    }
#endif

 private:
    NumericalIntegratorExplicitInterface::Ptr _integrator;
    int _n = 10;

    Eigen::VectorXd _x_sub;
};

FACTORY_REGISTER_EXPLICIT_INTEGRATOR(IntegratorMultiStageScaled)

}  // namespace corbo

#endif  // SRC_NUMERICS_INCLUDE_CORBO_NUMERICS_EXPLICIT_INTEGRATORS_H_
