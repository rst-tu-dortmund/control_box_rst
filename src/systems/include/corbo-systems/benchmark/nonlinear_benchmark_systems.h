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

#ifndef SRC_SYSTEMS_INCLUDE_CORBO_SYSTEMS_BENCHMARK_NONLINEAR_BENCHMARK_SYSTEMS_H_
#define SRC_SYSTEMS_INCLUDE_CORBO_SYSTEMS_BENCHMARK_NONLINEAR_BENCHMARK_SYSTEMS_H_

#include <corbo-systems/system_dynamics_interface.h>

namespace corbo {

class VanDerPolOscillator : public SystemDynamicsInterface
{
 public:
    //! Default constructor
    VanDerPolOscillator() {}

    // implements interface method
    Ptr getInstance() const override { return std::make_shared<VanDerPolOscillator>(); }

    // implements interface method
    bool isContinuousTime() const override { return true; }
    // implements interface method
    bool isLinear() const override { return false; }

    // implements interface method
    int getInputDimension() const override { return 1; }
    // implements interface method
    int getStateDimension() const override { return 2; }

    // implements interface method
    void dynamics(const Eigen::Ref<const StateVector>& x, const Eigen::Ref<const ControlVector>& u, Eigen::Ref<StateVector> f) const override
    {
        assert(x.size() == getStateDimension());
        assert(u.size() == getInputDimension());
        assert(x.size() == f.size() && "VanDerPolOscillator::dynamics(): x and f are not of the same size, do not forget to pre-allocate f.");

        f[0] = x[1];  // xdot = xdot
        f[1] = -_a * (x[0] * x[0] - 1) * x[1] - x[0] + u[0];
    }

    // access parameters
    const double& getDampingCoefficient() const { return _a; }
    void setDampingCoefficient(double a) { _a = a; }

#ifdef MESSAGE_SUPPORT
    // implements interface method
    void toMessage(messages::SystemDynamics& message) const override
    {
        SystemDynamicsInterface::toMessage(message);

        message.mutable_van_der_pol_oscillator()->set_a(_a);
    }
    // implements interface method
    void fromMessage(const messages::SystemDynamics& message, std::stringstream* issues) override
    {
        SystemDynamicsInterface::fromMessage(message, issues);

        _a = message.van_der_pol_oscillator().a();
    }
#endif

 private:
    double _a = 1;
};
FACTORY_REGISTER_SYSTEM_DYNAMICS(VanDerPolOscillator)

class DuffingOscillator : public SystemDynamicsInterface
{
 public:
    //! Default constructor
    DuffingOscillator() {}

    // implements interface method
    Ptr getInstance() const override { return std::make_shared<DuffingOscillator>(); }

    // implements interface method
    bool isContinuousTime() const override { return true; }
    // implements interface method
    bool isLinear() const override { return false; }

    // implements interface method
    int getInputDimension() const override { return 1; }
    // implements interface method
    int getStateDimension() const override { return 2; }

    // implements interface method
    void dynamics(const Eigen::Ref<const StateVector>& x, const Eigen::Ref<const ControlVector>& u, Eigen::Ref<StateVector> f) const override
    {
        assert(x.size() == getStateDimension());
        assert(u.size() == getInputDimension());
        assert(x.size() == f.size() && "DuffingOscillator::dynamics(): x and f are not of the same size, do not forget to pre-allocate f.");

        f[0] = x[1];  // xdot = xdot
        f[1] = -_damping * x[1] - _spring_alpha * x[0] - _spring_beta * x[0] * x[0] * x[0] + u[0];
    }

    // access parameters
    void setParameters(double damping, double spring_alpha, double spring_beta)
    {
        _damping      = damping;
        _spring_alpha = spring_alpha;
        _spring_beta  = spring_beta;
    }

#ifdef MESSAGE_SUPPORT
    // implements interface method
    void toMessage(messages::SystemDynamics& message) const override
    {
        SystemDynamicsInterface::toMessage(message);

        message.mutable_duffing_oscillator()->set_damping(_damping);
        message.mutable_duffing_oscillator()->set_spring_alpha(_spring_alpha);
        message.mutable_duffing_oscillator()->set_spring_beta(_spring_beta);
    }
    // implements interface method
    void fromMessage(const messages::SystemDynamics& message, std::stringstream* issues) override
    {
        SystemDynamicsInterface::fromMessage(message, issues);

        _damping      = message.duffing_oscillator().damping();
        _spring_alpha = message.duffing_oscillator().spring_alpha();
        _spring_beta  = message.duffing_oscillator().spring_beta();
    }
#endif

 private:
    double _damping      = 1;
    double _spring_alpha = 1;
    double _spring_beta  = 1;
};
FACTORY_REGISTER_SYSTEM_DYNAMICS(DuffingOscillator)

class FreeSpaceRocket : public SystemDynamicsInterface
{
 public:
    //! Default constructor
    FreeSpaceRocket() {}

    // implements interface method
    Ptr getInstance() const override { return std::make_shared<FreeSpaceRocket>(); }

    // implements interface method
    bool isContinuousTime() const override { return true; }
    // implements interface method
    bool isLinear() const override { return false; }

    // implements interface method
    int getInputDimension() const override { return 1; }
    // implements interface method
    int getStateDimension() const override { return 3; }

    // implements interface method
    void dynamics(const Eigen::Ref<const StateVector>& x, const Eigen::Ref<const ControlVector>& u, Eigen::Ref<StateVector> f) const override
    {
        assert(x.size() == getStateDimension());
        assert(u.size() == getInputDimension());
        assert(x.size() == f.size() && "FreeSpaceRocket::dynamics(): x and f are not of the same size, do not forget to pre-allocate f.");

        f[0] = x[1];                                // sdot = v
        f[1] = (u[0] - 0.02 * x[1] * x[1]) / x[2];  // (u-0.02*v^2)/m
        f[2] = -0.01 * u[0] * u[0];                 // -0.01 * u^2
    }
};
FACTORY_REGISTER_SYSTEM_DYNAMICS(FreeSpaceRocket)

class SimplePendulum : public SystemDynamicsInterface
{
 public:
    //! Default constructor
    SimplePendulum() {}

    // implements interface method
    Ptr getInstance() const override { return std::make_shared<SimplePendulum>(); }

    // implements interface method
    bool isContinuousTime() const override { return true; }
    // implements interface method
    bool isLinear() const override { return false; }

    // implements interface method
    int getInputDimension() const override { return 1; }
    // implements interface method
    int getStateDimension() const override { return 2; }

    // implements interface method
    void dynamics(const Eigen::Ref<const StateVector>& x, const Eigen::Ref<const ControlVector>& u, Eigen::Ref<StateVector> f) const override
    {
        assert(u.size() == getInputDimension());
        assert(x.size() == getStateDimension());
        assert(x.size() == f.size() && "SimplePendulum::dynamics(): x and f are not of the same size, do not forget to pre-allocate f.");

        f[0] = x[1];  // phidot = phidot
        f[1] = u[0] - _rho / (_m * _l * _l) * x[1] - _g / _l * std::sin(x[0]);
    }

    // access parameters
    void setParameters(double mass, double length, double gravitation, double friction)
    {
        _m   = mass;
        _l   = length;
        _g   = gravitation;
        _rho = friction;
    }

#ifdef MESSAGE_SUPPORT
    // implements interface method
    void toMessage(messages::SystemDynamics& message) const override
    {
        SystemDynamicsInterface::toMessage(message);

        messages::SimplePendulum* msg = message.mutable_simple_pendulum();

        msg->set_m(_m);
        msg->set_l(_l);
        msg->set_g(_g);
        msg->set_rho(_rho);
    }
    // implements interface method
    void fromMessage(const messages::SystemDynamics& message, std::stringstream* issues) override
    {
        SystemDynamicsInterface::fromMessage(message, issues);

        const messages::SimplePendulum& msg = message.simple_pendulum();

        _m   = msg.m();
        _l   = msg.l();
        _g   = msg.g();
        _rho = msg.rho();
    }
#endif

 private:
    double _m   = 0.205;
    double _l   = 0.34;
    double _g   = 9.81;
    double _rho = 0;
};
FACTORY_REGISTER_SYSTEM_DYNAMICS(SimplePendulum)

class MasslessPendulum : public SystemDynamicsInterface
{
 public:
    //! Default constructor
    MasslessPendulum() {}

    // implements interface method
    Ptr getInstance() const override { return std::make_shared<MasslessPendulum>(); }

    // implements interface method
    bool isContinuousTime() const override { return true; }
    // implements interface method
    bool isLinear() const override { return false; }

    // implements interface method
    int getInputDimension() const override { return 1; }
    // implements interface method
    int getStateDimension() const override { return 2; }

    // implements interface method
    void dynamics(const Eigen::Ref<const StateVector>& x, const Eigen::Ref<const ControlVector>& u, Eigen::Ref<StateVector> f) const override
    {
        assert(u.size() == getInputDimension());
        assert(x.size() == getStateDimension());
        assert(x.size() == f.size() && "MasslessPendulum::dynamics(): x and f are not of the same size, do not forget to pre-allocate f.");

        f[0] = x[1];  // phidot = phidot
        f[1] = u[0] - _omega0 * std::sin(x[0]);
    }

    // access parameters
    void setParameter(double omega0) { _omega0 = omega0; }

#ifdef MESSAGE_SUPPORT
    // implements interface method
    void toMessage(messages::SystemDynamics& message) const override
    {
        SystemDynamicsInterface::toMessage(message);
        messages::MasslessPendulum* msg = message.mutable_massless_pendulum();
        msg->set_omega0(_omega0);
    }
    // implements interface method
    void fromMessage(const messages::SystemDynamics& message, std::stringstream* issues) override
    {
        SystemDynamicsInterface::fromMessage(message, issues);

        const messages::MasslessPendulum& msg = message.massless_pendulum();
        _omega0                               = msg.omega0();
    }
#endif

 private:
    double _omega0 = 1.0;
};
FACTORY_REGISTER_SYSTEM_DYNAMICS(MasslessPendulum)

class CartPole : public SystemDynamicsInterface
{
 public:
    //! Default constructor
    CartPole() {}

    // implements interface method
    Ptr getInstance() const override { return std::make_shared<CartPole>(); }

    // implements interface method
    bool isContinuousTime() const override { return true; }
    // implements interface method
    bool isLinear() const override { return false; }

    // implements interface method
    int getInputDimension() const override { return 1; }
    // implements interface method
    int getStateDimension() const override { return 4; }

    // implements interface method
    void dynamics(const Eigen::Ref<const StateVector>& x, const Eigen::Ref<const ControlVector>& u, Eigen::Ref<StateVector> f) const override
    {
        assert(u.size() == getInputDimension());
        assert(x.size() == getStateDimension());
        assert(x.size() == f.size() && "CartPoleSystem::dynamics(): x and f are not of the same size, do not forget to pre-allocate f.");

        // [x phi xdot phidot]

        double sin_phi_phidot_sq = std::sin(x[1]) * x[3] * x[3];
        double denum             = _mc + _mp * (1 - std::pow(std::cos(x[1]), 2));

        f[0] = x[2];  // xdot = xdot
        f[1] = x[3];  // phidot = phidot
        f[2] = (_l * _mp * sin_phi_phidot_sq + u[0] + _mp * _g * std::cos(x[1]) * std::sin(x[1])) / denum;
        f[3] = -(_l * _mp * std::cos(x[1]) * sin_phi_phidot_sq + u[0] * std::cos(x[1]) + (_mp + _mc) * _g * std::sin(x[1])) / (_l * denum);
    }

// access parameters
//    void setParameters(double mass, double length, double gravitation, double friction)
//    {
//        _m   = mass;
//        _l   = length;
//        _g   = gravitation;
//        _rho = friction;
//    }

#ifdef MESSAGE_SUPPORT
    // implements interface method
    void toMessage(messages::SystemDynamics& message) const override
    {
        SystemDynamicsInterface::toMessage(message);

        messages::CartPole* msg = message.mutable_cart_pole();

        msg->set_mc(_mc);
        msg->set_mp(_mp);
        msg->set_l(_l);
        msg->set_g(_g);
    }
    // implements interface method
    void fromMessage(const messages::SystemDynamics& message, std::stringstream* issues) override
    {
        SystemDynamicsInterface::fromMessage(message, issues);

        const messages::CartPole& msg = message.cart_pole();

        _mc = msg.mc();
        _mp = msg.mp();
        _l  = msg.l();
        _g  = msg.g();
    }
#endif

 private:
    double _mc = 1.0;
    double _mp = 0.3;
    double _l  = 0.5;
    double _g  = 9.81;

    // suggested u_max: 20N
    // x_max = 2 m
};
FACTORY_REGISTER_SYSTEM_DYNAMICS(CartPole)

// Taken from:
// Verschueren et al., "A Stabilizing Nonlinear Model Predictive Control Scheme for Time-optimal Point-to-point Motions", CDC, 2017.
// Chen and Allgoewer, "A quasi-infinite horizon nonlinear model predictive control scheme with guaranteed stability", Automatica, 2002.
// Recommended bounds: -10 < u < 10
// N = 50
class ToyExample : public SystemDynamicsInterface
{
 public:
    //! Default constructor
    ToyExample() {}

    // implements interface method
    Ptr getInstance() const override { return std::make_shared<ToyExample>(); }

    // implements interface method
    bool isContinuousTime() const override { return true; }
    // implements interface method
    bool isLinear() const override { return false; }

    // implements interface method
    int getInputDimension() const override { return 1; }
    // implements interface method
    int getStateDimension() const override { return 2; }

    // implements interface method
    void dynamics(const Eigen::Ref<const StateVector>& x, const Eigen::Ref<const ControlVector>& u, Eigen::Ref<StateVector> f) const override
    {
        assert(x.size() == getStateDimension());
        assert(u.size() == getInputDimension());
        assert(x.size() == f.size() && "ToyExample::dynamics(): x and f are not of the same size, do not forget to pre-allocate f.");

        // x = [p; q]

        f[0] = x[1] + u[0] * (_mu + (1.0 - _mu) * x[0]);        // pdot = q + u(mu + (1-mu)p)
        f[1] = x[0] + u[0] * (_mu - 4.0 * (1.0 - _mu) * x[1]);  // qdot = p + u(mu-4(1-mu)q)
    }

    // access parameters
    void setParameters(double mu) { _mu = mu; }

#ifdef MESSAGE_SUPPORT
    // implements interface method
    void toMessage(messages::SystemDynamics& message) const override
    {
        SystemDynamicsInterface::toMessage(message);

        message.mutable_toy_example()->set_mu(_mu);
    }
    // implements interface method
    void fromMessage(const messages::SystemDynamics& message, std::stringstream* issues) override
    {
        SystemDynamicsInterface::fromMessage(message, issues);

        _mu = message.toy_example().mu();
    }
#endif

 private:
    double _mu = 0.5;
};
FACTORY_REGISTER_SYSTEM_DYNAMICS(ToyExample)

class ArtsteinsCircle : public SystemDynamicsInterface
{
 public:
    //! Default constructor
    ArtsteinsCircle() {}

    // implements interface method
    Ptr getInstance() const override { return std::make_shared<ArtsteinsCircle>(); }

    // implements interface method
    bool isContinuousTime() const override { return true; }
    // implements interface method
    bool isLinear() const override { return false; }

    // implements interface method
    int getInputDimension() const override { return 1; }
    // implements interface method
    int getStateDimension() const override { return 2; }

    // implements interface method
    void dynamics(const Eigen::Ref<const StateVector>& x, const Eigen::Ref<const ControlVector>& u, Eigen::Ref<StateVector> f) const override
    {
        assert(x.size() == getStateDimension());
        assert(u.size() == getInputDimension());
        assert(x.size() == f.size() && "ArtsteinsCircle::dynamics(): x and f are not of the same size, do not forget to pre-allocate f.");

        f[0] = (x[0] * x[0] - x[1] * x[1]) * u[0];
        f[1] = 2 * x[0] * x[1] * u[0];
    }

#ifdef MESSAGE_SUPPORT
    // implements interface method
    void toMessage(messages::SystemDynamics& message) const override
    {
        SystemDynamicsInterface::toMessage(message);

        message.mutable_artsteins_circle();
    }
    // implements interface method
    void fromMessage(const messages::SystemDynamics& message, std::stringstream* issues) override
    {
        SystemDynamicsInterface::fromMessage(message, issues);
    }
#endif
};
FACTORY_REGISTER_SYSTEM_DYNAMICS(ArtsteinsCircle)

}  // namespace corbo

#endif  // SRC_SYSTEMS_INCLUDE_CORBO_SYSTEMS_BENCHMARK_NONLINEAR_BENCHMARK_SYSTEMS_H_
