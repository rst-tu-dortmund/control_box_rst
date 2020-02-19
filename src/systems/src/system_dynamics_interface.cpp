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

#include <corbo-systems/system_dynamics_interface.h>

#include <corbo-numerics/finite_differences.h>

namespace corbo {

SystemDynamicsInterface::SystemDynamicsInterface() : _linearization_method(std::make_shared<ForwardDifferences>()) {}

void SystemDynamicsInterface::getLinearA(const StateVector& x0, const ControlVector& u0, Eigen::MatrixXd& A) const
{
    assert(getStateDimension() == x0.size());
    assert(getInputDimension() == u0.size());
    assert(A.rows() == x0.size());
    assert(A.cols() == x0.size());

    StateVector x(x0);

    auto inc  = [&x](int idx, double inc) { x[idx] += inc; };
    auto eval = [&, this](StateVector& values) { dynamics(x, u0, values); };
    _linearization_method->computeJacobian2(inc, eval, A);
}

void SystemDynamicsInterface::getLinearB(const StateVector& x0, const ControlVector& u0, Eigen::MatrixXd& B) const
{
    assert(getStateDimension() == x0.size());
    assert(getInputDimension() == u0.size());
    assert(B.rows() == x0.size());
    assert(B.cols() == u0.size());

    ControlVector u(u0);

    auto inc  = [&u](int idx, double inc) { u[idx] += inc; };
    auto eval = [&, this](StateVector& values) { dynamics(x0, u, values); };
    _linearization_method->computeJacobian2(inc, eval, B);
}

void SystemDynamicsInterface::setLinearizationMethod(std::shared_ptr<FiniteDifferencesInterface> lin_method) { _linearization_method = lin_method; }

#ifdef MESSAGE_SUPPORT
void SystemDynamicsInterface::toMessage(corbo::messages::SystemDynamics& message) const
{
    message.set_deadtime(_deadtime);
    _linearization_method->toMessage(*message.mutable_linearization_method());
}

void SystemDynamicsInterface::fromMessage(const corbo::messages::SystemDynamics& message, std::stringstream* issues)
{
    if (message.deadtime() >= 0)
    {
        _deadtime = message.deadtime();
    }
    else
    {
        if (issues) *issues << "Deadtime must be >= 0" << std::endl;
        _deadtime = 0;
        PRINT_ERROR("Deadtime must be >= 0");
    }
    _linearization_method->fromMessage(message.linearization_method(), issues);
}
#endif

}  // namespace corbo
