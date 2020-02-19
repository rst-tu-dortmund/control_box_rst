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

#include <corbo-communication/utilities.h>
#include <corbo-controllers/simple_state_controller.h>
#include <corbo-core/console.h>

#include <algorithm>

namespace corbo {

void SimpleStateController::setGainMatrixK(const Eigen::Ref<const Eigen::MatrixXd>& K) { _K = K; }

void SimpleStateController::setFilterMatrixV(const Eigen::Ref<const Eigen::MatrixXd>& V) { _V = V; }

bool SimpleStateController::initialize(const StateVector& x, ReferenceTrajectoryInterface& expected_xref, ReferenceTrajectoryInterface& expected_uref,
                                       const Duration& expected_dt, const Time& t, ReferenceTrajectoryInterface* sref)
{
    assert(x.rows() == expected_xref.getDimension() && "Dimension mismatch in controller: current state x and reference");
    assert(getStateDimension() == property::INHERITED || x.rows() == getStateDimension());

    return true;
}

bool SimpleStateController::step(const ControllerInterface::StateVector& x, ReferenceTrajectoryInterface& xref, ReferenceTrajectoryInterface& uref,
                                 const Duration& dt, const Time& t, TimeSeries::Ptr u_sequence, TimeSeries::Ptr x_sequence,
                                 SignalTargetInterface* signal_target, ReferenceTrajectoryInterface* sref, ReferenceTrajectoryInterface* xinit,
                                 ReferenceTrajectoryInterface* uinit, const std::string& ns)
{
    // const int q = uref.getDimension();
    ControlVector u;

    ReferenceTrajectoryInterface::OutputVector xref_vec;
    xref.getReference(t, xref_vec);
    ReferenceTrajectoryInterface::OutputVector uref_vec;
    uref.getReference(t, uref_vec);

    if (_V.rows() > 0 && _V.cols() > 0)
    {
        u = -_K * x + _V * xref_vec;

        if (signal_target && _publish_error)
        {
            signal_target->sendMeasurement(ns + "controller/error_norml2", t.toSec(), {x.norm()});
        }
    }
    else
    {
        StateVector state_error = xref_vec - x;

        u = _K * state_error + uref_vec;

        if (signal_target && _publish_error)
        {
            signal_target->sendMeasurement(ns + "controller/error_norml2", t.toSec(), {state_error.norm()});
        }
    }

    u_sequence->clear();
    x_sequence->clear();
    u_sequence->add(0.0, u);
    x_sequence->add(0.0, x);

    return true;
}

void SimpleStateController::getAvailableSignals(SignalTargetInterface& signal_target, const std::string& ns) const
{
    if (_publish_error)
    {
        signal_target.registerMeasurement(ns + "controller/error_norml2", 1);
    }
}

void SimpleStateController::reset() {}

#ifdef MESSAGE_SUPPORT
void SimpleStateController::toMessage(corbo::messages::Controller& message) const
{
    // feedback gain matrix K
    message.mutable_simple_state_controller()->mutable_k()->Resize(_K.rows() * _K.cols(), 0);
    Eigen::Map<Eigen::Matrix<double, -1, -1, Eigen::RowMajor>>(message.mutable_simple_state_controller()->mutable_k()->mutable_data(), _K.rows(),
                                                               _K.cols()) = _K;

    // filter matrix V
    message.mutable_simple_state_controller()->mutable_v()->Resize(_V.rows() * _V.cols(), 0);
    Eigen::Map<Eigen::Matrix<double, -1, -1, Eigen::RowMajor>>(message.mutable_simple_state_controller()->mutable_v()->mutable_data(), _V.rows(),
                                                               _V.cols()) = _V;

    message.mutable_simple_state_controller()->set_state_dim(getStateDimension());

    // publish error
    message.mutable_simple_state_controller()->set_publish_error(_publish_error);
}

void SimpleStateController::fromMessage(const corbo::messages::Controller& message, std::stringstream* issues)
{
    int state_dim   = 0;
    int control_dim = 0;
    int output_dim  = 0;

    if (message.simple_state_controller().v_size() == 0)
    {
        // no prefilter specified
        state_dim   = message.simple_state_controller().state_dim();
        output_dim  = state_dim;
        control_dim = message.simple_state_controller().k_size() / state_dim;
        if (message.simple_state_controller().output_dim() != output_dim && issues)
        {
            *issues << "SimpleStateController: output_dim must match state_dim since no prefilter matrix V is specified.\n";
        }
    }
    else
    {
        output_dim  = message.simple_state_controller().output_dim();
        control_dim = message.simple_state_controller().v_size() / output_dim;
        state_dim   = message.simple_state_controller().state_dim();
        if (control_dim != message.simple_state_controller().k_size() / state_dim && issues)
        {
            *issues << "SimpleStateController: dimension mismatch. cols(K) must be equal to cols(V).\n";
        }

        if (output_dim * control_dim == message.simple_state_controller().v_size())
        {
            setFilterMatrixV(Eigen::Map<const Eigen::Matrix<double, -1, -1, Eigen::RowMajor>>(message.simple_state_controller().v().data(),
                                                                                              control_dim, output_dim));
        }
        else if (issues)
            *issues << "SimpleStateController: invalid size of filter matrix V.\n";
    }
    if (state_dim * control_dim == message.simple_state_controller().k_size())
    {
        setGainMatrixK(
            Eigen::Map<const Eigen::Matrix<double, -1, -1, Eigen::RowMajor>>(message.simple_state_controller().k().data(), control_dim, state_dim));
    }
    else if (issues)
        *issues << "SimpleStateController: invalid size of feedback gain matrix K.\n";

    // publish error
    _publish_error = message.simple_state_controller().publish_error();
}
#endif

}  // namespace corbo
