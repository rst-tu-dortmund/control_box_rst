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

#include <corbo-controllers/dual_mode_controller.h>

#include <corbo-communication/utilities.h>
#include <corbo-core/console.h>
#include <corbo-numerics/matrix_utilities.h>

namespace corbo {

bool DualModeController::initialize(const StateVector& x, ReferenceTrajectoryInterface& expected_xref, ReferenceTrajectoryInterface& expected_uref,
                                    const Duration& expected_dt, const Time& t, ReferenceTrajectoryInterface* sref)
{
    if (!_pred_controller.initialize(x, expected_xref, expected_uref, expected_dt, t, sref))
    {
        PRINT_ERROR("DualModeController::initialize(): predictive controller initialization failed.");
        return false;
    }
    if (!_local_controller->initialize(x, expected_xref, expected_uref, expected_dt, t, sref))
    {
        PRINT_ERROR("DualModeController::initialize(): LQR controller initialization failed.");
        return false;
    }

    if (!_switch_dt && !_switch_terminal_ball)
    {
        PRINT_WARNING("DualModeController::initialize(): no switch condition specified.");
    }

    _initialized = true;
    return true;
}

bool DualModeController::step(const ControllerInterface::StateVector& x, ReferenceTrajectoryInterface& xref, ReferenceTrajectoryInterface& uref,
                              const Duration& dt, const Time& t, TimeSeries::Ptr u_sequence, TimeSeries::Ptr x_sequence,
                              SignalTargetInterface* signal_target, ReferenceTrajectoryInterface* sref, ReferenceTrajectoryInterface* xinit,
                              ReferenceTrajectoryInterface* uinit, const std::string& ns)
{
    if (!_initialized)
    {
        if (!initialize(x, xref, uref, dt, t, sref)) return false;
    }
    if (!_local_controller)
    {
        PRINT_ERROR_NAMED("No local controller provided.");
        return false;
    }

    Eigen::VectorXd xf(x.size());
    xref.getReference(t, xf);
    // for dualmode, assume we have only static references.... // xref.getNextSteadyState(t); // TODO(roesmann): this might not be good?!

    bool success = false;

    Time t_pre = Time::now();

    if (_switch_terminal_ball) _local_ctrl_active = isInsideInTerminalBall(x, xf);

    if (!_first_run && _switch_dt)
    {
        _local_ctrl_active = (_pred_controller.getControlDuration() <= _min_dt);  // this is from the last run!
    }

    if (_local_ctrl_active)
    {
        success = _local_controller->step(x, xref, uref, dt, t, u_sequence, x_sequence);
    }
    else
    {
        success = _pred_controller.step(x, xref, uref, dt, t, u_sequence, x_sequence);
    }

    _statistics.step_time = Time::now() - t_pre;

    _first_run = false;
    return success;
}

bool DualModeController::isInsideInTerminalBall(const Eigen::Ref<const Eigen::VectorXd>& x0, const Eigen::Ref<const Eigen::VectorXd>& xf) const
{
    Eigen::VectorXd xd = xf - x0;
    return (xd.transpose() * _S * xd <= _gamma);
}

void DualModeController::getAvailableSignals(SignalTargetInterface& signal_target, const std::string& ns) const
{
    _pred_controller.getAvailableSignals(signal_target, ns + "pred_ctrl/");
    _local_controller->getAvailableSignals(signal_target, ns + "local_ctrl/");
    signal_target.registerMeasurement(ns + "cpu_time", 1);
    signal_target.registerMeasurement(ns + "local_active", 1);
}

void DualModeController::reset()
{
    _pred_controller.reset();
    _local_controller.reset();
    _local_ctrl_active = false;
    _first_run         = true;
}

void DualModeController::sendSignals(double t, SignalTargetInterface& signal_target, const std::string& ns) const
{
    if (_local_ctrl_active)
        _local_controller->sendSignals(t, signal_target, ns + "local_ctrl/");
    else
        _pred_controller.sendSignals(t, signal_target, ns + "pred_ctrl/");

    signal_target.sendMeasurement(ns + "cpu_time", t, {_statistics.step_time.toSec()});
    signal_target.sendMeasurement(ns + "local_active", t, {(double)_local_ctrl_active});
}

bool DualModeController::setWeightS(const Eigen::Ref<const Eigen::MatrixXd>& S)
{
    _S = S;
    return true;
}

#ifdef MESSAGE_SUPPORT
void DualModeController::toMessage(corbo::messages::DualModeController& message) const
{
    _pred_controller.toMessage(*message.mutable_predictive_controller());
    _local_controller->toMessage(*message.mutable_local_controller());

    message.set_switch_dt(_switch_dt);
    message.set_min_dt(_min_dt);
    message.set_switch_terminal_ball(_switch_terminal_ball);

    // Terminal Ball
    // weight matrix S
    if (_S.isDiagonal())
    {
        Eigen::VectorXd Sdiag = _S.diagonal();
        message.mutable_ball_s()->Resize(Sdiag.size(), 0);
        Eigen::Map<Eigen::VectorXd>(message.mutable_ball_s()->mutable_data(), Sdiag.size()) = Sdiag;

        message.set_ball_s_diagonal_only(true);
    }
    else
    {
        message.mutable_ball_s()->Resize(_S.rows() * _S.cols(), 0);
        Eigen::Map<Eigen::Matrix<double, -1, -1, Eigen::RowMajor>>(message.mutable_ball_s()->mutable_data(), _S.rows(), _S.cols()) = _S;

        message.set_ball_s_diagonal_only(false);
    }

    // gamma
    message.set_ball_gamma(_gamma);
}
void DualModeController::fromMessage(const corbo::messages::DualModeController& message, std::stringstream* issues)
{
    _pred_controller.fromMessage(message.predictive_controller(), issues);

    // construct local controller
    std::string type;
    util::get_oneof_field_type_expand_isolated(message.local_controller(), "controller", type, false, 1);
    ControllerInterface::Ptr local_controller = ControllerFactory::instance().create(type);
    // import parameters
    if (local_controller)
    {
        local_controller->fromMessage(message.local_controller(), issues);
        setLocalController(local_controller);
    }
    else
    {
        if (issues) *issues << "Cannot load local controller.\n";
        return;
    }

    int dim_x_pred  = _pred_controller.getStateDimension();
    int dim_x_local = _local_controller->getStateDimension();

    if (dim_x_pred != dim_x_local)
    {
        if (issues) *issues << "DualModeController: state dimension mismatch (lqr: " << dim_x_local << ", pred: " << dim_x_pred << ")." << std::endl;
    }

    int dim_u_pred  = _pred_controller.getControlInputDimension();
    int dim_u_local = _local_controller->getControlInputDimension();

    if (dim_u_pred != dim_u_local)
    {
        if (issues)
            *issues << "DualModeController: control input dimension mismatch (lqr: " << dim_u_local << ", pred: " << dim_u_pred << ")." << std::endl;
    }

    _switch_dt            = message.switch_dt();
    _min_dt               = message.min_dt();
    _switch_terminal_ball = message.switch_terminal_ball();

    // terminal ball
    if (_switch_terminal_ball)
    {
        if (message.ball_s_diagonal_only())
        {
            Eigen::MatrixXd S = Eigen::Map<const Eigen::Matrix<double, -1, 1>>(message.ball_s().data(), message.ball_s_size()).asDiagonal();
            if (!setWeightS(S))
            {
                *issues << "DualModeController: cannot set diagonal weight matrix S.\n";
                return;
            }
        }
        else
        {
            if (!is_square(message.ball_s_size()))
            {
                *issues << "DualModeController: weight matrix S is not square.\n";
                return;
            }
            int p = std::sqrt(message.ball_s_size());

            if (p != dim_x_pred)
            {
                if (issues)
                    *issues << "DualModeController: state dimension mismatch for terminal ball (dim_x_pred " << dim_x_pred << ", matrix S: " << p
                            << ")." << std::endl;
            }

            if (!setWeightS(Eigen::Map<const Eigen::Matrix<double, -1, -1, Eigen::RowMajor>>(message.ball_s().data(), p, p)) && issues)
            {
                *issues << "DualModeController: cannot set weight matrix S.\n";
                return;
            }
        }
    }

    // gamma
    _gamma = message.ball_gamma();
}
#endif

}  // namespace corbo
