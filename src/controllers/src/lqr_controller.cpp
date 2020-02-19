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
#include <corbo-controllers/lqr_controller.h>
#include <corbo-core/console.h>
#include <corbo-numerics/algebraic_riccati_continuous.h>
#include <corbo-numerics/algebraic_riccati_discrete.h>

#include <Eigen/Cholesky>
#include <Eigen/Core>
#include <Eigen/Dense>

#include <algorithm>

namespace corbo {

LqrController::LqrController() {}

LqrController::LqrController(SystemDynamicsInterface::Ptr system_model) { setSystemModel(system_model); }

void LqrController::setSystemModel(SystemDynamicsInterface::Ptr system_model)
{
    _system_model = system_model;

    const int p = getStateDimension();
    const int q = getControlInputDimension();

    setWeights(Eigen::MatrixXd::Identity(p, p), Eigen::MatrixXd::Identity(q, q));
}

bool LqrController::setWeights(const Eigen::Ref<const Eigen::MatrixXd>& R, const Eigen::Ref<const Eigen::MatrixXd>& Q)
{
    bool ret_val = true;
    if (!setWeightR(R)) ret_val = false;
    if (!setWeightQ(Q)) ret_val = false;
    return ret_val;
}

bool LqrController::setWeightR(const Eigen::Ref<const Eigen::MatrixXd>& R)
{
    if (!_system_model || R.rows() != _system_model->getInputDimension() || R.cols() != _system_model->getInputDimension()) return false;
    _R = R;
    return true;
}

bool LqrController::setWeightQ(const Eigen::Ref<const Eigen::MatrixXd>& Q)
{
    if (!_system_model || Q.rows() != _system_model->getStateDimension() || Q.cols() != _system_model->getStateDimension()) return false;
    _Q = Q;
    return true;
}

bool LqrController::initialize(const StateVector& x, ReferenceTrajectoryInterface& expected_xref, ReferenceTrajectoryInterface& expected_uref,
                               const Duration& expected_dt, const Time& t, ReferenceTrajectoryInterface* sref)
{
    assert(x.rows() == expected_xref.getDimension() && "Dimension mismatch in controller: current state x and reference");
    assert(getStateDimension() == property::INHERITED || x.rows() == getStateDimension());

    if (!_system_model)
    {
        PRINT_ERROR("LqrController: no system model provided. Cancelling.");
        return false;
    }

    if (_system_model->isContinuousTime())
    {
        _discrete_time = false;
        PRINT_INFO("LqrController configured for continuous-time system dynamics.");
    }
    else
    {
        _discrete_time = true;
        PRINT_INFO("LqrController configured for discrete-time system dynamics.");
    }

    ReferenceTrajectoryInterface::OutputVector xref = expected_xref.getNextSteadyState(t);

    if (xref.rows() != x.rows())
    {
        PRINT_WARNING("LqrController currently supports only full state reference trajectories");
        return false;
    }

    ReferenceTrajectoryInterface::OutputVector uref = expected_uref.getNextSteadyState(t);

    PRINT_INFO_COND(!_system_model->isLinear(), "LQR is provided with a nonlinear model: linearizing model at target state...");

    const int p = xref.rows();
    const int q = uref.rows();

    _A.resize(p, p);
    _system_model->getLinearA(xref, uref, _A);
    _B.resize(p, q);
    _system_model->getLinearB(xref, uref, _B);

    bool care_success;
    if (_discrete_time)
        care_success = AlgebraicRiccatiDiscrete::solve(_A, _B, _Q, _R, _S, &_K);
    else
    {
        PRINT_WARNING_COND_NAMED(!AlgebraicRiccatiContinuous::isNumericallyStable(_A, _B, _Q, _R),
                                 "Eigenvalues of the Hamiltonian are close to the imaginary axis. Numerically unstable results are expected.");
        care_success = AlgebraicRiccatiContinuous::solve(_A, _B, _Q, _R, _S, &_K);
    }

    PRINT_WARNING_COND(!care_success || !_S.allFinite(), "Riccati solving failed or NaN coefficients found in solution.");

    _initialized = true;

    return true;
}

bool LqrController::step(const ControllerInterface::StateVector& x, ReferenceTrajectoryInterface& xref, ReferenceTrajectoryInterface& uref,
                         const Duration& dt, const Time& t, TimeSeries::Ptr u_sequence, TimeSeries::Ptr x_sequence,
                         SignalTargetInterface* signal_target, ReferenceTrajectoryInterface* sref, ReferenceTrajectoryInterface* xinit,
                         ReferenceTrajectoryInterface* uinit, const std::string& ns)
{
    const int q = uref.getDimension();
    ControlVector u;

    _dt = dt.toSec();

    if (!_initialized)
    {
        if (!initialize(x, xref, uref, dt, t, sref))
        {
            u.setZero(q);
            return false;
        }
    }

    ReferenceTrajectoryInterface::OutputVector xref_vec;
    xref.getReference(t, xref_vec);
    ReferenceTrajectoryInterface::OutputVector uref_vec;
    uref.getReference(t, uref_vec);

    StateVector state_error = xref_vec - x;

    u = _K * state_error + uref_vec;

    if (signal_target && _publish_error)
    {
        signal_target->sendMeasurement(ns + "controller/error_norml2", t.toSec(), {state_error.norm()});
    }

    u_sequence->clear();
    x_sequence->clear();
    u_sequence->add(0.0, u);
    x_sequence->add(0.0, x);

    return true;
}

void LqrController::getAvailableSignals(SignalTargetInterface& signal_target, const std::string& ns) const
{
    if (_publish_error)
    {
        signal_target.registerMeasurement(ns + "controller/error_norml2", 1);
    }
}

void LqrController::reset() { _initialized = false; }

#ifdef MESSAGE_SUPPORT
void LqrController::toMessage(corbo::messages::LqrController& message) const
{
    if (_system_model) _system_model->toMessage(*message.mutable_system_model());

    // weight matrix Q
    if (_Q.isDiagonal())
    {
        Eigen::VectorXd Qdiag = _Q.diagonal();
        message.mutable_q()->Resize(Qdiag.size(), 0);
        Eigen::Map<Eigen::VectorXd>(message.mutable_q()->mutable_data(), Qdiag.size()) = Qdiag;

        message.set_q_diagonal_only(true);
    }
    else
    {
        message.mutable_q()->Resize(_Q.rows() * _Q.cols(), 0);
        Eigen::Map<Eigen::Matrix<double, -1, -1, Eigen::RowMajor>>(message.mutable_q()->mutable_data(), _Q.rows(), _Q.cols()) = _Q;

        message.set_q_diagonal_only(false);
    }

    // weight matrix R
    if (_R.isDiagonal())
    {
        Eigen::VectorXd Rdiag = _R.diagonal();
        message.mutable_r()->Resize(Rdiag.size(), 0);
        Eigen::Map<Eigen::VectorXd>(message.mutable_r()->mutable_data(), Rdiag.size()) = Rdiag;

        message.set_r_diagonal_only(true);
    }
    else
    {
        message.mutable_r()->Resize(_R.rows() * _R.cols(), 0);
        Eigen::Map<Eigen::Matrix<double, -1, -1, Eigen::RowMajor>>(message.mutable_r()->mutable_data(), _R.rows(), _R.cols()) = _R;

        message.set_r_diagonal_only(false);
    }

    // publish error
    message.set_publish_error(_publish_error);
}

void LqrController::fromMessage(const corbo::messages::LqrController& message, std::stringstream* issues)
{
    _system_model.reset();

    // system model
    if (!message.has_system_model())
    {
        if (issues) *issues << "LqrController: no system model specified.\n";
        return;
    }

    // construct object
    std::string type;
    util::get_oneof_field_type_expand_isolated(message.system_model(), "system_dynamics", type, false, 1);
    SystemDynamicsInterface::Ptr dynamics = SystemDynamicsFactory::instance().create(type);
    // import parameters
    if (dynamics)
    {
        dynamics->fromMessage(message.system_model(), issues);
        setSystemModel(dynamics);
    }
    else
    {
        if (issues) *issues << "LqrController: unknown system model specified.\n";
        return;
    }

    const int p = dynamics->getStateDimension();
    const int q = dynamics->getInputDimension();

    // weight matrix Q
    if (message.q_diagonal_only())
    {
        if (p == message.q_size())
        {
            Eigen::MatrixXd Q = Eigen::Map<const Eigen::Matrix<double, -1, 1>>(message.q().data(), message.q_size()).asDiagonal();
            if (!setWeightQ(Q) && issues) *issues << "LqrController: invalid size of weight matrix Q: it shout be [" << p << " x " << p << "].\n";
        }
        else if (issues)
            *issues << "LqrController: invalid size of weight matrix diag(Q): it shout be [ 1 x " << p << "].\n";
    }
    else
    {
        if (p * p == message.q_size())
        {
            if (!setWeightQ(Eigen::Map<const Eigen::Matrix<double, -1, -1, Eigen::RowMajor>>(message.q().data(), p, p)) && issues)
                *issues << "LqrController: invalid size of weight matrix Q: it shout be [" << p << " x " << p << "].\n";
        }
        else if (issues)
            *issues << "LqrController: invalid size of weight matrix Q: it shout be [" << p << " x " << p << "].\n";
    }
    // weight matrix R
    if (message.r_diagonal_only())
    {
        if (q == message.r_size())
        {
            Eigen::MatrixXd R = Eigen::Map<const Eigen::Matrix<double, -1, 1>>(message.r().data(), message.r_size()).asDiagonal();
            if (!setWeightR(R) && issues) *issues << "LqrController: invalid size of weight matrix R: it shout be [" << q << " x " << q << "].\n";
        }
        else if (issues)
            *issues << "LqrController: invalid size of weight matrix diag(R): it shout be [ 1 x " << q << "].\n";
    }
    else
    {
        if (q * q == message.r_size())
        {
            if (!setWeightR(Eigen::Map<const Eigen::Matrix<double, -1, -1, Eigen::RowMajor>>(message.r().data(), q, q)) && issues)
                *issues << "LqrController: invalid size of weight matrix R: it shout be [" << q << " x " << q << "].\n";
        }
        else if (issues)
            *issues << "LqrController: invalid size of weight matrix R: it shout be [" << q << " x " << q << "].\n";
    }

    // publish error
    _publish_error = message.publish_error();
}
#endif

}  // namespace corbo
