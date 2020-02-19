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

#include <corbo-optimal-control/functions/final_state_cost.h>
#include <corbo-optimal-control/structured_ocp/discretization_grids/discretization_grid_interface.h>

#include <corbo-numerics/algebraic_riccati_continuous.h>
#include <corbo-numerics/algebraic_riccati_discrete.h>
#include <corbo-numerics/controllability.h>
#include <corbo-numerics/matrix_utilities.h>

#include <corbo-communication/utilities.h>

#include <cmath>

namespace corbo {

bool QuadraticFinalStateCost::setWeightQf(const Eigen::Ref<const Eigen::MatrixXd>& Qf)
{
    _diagonal_mode_intentionally = false;

    if (!is_square(Qf)) return false;

    // check if we have a diagonal matrix
    if (Qf.isDiagonal(1e-10)) return setWeightQf(Qf.diagonal().asDiagonal());

    _diagonal_mode = false;

    _Qf = Qf;  // also store original Qf
    if (Qf.isZero())
    {
        _Qf_sqrt.setZero();
        return true;
    }
    Eigen::LLT<Eigen::MatrixXd, Eigen::Upper> cholesky(Qf);
    if (cholesky.info() == Eigen::NumericalIssue) return false;
    _Qf_sqrt = cholesky.matrixU();
    return true;
}

bool QuadraticFinalStateCost::setWeightQf(const Eigen::DiagonalMatrix<double, -1>& Qf)
{
    _diagonal_mode_intentionally = true;
    _diagonal_mode               = true;
    _Qf_diag                     = Qf;
    _Qf_diag_sqrt                = Qf.diagonal().cwiseSqrt().asDiagonal();
    _Qf                          = _Qf_diag.toDenseMatrix();
    return true;
}

void QuadraticFinalStateCost::computeNonIntegralStateTerm(int k, const Eigen::Ref<const Eigen::VectorXd>& x_k, Eigen::Ref<Eigen::VectorXd> cost) const
{
    assert(cost.size() == getNonIntegralStateTermDimension(k));

    if (_lsq_form)
    {
        if (_zero_x_ref)
        {
            if (_diagonal_mode)
                cost.noalias() = _Qf_diag_sqrt * x_k;
            else
                cost.noalias() = _Qf_sqrt * x_k;
        }
        else
        {
            Eigen::VectorXd xd = x_k - _x_ref->getReferenceCached(k);
            if (_diagonal_mode)
                cost.noalias() = _Qf_diag_sqrt * xd;
            else
                cost.noalias() = _Qf_sqrt * xd;
        }
    }
    else
    {
        if (_zero_x_ref)
        {
            if (_diagonal_mode)
                cost.noalias() = x_k.transpose() * _Qf_diag * x_k;
            else
                cost.noalias() = x_k.transpose() * _Qf * x_k;
        }
        else
        {
            Eigen::VectorXd xd = x_k - _x_ref->getReferenceCached(k);
            if (_diagonal_mode)
                cost.noalias() = xd.transpose() * _Qf_diag * xd;
            else
                cost.noalias() = xd.transpose() * _Qf * xd;
        }
    }
}

bool QuadraticFinalStateCost::checkParameters(int state_dim, int control_dim, std::stringstream* issues) const
{
    bool success = true;

    if (_diagonal_mode_intentionally && _diagonal_mode)
    {
        if (_Qf_diag.diagonal().size() != state_dim)
        {
            if (issues)
                *issues << "QuadraticFinalStateCost: Diagonal matrix dimension of Qf (" << _Qf_diag.diagonal().size()
                        << ") does not match state vector dimension (" << state_dim << "); Please specify diagonal elements only." << std::endl;
            success = false;
        }
    }
    else
    {
        if (_Qf.rows() != state_dim || _Qf.cols() != state_dim)
        {
            if (issues)
                *issues << "QuadraticFinalStateCost: Matrix dimension of Qf (" << _Qf.rows() << "x" << _Qf.cols()
                        << ") does not match state vector dimension (" << state_dim << "); Please specify " << (state_dim * state_dim)
                        << " elements (Row-Major)." << std::endl;
            success = false;
        }
    }

    return success;
}

#ifdef MESSAGE_SUPPORT
bool QuadraticFinalStateCost::fromMessage(const messages::FinalStageCost& message, std::stringstream* issues)
{
    _diagonal_mode = message.quadratic_state_cost().diagonal_only();

    if (_diagonal_mode)
    {
        if (!setWeightQf(
                Eigen::Map<const Eigen::Matrix<double, -1, 1>>(message.quadratic_state_cost().qf().data(), message.quadratic_state_cost().qf_size())
                    .asDiagonal()))
        {
            *issues << "QuadraticFinalStateCost: cannot set diagonal weight matrix Qf.\n";
            return false;
        }
    }
    else
    {
        int p = std::sqrt(message.quadratic_state_cost().qf_size());

        if (p * p != message.quadratic_state_cost().qf_size())
        {
            *issues << "QuadraticFinalStateCost: weight matrix Qf is not square.\n";
            return false;
        }

        if (!setWeightQf(Eigen::Map<const Eigen::Matrix<double, -1, -1, Eigen::RowMajor>>(message.quadratic_state_cost().qf().data(), p, p)) &&
            issues)
        {
            *issues << "QuadraticFinalStateCost: weight matrix Qf is not positive definite.\n";
            return false;
        }
    }

    // lsq form
    _lsq_form = message.quadratic_state_cost().lsq_form();

    return true;
}

void QuadraticFinalStateCost::toMessage(messages::FinalStageCost& message) const
{
    // weight matrix Q
    if (_diagonal_mode_intentionally && _diagonal_mode)
    {
        Eigen::VectorXd Qfdiag = _Qf_diag.diagonal();
        message.mutable_quadratic_state_cost()->mutable_qf()->Resize(Qfdiag.size(), 0);
        Eigen::Map<Eigen::VectorXd>(message.mutable_quadratic_state_cost()->mutable_qf()->mutable_data(), Qfdiag.size()) = Qfdiag;

        message.mutable_quadratic_state_cost()->set_diagonal_only(true);
    }
    else
    {
        Eigen::MatrixXd Qf = _Qf_sqrt.adjoint() * _Qf_sqrt;
        message.mutable_quadratic_state_cost()->mutable_qf()->Resize(Qf.rows() * Qf.cols(), 0);
        Eigen::Map<Eigen::Matrix<double, -1, -1, Eigen::RowMajor>>(message.mutable_quadratic_state_cost()->mutable_qf()->mutable_data(), Qf.rows(),
                                                                   Qf.cols()) = Qf;

        message.mutable_quadratic_state_cost()->set_diagonal_only(false);
    }

    // lsq_form
    message.mutable_quadratic_state_cost()->set_lsq_form(_lsq_form);
}
#endif

bool QuadraticFinalStateCostRiccati::setWeightQ(const Eigen::Ref<const Eigen::MatrixXd>& Q)
{
    if (!is_square(Q)) return false;
    if (!is_positive_definite(Q)) return false;
    _Q_diagonal_mode_intentionally = false;

    _Q = Q;
    return true;
}

bool QuadraticFinalStateCostRiccati::setWeightQ(const Eigen::DiagonalMatrix<double, -1>& Q)
{
    _Q_diagonal_mode_intentionally = true;

    _Q = Q.toDenseMatrix();
    if (!is_positive_definite(_Q)) return false;
    return true;
}

bool QuadraticFinalStateCostRiccati::setWeightR(const Eigen::Ref<const Eigen::MatrixXd>& R)
{
    if (!is_square(R)) return false;
    if (!is_positive_definite(R)) return false;
    _R_diagonal_mode_intentionally = false;

    _R = R;
    return true;
}

bool QuadraticFinalStateCostRiccati::setWeightR(const Eigen::DiagonalMatrix<double, -1>& R)
{
    _R_diagonal_mode_intentionally = true;

    _R = R.toDenseMatrix();
    if (!is_positive_definite(_R)) return false;
    return true;
}

bool QuadraticFinalStateCostRiccati::update(int n, double t, ReferenceTrajectoryInterface& xref, ReferenceTrajectoryInterface& uref,
                                            ReferenceTrajectoryInterface* sref, bool single_dt, const Eigen::VectorXd& x0,
                                            StagePreprocessor::Ptr stage_preprocessor, const std::vector<double>& dts,
                                            const DiscretizationGridInterface* grid)
{
    bool dimension_modified = BaseQuadraticFinalStateCost::update(n, t, xref, uref, sref, single_dt, x0, stage_preprocessor, dts, grid);

    // store reference trajectory
    _x_ref = &xref;

    _zero_x_ref = _x_ref->isZero();

    assert(_dynamics);

    Time t_ref(t);
    Eigen::VectorXd steady_state_x = xref.getNextSteadyState(t_ref);
    Eigen::VectorXd steady_state_u = uref.getNextSteadyState(t_ref);

    if (!_are_solved_before || !_dynamics->isLinear())
    {
        // don't solve if steady steade has not changed
        if (_steady_state_x.size() == 0 || _steady_state_u.size() == 0 || _steady_state_x != steady_state_x || _steady_state_u != steady_state_u)
        {
            Eigen::MatrixXd A(_dynamics->getStateDimension(), _dynamics->getStateDimension());
            Eigen::MatrixXd B(_dynamics->getStateDimension(), _dynamics->getInputDimension());
            _dynamics->getLinearA(steady_state_x, steady_state_u, A);
            _dynamics->getLinearB(steady_state_x, steady_state_u, B);

            assert(have_equal_size(_Q, A));
            assert(_R.rows() == _dynamics->getInputDimension());
            assert(_R.cols() == _dynamics->getInputDimension());

            _Qf.resize(_Q.rows(), _Q.cols());

#ifndef NDEBUG
            if (!Controllability::checkLinearTimeInvariantSystem(A, B))
            {
                PRINT_WARNING("QuadraticFinalStateCostRiccati::update(): the linearized model is not fully controllable.");
            }
#endif

            if (_dynamics->isContinuousTime())
            {
                if (!AlgebraicRiccatiContinuous::solve(A, B, _Q, _R, _Qf))
                {
                    PRINT_ERROR("QuadraticFinalStateCostRiccati::update(): continuous-time algebraic riccati solver failed. Setting Qf = Q.");
                    _Qf = _Q;
                }
            }
            else
            {
                if (!AlgebraicRiccatiDiscrete::solve(A, B, _Q, _R, _Qf))
                {
                    PRINT_ERROR("QuadraticFinalStateCostRiccati::update(): discrete-time algebraic riccati solver failed. Setting Qf = Q.");
                    _Qf = _Q;
                }
            }
            // update internals
            if (!computeWeightQfSqrt() && _lsq_form)
            {
                PRINT_ERROR(
                    "QuadraticFinalStateCostRiccati::update(): Cholesky solution on Qf failed. Since lsq_mode is on, setting Qf_sqrt = Q_sqrt.");
                _Qf = _Q;
                computeWeightQfSqrt();
            }
            _steady_state_x    = steady_state_x;
            _steady_state_u    = steady_state_u;
            _are_solved_before = true;
        }
    }
    return dimension_modified;
}

bool QuadraticFinalStateCostRiccati::computeWeightQfSqrt()
{
    if (_Qf.isZero())
    {
        _Qf_sqrt.setZero();
        return true;
    }
    Eigen::LLT<Eigen::MatrixXd, Eigen::Upper> cholesky(_Qf);
    if (cholesky.info() == Eigen::NumericalIssue) return false;
    _Qf_sqrt = cholesky.matrixU();
    return true;
}

void QuadraticFinalStateCostRiccati::computeNonIntegralStateTerm(int k, const Eigen::Ref<const Eigen::VectorXd>& x_k,
                                                                 Eigen::Ref<Eigen::VectorXd> cost) const
{
    assert(cost.size() == getNonIntegralStateTermDimension(k));

    if (_zero_x_ref)
    {
        cost.noalias() = x_k.transpose() * _Qf * x_k;
    }
    else
    {
        Eigen::VectorXd xd = x_k - _x_ref->getReferenceCached(k);
        cost.noalias()     = xd.transpose() * _Qf * xd;
    }
}

bool QuadraticFinalStateCostRiccati::checkParameters(int state_dim, int control_dim, std::stringstream* issues) const
{
    bool success = true;

    if (!_dynamics)
    {
        if (issues) *issues << "QuadraticFinalStateCostRiccati: No system model specified. Cannot solve algebraic riccati equation." << std::endl;
        return false;
    }

    if (_dynamics->getStateDimension() != state_dim)
    {
        if (issues)
            *issues << "QuadraticFinalStateCostRiccati: State dimension of the specified sytem model (" << _dynamics->getStateDimension()
                    << ") does not match state vector dimension (" << state_dim << ")." << std::endl;
        success = false;
    }

    if (_dynamics->getInputDimension() != control_dim)
    {
        if (issues)
            *issues << "QuadraticFinalStateCostRiccati: Control input dimension of the specified sytem model (" << _dynamics->getStateDimension()
                    << ") does not match control input vector dimension (" << control_dim << ")." << std::endl;
        success = false;
    }

    if (_Q.rows() != state_dim || _Q.cols() != state_dim)
    {
        if (issues)
            *issues << "QuadraticFinalStateCostRiccati: Matrix dimension of Q (" << _Q.rows() << "x" << _Q.cols()
                    << ") does not match state vector dimension (" << state_dim << "); Please specify " << (state_dim * state_dim)
                    << " elements (Row-Major)." << std::endl;
        success = false;
    }

    if (_R.rows() != control_dim || _R.cols() != control_dim)
    {
        if (issues)
            *issues << "QuadraticFinalStateCostRiccati: Matrix dimension of R (" << _R.rows() << "x" << _R.cols()
                    << ") does not match control input vector dimension (" << control_dim << "); Please specify " << (control_dim * control_dim)
                    << " elements (Row-Major)." << std::endl;
        success = false;
    }

    return success;
}

#ifdef MESSAGE_SUPPORT
bool QuadraticFinalStateCostRiccati::fromMessage(const messages::FinalStageCost& message, std::stringstream* issues)
{
    const messages::QuadraticFinalStateCostRiccati& msg = message.quadratic_state_cost_riccati();

    _Q_diagonal_mode_intentionally = msg.q_diagonal_only();
    _R_diagonal_mode_intentionally = msg.r_diagonal_only();

    // system dynamics
    if (!msg.has_system_dynamics())
    {
        if (issues) *issues << "QuadraticFinalStateCostRiccati: no system model specified.\n";
        return false;
    }
    // construct object
    std::string type;
    util::get_oneof_field_type_expand_isolated(msg.system_dynamics(), "system_dynamics", type, false, 1);
    SystemDynamicsInterface::Ptr dynamics = SystemDynamicsFactory::instance().create(type);
    // import parameters
    if (dynamics)
    {
        dynamics->fromMessage(msg.system_dynamics(), issues);
        setSystemDynamics(dynamics);
    }
    else
    {
        if (issues) *issues << "QuadraticFinalStateCostRiccati: unknown system model specified.\n";
        return false;
    }

    // Weight matrix Q
    if (_Q_diagonal_mode_intentionally)
    {
        if (!setWeightQ(Eigen::Map<const Eigen::Matrix<double, -1, 1>>(msg.q().data(), msg.q_size()).asDiagonal()))
        {
            *issues << "QuadraticFinalStateCostRiccati: cannot set diagonal weight matrix Q.\n";
            return false;
        }
    }
    else
    {
        int p = std::sqrt(msg.q_size());

        if (p * p != msg.q_size())
        {
            *issues << "QuadraticFinalStateCostRiccati: weight matrix Q is not square.\n";
            return false;
        }

        if (!setWeightQ(Eigen::Map<const Eigen::Matrix<double, -1, -1, Eigen::RowMajor>>(msg.q().data(), p, p)) && issues)
        {
            *issues << "QuadraticFinalStateCostRiccati: weight matrix Q is not positive definite.\n";
            return false;
        }
    }

    // Weight matrix R
    if (_R_diagonal_mode_intentionally)
    {
        if (!setWeightR(Eigen::Map<const Eigen::Matrix<double, -1, 1>>(msg.r().data(), msg.r_size()).asDiagonal()))
        {
            *issues << "QuadraticFinalStateCostRiccati: cannot set diagonal weight matrix R.\n";
            return false;
        }
    }
    else
    {
        double q = std::sqrt(msg.r_size());

        if (q * q != msg.r_size())
        {
            *issues << "QuadraticFinalStateCostRiccati: weight matrix R is not square.\n";
            return false;
        }

        if (!setWeightR(Eigen::Map<const Eigen::Matrix<double, -1, -1, Eigen::RowMajor>>(msg.r().data(), q, q)) && issues)
        {
            *issues << "QuadraticFinalStateCostRiccati: weight matrix R is not positive definite.\n";
            return false;
        }
    }

    // lsq form
    _lsq_form = msg.lsq_form();

    return true;
}

void QuadraticFinalStateCostRiccati::toMessage(messages::FinalStageCost& message) const
{
    messages::QuadraticFinalStateCostRiccati* msg = message.mutable_quadratic_state_cost_riccati();

    // system model
    if (_dynamics) _dynamics->toMessage(*msg->mutable_system_dynamics());

    // weight matrix Q
    if (_Q_diagonal_mode_intentionally)
    {
        Eigen::VectorXd Qdiag = _Q.diagonal();
        msg->mutable_q()->Resize(Qdiag.size(), 0);
        Eigen::Map<Eigen::VectorXd>(msg->mutable_q()->mutable_data(), Qdiag.size()) = Qdiag;

        msg->set_q_diagonal_only(true);
    }
    else
    {
        Eigen::MatrixXd Q = _Q;  // Q = _Q_sqrt.adjoint() * _Q_sqrt;
        msg->mutable_q()->Resize(Q.rows() * Q.cols(), 0);
        Eigen::Map<Eigen::Matrix<double, -1, -1, Eigen::RowMajor>>(msg->mutable_q()->mutable_data(), Q.rows(), Q.cols()) = Q;

        msg->set_q_diagonal_only(false);
    }

    // weight matrix R
    if (_R_diagonal_mode_intentionally)
    {
        Eigen::VectorXd Rdiag = _R.diagonal();
        msg->mutable_r()->Resize(Rdiag.size(), 0);
        Eigen::Map<Eigen::VectorXd>(msg->mutable_r()->mutable_data(), Rdiag.size()) = Rdiag;

        msg->set_r_diagonal_only(true);
    }
    else
    {
        Eigen::MatrixXd R = _R;  // R = _R_sqrt.adjoint() * _R_sqrt;
        msg->mutable_r()->Resize(R.rows() * R.cols(), 0);
        Eigen::Map<Eigen::Matrix<double, -1, -1, Eigen::RowMajor>>(msg->mutable_r()->mutable_data(), R.rows(), R.cols()) = R;

        msg->set_r_diagonal_only(false);
    }

    // lsq_form
    message.mutable_quadratic_state_cost()->set_lsq_form(_lsq_form);
}
#endif

}  // namespace corbo
