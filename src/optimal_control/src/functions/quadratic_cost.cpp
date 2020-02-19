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

#include <corbo-optimal-control/functions/quadratic_cost.h>

#include <corbo-numerics/matrix_utilities.h>

#include <cmath>

namespace corbo {

bool QuadraticFormCost::setWeightQ(const Eigen::Ref<const Eigen::MatrixXd>& Q)
{
    _Q_diagonal_mode_intentionally = false;

    if (!is_square(Q)) return false;

    // check if we have a diagonal matrix
    if (Q.isDiagonal(1e-10)) return setWeightQ(Q.diagonal().asDiagonal());

    _Q_diagonal_mode               = false;
    _Q_diagonal_mode_intentionally = false;

    _Q = Q;  // also store original matrix
    if (Q.isZero())
    {
        _Q_sqrt.setZero();
        return true;
    }
    Eigen::LLT<Eigen::MatrixXd, Eigen::Upper> cholesky(Q);
    if (cholesky.info() == Eigen::NumericalIssue) return false;
    _Q_sqrt = cholesky.matrixU();
    return true;
}

bool QuadraticFormCost::setWeightQ(const Eigen::DiagonalMatrix<double, -1>& Q)
{
    _Q_diagonal_mode_intentionally = true;
    _Q_diagonal_mode               = true;
    _Q_diag                        = Q;
    _Q_diag_sqrt                   = Q.diagonal().cwiseSqrt().asDiagonal();
    _Q                             = Q.toDenseMatrix();
    return true;
}

bool QuadraticFormCost::setWeightR(const Eigen::Ref<const Eigen::MatrixXd>& R)
{
    _R_diagonal_mode_intentionally = false;

    if (!is_square(R)) return false;

    // check if we have a diagonal matrix
    if (R.isDiagonal(1e-10)) return setWeightR(R.diagonal().asDiagonal());

    _R_diagonal_mode = false;

    _R = R;  // also store original matrix
    if (R.isZero())
    {
        _R_sqrt.setZero();
        return true;
    }
    Eigen::LLT<Eigen::MatrixXd, Eigen::Upper> cholesky(R);
    if (cholesky.info() == Eigen::NumericalIssue) return false;
    _R_sqrt = cholesky.matrixU();
    return true;
}

bool QuadraticFormCost::setWeightR(const Eigen::DiagonalMatrix<double, -1>& R)
{
    _R_diagonal_mode_intentionally = true;
    _R_diagonal_mode               = true;
    _R_diag                        = R;
    _R_diag_sqrt                   = R.diagonal().cwiseSqrt().asDiagonal();
    _R                             = R.toDenseMatrix();
    return true;
}

void QuadraticFormCost::computeNonIntegralStateTerm(int k, const Eigen::Ref<const Eigen::VectorXd>& x_k, Eigen::Ref<Eigen::VectorXd> cost) const
{
    assert(!_integral_form);
    assert(cost.size() == getNonIntegralStateTermDimension(k));

    if (_lsq_form)
    {
        if (_zero_x_ref)
        {
            if (_Q_diagonal_mode)
                cost.noalias() = _Q_diag_sqrt * x_k;
            else
                cost.noalias() = x_k.transpose() * _Q_sqrt * x_k;
        }
        else
        {
            Eigen::VectorXd xd = x_k - _x_ref->getReferenceCached(k);
            if (_Q_diagonal_mode)
                cost.noalias() = _Q_diag_sqrt * xd;
            else
                cost.noalias() = _Q_sqrt * xd;
        }
    }
    else
    {
        if (_zero_x_ref)
        {
            if (_Q_diagonal_mode)
                cost.noalias() = x_k.transpose() * _Q_diag * x_k;
            else
                cost.noalias() = x_k.transpose() * _Q * x_k;
        }
        else
        {
            Eigen::VectorXd xd = x_k - _x_ref->getReferenceCached(k);
            if (_Q_diagonal_mode)
                cost.noalias() = xd.transpose() * _Q_diag * xd;
            else
                cost.noalias() = xd.transpose() * _Q * xd;
        }
    }
}

void QuadraticFormCost::computeNonIntegralControlTerm(int k, const Eigen::Ref<const Eigen::VectorXd>& u_k, Eigen::Ref<Eigen::VectorXd> cost) const
{
    assert(!_integral_form);
    assert(cost.size() == getNonIntegralControlTermDimension(k));

    if (_lsq_form)
    {
        if (_zero_u_ref)
        {
            if (_R_diagonal_mode)
                cost.noalias() = _R_diag_sqrt * u_k;
            else
                cost.noalias() = _R_sqrt * u_k;
        }
        else
        {
            Eigen::VectorXd ud = u_k - _u_ref->getReferenceCached(k);
            if (_R_diagonal_mode)
                cost.noalias() = ud.transpose() * _R_diag_sqrt * ud;
            else
                cost.noalias() = ud.transpose() * _R_sqrt * ud;
        }
    }
    else
    {
        if (_zero_u_ref)
        {
            if (_R_diagonal_mode)
                cost.noalias() = u_k.transpose() * _R_diag * u_k;
            else
                cost.noalias() = u_k.transpose() * _R * u_k;
        }
        else
        {
            Eigen::VectorXd ud = u_k - _u_ref->getReferenceCached(k);
            if (_R_diagonal_mode)
                cost.noalias() = ud.transpose() * _R_diag * ud;
            else
                cost.noalias() = ud.transpose() * _R * ud;
        }
    }
}

void QuadraticFormCost::computeIntegralStateControlTerm(int k, const Eigen::Ref<const Eigen::VectorXd>& x_k,
                                                        const Eigen::Ref<const Eigen::VectorXd>& u_k, Eigen::Ref<Eigen::VectorXd> cost) const
{
    assert(_integral_form);
    assert(cost.size() == 1);

    cost[0] = 0;

    if (_zero_x_ref)
    {
        if (_Q_diagonal_mode)
            cost[0] += x_k.transpose() * _Q_diag * x_k;
        else
            cost[0] += x_k.transpose() * _Q * x_k;
    }
    else
    {
        Eigen::VectorXd xd = x_k - _x_ref->getReferenceCached(k);
        if (_Q_diagonal_mode)
            cost[0] += xd.transpose() * _Q_diag * xd;
        else
            cost[0] += xd.transpose() * _Q * xd;
    }
    if (_zero_u_ref)
    {
        if (_R_diagonal_mode)
            cost[0] += u_k.transpose() * _R_diag * u_k;
        else
            cost[0] += u_k.transpose() * _R * u_k;
    }
    else
    {
        Eigen::VectorXd ud = u_k - _u_ref->getReferenceCached(k);
        if (_R_diagonal_mode)
            cost[0] += ud.transpose() * _R_diag * ud;
        else
            cost[0] += ud.transpose() * _R * ud;
    }
}

bool QuadraticFormCost::checkParameters(int state_dim, int control_dim, std::stringstream* issues) const
{
    bool success = true;

    if (_Q_diagonal_mode_intentionally && _Q_diagonal_mode)
    {
        if (_Q_diag.diagonal().size() != state_dim)
        {
            if (issues)
                *issues << "QuadraticFormCost: Diagonal matrix dimension of Q (" << _Q_diag.diagonal().size()
                        << ") does not match state vector dimension (" << state_dim << "); Please specify diagonal elements only." << std::endl;
            success = false;
        }
    }
    else
    {
        if (_Q.rows() != state_dim || _Q.cols() != state_dim)
        {
            if (issues)
                *issues << "QuadraticFormCost: Matrix dimension of Q (" << _Q.rows() << "x" << _Q.cols()
                        << ") does not match state vector dimension (" << state_dim << "); Please specify " << (state_dim * state_dim)
                        << " elements (Row-Major)." << std::endl;
            success = false;
        }
    }

    if (_R_diagonal_mode_intentionally && _R_diagonal_mode)
    {
        if (_R_diag.diagonal().size() != control_dim)
        {
            if (issues)
                *issues << "QuadraticFormCost: diagonal matrix dimension of R (" << _R_diag.diagonal().size()
                        << ") does not match control input vector dimension (" << control_dim << "); Please specify diagonal elements only."
                        << std::endl;
            success = false;
        }
    }
    else
    {
        if (_R.rows() != control_dim || _R.cols() != control_dim)
        {
            if (issues)
                *issues << "QuadraticFormCost: Matrix dimension of R (" << _R.rows() << "x" << _R.cols()
                        << ") does not match control input vector dimension (" << control_dim << "); Please specify " << (control_dim * control_dim)
                        << " elements (Row-Major)." << std::endl;
            success = false;
        }
    }

    return success;
}

void QuadraticFormCost::scaleCurrentWeightQ(double scale)
{
    _Q *= scale;
    _Q_diag.diagonal() *= scale;
    _Q_sqrt *= scale;
    _Q_diag_sqrt.diagonal() *= scale;
}

void QuadraticFormCost::scaleCurrentWeightR(double scale)
{
    _R *= scale;
    _R_diag.diagonal() *= scale;
    _R_sqrt *= scale;
    _R_diag_sqrt.diagonal() *= scale;
}

#ifdef MESSAGE_SUPPORT
bool QuadraticFormCost::fromMessage(const messages::QuadraticFormCost& message, std::stringstream* issues)
{
    const messages::QuadraticFormCost& msg = message;

    _Q_diagonal_mode = msg.q_diagonal_only();
    _R_diagonal_mode = msg.r_diagonal_only();

    // Weight matrix Q
    if (_Q_diagonal_mode)
    {
        if (!setWeightQ(Eigen::Map<const Eigen::Matrix<double, -1, 1>>(msg.q().data(), msg.q_size()).asDiagonal()))
        {
            *issues << "QuadraticStateCost: cannot set diagonal weight matrix Q.\n";
            return false;
        }
    }
    else
    {
        int p = std::sqrt(msg.q_size());

        if (p * p != msg.q_size())
        {
            *issues << "QuadraticStateCost: weight matrix Q is not square.\n";
            return false;
        }

        // weight matrix Q
        // if (p * p == message.full_discretization_ocp().q_size())
        //{
        if (!setWeightQ(Eigen::Map<const Eigen::Matrix<double, -1, -1, Eigen::RowMajor>>(msg.q().data(), p, p)) && issues)
        {
            *issues << "QuadraticStateCost: weight matrix Q is not positive definite.\n";
            return false;
        }
        //}
        // else if (issues)
        //    *issues << "QuadraticStateCost: invalid size of weight matrix Q: it shout be [" << p << " x " << p << "].\n";
    }

    // Weight matrix R
    if (_R_diagonal_mode)
    {
        if (!setWeightR(Eigen::Map<const Eigen::Matrix<double, -1, 1>>(msg.r().data(), msg.r_size()).asDiagonal()))
        {
            *issues << "QuadraticFormCost: cannot set diagonal weight matrix R.\n";
            return false;
        }
    }
    else
    {
        double q = std::sqrt(msg.r_size());

        if (q * q != msg.r_size())
        {
            *issues << "QuadraticFormCost: weight matrix R is not square.\n";
            return false;
        }

        // weight matrix Q
        // if (p * p == message.full_discretization_ocp().q_size())
        //{
        if (!setWeightR(Eigen::Map<const Eigen::Matrix<double, -1, -1, Eigen::RowMajor>>(msg.r().data(), q, q)) && issues)
        {
            *issues << "QuadraticFormCost: weight matrix R is not positive definite.\n";
            return false;
        }
        //}
        // else if (issues)
        //    *issues << "QuadraticControlCost: invalid size of weight matrix Q: it shout be [" << p << " x " << p << "].\n";
    }

    // others
    _lsq_form      = msg.lsq_form();
    _integral_form = msg.integral_form();

    return true;
}

void QuadraticFormCost::toMessage(messages::QuadraticFormCost& message) const
{
    messages::QuadraticFormCost* msg = &message;

    // weight matrix Q
    if (_Q_diagonal_mode_intentionally && _Q_diagonal_mode)
    {
        Eigen::VectorXd Qdiag = _Q_diag.diagonal();
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
    if (_R_diagonal_mode_intentionally && _R_diagonal_mode)
    {
        Eigen::VectorXd Rdiag = _R_diag.diagonal();
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

    // others
    msg->set_lsq_form(_lsq_form);
    msg->set_integral_form(_integral_form);
}
#endif

}  // namespace corbo
