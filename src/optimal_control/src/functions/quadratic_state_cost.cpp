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

#include <corbo-optimal-control/functions/quadratic_state_cost.h>

#include <corbo-numerics/matrix_utilities.h>

#include <cmath>

namespace corbo {

bool QuadraticStateCost::setWeightQ(const Eigen::Ref<const Eigen::MatrixXd>& Q)
{
    _diagonal_mode_intentionally = false;

    if (!is_square(Q)) return false;

    // check if we have a diagonal matrix
    if (Q.isDiagonal(1e-10)) return setWeightQ(Q.diagonal().asDiagonal());

    _diagonal_mode = false;

    _Q = Q;  // also store original Q
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

bool QuadraticStateCost::setWeightQ(const Eigen::DiagonalMatrix<double, -1>& Q)
{
    _diagonal_mode_intentionally = true;
    _diagonal_mode               = true;
    _Q_diag                      = Q;
    _Q_diag_sqrt                 = Q.diagonal().cwiseSqrt().asDiagonal();
    return true;
}

void QuadraticStateCost::computeNonIntegralStateTerm(int k, const Eigen::Ref<const Eigen::VectorXd>& x_k, Eigen::Ref<Eigen::VectorXd> cost) const
{
    assert(!_integral_form);
    assert(cost.size() == getNonIntegralStateTermDimension(k));

    if (_lsq_form)
    {
        if (_zero_x_ref)
        {
            if (_diagonal_mode)
                cost.noalias() = _Q_diag_sqrt * x_k;
            else
                cost.noalias() = _Q_sqrt * x_k;
        }
        else
        {
            Eigen::VectorXd xd = x_k - _x_ref->getReferenceCached(k);
            if (_diagonal_mode)
                cost.noalias() = _Q_diag_sqrt * xd;
            else
                cost.noalias() = _Q_sqrt * xd;
        }
    }
    else
    {
        if (_zero_x_ref)
        {
            if (_diagonal_mode)
                cost.noalias() = x_k.transpose() * _Q_diag * x_k;
            else
                cost.noalias() = x_k.transpose() * _Q * x_k;
        }
        else
        {
            Eigen::VectorXd xd = x_k - _x_ref->getReferenceCached(k);
            if (_diagonal_mode)
                cost.noalias() = xd.transpose() * _Q_diag * xd;
            else
                cost.noalias() = xd.transpose() * _Q * xd;
        }
    }
}

void QuadraticStateCost::computeIntegralStateControlTerm(int k, const Eigen::Ref<const Eigen::VectorXd>& x_k,
                                                         const Eigen::Ref<const Eigen::VectorXd>& u_k, Eigen::Ref<Eigen::VectorXd> cost) const
{
    assert(_integral_form);
    assert(cost.size() == 1);

    cost[0] = 0;

    if (_zero_x_ref)
    {
        if (_diagonal_mode)
            cost[0] += x_k.transpose() * _Q_diag * x_k;
        else
            cost[0] += x_k.transpose() * _Q * x_k;
    }
    else
    {
        Eigen::VectorXd xd = x_k - _x_ref->getReferenceCached(k);
        if (_diagonal_mode)
            cost[0] += xd.transpose() * _Q_diag * xd;
        else
            cost[0] += xd.transpose() * _Q * xd;
    }
}

bool QuadraticStateCost::checkParameters(int state_dim, int control_dim, std::stringstream* issues) const
{
    bool success = true;

    if (_diagonal_mode_intentionally && _diagonal_mode)
    {
        if (_Q_diag.diagonal().size() != state_dim)
        {
            if (issues)
                *issues << "QuadraticStateCost: Diagonal matrix dimension of Q (" << _Q_diag.diagonal().size()
                        << ") does not match state vector dimension (" << state_dim << "); Please specify diagonal elements only." << std::endl;
            success = false;
        }
    }
    else
    {
        if (_Q.rows() != state_dim || _Q.cols() != state_dim)
        {
            if (issues)
                *issues << "QuadraticStateCost: Matrix dimension of Q (" << _Q.rows() << "x" << _Q.cols()
                        << ") does not match state vector dimension (" << state_dim << "); Please specify " << (state_dim * state_dim)
                        << " elements (Row-Major)." << std::endl;
            success = false;
        }
    }

    return success;
}

#ifdef MESSAGE_SUPPORT
bool QuadraticStateCost::fromMessage(const messages::QuadraticStateCost& message, std::stringstream* issues)
{
    _diagonal_mode = message.diagonal_only();

    if (_diagonal_mode)
    {
        if (!setWeightQ(
                Eigen::Map<const Eigen::Matrix<double, -1, 1>>(message.q().data(), message.q_size()).asDiagonal()))
        {
            *issues << "QuadraticStateCost: cannot set diagonal weight matrix Q.\n";
            return false;
        }
    }
    else
    {
        int p = std::sqrt(message.q_size());

        if (p * p != message.q_size())
        {
            *issues << "QuadraticStateCost: weight matrix Q is not square.\n";
            return false;
        }

        // weight matrix Q
        // if (p * p == message.full_discretization_ocp().q_size())
        //{
        if (!setWeightQ(Eigen::Map<const Eigen::Matrix<double, -1, -1, Eigen::RowMajor>>(message.q().data(), p, p)) && issues)
        {
            *issues << "QuadraticStateCost: weight matrix Q is not positive definite.\n";
            return false;
        }
        //}
        // else if (issues)
        //    *issues << "QuadraticStateCost: invalid size of weight matrix Q: it shout be [" << p << " x " << p << "].\n";
    }

    // others
    _lsq_form      = message.lsq_form();
    _integral_form = message.integral_form();
    return true;
}

void QuadraticStateCost::toMessage(messages::QuadraticStateCost& message) const
{
    // weight matrix Q
    if (_diagonal_mode_intentionally && _diagonal_mode)
    {
        Eigen::VectorXd Qdiag = _Q_diag.diagonal();
        message.mutable_q()->Resize(Qdiag.size(), 0);
        Eigen::Map<Eigen::VectorXd>(message.mutable_q()->mutable_data(), Qdiag.size()) = Qdiag;

        message.set_diagonal_only(true);
    }
    else
    {
        Eigen::MatrixXd Q = _Q;  // Q = _Q_sqrt.adjoint() * _Q_sqrt;
        message.mutable_q()->Resize(Q.rows() * Q.cols(), 0);
        Eigen::Map<Eigen::Matrix<double, -1, -1, Eigen::RowMajor>>(message.mutable_q()->mutable_data(), Q.rows(), Q.cols()) = Q;

        message.set_diagonal_only(false);
    }
    // others
    message.set_lsq_form(_lsq_form);
    message.set_integral_form(_integral_form);
}
#endif

}  // namespace corbo
