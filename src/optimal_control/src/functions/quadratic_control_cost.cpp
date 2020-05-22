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

#include <corbo-optimal-control/functions/quadratic_control_cost.h>

#include <corbo-numerics/matrix_utilities.h>

#include <cmath>

namespace corbo {

bool QuadraticControlCost::setWeightR(const Eigen::Ref<const Eigen::MatrixXd>& R)
{
    _diagonal_mode_intentionally = false;

    if (!is_square(R)) return false;

    // check if we have a diagonal matrix
    if (R.isDiagonal(1e-10)) return setWeightR(R.diagonal().asDiagonal());

    _diagonal_mode = false;

    _R = R;  // also store original R
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

bool QuadraticControlCost::setWeightR(const Eigen::DiagonalMatrix<double, -1>& R)
{
    _diagonal_mode_intentionally = true;
    _diagonal_mode               = true;
    _R_diag                      = R;
    _R_diag_sqrt                 = R.diagonal().cwiseSqrt().asDiagonal();
    return true;
}

void QuadraticControlCost::computeNonIntegralControlTerm(int k, const Eigen::Ref<const Eigen::VectorXd>& u_k, Eigen::Ref<Eigen::VectorXd> cost) const
{
    assert(!_integral_form);
    assert(cost.size() == getNonIntegralControlTermDimension(k));

    if (_lsq_form)
    {
        if (_zero_u_ref)
        {
            if (_diagonal_mode)
                cost.noalias() = _R_diag_sqrt * u_k;
            else
                cost.noalias() = _R_sqrt * u_k;
        }
        else
        {
            Eigen::VectorXd ud = u_k - _u_ref->getReferenceCached(k);
            if (_diagonal_mode)
                cost.noalias() = _R_diag_sqrt * ud;
            else
                cost.noalias() = _R_sqrt * ud;
        }
    }
    else
    {
        if (_zero_u_ref)
        {
            if (_diagonal_mode)
                cost.noalias() = u_k.transpose() * _R_diag * u_k;
            else
                cost.noalias() = u_k.transpose() * _R * u_k;
        }
        else
        {
            Eigen::VectorXd ud = u_k - _u_ref->getReferenceCached(k);
            if (_diagonal_mode)
                cost.noalias() = ud.transpose() * _R_diag * ud;
            else
                cost.noalias() = ud.transpose() * _R * ud;
        }
    }
}

void QuadraticControlCost::computeIntegralStateControlTerm(int k, const Eigen::Ref<const Eigen::VectorXd>& x_k,
                                                           const Eigen::Ref<const Eigen::VectorXd>& u_k, Eigen::Ref<Eigen::VectorXd> cost) const
{
    assert(_integral_form);
    assert(cost.size() == 1);

    cost[0] = 0;

    if (_zero_u_ref)
    {
        if (_diagonal_mode)
            cost[0] += u_k.transpose() * _R_diag * u_k;
        else
            cost[0] += u_k.transpose() * _R * u_k;
    }
    else
    {
        Eigen::VectorXd ud = u_k - _u_ref->getReferenceCached(k);
        if (_diagonal_mode)
            cost[0] += ud.transpose() * _R_diag * ud;
        else
            cost[0] += ud.transpose() * _R * ud;
    }
}

bool QuadraticControlCost::checkParameters(int state_dim, int control_dim, std::stringstream* issues) const
{
    bool success = true;

    if (_diagonal_mode_intentionally && _diagonal_mode)
    {
        if (_R_diag.diagonal().size() != control_dim)
        {
            if (issues)
                *issues << "QuadraticControlCost: diagonal matrix dimension of R (" << _R_diag.diagonal().size()
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
                *issues << "QuadraticControlCost: Matrix dimension of R (" << _R.rows() << "x" << _R.cols()
                        << ") does not match control input vector dimension (" << control_dim << "); Please specify " << (control_dim * control_dim)
                        << " elements (Row-Major)." << std::endl;
            success = false;
        }
    }

    return success;
}

#ifdef MESSAGE_SUPPORT
bool QuadraticControlCost::fromMessage(const messages::QuadraticControlCost& message, std::stringstream* issues)
{
    _diagonal_mode = message.diagonal_only();

    if (_diagonal_mode)
    {
        if (!setWeightR(Eigen::Map<const Eigen::Matrix<double, -1, 1>>(message.r().data(), message.r_size()).asDiagonal()))
        {
            *issues << "QuadraticControlCost: cannot set diagonal weight matrix R.\n";
            return false;
        }
    }
    else
    {
        int q = std::sqrt(message.r_size());

        if (q * q != message.r_size())
        {
            *issues << "QuadraticControlCost: weight matrix R is not square.\n";
            return false;
        }

        // weight matrix Q
        // if (p * p == message.full_discretization_ocp().q_size())
        //{
        if (!setWeightR(Eigen::Map<const Eigen::Matrix<double, -1, -1, Eigen::RowMajor>>(message.r().data(), q, q)) && issues)
        {
            *issues << "QuadraticControlCost: weight matrix R is not positive definite.\n";
            return false;
        }
        //}
        // else if (issues)
        //    *issues << "QuadraticControlCost: invalid size of weight matrix Q: it shout be [" << p << " x " << p << "].\n";
    }

    // others
    _lsq_form      = message.lsq_form();
    _integral_form = message.integral_form();

    return true;
}

void QuadraticControlCost::toMessage(messages::QuadraticControlCost& message) const
{
    // weight matrix R
    if (_diagonal_mode_intentionally && _diagonal_mode)
    {
        Eigen::VectorXd Rdiag = _R_diag.diagonal();
        message.mutable_r()->Resize(Rdiag.size(), 0);
        Eigen::Map<Eigen::VectorXd>(message.mutable_r()->mutable_data(), Rdiag.size()) = Rdiag;

        message.set_diagonal_only(true);
    }
    else
    {
        Eigen::MatrixXd R = _R;  // R = _R_sqrt.adjoint() * _R_sqrt;
        message.mutable_r()->Resize(R.rows() * R.cols(), 0);
        Eigen::Map<Eigen::Matrix<double, -1, -1, Eigen::RowMajor>>(message.mutable_r()->mutable_data(), R.rows(),
                                                                   R.cols()) = R;

        message.set_diagonal_only(false);
    }

    // others
    message.set_lsq_form(_lsq_form);
    message.set_integral_form(_integral_form);
}
#endif

}  // namespace corbo
