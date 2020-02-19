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

#include <corbo-optimal-control/functions/final_state_constraints.h>

#include <corbo-numerics/matrix_utilities.h>

#include <corbo-communication/utilities.h>

#include <cmath>

namespace corbo {

bool TerminalBall::setWeightS(const Eigen::Ref<const Eigen::MatrixXd>& S)
{
    if (S.size() == 0) return false;
    if (!is_square(S)) return false;

    // check if we have a diagonal matrix
    if (S.isDiagonal(1e-10)) return setWeightS(S.diagonal().asDiagonal());

    _diagonal_mode               = false;
    _diagonal_mode_intentionally = false;
    _S                           = S;
    return true;
}

bool TerminalBall::setWeightS(const Eigen::DiagonalMatrix<double, -1>& S)
{
    if (S.size() == 0) return false;

    _diagonal_mode_intentionally = true;
    _diagonal_mode               = true;
    _S_diag                      = S;
    _S                           = _S_diag.toDenseMatrix();
    return true;
}

void TerminalBall::computeNonIntegralStateTerm(int k, const Eigen::Ref<const Eigen::VectorXd>& x_k, Eigen::Ref<Eigen::VectorXd> cost) const
{
    assert(cost.size() == getNonIntegralStateTermDimension(k));

    if (_zero_x_ref)
    {
        if (_diagonal_mode)
            cost[0] = x_k.transpose() * _S_diag * x_k - _gamma;
        else
            cost[0] = x_k.transpose() * _S * x_k - _gamma;
    }
    else
    {
        Eigen::VectorXd xd = x_k - _x_ref->getReferenceCached(k);
        if (_diagonal_mode)
            cost[0] = xd.transpose() * _S_diag * xd - _gamma;
        else
            cost[0] = xd.transpose() * _S * xd - _gamma;
    }
}

bool TerminalBall::checkParameters(int state_dim, int control_dim, FinalStageCost::ConstPtr final_stage_cost, std::stringstream* issues) const
{
    bool success = true;

    if (_diagonal_mode_intentionally && _diagonal_mode)
    {
        if (_S.diagonal().size() != state_dim)
        {
            if (issues)
                *issues << "TerminalBall: Diagonal matrix dimension of S (" << _S_diag.diagonal().size()
                        << ") does not match state vector dimension (" << state_dim << "); Please specify diagonal elements only." << std::endl;
            success = false;
        }
    }
    else
    {
        if (_S.rows() != state_dim || _S.cols() != state_dim)
        {
            if (issues)
                *issues << "TerminalBall: Matrix dimension of S (" << _S.rows() << "x" << _S.cols() << ") does not match state vector dimension ("
                        << state_dim << "); Please specify " << (state_dim * state_dim) << " elements (Row-Major)." << std::endl;
            success = false;
        }
    }

    return success;
}

#ifdef MESSAGE_SUPPORT
bool TerminalBall::fromMessage(const messages::FinalStageConstraints& message, std::stringstream* issues)
{
    const messages::TerminalBall& msg = message.terminal_ball();

    _diagonal_mode = msg.s_diagonal_only();

    if (_diagonal_mode)
    {
        if (!setWeightS(Eigen::Map<const Eigen::Matrix<double, -1, 1>>(msg.s().data(), msg.s_size()).asDiagonal()))
        {
            *issues << "TerminalBall: cannot set diagonal weight matrix S.\n";
            return false;
        }
    }
    else
    {
        if (!is_square(msg.s_size()))
        {
            *issues << "TerminalBall: weight matrix S is not square.\n";
            return false;
        }
        int p = std::sqrt(msg.s_size());

        if (!setWeightS(Eigen::Map<const Eigen::Matrix<double, -1, -1, Eigen::RowMajor>>(msg.s().data(), p, p)) && issues)
        {
            *issues << "TerminalBall: cannot set weight matrix S.\n";
            return false;
        }
    }

    // gamma
    _gamma = msg.gamma();
    return true;
}

void TerminalBall::toMessage(messages::FinalStageConstraints& message) const
{
    messages::TerminalBall* msg = message.mutable_terminal_ball();

    // weight matrix S
    if (_diagonal_mode_intentionally && _diagonal_mode)
    {
        Eigen::VectorXd Sdiag = _S_diag.diagonal();
        msg->mutable_s()->Resize(Sdiag.size(), 0);
        Eigen::Map<Eigen::VectorXd>(msg->mutable_s()->mutable_data(), Sdiag.size()) = Sdiag;

        msg->set_s_diagonal_only(true);
    }
    else
    {
        msg->mutable_s()->Resize(_S.rows() * _S.cols(), 0);
        Eigen::Map<Eigen::Matrix<double, -1, -1, Eigen::RowMajor>>(msg->mutable_s()->mutable_data(), _S.rows(), _S.cols()) = _S;

        msg->set_s_diagonal_only(false);
    }

    // gamma
    msg->set_gamma(_gamma);
}
#endif

bool TerminalBallInheritFromCost::update(int n, double t, ReferenceTrajectoryInterface& xref, ReferenceTrajectoryInterface& uref,
                                         ReferenceTrajectoryInterface* sref, bool single_dt, const Eigen::VectorXd& x0,
                                         FinalStageCost::ConstPtr final_state_cost, StagePreprocessor::Ptr stage_preprocessor,
                                         const std::vector<double>& dts, const DiscretizationGridInterface* grid)
{
    // call parent update method
    bool dimension_modified = TerminalBall::update(n, t, xref, uref, sref, single_dt, x0, final_state_cost, stage_preprocessor, dts, grid);

    BaseQuadraticFinalStateCost::ConstPtr quadratic_cost = std::dynamic_pointer_cast<const BaseQuadraticFinalStateCost>(final_state_cost);
    if (quadratic_cost)
    {
        if (quadratic_cost->getWeightQf().rows() == 0)
        {
            PRINT_ERROR("TerminalBallInheritFromCost::update(): weight matrix obtained from terminal cost function has zero size!");
            TerminalBall::setWeightS(Eigen::MatrixXd::Zero(xref.getDimension(), xref.getDimension()));
        }
        TerminalBall::setWeightS(quadratic_cost->getWeightQf());
    }
    else
    {
        PRINT_ERROR("TerminalBallInheritFromCost::update(): this constraint requires quadratic final cost! Setting final weight matrix to zero!");
        TerminalBall::setWeightS(Eigen::MatrixXd::Zero(xref.getDimension(), xref.getDimension()));
    }
    return dimension_modified;
}

void TerminalBallInheritFromCost::computeNonIntegralStateTerm(int k, const Eigen::Ref<const Eigen::VectorXd>& x_k,
                                                              Eigen::Ref<Eigen::VectorXd> cost) const
{
    assert(cost.size() == getNonIntegralStateTermDimension(k));
    assert(_S_diag.diagonal().size() == x_k.size());
    assert(_S.rows() == x_k.size());

    if (_zero_x_ref)
    {
        if (_diagonal_mode)
            cost[0] = x_k.transpose() * _S_diag * x_k - _gamma;
        else
            cost[0] = x_k.transpose() * _S * x_k - _gamma;
    }
    else
    {
        Eigen::VectorXd xd = x_k - _x_ref->getReferenceCached(k);
        if (_diagonal_mode)
            cost[0] = xd.transpose() * _S_diag * xd - _gamma;
        else
            cost[0] = xd.transpose() * _S * xd - _gamma;
    }
}

bool TerminalBallInheritFromCost::checkParameters(int state_dim, int control_dim, FinalStageCost::ConstPtr final_stage_cost,
                                                  std::stringstream* issues) const
{
    bool success = true;

    BaseQuadraticFinalStateCost::ConstPtr quadratic_cost = std::dynamic_pointer_cast<const BaseQuadraticFinalStateCost>(final_stage_cost);
    if (!quadratic_cost)
    {
        if (issues)
            *issues << "TerminalBallInheritFromCost: This final cost requires a complementary quadratic final/terminal cost function." << std::endl;
        success = false;
    }

    // success = success && TerminalBall::checkParameters(state_dim, control_dim, final_stage_cost, issues);

    return success;
}

#ifdef MESSAGE_SUPPORT
bool TerminalBallInheritFromCost::fromMessage(const messages::FinalStageConstraints& message, std::stringstream* issues)
{
    const messages::TerminalBallInheritFromCost& msg = message.terminal_ball_from_cost();

    _gamma = msg.gamma();
    return true;
}

void TerminalBallInheritFromCost::toMessage(messages::FinalStageConstraints& message) const
{
    messages::TerminalBallInheritFromCost* msg = message.mutable_terminal_ball_from_cost();
    msg->set_gamma(msg->gamma());
}
#endif

}  // namespace corbo
