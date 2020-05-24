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

#ifndef SRC_OPTIMAL_CONTROL_INCLUDE_CORBO_OPTIMAL_CONTROL_FUNCTIONS_QUADRATIC_CONTROL_COST_H_
#define SRC_OPTIMAL_CONTROL_INCLUDE_CORBO_OPTIMAL_CONTROL_FUNCTIONS_QUADRATIC_CONTROL_COST_H_

#include <corbo-core/reference_trajectory.h>
#include <corbo-optimal-control/functions/stage_functions.h>

#include <memory>

namespace corbo {

class QuadraticControlCost : public StageCost
{
 public:
    QuadraticControlCost() { _R_sqrt = Eigen::MatrixXd::Constant(1, 1, 1); }

    QuadraticControlCost(const Eigen::Ref<const Eigen::MatrixXd>& R, bool integral_form = false, bool lsq_form = false)
        : _integral_form(integral_form), _lsq_form(lsq_form)
    {
        setWeightR(R);
    }

    StageCost::Ptr getInstance() const override { return std::make_shared<QuadraticControlCost>(); }

    bool setWeightR(const Eigen::Ref<const Eigen::MatrixXd>& R);
    bool setWeightR(const Eigen::DiagonalMatrix<double, -1>& R);
    void setLsqForm(bool lsq_form) { _lsq_form = lsq_form; }
    void setIntegralForm(bool integral_form) { _integral_form = integral_form; }

    bool hasIntegralTerms(int k) const override { return _integral_form; }
    bool hasNonIntegralTerms(int k) const override { return !_integral_form; }

    int getNonIntegralControlTermDimension(int k) const override { return !_integral_form ? (_lsq_form ? _R.rows() : 1) : 0; }
    int getIntegralStateControlTermDimension(int k) const override { return _integral_form ? 1 : 0; }

    bool isLsqFormNonIntegralControlTerm(int k) const override { return !_integral_form && _lsq_form; }

    void computeNonIntegralControlTerm(int k, const Eigen::Ref<const Eigen::VectorXd>& u_k, Eigen::Ref<Eigen::VectorXd> cost) const override;
    void computeIntegralStateControlTerm(int k, const Eigen::Ref<const Eigen::VectorXd>& x_k, const Eigen::Ref<const Eigen::VectorXd>& u_k,
                                         Eigen::Ref<Eigen::VectorXd> cost) const override;

    bool update(int n, double t, ReferenceTrajectoryInterface& xref, ReferenceTrajectoryInterface& uref, ReferenceTrajectoryInterface* sref,
                bool single_dt, const Eigen::VectorXd& x0, StagePreprocessor::Ptr stage_preprocessor, const std::vector<double>& dts,
                const DiscretizationGridInterface* /*grid*/) override
    {
        _u_ref = &uref;

        _zero_u_ref = _u_ref->isZero();

        return false;
    }

    bool checkParameters(int state_dim, int control_dim, std::stringstream* issues) const override;

#ifdef MESSAGE_SUPPORT
    virtual bool fromMessage(const messages::QuadraticControlCost& message, std::stringstream* issues);
    virtual void toMessage(messages::QuadraticControlCost& message) const;

    bool fromMessage(const messages::StageCost& message, std::stringstream* issues) override
    {
        return fromMessage(message.quadratic_control_cost(), issues);
    }
    void toMessage(messages::StageCost& message) const override { toMessage(*message.mutable_quadratic_control_cost()); }
#endif

 protected:
    Eigen::MatrixXd _R_sqrt;
    Eigen::MatrixXd _R;
    Eigen::DiagonalMatrix<double, -1> _R_diag_sqrt;
    Eigen::DiagonalMatrix<double, -1> _R_diag;
    bool _diagonal_mode               = false;
    bool _diagonal_mode_intentionally = false;
    bool _lsq_form                    = true;
    bool _integral_form               = false;

    const ReferenceTrajectoryInterface* _u_ref = nullptr;
    bool _zero_u_ref                           = false;
};
FACTORY_REGISTER_STAGE_COST(QuadraticControlCost)

}  // namespace corbo

#endif  // SRC_OPTIMAL_CONTROL_INCLUDE_CORBO_OPTIMAL_CONTROL_FUNCTIONS_QUADRATIC_CONTROL_COST_H_
