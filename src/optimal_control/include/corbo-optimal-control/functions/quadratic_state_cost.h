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

#ifndef SRC_OPTIMAL_CONTROL_INCLUDE_CORBO_OPTIMAL_CONTROL_FUNCTIONS_QUADRATIC_STATE_COST_H_
#define SRC_OPTIMAL_CONTROL_INCLUDE_CORBO_OPTIMAL_CONTROL_FUNCTIONS_QUADRATIC_STATE_COST_H_

#include <corbo-core/reference_trajectory.h>
#include <corbo-optimal-control/functions/stage_functions.h>

#include <memory>

namespace corbo {

class QuadraticStateCost : public StageCost
{
 public:
    QuadraticStateCost() { _Q_sqrt = Eigen::MatrixXd::Constant(1, 1, 1); }

    QuadraticStateCost(const Eigen::Ref<const Eigen::MatrixXd>& Q, bool integral_form, bool lsq_form = false)
        : _integral_form(integral_form), _lsq_form(lsq_form)
    {
        setWeightQ(Q);
    }

    StageCost::Ptr getInstance() const override { return std::make_shared<QuadraticStateCost>(); }

    bool hasIntegralTerms(int k) const override { return _integral_form; }
    bool hasNonIntegralTerms(int k) const override { return !_integral_form; }

    int getNonIntegralStateTermDimension(int k) const override { return !_integral_form ? (_lsq_form ? _Q.rows() : 1) : 0; }
    int getIntegralStateControlTermDimension(int k) const override { return _integral_form ? 1 : 0; }

    bool isLsqFormNonIntegralStateTerm(int k) const override { return !_integral_form && _lsq_form; }

    bool setWeightQ(const Eigen::Ref<const Eigen::MatrixXd>& Q);
    bool setWeightQ(const Eigen::DiagonalMatrix<double, -1>& Q);
    void setLsqForm(bool lsq_form) { _lsq_form = lsq_form; }
    void setIntegralForm(bool integral_form) { _integral_form = integral_form; }

    void computeNonIntegralStateTerm(int k, const Eigen::Ref<const Eigen::VectorXd>& x_k, Eigen::Ref<Eigen::VectorXd> cost) const override;
    void computeIntegralStateControlTerm(int k, const Eigen::Ref<const Eigen::VectorXd>& x_k, const Eigen::Ref<const Eigen::VectorXd>& u_k,
                                         Eigen::Ref<Eigen::VectorXd> cost) const override;

    bool update(int n, double t, ReferenceTrajectoryInterface& xref, ReferenceTrajectoryInterface& uref, ReferenceTrajectoryInterface* sref,
                bool single_dt, const Eigen::VectorXd& x0, StagePreprocessor::Ptr stage_preprocessor, const std::vector<double>& dts,
                const DiscretizationGridInterface* /*grid*/) override
    {
        _x_ref = &xref;

        _zero_x_ref = _x_ref->isZero();

        return false;
    }

    bool checkParameters(int state_dim, int control_dim, std::stringstream* issues) const override;

#ifdef MESSAGE_SUPPORT
    virtual bool fromMessage(const messages::QuadraticStateCost& message, std::stringstream* issues);
    virtual void toMessage(messages::QuadraticStateCost& message) const;

    bool fromMessage(const messages::StageCost& message, std::stringstream* issues) override
    {
        return fromMessage(message.quadratic_state_cost(), issues);
    }
    void toMessage(messages::StageCost& message) const override { toMessage(*message.mutable_quadratic_state_cost()); }
#endif

 protected:
    Eigen::MatrixXd _Q_sqrt;
    Eigen::MatrixXd _Q;
    Eigen::DiagonalMatrix<double, -1> _Q_diag_sqrt;
    Eigen::DiagonalMatrix<double, -1> _Q_diag;
    bool _diagonal_mode               = false;
    bool _diagonal_mode_intentionally = false;
    bool _lsq_form                    = true;
    bool _integral_form               = false;

    const ReferenceTrajectoryInterface* _x_ref = nullptr;
    bool _zero_x_ref                           = false;
};
FACTORY_REGISTER_STAGE_COST(QuadraticStateCost)

}  // namespace corbo

#endif  // SRC_OPTIMAL_CONTROL_INCLUDE_CORBO_OPTIMAL_CONTROL_FUNCTIONS_QUADRATIC_STATE_COST_H_
