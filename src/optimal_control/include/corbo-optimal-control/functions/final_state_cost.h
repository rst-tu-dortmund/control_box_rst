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

#ifndef SRC_OPTIMAL_CONTROL_INCLUDE_CORBO_OPTIMAL_CONTROL_FUNCTIONS_FINAL_STATE_COST_H_
#define SRC_OPTIMAL_CONTROL_INCLUDE_CORBO_OPTIMAL_CONTROL_FUNCTIONS_FINAL_STATE_COST_H_

#include <corbo-optimal-control/functions/quadratic_state_cost.h>

#include <corbo-core/reference_trajectory.h>
#include <corbo-optimal-control/functions/stage_functions.h>
#include <corbo-systems/system_dynamics_interface.h>

#include <memory>

namespace corbo {

class BaseQuadraticFinalStateCost : public FinalStageCost
{
 public:
    using Ptr      = std::shared_ptr<BaseQuadraticFinalStateCost>;
    using ConstPtr = std::shared_ptr<const BaseQuadraticFinalStateCost>;

    virtual const Eigen::MatrixXd& getWeightQf() const = 0;
};

class QuadraticFinalStateCost : public BaseQuadraticFinalStateCost
{
 public:
    using Ptr = std::shared_ptr<QuadraticFinalStateCost>;

    QuadraticFinalStateCost() { _Qf_sqrt = Eigen::MatrixXd::Constant(1, 1, 1); }

    QuadraticFinalStateCost(const Eigen::Ref<const Eigen::MatrixXd>& Qf, bool lsq_form) : _lsq_form(lsq_form) { setWeightQf(Qf); }

    FinalStageCost::Ptr getInstance() const override { return std::make_shared<QuadraticFinalStateCost>(); }

    int getNonIntegralStateTermDimension(int k) const override { return _lsq_form ? _Qf.rows() : 1; }
    bool isLsqFormNonIntegralStateTerm(int k) const override { return _lsq_form; }

    bool setWeightQf(const Eigen::Ref<const Eigen::MatrixXd>& Qf);
    bool setWeightQf(const Eigen::DiagonalMatrix<double, -1>& Qf);

    const Eigen::MatrixXd& getWeightQf() const override { return _Qf; }

    void setLsqForm(bool lsq_form) { _lsq_form = lsq_form; }

    void computeNonIntegralStateTerm(int k, const Eigen::Ref<const Eigen::VectorXd>& x_k, Eigen::Ref<Eigen::VectorXd> cost) const override;

    bool update(int n, double t, ReferenceTrajectoryInterface& xref, ReferenceTrajectoryInterface& uref, ReferenceTrajectoryInterface* sref,
                bool single_dt, const Eigen::VectorXd& x0, StagePreprocessor::Ptr stage_preprocessor, const std::vector<double>& dts,
                const DiscretizationGridInterface* grid) override
    {
        bool dimension_modified = BaseQuadraticFinalStateCost::update(n, t, xref, uref, sref, single_dt, x0, stage_preprocessor, dts, grid);

        _x_ref      = &xref;
        _zero_x_ref = _x_ref->isZero();

        return dimension_modified;
    }

    bool checkParameters(int state_dim, int control_dim, std::stringstream* issues) const override;

#ifdef MESSAGE_SUPPORT
    bool fromMessage(const messages::FinalStageCost& message, std::stringstream* issues) override;
    void toMessage(messages::FinalStageCost& message) const override;
#endif

 protected:
    Eigen::MatrixXd _Qf_sqrt;
    Eigen::MatrixXd _Qf;
    Eigen::DiagonalMatrix<double, -1> _Qf_diag_sqrt;
    Eigen::DiagonalMatrix<double, -1> _Qf_diag;
    bool _diagonal_mode               = false;
    bool _diagonal_mode_intentionally = false;
    bool _lsq_form                    = true;

    const ReferenceTrajectoryInterface* _x_ref = nullptr;
    bool _zero_x_ref                           = false;
};
FACTORY_REGISTER_FINAL_STAGE_COST(QuadraticFinalStateCost)

class QuadraticFinalStateCostRiccati : public BaseQuadraticFinalStateCost
{
 public:
    QuadraticFinalStateCostRiccati() { _Qf_sqrt = Eigen::MatrixXd::Constant(1, 1, 1); }

    QuadraticFinalStateCostRiccati(SystemDynamicsInterface::Ptr dynamics, const Eigen::Ref<const Eigen::MatrixXd>& Q,
                                   const Eigen::Ref<const Eigen::MatrixXd>& R)
    {
        setSystemDynamics(dynamics);
        setWeightQ(Q);
        setWeightR(R);
    }

    FinalStageCost::Ptr getInstance() const override { return std::make_shared<QuadraticFinalStateCostRiccati>(); }

    int getNonIntegralStateTermDimension(int k) const override { return 1; }

    void setSystemDynamics(SystemDynamicsInterface::Ptr dynamics)
    {
        _are_solved_before = false;
        _dynamics          = dynamics;
    }

    bool setWeightQ(const Eigen::Ref<const Eigen::MatrixXd>& Q);
    bool setWeightQ(const Eigen::DiagonalMatrix<double, -1>& Q);
    bool setWeightR(const Eigen::Ref<const Eigen::MatrixXd>& R);
    bool setWeightR(const Eigen::DiagonalMatrix<double, -1>& R);

    void setLsqForm(bool lsq_form) { _lsq_form = lsq_form; }

    const Eigen::MatrixXd& getWeightQf() const override { return _Qf; }

    void computeNonIntegralStateTerm(int k, const Eigen::Ref<const Eigen::VectorXd>& x_k, Eigen::Ref<Eigen::VectorXd> cost) const override;

    bool update(int n, double t, ReferenceTrajectoryInterface& xref, ReferenceTrajectoryInterface& uref, ReferenceTrajectoryInterface* sref,
                bool single_dt, const Eigen::VectorXd& x0, StagePreprocessor::Ptr stage_preprocessor, const std::vector<double>& dts,
                const DiscretizationGridInterface* /*grid*/) override;

    bool checkParameters(int state_dim, int control_dim, std::stringstream* issues) const override;

#ifdef MESSAGE_SUPPORT
    bool fromMessage(const messages::FinalStageCost& message, std::stringstream* issues) override;
    void toMessage(messages::FinalStageCost& message) const override;
#endif

 protected:
    bool computeWeightQfSqrt();

    Eigen::MatrixXd _Q;
    Eigen::MatrixXd _R;
    bool _Q_diagonal_mode_intentionally = false;
    bool _R_diagonal_mode_intentionally = false;

    Eigen::MatrixXd _Qf_sqrt;
    Eigen::MatrixXd _Qf;
    SystemDynamicsInterface::Ptr _dynamics;
    Eigen::VectorXd _steady_state_x;
    Eigen::VectorXd _steady_state_u;
    bool _lsq_form          = true;
    bool _are_solved_before = false;  // algebraic riccati equation solved before (e.g., for linear systems solving once is sufficient)

    const ReferenceTrajectoryInterface* _x_ref = nullptr;
    bool _zero_x_ref                           = false;
};
FACTORY_REGISTER_FINAL_STAGE_COST(QuadraticFinalStateCostRiccati)

}  // namespace corbo

#endif  // SRC_OPTIMAL_CONTROL_INCLUDE_CORBO_OPTIMAL_CONTROL_FUNCTIONS_FINAL_STATE_COST_H_
