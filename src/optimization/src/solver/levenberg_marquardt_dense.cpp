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

#include <corbo-optimization/solver/levenberg_marquardt_dense.h>

#include <corbo-core/console.h>

#include <Eigen/Cholesky>

namespace corbo {

bool LevenbergMarquardtDense::initialize(OptimizationProblemInterface* problem)
{
    if (problem && !problem->isLeastSquaresProblem())
    {
        PRINT_ERROR("LevenbergMarquardtDense(): cannot handle non-least-squares objectives or LS objectives in non-LS form.");
        return false;
    }

    return true;
}

SolverStatus LevenbergMarquardtDense::solve(OptimizationProblemInterface& problem, bool new_structure, bool new_run, double* obj_value)
{
    if (obj_value) *obj_value = -1;  // set to invalid first

    if (new_structure)
    {
        if (!problem.isLeastSquaresProblem())
        {
            PRINT_ERROR("LevenbergMarquardtDense(): cannot handle non-least-squares objectives or LS objectives in non-LS form.");
            return SolverStatus::Error;
        }

        // get dimension for the value/cost vector
        _obj_dim           = problem.getLsqObjectiveDimension();
        _eq_dim            = problem.getEqualityDimension();
        _ineq_dim          = problem.getInequalityDimension();
        _finite_bounds_dim = problem.finiteCombinedBoundsDimension();
        _val_dim           = _obj_dim + _eq_dim + _ineq_dim + _finite_bounds_dim;
        _param_dim         = problem.getParameterDimension();

        if (_param_dim == 0)
        {
            PRINT_WARNING("LevenbergMarquardtDense: problem has zero dimension");
            return SolverStatus::Error;
        }

        _values.resize(_val_dim);
        _jacobian.resize(_val_dim, _param_dim);
        _jacobian.setZero();
        _hessian.resize(_param_dim, _param_dim);
        _delta.resize(_param_dim);
        _rhs.resize(_param_dim);
    }

    // adapt weights
    if (new_run)
        resetWeights();
    else
        adaptWeights();

    // compute complete value vector
    computeValues(problem);

    // compute jacobian
    computeJacobian(problem);

    // construct quasi-newton hessian approximation
    _hessian.noalias() = _jacobian.transpose() * _jacobian;

    // construct right-hand-side of the approxiamted linear system
    _rhs.noalias() = _jacobian.transpose() * -_values;

    // convergence consts
    constexpr const double eps1 = 1e-5;
    constexpr const double eps2 = 1e-5;
    constexpr const double eps3 = 1e-5;
    constexpr const double eps4 = 0;

    // lm variables
    unsigned int v = 2;
    double tau     = 1e-5;

    constexpr const double goodStepUpperScale = 2. / 3.;
    constexpr const double goodStepLowerScale = 1. / 3.;

    bool stop = (_rhs.lpNorm<Eigen::Infinity>() <= eps1);

    double mu = tau * _hessian.diagonal().maxCoeff();
    if (mu < 0) mu = 0;

    double rho = 0;

    // get old chi2
    double chi2_old = _values.squaredNorm();
    if (obj_value) *obj_value = chi2_old;

    // start levenberg marquardt optimization loop
    for (int k = 0; k < _iterations; ++k)
    {
        do
        {
            // augment hessian diagonal with damping factor
            _hessian.diagonal().array() += mu;

            // solve linear system
            _delta = _hessian.ldlt().solve(_rhs);

            // if (delta.norm() <= eps2 * (x.norm() + eps2) stop = true; // x -> values of vertices, not constructed until now. modify check
            if (_delta.norm() <= eps2)
                stop = true;
            else
            {
                // backup current parameter values
                problem.backupParameters();

                // apply new delta/increment to vertices
                problem.applyIncrement(_delta);

                // calculate new costs/values
                computeValues(problem);

                // get new chi2
                double chi2_new = _values.squaredNorm();

                rho = (chi2_old - chi2_new) / (_delta.transpose() * (mu * _delta + _rhs));

                if (rho > 0 && !std::isnan(chi2_new) && !std::isinf(chi2_new))  // adaptation of mu to find a sufficient descent step
                {
                    stop = (std::sqrt(chi2_old) - std::sqrt(chi2_new) < eps4 * std::sqrt(chi2_old));

                    // accept update and discard backup
                    problem.discardBackupParameters();

                    if (!stop && k < _iterations - 1)  // do not recalculate jacobian and hessian in the last iteration
                    {
                        // get new cost/value vector
                        // computeValues(problem); // redundant

                        // calculate new jacobian
                        computeJacobian(problem);

                        // get new hessian
                        _hessian.noalias() = _jacobian.transpose() * _jacobian;

                        // calculate new right hand side
                        _rhs.noalias() = _jacobian.transpose() * -_values;

                        stop = stop || (_rhs.lpNorm<Eigen::Infinity>() <= eps1);

                        double alpha       = std::min(goodStepUpperScale, 1 - std::pow((2 * rho - 1), 3));
                        double scaleFactor = std::max(goodStepLowerScale, alpha);
                        mu *= scaleFactor;
                        v = 2;
                    }

                    chi2_old = chi2_new;
                    if (obj_value) *obj_value = chi2_old;
                }
                else
                {
                    // restore from backup
                    problem.restoreBackupParameters(false);
                    // hessian.diagonal().array() -= mu; // this should be done here, but it works in the matlab version without for some
                    // time now.

                    mu = mu * v;
                    v  = 2 * v;
                }
            }
        } while (rho <= 0 && !stop);
        stop = (_values.norm() <= eps3);
    }
    return (stop || rho <= 0) ? SolverStatus::Converged
                              : SolverStatus::EarlyTerminated;  // TODO(roesmann): proper solver status (especially w.r.t. rho)
}

void LevenbergMarquardtDense::computeValues(OptimizationProblemInterface& problem)
{
    int idx = 0;
    if (_obj_dim > 0)
    {
        problem.computeValuesLsqObjective(_values.segment(idx, _obj_dim));
        idx += _obj_dim;
    }
    if (_eq_dim > 0)
    {
        problem.computeValuesEquality(_values.segment(idx, _eq_dim));
        _values.segment(idx, _eq_dim) *= _weight_eq;
        idx += _eq_dim;
    }
    if (_ineq_dim > 0)
    {
        problem.computeValuesActiveInequality(_values.segment(idx, _ineq_dim), _weight_ineq);
        idx += _ineq_dim;
    }
    if (_finite_bounds_dim > 0)
    {
        problem.computeDistanceFiniteCombinedBounds(_values.segment(idx, _finite_bounds_dim));
        _values.segment(idx, _finite_bounds_dim) *= _weight_bounds;
    }
}

void LevenbergMarquardtDense::computeJacobian(OptimizationProblemInterface& problem)
{
    int idx = 0;
    if (_obj_dim > 0)
    {
        problem.computeDenseJacobianLsqObjective(_jacobian.block(idx, 0, _obj_dim, _param_dim));
        idx += _obj_dim;
    }
    if (_eq_dim > 0)
    {
        problem.computeDenseJacobianEqualities(_jacobian.block(idx, 0, _eq_dim, _param_dim));
        _jacobian.block(idx, 0, _eq_dim, _param_dim) *= _weight_eq;  // TODO(roesmann) include scalar multiplier to OptimizationProblemInterface
        idx += _eq_dim;
    }
    if (_ineq_dim > 0)
    {
        // compute only for active inequalities in order to apply quadratic penalties
        problem.computeDenseJacobianActiveInequalities(_jacobian.block(idx, 0, _ineq_dim, _param_dim), _weight_ineq);
        idx += _ineq_dim;
    }
    if (_finite_bounds_dim > 0)
    {
        problem.computeDenseJacobianFiniteCombinedBounds(_jacobian.block(idx, 0, _finite_bounds_dim, _param_dim), _weight_bounds);
    }
}

void LevenbergMarquardtDense::setPenaltyWeights(double weight_eq, double weight_ineq, double weight_bounds)
{
    _weight_init_eq     = weight_eq;
    _weight_init_ineq   = weight_ineq;
    _weight_init_bounds = weight_bounds;

    PRINT_WARNING_COND(_weight_init_eq > _weight_adapt_max_eq, "LevenbergMarquardtDense(): weight_eq > max_eq");
    PRINT_WARNING_COND(_weight_init_ineq > _weight_adapt_max_ineq, "LevenbergMarquardtDense(): weight_ineq > max_ineq");
    PRINT_WARNING_COND(_weight_init_bounds > _weight_adapt_max_bounds, "LevenbergMarquardtDense(): weight_bounds > max_bounds");
}

void LevenbergMarquardtDense::setWeightAdapation(double factor_eq, double factor_ineq, double factor_bounds, double max_eq, double max_ineq,
                                                 double max_bounds)
{
    _weight_adapt_factor_eq     = factor_eq;
    _weight_adapt_factor_ineq   = factor_ineq;
    _weight_adapt_factor_bounds = factor_bounds;
    _weight_adapt_max_eq        = max_eq;
    _weight_adapt_max_ineq      = max_ineq;
    _weight_adapt_max_bounds    = max_bounds;
}

void LevenbergMarquardtDense::resetWeights()
{
    _weight_eq     = _weight_init_eq;
    _weight_ineq   = _weight_init_ineq;
    _weight_bounds = _weight_init_bounds;
}

void LevenbergMarquardtDense::adaptWeights()
{
    _weight_eq *= _weight_adapt_factor_eq;
    if (_weight_eq > _weight_adapt_max_eq) _weight_eq = _weight_adapt_max_eq;

    _weight_ineq *= _weight_adapt_factor_ineq;
    if (_weight_ineq > _weight_adapt_max_ineq) _weight_ineq = _weight_adapt_max_ineq;

    _weight_bounds *= _weight_adapt_factor_bounds;
    if (_weight_bounds > _weight_adapt_max_bounds) _weight_bounds = _weight_adapt_max_bounds;
}

void LevenbergMarquardtDense::clear()
{
    _param_dim         = 0;
    _obj_dim           = 0;
    _eq_dim            = 0;
    _ineq_dim          = 0;
    _finite_bounds_dim = 0;
    _val_dim           = 0;

    _weight_eq     = _weight_init_eq;
    _weight_ineq   = _weight_init_ineq;
    _weight_bounds = _weight_init_bounds;
}

#ifdef MESSAGE_SUPPORT
void LevenbergMarquardtDense::toMessage(corbo::messages::NlpSolver& message) const
{
    message.mutable_levenberg_marquardt_dense()->set_iterations(_iterations);
    message.mutable_levenberg_marquardt_dense()->set_weight_eq(_weight_init_eq);
    message.mutable_levenberg_marquardt_dense()->set_weight_ineq(_weight_init_ineq);
    message.mutable_levenberg_marquardt_dense()->set_weight_bounds(_weight_init_bounds);
    message.mutable_levenberg_marquardt_dense()->set_weight_eq_factor(_weight_adapt_factor_eq);
    message.mutable_levenberg_marquardt_dense()->set_weight_ineq_factor(_weight_adapt_factor_ineq);
    message.mutable_levenberg_marquardt_dense()->set_weight_bounds_factor(_weight_adapt_factor_bounds);
    message.mutable_levenberg_marquardt_dense()->set_weight_eq_max(_weight_adapt_max_eq);
    message.mutable_levenberg_marquardt_dense()->set_weight_ineq_max(_weight_adapt_max_ineq);
    message.mutable_levenberg_marquardt_dense()->set_weight_bounds_max(_weight_adapt_max_bounds);
}

void LevenbergMarquardtDense::fromMessage(const corbo::messages::NlpSolver& message, std::stringstream* issues)
{
    _iterations                 = message.levenberg_marquardt_dense().iterations();
    _weight_init_eq             = message.levenberg_marquardt_dense().weight_eq();
    _weight_init_ineq           = message.levenberg_marquardt_dense().weight_ineq();
    _weight_init_bounds         = message.levenberg_marquardt_dense().weight_bounds();
    _weight_adapt_factor_eq     = message.levenberg_marquardt_dense().weight_eq_factor();
    _weight_adapt_factor_ineq   = message.levenberg_marquardt_dense().weight_ineq_factor();
    _weight_adapt_factor_bounds = message.levenberg_marquardt_dense().weight_bounds_factor();
    _weight_adapt_max_eq        = message.levenberg_marquardt_dense().weight_eq_max();
    _weight_adapt_max_ineq      = message.levenberg_marquardt_dense().weight_ineq_max();
    _weight_adapt_max_bounds    = message.levenberg_marquardt_dense().weight_bounds_max();
}
#endif

}  // namespace corbo
