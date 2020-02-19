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

#include <corbo-optimization/simple_optimization_problem.h>

namespace corbo {

void SimpleOptimizationProblem::resizeParameterVector(int parameter_dim)
{
    _x.setZero(parameter_dim);
    _lb.setConstant(parameter_dim, -CORBO_INF_DBL);
    _ub.setConstant(parameter_dim, CORBO_INF_DBL);
    _x_backup.clear();
}

void SimpleOptimizationProblem::restoreBackupParameters(bool keep_backup)
{
    assert(!_x_backup.empty());
    _x = _x_backup.back();
    if (!keep_backup) _x_backup.pop_back();
}

void SimpleOptimizationProblem::getBounds(Eigen::Ref<Eigen::VectorXd> lb, Eigen::Ref<Eigen::VectorXd> ub)
{
    lb = _lb;
    ub = _ub;
}

void SimpleOptimizationProblem::setBounds(const Eigen::Ref<const Eigen::VectorXd>& lb, const Eigen::Ref<const Eigen::VectorXd>& ub)
{
    _lb = lb;
    _ub = ub;
}

void SimpleOptimizationProblemWithCallbacks::setObjectiveFunction(std::function<void(const Eigen::VectorXd&, Eigen::Ref<Eigen::VectorXd>)> obj_fun,
                                                                  int obj_dim, bool lsq_form)
{
    _obj_fun  = obj_fun;
    _obj_dim  = obj_dim;
    _lsq_form = lsq_form;
}

void SimpleOptimizationProblemWithCallbacks::setEqualityConstraint(std::function<void(const Eigen::VectorXd&, Eigen::Ref<Eigen::VectorXd>)> eq_fun,
                                                                   int eq_dim)
{
    _eq_fun = eq_fun;
    _eq_dim = eq_dim;
}

void SimpleOptimizationProblemWithCallbacks::setInequalityConstraint(
    std::function<void(const Eigen::VectorXd&, Eigen::Ref<Eigen::VectorXd>)> ineq_fun, int ineq_dim)
{
    _ineq_fun = ineq_fun;
    _ineq_dim = ineq_dim;
}

}  // namespace corbo
