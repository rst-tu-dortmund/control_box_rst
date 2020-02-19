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

#include <corbo-numerics/finite_differences.h>

namespace corbo {

void ForwardDifferences::computeJacobian(std::function<void(int, const double&)> inc_fun, std::function<void(Eigen::Ref<Eigen::VectorXd>)> eval_fun,
                                         Eigen::Ref<Eigen::MatrixXd> jacobian)
{
    ForwardDifferences::jacobian(inc_fun, eval_fun, jacobian);
}

void ForwardDifferences::computeJacobian2(std::function<void(int, const double&)> inc_fun, std::function<void(Eigen::VectorXd&)> eval_fun,
                                          Eigen::Ref<Eigen::MatrixXd> jacobian)
{
    ForwardDifferences::jacobian(inc_fun, eval_fun, jacobian);
}

void ForwardDifferences::computeHessian(std::function<void(int, const double&)> inc_fun, std::function<void(Eigen::Ref<Eigen::VectorXd>)> eval_fun,
                                        int dim_f, Eigen::Ref<Eigen::MatrixXd> hessian, const double* multipliers)
{
    ForwardDifferences::hessian(inc_fun, eval_fun, dim_f, hessian, multipliers);
}

void ForwardDifferences::computeHessian2(std::function<void(int, const double&)> inc_fun, std::function<void(Eigen::VectorXd&)> eval_fun, int dim_f,
                                         Eigen::Ref<Eigen::MatrixXd> hessian, const double* multipliers)
{
    ForwardDifferences::hessian(inc_fun, eval_fun, dim_f, hessian, multipliers);
}

void ForwardDifferences::computeJacobianAndHessian(std::function<void(int, const double&)> inc_fun,
                                                   std::function<void(Eigen::Ref<Eigen::VectorXd>)> eval_fun, Eigen::Ref<Eigen::MatrixXd> jacobian,
                                                   Eigen::Ref<Eigen::MatrixXd> hessian, const double* multipliers)
{
    ForwardDifferences::jacobianHessian(inc_fun, eval_fun, jacobian, hessian, multipliers);
}

void ForwardDifferences::computeJacobianAndHessian2(std::function<void(int, const double&)> inc_fun, std::function<void(Eigen::VectorXd&)> eval_fun,
                                                    Eigen::Ref<Eigen::MatrixXd> jacobian, Eigen::Ref<Eigen::MatrixXd> hessian,
                                                    const double* multipliers)
{
    ForwardDifferences::jacobianHessian(inc_fun, eval_fun, jacobian, hessian, multipliers);
}

void CentralDifferences::computeJacobian(std::function<void(int, const double&)> inc_fun, std::function<void(Eigen::Ref<Eigen::VectorXd>)> eval_fun,
                                         Eigen::Ref<Eigen::MatrixXd> jacobian)
{
    CentralDifferences::jacobian(inc_fun, eval_fun, jacobian);
}

void CentralDifferences::computeJacobian2(std::function<void(int, const double&)> inc_fun, std::function<void(Eigen::VectorXd&)> eval_fun,
                                          Eigen::Ref<Eigen::MatrixXd> jacobian)
{
    CentralDifferences::jacobian(inc_fun, eval_fun, jacobian);
}

void CentralDifferences::computeHessian(std::function<void(int, const double&)> inc_fun, std::function<void(Eigen::Ref<Eigen::VectorXd>)> eval_fun,
                                        int dim_f, Eigen::Ref<Eigen::MatrixXd> hessian, const double* multipliers)
{
    CentralDifferences::hessian(inc_fun, eval_fun, dim_f, hessian, multipliers);
}

void CentralDifferences::computeHessian2(std::function<void(int, const double&)> inc_fun, std::function<void(Eigen::VectorXd&)> eval_fun, int dim_f,
                                         Eigen::Ref<Eigen::MatrixXd> hessian, const double* multipliers)
{
    CentralDifferences::hessian(inc_fun, eval_fun, dim_f, hessian, multipliers);
}

void CentralDifferences::computeJacobianAndHessian(std::function<void(int, const double&)> inc_fun,
                                                   std::function<void(Eigen::Ref<Eigen::VectorXd>)> eval_fun, Eigen::Ref<Eigen::MatrixXd> jacobian,
                                                   Eigen::Ref<Eigen::MatrixXd> hessian, const double* multipliers)
{
    CentralDifferences::jacobianHessian(inc_fun, eval_fun, jacobian, hessian, multipliers);
}

void CentralDifferences::computeJacobianAndHessian2(std::function<void(int, const double&)> inc_fun, std::function<void(Eigen::VectorXd&)> eval_fun,
                                                    Eigen::Ref<Eigen::MatrixXd> jacobian, Eigen::Ref<Eigen::MatrixXd> hessian,
                                                    const double* multipliers)
{
    CentralDifferences::jacobianHessian(inc_fun, eval_fun, jacobian, hessian, multipliers);
}

}  // namespace corbo
