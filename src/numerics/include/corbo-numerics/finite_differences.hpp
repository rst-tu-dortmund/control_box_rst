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

template <typename IncFun, typename EvalFun>
void ForwardDifferences::jacobian(IncFun inc_fun, EvalFun eval_fun, Eigen::Ref<Eigen::MatrixXd> jacobian)
{
    constexpr const double delta  = 1e-9;         // discretization width
    constexpr const double scalar = 1.0 / delta;  // auxillary variable to avoid some repeated flops.

    int dim_x   = jacobian.cols();
    int dim_val = jacobian.rows();

    Eigen::VectorXd f0(dim_val);
    Eigen::VectorXd f1(dim_val);

    eval_fun(f0);  // TODO(roesmann) should we support caching?
    for (int i = 0; i < dim_x; ++i)
    {
        inc_fun(i, delta);
        eval_fun(f1);
        inc_fun(i, -delta);
        jacobian.col(i).noalias() = scalar * (f1 - f0);
    }
}

template <typename IncFun, typename EvalFun>
void ForwardDifferences::hessian(IncFun inc_fun, EvalFun eval_fun, int dim_f, Eigen::Ref<Eigen::MatrixXd> hessian, const double* multipliers)
{
    constexpr const double delta  = 1e-5;                 // discretization width
    constexpr const double scalar = 1 / (delta * delta);  // auxillary variable to avoid some repeated flops.

    assert(hessian.rows() == hessian.cols());

    int dim_x = hessian.cols();

    // http://www.okstate.edu/sas/v8/sashtml/ormp/chap5/sect28.htm
    // http://www.mathematik.uni-dortmund.de/~kuzmin/cfdintro/lecture4.pdf <- higher order approximations

    Eigen::VectorXd f0(dim_f);  // f(x, y)
    Eigen::VectorXd f1(dim_f);  // f(x+h, y)
    Eigen::VectorXd f2(dim_f);  // f(x, y+h)
    Eigen::VectorXd f3(dim_f);  // f(x+h, y+h)

    // df(x,y)/(dxdy) = 1/(h^2) ( f3 - f1 - f2 + f0 )

    // TODO(roesmann) symmetry?
    for (int i = 0; i < dim_x; ++i)
    {
        for (int j = 0; j < dim_x; ++j)
        {
            inc_fun(i, delta);
            eval_fun(f1);
            inc_fun(j, delta);
            eval_fun(f3);
            inc_fun(i, -delta);
            eval_fun(f2);
            inc_fun(j, -delta);
            eval_fun(f0);

            if (multipliers)
            {
                hessian(i, j) = scalar * (f3[0] - f1[0] - f2[0] + f0[0]) * multipliers[0];
                for (int v = 1; v < dim_f; ++v)
                {
                    hessian(i, j) += scalar * (f3[v] - f1[v] - f2[v] + f0[v]) * multipliers[v];
                }
            }
            else
            {
                hessian(i, j) = scalar * (f3[0] - f1[0] - f2[0] + f0[0]);
                for (int v = 1; v < dim_f; ++v)
                {
                    hessian(i, j) += scalar * (f3[v] - f1[v] - f2[v] + f0[v]);
                }
            }
        }
    }
}

template <typename IncFun, typename EvalFun>
void ForwardDifferences::jacobianHessian(IncFun inc_fun, EvalFun eval_fun, Eigen::Ref<Eigen::MatrixXd> jacobian, Eigen::Ref<Eigen::MatrixXd> hessian,
                                         const double* multipliers)
{
    constexpr const double delta     = 1e-5;                 // discretization width
    constexpr const double scalar_x  = 1 / (delta);          // auxillary variable to avoid some repeated flops.
    constexpr const double scalar_xy = 1 / (delta * delta);  // auxillary variable to avoid some repeated flops.

    assert(hessian.rows() == hessian.cols());
    assert(jacobian.cols() == hessian.cols());

    int dim_x = hessian.cols();
    int dim_f = jacobian.rows();

    Eigen::Map<const Eigen::VectorXd> multipliers_vec(multipliers, dim_f);

    // http://www.okstate.edu/sas/v8/sashtml/ormp/chap5/sect28.htm

    Eigen::VectorXd f0(dim_f);  // f(x, y)
    Eigen::VectorXd f1(dim_f);  // f(x+h, y)
    Eigen::VectorXd f2(dim_f);  // f(x, y+h)
    Eigen::VectorXd f3(dim_f);  // f(x+h, y+h)

    // df(x,y)/(dxdy) = 1/(h^2) ( f3 - f1 - f2 + f0 )

    // TODO(roesmann) symmetry?
    for (int i = 0; i < dim_x; ++i)
    {
        for (int j = 0; j < dim_x; ++j)
        {
            inc_fun(i, delta);
            eval_fun(f1);
            inc_fun(j, delta);
            eval_fun(f3);
            inc_fun(i, -delta);
            eval_fun(f2);
            inc_fun(j, -delta);
            eval_fun(f0);

            if (multipliers)
            {
                hessian(i, j) = scalar_xy * (f3[0] - f1[0] - f2[0] + f0[0]) * multipliers[0];
                for (int v = 1; v < dim_f; ++v)
                {
                    hessian(i, j) += scalar_xy * (f3[v] - f1[v] - f2[v] + f0[v]) * multipliers[v];
                }
                if (i == j) jacobian.col(i).noalias() = scalar_x * (f1 - f0).cwiseProduct(multipliers_vec);
            }
            else
            {
                hessian(i, j) = scalar_xy * (f3[0] - f1[0] - f2[0] + f0[0]);
                for (int v = 1; v < dim_f; ++v)
                {
                    hessian(i, j) += scalar_xy * (f3[v] - f1[v] - f2[v] + f0[v]);
                }
                if (i == j) jacobian.col(i).noalias() = scalar_x * (f1 - f0);
            }
        }
    }
}

template <typename IncFun, typename EvalFun>
void CentralDifferences::jacobian(IncFun inc_fun, EvalFun eval_fun, Eigen::Ref<Eigen::MatrixXd> jacobian)
{
    constexpr const double delta  = 1e-9;          // discretization width
    constexpr const double ddelta = 2 * delta;     // two times the discretization width
    constexpr const double scalar = 1.0 / ddelta;  // auxillary variable to avoid some repeated flops.

    int dim_x   = jacobian.cols();
    int dim_val = jacobian.rows();

    Eigen::VectorXd f0(dim_val);
    Eigen::VectorXd f1(dim_val);

    for (int i = 0; i < dim_x; ++i)
    {
        inc_fun(i, delta);
        eval_fun(f1);
        inc_fun(i, -ddelta);
        eval_fun(f0);
        jacobian.col(i).noalias() = scalar * (f1 - f0);
        inc_fun(i, delta);
    }
}

template <typename IncFun, typename EvalFun>
void CentralDifferences::hessian(IncFun inc_fun, EvalFun eval_fun, int dim_f, Eigen::Ref<Eigen::MatrixXd> hessian, const double* multipliers)
{
    constexpr const double delta     = 1e-5;                       // discretization width
    constexpr const double ddelta    = 2 * delta;                  // two times the discretization width
    constexpr const double scalar_xx = 1 / (delta * delta);        // auxillary variable to avoid some repeated flops.
    constexpr const double scalar_xy = 1 / (4.0 * delta * delta);  // auxillary variable to avoid some repeated flops.

    assert(hessian.rows() == hessian.cols());

    int dim_x = hessian.cols();

    Eigen::VectorXd f1(dim_f);  // f(x+h, y+h)
    Eigen::VectorXd f2(dim_f);  // f(x+h, y-h)
    Eigen::VectorXd f3(dim_f);  // f(x-h, y+h)
    Eigen::VectorXd f4(dim_f);  // f(x-h, y-h)

    // df(x,y)/(dxdx) = 1/(h^2) ( f(x+h,y) - 2f(x,y) + f(x-h,y)
    // df(x,y)/(dxdy) = 1/(4h^2) ( f1 - f2 - f3 + f4 )

    // TODO(roesmann) symmetry?
    for (int i = 0; i < dim_x; ++i)
    {
        for (int j = 0; j < dim_x; ++j)
        {
            if (i == j)
            {
                inc_fun(i, delta);
                eval_fun(f1);
                inc_fun(i, -ddelta);
                eval_fun(f3);
                inc_fun(i, delta);
                eval_fun(f2);

                if (multipliers)
                {
                    hessian(i, j) = scalar_xx * (f1[0] - 2 * f2[0] + f3[0]) * multipliers[0];
                    for (int v = 1; v < dim_f; ++v)
                    {
                        hessian(i, j) += scalar_xx * (f1[v] - 2 * f2[v] + f3[v]) * multipliers[v];
                    }
                }
                else
                {
                    hessian(i, j) = scalar_xx * (f1[0] - 2 * f2[0] + f3[0]);
                    for (int v = 1; v < dim_f; ++v)
                    {
                        hessian(i, j) += scalar_xx * (f1[v] - 2 * f2[v] + f3[v]);
                    }
                }
            }
            else
            {
                inc_fun(i, delta);
                inc_fun(j, delta);
                eval_fun(f1);
                inc_fun(j, -ddelta);
                eval_fun(f2);
                inc_fun(i, -ddelta);
                eval_fun(f4);
                inc_fun(j, ddelta);
                eval_fun(f3);

                // reset values
                inc_fun(i, delta);
                inc_fun(j, -delta);

                if (multipliers)
                {
                    hessian(i, j) = scalar_xy * (f1[0] - f2[0] - f3[0] + f4[0]) * multipliers[0];
                    for (int v = 1; v < dim_f; ++v)
                    {
                        hessian(i, j) += scalar_xy * (f1[v] - f2[v] - f3[v] + f4[v]) * multipliers[v];
                    }
                }
                else
                {
                    hessian(i, j) = scalar_xy * (f1[0] - f2[0] - f3[0] + f4[0]);
                    for (int v = 1; v < dim_f; ++v)
                    {
                        hessian(i, j) += scalar_xy * (f1[v] - f2[v] - f3[v] + f4[v]);
                    }
                }
            }
        }
    }
}

template <typename IncFun, typename EvalFun>
void CentralDifferences::jacobianHessian(IncFun inc_fun, EvalFun eval_fun, Eigen::Ref<Eigen::MatrixXd> jacobian, Eigen::Ref<Eigen::MatrixXd> hessian,
                                         const double* multipliers)
{
    constexpr const double delta     = 1e-5;                       // discretization width
    constexpr const double ddelta    = 2 * delta;                  // two times the discretization width
    constexpr const double scalar_x  = 1 / ddelta;                 // auxillary variable to avoid some repeated flops.
    constexpr const double scalar_xx = 1 / (delta * delta);        // auxillary variable to avoid some repeated flops.
    constexpr const double scalar_xy = 1 / (2.0 * delta * delta);  // auxillary variable to avoid some repeated flops.

    assert(hessian.rows() == hessian.cols());
    assert(jacobian.cols() == hessian.cols());

    int dim_x = hessian.cols();
    int dim_f = jacobian.rows();

    Eigen::Map<const Eigen::VectorXd> multipliers_vec(multipliers, dim_f);

    // https://en.wikipedia.org/wiki/Finite_difference

    Eigen::VectorXd f0(dim_f);  // f(x+h,y+h)
    Eigen::VectorXd f1(dim_f);  // f(x+h,y)
    Eigen::VectorXd f2(dim_f);  // f(x, y+h)
    Eigen::VectorXd f3(dim_f);  // f(x,y)
    Eigen::VectorXd f4(dim_f);  // f(x-h,y)
    Eigen::VectorXd f5(dim_f);  // f(x, y-h)
    Eigen::VectorXd f6(dim_f);  // f(x-h,y-h)

    // TODO(roesmann) symmetry?

    for (int i = 0; i < dim_x; ++i)
    {
        for (int j = 0; j < dim_x; ++j)
        {
            if (i == j)
            {
                inc_fun(i, delta);
                eval_fun(f1);
                inc_fun(i, -ddelta);
                eval_fun(f4);
                inc_fun(i, delta);
                eval_fun(f3);

                if (multipliers)
                {
                    hessian(i, j) = scalar_xx * (f1[0] - 2 * f3[0] + f4[0]) * multipliers[0];
                    for (int v = 1; v < dim_f; ++v)
                    {
                        hessian(i, j) += scalar_xx * (f1[v] - 2 * f3[v] + f4[v]) * multipliers[v];
                    }
                    jacobian.col(i).noalias() = scalar_x * (f1 - f4).cwiseProduct(multipliers_vec);  // only for i==j
                }
                else
                {
                    hessian(i, j) = scalar_xx * (f1[0] - 2 * f3[0] + f4[0]);
                    for (int v = 1; v < dim_f; ++v)
                    {
                        hessian(i, j) += scalar_xx * (f1[v] - 2 * f3[v] + f4[v]);
                    }
                    jacobian.col(i).noalias() = scalar_x * (f1 - f4);  // only for i==j
                }
            }
            else
            {
                inc_fun(i, delta);
                eval_fun(f1);
                inc_fun(j, delta);
                eval_fun(f0);
                inc_fun(i, -delta);
                eval_fun(f2);
                inc_fun(j, -ddelta);
                eval_fun(f5);
                inc_fun(i, -delta);
                eval_fun(f6);
                inc_fun(j, delta);
                eval_fun(f4);
                inc_fun(i, delta);
                eval_fun(f3);

                if (multipliers)
                {
                    hessian(i, j) = scalar_xy * (f0[0] - f1[0] - f2[0] + 2 * f3[0] - f4[0] - f5[0] + f6[0]) * multipliers[0];
                    for (int v = 1; v < dim_f; ++v)
                    {
                        hessian(i, j) += scalar_xy * (f0[v] - f1[v] - f2[v] + 2 * f3[v] - f4[v] - f5[v] + f6[v]) * multipliers[v];
                    }
                }
                else
                {
                    hessian(i, j) = scalar_xy * (f0[0] - f1[0] - f2[0] + 2 * f3[0] - f4[0] - f5[0] + f6[0]);
                    for (int v = 1; v < dim_f; ++v)
                    {
                        hessian(i, j) += scalar_xy * (f0[v] - f1[v] - f2[v] + 2 * f3[v] - f4[v] - f5[v] + f6[v]);
                    }
                }
            }
        }
    }
}

}  // namespace corbo
