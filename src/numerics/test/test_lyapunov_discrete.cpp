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

#include <corbo-numerics/lyapunov_discrete.h>

#include <corbo-core/macros.h>

#include "gtest/gtest.h"

using corbo::LyapunovDiscrete;

class TestLyapunovEquationDiscreteTime : public testing::Test
{
 protected:
    // You can do set-up work for each test here.
    TestLyapunovEquationDiscreteTime() {}
    // You can do clean-up work that doesn't throw exceptions here.
    virtual ~TestLyapunovEquationDiscreteTime() {}
    // If the constructor and destructor are not enough for setting up
    // and cleaning up each test, you can define the following methods:

    // Code here will be called immediately after the constructor (right
    // before each test).
    // virtual void SetUp() {}
    // Code here will be called immediately after each test (right
    // before the destructor).
    // virtual void TearDown();
};

TEST_F(TestLyapunovEquationDiscreteTime, solve_feasible1)
{
    Eigen::Matrix2d A;
    A << 0.2, 1, 0, -0.25;

    Eigen::Matrix2d Q;
    Q << 2, 0, 0, 2;

    EXPECT_TRUE(LyapunovDiscrete::hasUniqueSolution(A));

    Eigen::MatrixXd X;
    bool solve_success = LyapunovDiscrete::solve(A, Q, X);
    EXPECT_TRUE(solve_success);

    Eigen::Matrix2d X_sol;
    X_sol << 4.0939153, -0.5079365, -0.5079365, 2.133333;

    EXPECT_EQ_MATRIX(X, X_sol, 1e-5);
}

TEST_F(TestLyapunovEquationDiscreteTime, solve_feasible2)
{
    // G(z) = 1/(z+0.9)/(z+0.7)/(z+0.5)/(z+0.3)/(z+0.1)
    Eigen::MatrixXd A(5, 5);
    A << -2.5, -2.3, -0.95, -0.1689, -0.0095, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0;

    Eigen::MatrixXd Q(5, 5);
    // Q = magic(5)/100
    Q << 0.17, 0.24, 0.01, 0.08, 0.15, 0.23, 0.05, 0.07, 0.14, 0.16, 0.04, 0.06, 0.13, 0.20, 0.22, 0.10, 0.12, 0.19, 0.21, 0.03, 0.11, 0.18, 0.25,
        0.02, 0.09;

    EXPECT_TRUE(LyapunovDiscrete::hasUniqueSolution(A));

    Eigen::MatrixXd X;
    bool solve_success = LyapunovDiscrete::solve(A, Q, X);
    EXPECT_TRUE(solve_success);

    Eigen::MatrixXd X_sol(5, 5);
    X_sol << 617.5979, -611.0219, 592.5646, -565.5248, 532.3784, -610.8914, 617.6479, -610.9519, 592.7046, -565.3648, 592.3457, -610.8314, 617.7779,
        -610.7519, 592.9246, -565.2179, 592.4657, -610.6414, 617.9879, -610.7219, 532.0166, -565.0379, 592.7157, -610.6214, 618.0779;

    EXPECT_EQ_MATRIX(X, X_sol, 1e-5);
}

TEST_F(TestLyapunovEquationDiscreteTime, solve_feasible3)
{
    // this test differs from solve_feasible 1-2 such that
    // the underlying Schur solver in the current implementation must be
    // of complex type (rather than real) in order to find the correct solution.
    Eigen::MatrixXd A(5, 5);
    A << -2.5, -2.3, -0.95, -0.1689, -0.095, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0;

    Eigen::MatrixXd Q = Eigen::MatrixXd::Identity(5, 5);

    EXPECT_TRUE(LyapunovDiscrete::hasUniqueSolution(A));

    Eigen::MatrixXd X;
    bool solve_success = LyapunovDiscrete::solve(A, Q, X);
    EXPECT_TRUE(solve_success);

    Eigen::MatrixXd X_sol(5, 5);
    X_sol << -0.9387, 11.8081, -41.2126, 76.7855, -109.8624, 11.8081, 0.0613, 11.8081, -41.2126, 76.7855, -41.2126, 11.8081, 1.0613, 11.8081,
        -41.2126, 76.7855, -41.2126, 11.8081, 2.0613, 11.8081, -109.8624, 76.7855, -41.2126, 11.8081, 3.0613;

    EXPECT_EQ_MATRIX(X, X_sol, 1e-5);
}

TEST_F(TestLyapunovEquationDiscreteTime, try_solve_A_nonsquare)
{
    Eigen::MatrixXd A = Eigen::MatrixXd::Random(5, 4);

    Eigen::Matrix2d Q;
    Q << 3, 1, 1, 1;

    Eigen::MatrixXd X;
    EXPECT_DEBUG_DEATH(LyapunovDiscrete::solve(A, Q, X), "");
}

TEST_F(TestLyapunovEquationDiscreteTime, try_solve_Q_nonsquare)
{
    Eigen::MatrixXd A = Eigen::MatrixXd::Random(4, 4);
    Eigen::MatrixXd Q = Eigen::MatrixXd::Random(5, 4);

    Eigen::MatrixXd X;
    EXPECT_DEBUG_DEATH(LyapunovDiscrete::solve(A, Q, X), "");
}

TEST_F(TestLyapunovEquationDiscreteTime, try_solve_A_Q_dim_mismatch)
{
    Eigen::MatrixXd A = Eigen::MatrixXd::Random(4, 4);
    Eigen::MatrixXd Q = Eigen::MatrixXd::Random(5, 5);

    Eigen::MatrixXd X;
    EXPECT_DEBUG_DEATH(LyapunovDiscrete::solve(A, Q, X), "");
}

TEST_F(TestLyapunovEquationDiscreteTime, has_unique_solution_non_square)
{
    Eigen::MatrixXd A = Eigen::MatrixXd::Random(5, 4);

    Eigen::MatrixXd X;
    EXPECT_FALSE(LyapunovDiscrete::hasUniqueSolution(A));
}

TEST_F(TestLyapunovEquationDiscreteTime, has_unique_solution_infeasible1)
{
    Eigen::Matrix3d A = Eigen::Matrix3d::Identity();

    EXPECT_FALSE(LyapunovDiscrete::hasUniqueSolution(A));
}
