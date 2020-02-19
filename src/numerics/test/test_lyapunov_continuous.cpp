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

#include <corbo-numerics/lyapunov_continuous.h>

#include <corbo-core/macros.h>

#include "gtest/gtest.h"

using corbo::LyapunovContinuous;

class TestLyapunovEquationContinuousTime : public testing::Test
{
 protected:
    // You can do set-up work for each test here.
    TestLyapunovEquationContinuousTime() {}
    // You can do clean-up work that doesn't throw exceptions here.
    virtual ~TestLyapunovEquationContinuousTime() {}
    // If the constructor and destructor are not enough for setting up
    // and cleaning up each test, you can define the following methods:

    // Code here will be called immediately after the constructor (right
    // before each test).
    // virtual void SetUp() {}
    // Code here will be called immediately after each test (right
    // before the destructor).
    // virtual void TearDown();
};

TEST_F(TestLyapunovEquationContinuousTime, solve_feasible1)
{
    Eigen::Matrix2d A;
    A << 1, 2, -3, -4;

    Eigen::Matrix2d Q;
    Q << 3, 1, 1, 1;

    EXPECT_TRUE(LyapunovContinuous::hasUniqueSolution(A));

    Eigen::MatrixXd X;
    bool solve_success = LyapunovContinuous::solve(A, Q, X);
    EXPECT_TRUE(solve_success);

    Eigen::Matrix2d X_sol;
    X_sol << 6.166666, -3.8333333, -3.8333333, 3;

    EXPECT_EQ_MATRIX(X, X_sol, 1e-5);
}

TEST_F(TestLyapunovEquationContinuousTime, solve_feasible2)
{
    Eigen::Matrix2d A;
    A << 0, 1, -6, -5;

    A.transposeInPlace();

    Eigen::Matrix2d Q = Eigen::Matrix2d::Identity();

    EXPECT_TRUE(LyapunovContinuous::hasUniqueSolution(A));

    Eigen::MatrixXd X;
    bool solve_success = LyapunovContinuous::solve(A, Q, X);
    EXPECT_TRUE(solve_success);

    Eigen::Matrix2d X_sol;
    X_sol << 67.0 / 60.0, 1.0 / 12.0, 1.0 / 12.0, 7.0 / 60.0;

    EXPECT_EQ_MATRIX(X, X_sol, 1e-5);
}

TEST_F(TestLyapunovEquationContinuousTime, solve_feasible3)
{
    Eigen::Matrix2d A;
    A << 0, 1, -6, -5;

    Eigen::Matrix2d Q = Eigen::Matrix2d::Identity();

    EXPECT_TRUE(LyapunovContinuous::hasUniqueSolution(A));

    Eigen::MatrixXd X;
    bool solve_success = LyapunovContinuous::solve(A, Q, X);
    EXPECT_TRUE(solve_success);

    Eigen::Matrix2d X_sol;
    X_sol << 0.533333, -0.5, -0.5, 0.7;

    EXPECT_EQ_MATRIX(X, X_sol, 1e-5);
}

TEST_F(TestLyapunovEquationContinuousTime, solve_feasible4)
{
    // G(s) =  1/(s+1)/(s+2)/(s+3)/(s+4)/(s+5)
    Eigen::MatrixXd A(5, 5);
    A << -15, -85, -225, -274, -120, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0;

    Eigen::MatrixXd Q(5, 5);
    // Q = magic(5)/10
    Q << 1.7, 2.4, 0.1, 0.8, 1.5, 2.3, 0.5, 0.7, 1.4, 1.6, 0.4, 0.6, 1.3, 2.0, 2.2, 1.0, 1.2, 1.9, 2.1, 0.3, 1.1, 1.8, 2.5, 0.2, 0.9;

    EXPECT_TRUE(LyapunovContinuous::hasUniqueSolution(A));

    Eigen::MatrixXd X;
    bool solve_success = LyapunovContinuous::solve(A, Q, X);
    EXPECT_TRUE(solve_success);

    Eigen::MatrixXd X_sol(5, 5);
    X_sol << 489.7681, 0.0322, -35.0314, -0.7613, 6.1874, -0.5322, 34.3314, -0.6387, -7.7874, -1.1722, -34.9314, -0.6613, 5.7874, -1.0278, -3.6708,
        -0.5387, -7.6874, -1.0722, 3.3708, -0.4565, 5.8874, -1.4278, -3.5708, -0.4435, 7.9881;

    EXPECT_EQ_MATRIX(X, X_sol, 1e-5);
}

TEST_F(TestLyapunovEquationContinuousTime, solve_feasible5)
{
    // this test differs from solve_feasible 1-4 such that
    // the underlying Schur solver in the current implementation must be
    // of complex type (rather than real) in order to find the correct solution.
    Eigen::MatrixXd A(5, 5);
    A << -2.5, -2.3, -0.95, -0.1689, -0.095, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0;

    Eigen::MatrixXd Q = Eigen::MatrixXd::Identity(5, 5);

    EXPECT_TRUE(LyapunovContinuous::hasUniqueSolution(A));

    Eigen::MatrixXd X;
    bool solve_success = LyapunovContinuous::solve(A, Q, X);
    EXPECT_TRUE(solve_success);

    Eigen::MatrixXd X_sol(5, 5);
    X_sol << 1.1495, -0.5000, -0.8201, 0.5000, -5.5701, -0.5000, 0.8201, -0.5000, 5.5701, 0.5000, -0.8201, -0.5000, -5.5701, -0.5000, 85.0134, 0.5000,
        5.5701, -0.5000, -85.0134, -0.5000, -5.5701, 0.5000, 85.0134, -0.5000, -709.5066;

    EXPECT_EQ_MATRIX(X, X_sol, 1e-5);
}

TEST_F(TestLyapunovEquationContinuousTime, try_solve_A_nonsquare)
{
    Eigen::MatrixXd A = Eigen::MatrixXd::Random(5, 4);

    Eigen::Matrix2d Q;
    Q << 3, 1, 1, 1;

    Eigen::MatrixXd X;

    // EXPECT_DEBUG_DEATH(LyapunovEquationContinuousTime::solve(A, Q, X), "[\\s\\S]*[is_square|have_equal_size][\\s\\S]*");
    EXPECT_DEBUG_DEATH(LyapunovContinuous::solve(A, Q, X), "");
}

TEST_F(TestLyapunovEquationContinuousTime, try_solve_Q_nonsquare)
{
    Eigen::MatrixXd A = Eigen::MatrixXd::Random(4, 4);
    Eigen::MatrixXd Q = Eigen::MatrixXd::Random(5, 4);

    Eigen::MatrixXd X;
    EXPECT_DEBUG_DEATH(LyapunovContinuous::solve(A, Q, X), "");
}

TEST_F(TestLyapunovEquationContinuousTime, try_solve_A_Q_dim_mismatch)
{
    Eigen::MatrixXd A = Eigen::MatrixXd::Random(4, 4);
    Eigen::MatrixXd Q = Eigen::MatrixXd::Random(5, 5);

    Eigen::MatrixXd X;
    EXPECT_DEBUG_DEATH(LyapunovContinuous::solve(A, Q, X), "");
}

TEST_F(TestLyapunovEquationContinuousTime, has_unique_solution_non_square)
{
    Eigen::MatrixXd A = Eigen::MatrixXd::Random(5, 4);

    Eigen::MatrixXd X;
    EXPECT_FALSE(LyapunovContinuous::hasUniqueSolution(A));
}

TEST_F(TestLyapunovEquationContinuousTime, has_unique_solution_infeasible1)
{
    Eigen::Matrix3d A = Eigen::Matrix3d::Ones();

    EXPECT_FALSE(LyapunovContinuous::hasUniqueSolution(A));
}
