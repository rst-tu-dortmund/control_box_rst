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

#include <corbo-numerics/sylvester_continuous.h>

#include <corbo-core/macros.h>

#include "gtest/gtest.h"

using corbo::SylvesterContinuous;

class TestSylvesterEquationContinuousTime : public testing::Test
{
 protected:
    // You can do set-up work for each test here.
    TestSylvesterEquationContinuousTime() {}
    // You can do clean-up work that doesn't throw exceptions here.
    virtual ~TestSylvesterEquationContinuousTime() {}
    // If the constructor and destructor are not enough for setting up
    // and cleaning up each test, you can define the following methods:

    // Code here will be called immediately after the constructor (right
    // before each test).
    // virtual void SetUp() {}
    // Code here will be called immediately after each test (right
    // before the destructor).
    // virtual void TearDown();
};

TEST_F(TestSylvesterEquationContinuousTime, solve_feasible1)
{
    Eigen::Matrix<double, 1, 1> A;
    A << 5;
    Eigen::Matrix2d B;
    B << 4, 3, 4, 3;

    Eigen::Matrix<double, 1, 2> C;
    C << 2, 1;

    EXPECT_TRUE(SylvesterContinuous::hasUniqueSolution(A, B));

    Eigen::MatrixXd X;
    bool solve_success = SylvesterContinuous::solve(A, B, C, X);
    EXPECT_TRUE(solve_success);

    Eigen::Matrix<double, 1, 2> X_sol;
    X_sol << -0.2, -0.05;

    EXPECT_EQ_MATRIX(X, X_sol, 1e-5);
}

TEST_F(TestSylvesterEquationContinuousTime, solve_feasible2)
{
    Eigen::Matrix2d A;
    A << 1, 2, -3, -4;

    Eigen::Matrix2d B = A.transpose();

    Eigen::Matrix2d C;
    C << 3, 1, 1, 1;

    EXPECT_TRUE(SylvesterContinuous::hasUniqueSolution(A, B));

    Eigen::MatrixXd X;
    bool solve_success = SylvesterContinuous::solve(A, B, C, X);
    EXPECT_TRUE(solve_success);

    Eigen::Matrix2d X_sol;
    X_sol << 6.166666, -3.8333333, -3.8333333, 3;

    EXPECT_EQ_MATRIX(X, X_sol, 1e-5);
}

TEST_F(TestSylvesterEquationContinuousTime, solve_feasible3)
{
    // this test differs from solve_feasible 1-2 such that
    // the underlying Schur solver in the current implementation must be
    // of complex type (rather than real) in order to find the correct solution.
    Eigen::MatrixXd A(5, 5);
    A << -2.5, -2.3, -0.95, -0.1689, -0.095, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0;

    Eigen::MatrixXd B(5, 5);
    B << -2, -1, -0.5, -0.1, -0.005, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0;

    Eigen::MatrixXd C = Eigen::MatrixXd::Identity(5, 5);

    EXPECT_TRUE(SylvesterContinuous::hasUniqueSolution(A, B));

    Eigen::MatrixXd X;
    bool solve_success = SylvesterContinuous::solve(A, B, C, X);
    EXPECT_TRUE(solve_success);

    Eigen::MatrixXd X_sol(5, 5);
    X_sol << -0.2099, -0.8702, -0.0400, 0.0221, 0.0016, 0.3106, 0.8311, 0.1808, 0.1953, 0.0089, 1.7890, 3.2674, 0.9579, -0.2863, -0.0164, -3.2772,
        -8.3434, -6.5446, -2.5965, -1.0414, -8.2753, -13.2734, 0.0681, 2.4069, 1.7689;

    EXPECT_EQ_MATRIX(X, X_sol, 1e-5);
}

TEST_F(TestSylvesterEquationContinuousTime, try_solve_dim_mismatch)
{
    Eigen::MatrixXd A = Eigen::MatrixXd::Random(4, 4);
    Eigen::MatrixXd B = Eigen::MatrixXd::Random(5, 5);
    Eigen::MatrixXd C = Eigen::MatrixXd::Random(5, 4);  // dimensions switched here

    Eigen::MatrixXd X;
    EXPECT_DEBUG_DEATH(SylvesterContinuous::solve(A, B, C, X), "");
}

TEST_F(TestSylvesterEquationContinuousTime, has_unique_solution_false1)
{
    Eigen::Matrix3d A = Eigen::Matrix3d::Ones();
    Eigen::Matrix3d B = A.transpose();

    EXPECT_FALSE(SylvesterContinuous::hasUniqueSolution(A, B));
}

TEST_F(TestSylvesterEquationContinuousTime, has_unique_solution_false2)
{
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(10, 10);
    Eigen::MatrixXd B = Eigen::MatrixXd::Zero(10, 10);

    EXPECT_FALSE(SylvesterContinuous::hasUniqueSolution(A, B));
}

TEST_F(TestSylvesterEquationContinuousTime, has_unique_solution_non_square)
{
    Eigen::MatrixXd A = Eigen::MatrixXd::Random(5, 4);
    Eigen::MatrixXd B = Eigen::MatrixXd::Random(5, 5);

    Eigen::MatrixXd X;
    EXPECT_FALSE(SylvesterContinuous::hasUniqueSolution(A, B));
    EXPECT_FALSE(SylvesterContinuous::hasUniqueSolution(B, A));
}
