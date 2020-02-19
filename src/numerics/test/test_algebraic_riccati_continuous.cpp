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

#include <corbo-core/console.h>
#include <corbo-core/macros.h>
#include <corbo-numerics/algebraic_riccati_continuous.h>

#include "gtest/gtest.h"

using corbo::AlgebraicRiccatiContinuous;

class TestAlgebraicRiccatiContinuous : public testing::Test
{
 protected:
    // You can do set-up work for each test here.
    TestAlgebraicRiccatiContinuous() {}
    // You can do clean-up work that doesn't throw exceptions here.
    virtual ~TestAlgebraicRiccatiContinuous() {}
    // If the constructor and destructor are not enough for setting up
    // and cleaning up each test, you can define the following methods:

    // Code here will be called immediately after the constructor (right
    // before each test).
    // virtual void SetUp() {}
    // Code here will be called immediately after each test (right
    // before the destructor).
    // virtual void TearDown();
};

TEST_F(TestAlgebraicRiccatiContinuous, solve_feasible_2d_stable_real)
{
    Eigen::Matrix2d A;
    A << -3, 2, 1, 1;  // only real eigenvalues (all stable)

    Eigen::Vector2d B;
    B << 0, 1;

    Eigen::Matrix<double, 1, 2> C;
    C << 1, -1;

    Eigen::MatrixXd Q = C.transpose() * C;
    Eigen::MatrixXd R = Eigen::MatrixXd::Identity(1, 1);

    EXPECT_TRUE(AlgebraicRiccatiContinuous::isNumericallyStable(A, B, Q, R));

    Eigen::MatrixXd X;
    Eigen::MatrixXd G;
    bool solve_success = AlgebraicRiccatiContinuous::solve(A, B, Q, R, X, &G);
    EXPECT_TRUE(solve_success);

    Eigen::Matrix2d X_sol;
    X_sol << 0.2949, 0.5199, 0.5199, 3.0198;
    EXPECT_EQ_MATRIX(X, X_sol, 1e-4);

    // check gain matrix
    Eigen::MatrixXd G_sol(1, 2);
    G_sol << 0.519862500040767, 3.019764837837087;
    EXPECT_EQ_MATRIX(G, G_sol, 1e-4);

    // check closed-loop stability
    EXPECT_TRUE(AlgebraicRiccatiContinuous::isClosedLoopStable(A, B, G));
}

TEST_F(TestAlgebraicRiccatiContinuous, solve_feasible_4d_stableA_complex)
{
    Eigen::Matrix4d A;
    A << -8, -19, -22, -10, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0;  // a complex pair of eigenvalues (all eigenvalues stable)

    Eigen::Vector4d B;
    B << 1, 0, 0, 0;

    Eigen::Matrix<double, 1, 4> C;
    C << 0, 0, 0, 1;

    Eigen::MatrixXd Q = C.transpose() * C;

    Eigen::MatrixXd R(1, 1);
    R << 2;

    EXPECT_TRUE(AlgebraicRiccatiContinuous::isNumericallyStable(A, B, Q, R));

    Eigen::MatrixXd X;
    Eigen::MatrixXd G;
    bool solve_success = AlgebraicRiccatiContinuous::solve(A, B, Q, R, X, &G);
    EXPECT_TRUE(solve_success);

    Eigen::Matrix4d X_sol;
    X_sol << 0.0029, 0.0234, 0.0538, 0.0499, 0.0234, 0.1889, 0.4445, 0.4288, 0.0538, 0.4445, 1.1076, 1.1833, 0.0499, 0.4288, 1.1833, 1.6375;
    EXPECT_EQ_MATRIX(X, X_sol, 1e-4);

    // check gain matrix
    Eigen::MatrixXd G_sol(1, 4);
    G_sol << 0.001461438418027, 0.011692575245343, 0.026875649174793, 0.024968827881711;
    EXPECT_EQ_MATRIX(G, G_sol, 1e-4);

    // check closed-loop stability
    EXPECT_TRUE(AlgebraicRiccatiContinuous::isClosedLoopStable(A, B, G));
}

TEST_F(TestAlgebraicRiccatiContinuous, solve_feasible_4d_unstableA_complex)
{
    Eigen::Matrix4d A;
    A << 6, -5, -2, 10, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0;  // a complex pair of eigenvalues (some eigenvalues unstable)

    Eigen::Vector4d B;
    B << 1, 0, 0, 0;

    Eigen::Matrix<double, 1, 4> C;
    C << 0, 0, 0, 1;

    Eigen::MatrixXd Q = C.transpose() * C;

    Eigen::MatrixXd R(1, 1);
    R << 2;

    EXPECT_TRUE(AlgebraicRiccatiContinuous::isNumericallyStable(A, B, Q, R));

    Eigen::MatrixXd X;
    Eigen::MatrixXd G;
    bool solve_success = AlgebraicRiccatiContinuous::solve(A, B, Q, R, X, &G);
    EXPECT_TRUE(solve_success);

    Eigen::Matrix4d X_sol;
    X_sol << 28.0029, 28.0234, 40.0538, 40.0499, 28.0234, 324.1889, 336.4445, 40.4288, 40.0538, 336.4445, 777.1076, 481.1833, 40.0499, 40.4288,
        481.1833, 481.6375;
    EXPECT_EQ_MATRIX(X, X_sol, 1e-4);

    // check gain matrix
    Eigen::MatrixXd G_sol(1, 4);
    G_sol << 14.001461438418874, 14.011692575248027, 20.026875649177440, 20.024968827882546;
    EXPECT_EQ_MATRIX(G, G_sol, 1e-4);

    // check closed-loop stability
    EXPECT_TRUE(AlgebraicRiccatiContinuous::isClosedLoopStable(A, B, G));
}

TEST_F(TestAlgebraicRiccatiContinuous, solve_linearized_pendulum_model)
{
    Eigen::Matrix4d A;
    A << 0, 1.0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.0, 0, 0, 21.9524, -0.1348;

    Eigen::Vector4d B;
    B << 0, 1.0, 0, 2.2385;

    Eigen::Matrix<double, 4, 4> C;
    C << 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1;

    Eigen::Matrix<double, 4, 4> Q;
    Q << 35, 0, 0, 0, 0, 1, 0, 0, 0, 0, 5, 0, 0, 0, 0, 35;

    Eigen::MatrixXd R(1, 1);
    R << 0.1;

    Eigen::MatrixXd X;
    Eigen::MatrixXd G;
    bool solve_success = AlgebraicRiccatiContinuous::solve(A, B, Q, R, X, &G);
    EXPECT_TRUE(solve_success);

    Eigen::Matrix4d X_sol;
    X_sol << 37.1092, 19.1727, -57.065, -9.4007, 19.1727, 14.5091, -50.7528, -7.3677, -57.065, -50.7528, 231.2142, 27.6002, -9.4007, -7.3677, 27.6002,
        4.6237;

    EXPECT_EQ_MATRIX(X, X_sol, 1e-4);

    // check gain matrix
    Eigen::MatrixXd G_sol(1, 4);
    G_sol << -18.7083, -19.8357, 110.3077, 29.8251;
    EXPECT_EQ_MATRIX(G, G_sol, 1e-4);

    // check closed-loop stability
    EXPECT_TRUE(AlgebraicRiccatiContinuous::isClosedLoopStable(A, B, G));
}

TEST_F(TestAlgebraicRiccatiContinuous, numerically_unstable)
{
    Eigen::Matrix4d A;
    A << 0, 0, 0, 0, 1.0, 0, 0, 0, 0, 0, 0, 21.0524, 0, 0, 1, -0.1348;

    Eigen::Vector4d B;
    B << 0, 1.0, 0, 2.2385;

    Eigen::Matrix<double, 4, 4> Q;
    Q << 35, 0, 0, 0, 0, 1, 0, 0, 0, 0, 5, 0, 0, 0, 0, 35;

    Eigen::MatrixXd R(1, 1);
    R << 0.1;

    EXPECT_FALSE(AlgebraicRiccatiContinuous::isNumericallyStable(A, B, Q, R));
}
