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

#include <corbo-numerics/algebraic_riccati_discrete.h>

#include <corbo-core/macros.h>

#include "gtest/gtest.h"

using corbo::AlgebraicRiccatiDiscrete;

class TestAlgebraicRiccatiDiscrete : public testing::Test
{
 protected:
    // You can do set-up work for each test here.
    TestAlgebraicRiccatiDiscrete() {}
    // You can do clean-up work that doesn't throw exceptions here.
    virtual ~TestAlgebraicRiccatiDiscrete() {}
    // If the constructor and destructor are not enough for setting up
    // and cleaning up each test, you can define the following methods:

    // Code here will be called immediately after the constructor (right
    // before each test).
    // virtual void SetUp() {}
    // Code here will be called immediately after each test (right
    // before the destructor).
    // virtual void TearDown();
};

TEST_F(TestAlgebraicRiccatiDiscrete, solve_feasible_2d_stable_real)
{
    Eigen::Matrix2d A;
    A << -0.05, 0.05, 1, 0;  // only real eigenvalues (all stable)

    Eigen::Vector2d B;
    B << 1, 0;

    Eigen::Matrix<double, 1, 2> C;
    C << 0, 1;

    Eigen::MatrixXd Q = C.transpose() * C;
    Eigen::MatrixXd R(1, 1);
    R << 1;

    Eigen::MatrixXd X;
    Eigen::MatrixXd G;
    bool solve_success = AlgebraicRiccatiDiscrete::solve(A, B, Q, R, X, &G);
    EXPECT_TRUE(solve_success);

    Eigen::Matrix2d X_sol;
    X_sol << 1.0026, -0.0013, -0.0013, 1.0013;
    EXPECT_EQ_MATRIX(X, X_sol, 1e-4);

    // check gain matrix
    Eigen::MatrixXd G_sol(1, 2);
    G_sol << -0.025673043427904, 0.025032039904055;
    EXPECT_EQ_MATRIX(G, G_sol, 1e-4);

    // check closed-loop stability
    EXPECT_TRUE(AlgebraicRiccatiDiscrete::isClosedLoopStable(A, B, G));
}

TEST_F(TestAlgebraicRiccatiDiscrete, solve_feasible_4d_stableA_complex)
{
    Eigen::Matrix4d A;
    A << -0.05, -0.2, -0.0125, 0.0125, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0;  // a complex pair of eigenvalues (all eigenvalues stable)

    Eigen::Vector4d B;
    B << 1, 0, 0, 0;

    Eigen::Matrix<double, 1, 4> C;
    C << 0, 0, 0, 1;

    Eigen::MatrixXd Q = C.transpose() * C;

    Eigen::MatrixXd R(1, 1);
    R << 2;

    Eigen::MatrixXd X;
    Eigen::MatrixXd G;
    bool solve_success = AlgebraicRiccatiDiscrete::solve(A, B, Q, R, X, &G);
    EXPECT_TRUE(solve_success);

    Eigen::Matrix4d X_sol;
    X_sol << 1.0289, 0.0075, -0.0013, -0.0004, 0.0075, 1.0277, 0.0016, -0.0017, -0.0013, 0.0016, 1.0002, -0.0001, -0.0004, -0.0017, -0.0001, 1.0001;
    EXPECT_EQ_MATRIX(X, X_sol, 1e-4);

    // check gain matrix
    Eigen::MatrixXd G_sol(1, 4);
    G_sol << -0.014509214644181, -0.068380590003095, -0.004366017750237, 0.004246262501251;
    EXPECT_EQ_MATRIX(G, G_sol, 1e-4);

    // check closed-loop stability
    EXPECT_TRUE(AlgebraicRiccatiDiscrete::isClosedLoopStable(A, B, G));
}

TEST_F(TestAlgebraicRiccatiDiscrete, solve_feasible_4d_unstableA_complex)
{
    Eigen::Matrix4d A;
    A << -1.5, -1.32, 3.742, 2.076, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0;  // a complex pair of eigenvalues (some eigenvalues unstable)

    Eigen::Vector4d B;
    B << 1, 0, 0, 0;

    Eigen::Matrix<double, 1, 4> C;
    C << 0, 0, 0, 1;

    Eigen::MatrixXd Q = C.transpose() * C;

    Eigen::MatrixXd R(1, 1);
    R << 2;

    Eigen::MatrixXd X;
    Eigen::MatrixXd G;
    bool solve_success = AlgebraicRiccatiDiscrete::solve(A, B, Q, R, X, &G);
    EXPECT_TRUE(solve_success);

    Eigen::Matrix4d X_sol;
    X_sol << 33.8189, 7.5236, -14.3996, -5.0081, 7.5236, 32.4107, 0.1086, -6.8438, -14.3996, 0.1086, 32.7867, 14.0887, -5.0081, -6.8438, 14.0887,
        9.1383;
    EXPECT_EQ_MATRIX(X, X_sol, 1e-4);

    // check gain matrix
    Eigen::MatrixXd G_sol(1, 4);
    G_sol << -1.206200888594686, -1.648306965169984, 3.393241337036358, 1.960083490678313;
    EXPECT_EQ_MATRIX(G, G_sol, 1e-4);

    // check closed-loop stability
    EXPECT_TRUE(AlgebraicRiccatiDiscrete::isClosedLoopStable(A, B, G));
}
