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

#include <corbo-core/console.h>
#include <corbo-core/macros.h>
#include <corbo-core/value_comparison.h>

#include <functional>

#include "gtest/gtest.h"

using corbo::FiniteDifferencesInterface;
using corbo::ForwardDifferences;
using corbo::CentralDifferences;

class TestFiniteDifferences : public testing::Test
{
 protected:
    // You can do set-up work for each test here.
    TestFiniteDifferences() {}
    // You can do clean-up work that doesn't throw exceptions here.
    virtual ~TestFiniteDifferences() {}
    // If the constructor and destructor are not enough for setting up
    // and cleaning up each test, you can define the following methods:

    // Code here will be called immediately after the constructor (right
    // before each test).
    virtual void SetUp()
    {
        // problem 1
        // x + 2y + 3z
        // 4y
        // 5y + 6z
        x1 << 1, 2, 3;
        jacobian1_sol << 1, 2, 3, 0, 4, 0, 0, 5, 6;

        // problem 2
        // x^2 + y^3 + x^2*y
        x2 << 2, 3;
        jacobian2_sol << 2 * x2[0] + 2 * x2[0] * x2[1], 3 * x2[1] * x2[1] + x2[0] * x2[0];
        hessian2_sol << 2 + 2 * x2[1], 2 * x2[0], 2 * x2[0], 6 * x2[1];

        // problem 3
        // x^2 + 2*y^2 + 3*z^2
        // x*y + 2*x*z
        x3 << 1, 1, 1;
        jacobian3_sol << 2 * x3[0], 4 * x3[1], 6 * x3[2], x3[1] + 2 * x3[2], x3[0], 2 * x3[0];
        hessian3_sol_a << 2, 0, 0, 0, 4, 0, 0, 0, 6;
        hessian3_sol_b << 0, 1, 2, 1, 0, 0, 2, 0, 0;
    }
    // Code here will be called immediately after each test (right
    // before the destructor).
    // virtual void TearDown();

    Eigen::Vector3d x1;
    Eigen::Matrix3d jacobian1;
    Eigen::Matrix3d jacobian1_sol;
    Eigen::Matrix3d hessian1;
    void incX1(int idx, double inc) { x1[idx] += inc; }
    void fun1(Eigen::VectorXd& values)
    {
        values[0] = x1[0] + 2 * x1[1] + 3 * x1[2];
        values[1] = 4 * x1[1];
        values[2] = 5 * x1[1] + 6 * x1[2];
    }
    double multipliers1[3]                             = {2, 3, 4};
    Eigen::Map<const Eigen::Vector3d> multipliers1_vec = Eigen::Map<const Eigen::Vector3d>(multipliers1);

    Eigen::Vector2d x2;
    Eigen::Matrix<double, 1, 2> jacobian2;
    Eigen::Matrix<double, 1, 2> jacobian2_sol;
    Eigen::Matrix2d hessian2;
    Eigen::Matrix2d hessian2_sol;
    void incX2(int idx, double inc) { x2[idx] += inc; }
    void fun2(Eigen::VectorXd& values) { values[0] = x2[0] * x2[0] + x2[1] * x2[1] * x2[1] + x2[0] * x2[0] * x2[1]; }
    double multipliers2[1]                         = {10};
    Eigen::Map<const Eigen::Matrix<double, 1, 1>> multipliers2_vec = Eigen::Map<const Eigen::Matrix<double, 1, 1>>(multipliers2);

    Eigen::Vector3d x3;
    Eigen::Matrix<double, 2, 3> jacobian3;
    Eigen::Matrix<double, 2, 3> jacobian3_sol;
    Eigen::Matrix3d hessian3;
    Eigen::Matrix3d hessian3_sol_a;
    Eigen::Matrix3d hessian3_sol_b;
    void incX3(int idx, double inc) { x3[idx] += inc; }
    void fun3(Eigen::VectorXd& values)
    {
        values[0] = x3[0] * x3[0] + 2 * x3[1] * x3[1] + 3 * x3[2] * x3[2];  // x^2 + 2*y^2 + 3*z^2
        values[1] = x3[0] * x3[1] + 2 * x3[0] * x3[2];                      // x*y + 2*x*z
    }
    double multipliers3[2]                             = {2, 4};
    Eigen::Map<const Eigen::Vector2d> multipliers3_vec = Eigen::Map<const Eigen::Vector2d>(multipliers3);
};

TEST_F(TestFiniteDifferences, forward_differences_jacobian)
{
    ForwardDifferences fd;
    // problem 1
    auto linc1 = [this](int idx, double inc) { incX1(idx, inc); };
    auto lfun1 = [this](Eigen::VectorXd& values) { fun1(values); };
    fd.computeJacobian2(linc1, lfun1, jacobian1);
    EXPECT_EQ_MATRIX(jacobian1, jacobian1_sol, 1e-5);

    // problem 2
    auto linc2 = [this](int idx, double inc) { incX2(idx, inc); };
    auto lfun2 = [this](Eigen::VectorXd& values) { fun2(values); };
    fd.computeJacobian2(linc2, lfun2, jacobian2);
    EXPECT_EQ_MATRIX(jacobian2, jacobian2_sol, 1e-5);

    // problem 3
    auto linc3 = [this](int idx, double inc) { incX3(idx, inc); };
    auto lfun3 = [this](Eigen::VectorXd& values) { fun3(values); };
    fd.computeJacobian2(linc3, lfun3, jacobian3);
    EXPECT_EQ_MATRIX(jacobian3, jacobian3_sol, 1e-5);
}

TEST_F(TestFiniteDifferences, central_differences_jacobian)
{
    CentralDifferences fd;
    // problem 1
    auto linc1 = [this](int idx, double inc) { incX1(idx, inc); };
    auto lfun1 = [this](Eigen::VectorXd& values) { fun1(values); };
    fd.computeJacobian2(linc1, lfun1, jacobian1);
    EXPECT_EQ_MATRIX(jacobian1, jacobian1_sol, 1e-5);

    // problem 2
    auto linc2 = [this](int idx, double inc) { incX2(idx, inc); };
    auto lfun2 = [this](Eigen::VectorXd& values) { fun2(values); };
    fd.computeJacobian2(linc2, lfun2, jacobian2);
    EXPECT_EQ_MATRIX(jacobian2, jacobian2_sol, 1e-5);

    // problem 3
    auto linc3 = [this](int idx, double inc) { incX3(idx, inc); };
    auto lfun3 = [this](Eigen::VectorXd& values) { fun3(values); };
    fd.computeJacobian2(linc3, lfun3, jacobian3);
    EXPECT_EQ_MATRIX(jacobian3, jacobian3_sol, 1e-5);
}

TEST_F(TestFiniteDifferences, forward_differences_hessian)
{
    ForwardDifferences fd;
    // problem 1
    auto linc1 = [this](int idx, double inc) { incX1(idx, inc); };
    auto lfun1 = [this](Eigen::VectorXd& values) { fun1(values); };
    fd.computeHessian2(linc1, lfun1, 3, hessian1);
    EXPECT_ZERO_MATRIX(hessian1, 1e-3);

    // problem 2
    auto linc2 = [this](int idx, double inc) { incX2(idx, inc); };
    auto lfun2 = [this](Eigen::VectorXd& values) { fun2(values); };
    fd.computeHessian2(linc2, lfun2, 1, hessian2);
    EXPECT_EQ_MATRIX(hessian2, hessian2_sol, 1e-5);

    // problem 3
    auto linc3 = [this](int idx, double inc) { incX3(idx, inc); };
    auto lfun3 = [this](Eigen::VectorXd& values) { fun3(values); };
    fd.computeHessian2(linc3, lfun3, 2, hessian3);
    EXPECT_EQ_MATRIX(hessian3, hessian3_sol_a + hessian3_sol_b, 1e-5);
}

TEST_F(TestFiniteDifferences, central_differences_hessian)
{
    CentralDifferences fd;
    // problem 1
    auto linc1 = [this](int idx, double inc) { incX1(idx, inc); };
    auto lfun1 = [this](Eigen::VectorXd& values) { fun1(values); };
    fd.computeHessian2(linc1, lfun1, 3, hessian1);
    EXPECT_ZERO_MATRIX(hessian1, 1e-5);

    // problem 2
    auto linc2 = [this](int idx, double inc) { incX2(idx, inc); };
    auto lfun2 = [this](Eigen::VectorXd& values) { fun2(values); };
    fd.computeHessian2(linc2, lfun2, 1, hessian2);
    EXPECT_EQ_MATRIX(hessian2, hessian2_sol, 1e-5);

    // problem 3
    auto linc3 = [this](int idx, double inc) { incX3(idx, inc); };
    auto lfun3 = [this](Eigen::VectorXd& values) { fun3(values); };
    fd.computeHessian2(linc3, lfun3, 2, hessian3);
    EXPECT_EQ_MATRIX(hessian3, hessian3_sol_a + hessian3_sol_b, 1e-5);
}

TEST_F(TestFiniteDifferences, forward_differences_hessian_multipliers)
{
    ForwardDifferences fd;
    // problem 1
    auto linc1 = [this](int idx, double inc) { incX1(idx, inc); };
    auto lfun1 = [this](Eigen::VectorXd& values) { fun1(values); };
    fd.computeHessian2(linc1, lfun1, 3, hessian1, multipliers1);
    EXPECT_ZERO_MATRIX(hessian1, 1e-3);

    // problem 2
    auto linc2 = [this](int idx, double inc) { incX2(idx, inc); };
    auto lfun2 = [this](Eigen::VectorXd& values) { fun2(values); };
    fd.computeHessian2(linc2, lfun2, 1, hessian2, multipliers2);
    EXPECT_EQ_MATRIX(hessian2, multipliers2[0] * hessian2_sol, 1e-5);

    // problem 3
    auto linc3 = [this](int idx, double inc) { incX3(idx, inc); };
    auto lfun3 = [this](Eigen::VectorXd& values) { fun3(values); };
    fd.computeHessian2(linc3, lfun3, 2, hessian3, multipliers3);
    EXPECT_EQ_MATRIX(hessian3, multipliers3[0] * hessian3_sol_a + multipliers3[1] * hessian3_sol_b, 1e-5);
}

TEST_F(TestFiniteDifferences, central_differences_hessian_multipliers)
{
    CentralDifferences fd;
    // problem 1
    auto linc1 = [this](int idx, double inc) { incX1(idx, inc); };
    auto lfun1 = [this](Eigen::VectorXd& values) { fun1(values); };
    fd.computeHessian2(linc1, lfun1, 3, hessian1, multipliers1);
    EXPECT_ZERO_MATRIX(hessian1, 1e-4);

    // problem 2
    auto linc2 = [this](int idx, double inc) { incX2(idx, inc); };
    auto lfun2 = [this](Eigen::VectorXd& values) { fun2(values); };
    fd.computeHessian2(linc2, lfun2, 1, hessian2, multipliers2);
    EXPECT_EQ_MATRIX(hessian2, multipliers2[0] * hessian2_sol, 1e-5);

    // problem 3
    auto linc3 = [this](int idx, double inc) { incX3(idx, inc); };
    auto lfun3 = [this](Eigen::VectorXd& values) { fun3(values); };
    fd.computeHessian2(linc3, lfun3, 2, hessian3, multipliers3);
    EXPECT_EQ_MATRIX(hessian3, multipliers3[0] * hessian3_sol_a + multipliers3[1] * hessian3_sol_b, 1e-5);
}

TEST_F(TestFiniteDifferences, forward_differences_jacobian_hessian)
{
    ForwardDifferences fd;
    // problem 1
    auto linc1 = [this](int idx, double inc) { incX1(idx, inc); };
    auto lfun1 = [this](Eigen::VectorXd& values) { fun1(values); };
    fd.computeJacobianAndHessian2(linc1, lfun1, jacobian1, hessian1);
    EXPECT_EQ_MATRIX(jacobian1, jacobian1_sol, 1e-5);
    EXPECT_ZERO_MATRIX(hessian1, 1e-3);

    // problem 2
    auto linc2 = [this](int idx, double inc) { incX2(idx, inc); };
    auto lfun2 = [this](Eigen::VectorXd& values) { fun2(values); };
    fd.computeJacobianAndHessian2(linc2, lfun2, jacobian2, hessian2);
    EXPECT_EQ_MATRIX(jacobian2, jacobian2_sol, 1e-5);
    EXPECT_EQ_MATRIX(hessian2, hessian2_sol, 1e-5);

    // problem 3
    auto linc3 = [this](int idx, double inc) { incX3(idx, inc); };
    auto lfun3 = [this](Eigen::VectorXd& values) { fun3(values); };
    fd.computeJacobianAndHessian2(linc3, lfun3, jacobian3, hessian3);
    EXPECT_EQ_MATRIX(jacobian3, jacobian3_sol, 1e-5);
    EXPECT_EQ_MATRIX(hessian3, hessian3_sol_a + hessian3_sol_b, 1e-5);
}

TEST_F(TestFiniteDifferences, central_differences_jacobian_hessian)
{
    CentralDifferences fd;
    // problem 1
    auto linc1 = [this](int idx, double inc) { incX1(idx, inc); };
    auto lfun1 = [this](Eigen::VectorXd& values) { fun1(values); };
    fd.computeJacobianAndHessian2(linc1, lfun1, jacobian1, hessian1);
    EXPECT_EQ_MATRIX(jacobian1, jacobian1_sol, 1e-5);
    EXPECT_ZERO_MATRIX(hessian1, 1e-5);

    // problem 2
    auto linc2 = [this](int idx, double inc) { incX2(idx, inc); };
    auto lfun2 = [this](Eigen::VectorXd& values) { fun2(values); };
    fd.computeJacobianAndHessian2(linc2, lfun2, jacobian2, hessian2);
    EXPECT_EQ_MATRIX(jacobian2, jacobian2_sol, 1e-5);
    EXPECT_EQ_MATRIX(hessian2, hessian2_sol, 1e-5);

    // problem 3
    auto linc3 = [this](int idx, double inc) { incX3(idx, inc); };
    auto lfun3 = [this](Eigen::VectorXd& values) { fun3(values); };
    fd.computeJacobianAndHessian2(linc3, lfun3, jacobian3, hessian3);
    EXPECT_EQ_MATRIX(jacobian3, jacobian3_sol, 1e-5);
    EXPECT_EQ_MATRIX(hessian3, hessian3_sol_a + hessian3_sol_b, 1e-5);
}

TEST_F(TestFiniteDifferences, forward_differences_jacobian_hessian_multipliers)
{
    ForwardDifferences fd;
    // problem 1
    auto linc1 = [this](int idx, double inc) { incX1(idx, inc); };
    auto lfun1 = [this](Eigen::VectorXd& values) { fun1(values); };
    fd.computeJacobianAndHessian2(linc1, lfun1, jacobian1, hessian1, multipliers1);
    Eigen::MatrixXd jacob1_scaled_sol = jacobian1_sol.array().colwise() * multipliers1_vec.array();
    EXPECT_EQ_MATRIX(jacobian1, jacob1_scaled_sol, 1e-5);
    EXPECT_ZERO_MATRIX(hessian1, 1e-3);

    // problem 2
    auto linc2 = [this](int idx, double inc) { incX2(idx, inc); };
    auto lfun2 = [this](Eigen::VectorXd& values) { fun2(values); };
    fd.computeJacobianAndHessian2(linc2, lfun2, jacobian2, hessian2, multipliers2);
    Eigen::MatrixXd jacob2_scaled_sol = jacobian2_sol.array().colwise() * multipliers2_vec.array();
    EXPECT_EQ_MATRIX(jacobian2, jacob2_scaled_sol, 1e-5);
    EXPECT_EQ_MATRIX(hessian2, multipliers2[0] * hessian2_sol, 1e-5);

    // problem 3
    auto linc3 = [this](int idx, double inc) { incX3(idx, inc); };
    auto lfun3 = [this](Eigen::VectorXd& values) { fun3(values); };
    fd.computeJacobianAndHessian2(linc3, lfun3, jacobian3, hessian3, multipliers3);
    Eigen::MatrixXd jacob3_scaled_sol = jacobian3_sol.array().colwise() * multipliers3_vec.array();
    EXPECT_EQ_MATRIX(jacobian3, jacob3_scaled_sol, 1e-5);
    EXPECT_EQ_MATRIX(hessian3, multipliers3[0] * hessian3_sol_a + multipliers3[1] * hessian3_sol_b, 1e-5);
}

TEST_F(TestFiniteDifferences, central_differences_jacobian_hessian_multipliers)
{
    CentralDifferences fd;
    // problem 1
    auto linc1 = [this](int idx, double inc) { incX1(idx, inc); };
    auto lfun1 = [this](Eigen::VectorXd& values) { fun1(values); };
    fd.computeJacobianAndHessian2(linc1, lfun1, jacobian1, hessian1, multipliers1);
    Eigen::MatrixXd jacob1_scaled_sol = jacobian1_sol.array().colwise() * multipliers1_vec.array();
    EXPECT_EQ_MATRIX(jacobian1, jacob1_scaled_sol, 1e-5);
    EXPECT_ZERO_MATRIX(hessian1, 1e-4);

    // problem 2
    auto linc2 = [this](int idx, double inc) { incX2(idx, inc); };
    auto lfun2 = [this](Eigen::VectorXd& values) { fun2(values); };
    fd.computeJacobianAndHessian2(linc2, lfun2, jacobian2, hessian2, multipliers2);
    Eigen::MatrixXd jacob2_scaled_sol = jacobian2_sol.array().colwise() * multipliers2_vec.array();
    EXPECT_EQ_MATRIX(jacobian2, jacob2_scaled_sol, 1e-5);
    EXPECT_EQ_MATRIX(hessian2, multipliers2[0] * hessian2_sol, 1e-5);

    // problem 3
    auto linc3 = [this](int idx, double inc) { incX3(idx, inc); };
    auto lfun3 = [this](Eigen::VectorXd& values) { fun3(values); };
    fd.computeJacobianAndHessian2(linc3, lfun3, jacobian3, hessian3, multipliers3);
    Eigen::MatrixXd jacob3_scaled_sol = jacobian3_sol.array().colwise() * multipliers3_vec.array();
    EXPECT_EQ_MATRIX(jacobian3, jacob3_scaled_sol, 1e-5);
    EXPECT_EQ_MATRIX(hessian3, multipliers3[0] * hessian3_sol_a + multipliers3[1] * hessian3_sol_b, 1e-5);
}
