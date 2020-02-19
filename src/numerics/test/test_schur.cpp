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

#include <corbo-numerics/schur.h>

#include <corbo-core/console.h>
#include <corbo-core/macros.h>
#include <corbo-core/value_comparison.h>

#include <Eigen/Eigenvalues>

#include "gtest/gtest.h"

using corbo::schur_decomposition_2d;
using corbo::swap_schur_blocks;
using corbo::reorder_schur_blocks;
using corbo::essentially_equal;
using corbo::approx_zero;

class TestSchur : public testing::Test
{
 protected:
    // You can do set-up work for each test here.
    TestSchur() {}
    // You can do clean-up work that doesn't throw exceptions here.
    virtual ~TestSchur() {}
    // If the constructor and destructor are not enough for setting up
    // and cleaning up each test, you can define the following methods:

    // Code here will be called immediately after the constructor (right
    // before each test).
    // virtual void SetUp() {}
    // Code here will be called immediately after each test (right
    // before the destructor).
    // virtual void TearDown();
};

TEST_F(TestSchur, schur_decomposition_2d1)
{
    Eigen::Matrix2d T;
    T << -5, -6, 1, 0;  // only real eigenvalues

    Eigen::Matrix2d U;
    schur_decomposition_2d(T, U);

    Eigen::Matrix2d Tsol;
    Tsol << -3, -7, 0, -2;

    Eigen::Matrix2d Usol;
    Usol << -0.9487, -0.3162, 0.3162, -0.9487;

    EXPECT_EQ_MATRIX(T, Tsol, 1e-4);
    EXPECT_EQ_MATRIX(U, Usol, 1e-4);
}

TEST_F(TestSchur, schur_decomposition_2d2)
{
    Eigen::Matrix2d T;
    T << -3, -7, 0, -2;  // already schur form

    Eigen::Matrix2d U;
    schur_decomposition_2d(T, U);

    Eigen::Matrix2d Tsol;
    Tsol << -3, -7, 0, -2;

    Eigen::Matrix2d Usol = Eigen::Matrix2d::Identity();

    EXPECT_EQ_MATRIX(T, Tsol, 1e-4);
    EXPECT_EQ_MATRIX(U, Usol, 1e-4);
}

TEST_F(TestSchur, schur_decomposition_2d3)
{
    Eigen::Matrix2d T;
    T << -1, -1.25, 1, 0;  // complex eigenvalues

    Eigen::Matrix2d U;
    schur_decomposition_2d(T, U);

    Eigen::Matrix2d Tsol;
    Tsol << -0.5, -1.6404, 0.6096, -0.5;

    Eigen::Matrix2d Usol;
    Usol << 0.7882, 0.6154, -0.6154, 0.7882;

    EXPECT_EQ_MATRIX(T, Tsol, 1e-4);
    EXPECT_EQ_MATRIX(U, Usol, 1e-4);
}

TEST_F(TestSchur, schur_decomposition_2d4)
{
    Eigen::Matrix2d T;
    T << -0.393678, -0.73017, 0.78598, -0.606322;  // already schur form, but non-standard

    Eigen::Matrix2d U;
    schur_decomposition_2d(T, U);

    Eigen::Matrix2d Tsol;
    Tsol << -0.5000, -0.6482, 0.8680, -0.5000;

    Eigen::Matrix2d Usol;
    Usol << 0.7918, 0.6108, -0.6108, 0.7918;

    EXPECT_EQ_MATRIX(T, Tsol, 1e-4);
    EXPECT_EQ_MATRIX(U, Usol, 1e-4);
}

TEST_F(TestSchur, schur_decomposition_2d_b_zero)
{
    // element T(0,1) = 0 (special implementation)
    Eigen::Matrix2d T;
    T << -5, 0, 1, 0;  // only real eigenvalues

    Eigen::Matrix2d U;
    schur_decomposition_2d(T, U);

    Eigen::Matrix2d Tsol;
    Tsol << 0, -1, 0, -5;

    Eigen::Matrix2d Usol;
    Usol << 0, -1, 1, 0;

    EXPECT_EQ_MATRIX(T, Tsol, 1e-4);
    EXPECT_EQ_MATRIX(U, Usol, 1e-4);
}

TEST_F(TestSchur, schur_decomposition_2d_a_minus_d_zero1)
{
    // special case: T(0,0) - T(1,1) = 0, but sign(1,b) != sign(1,c)
    Eigen::Matrix2d T;
    T << -5, -6, 1, -5;  // only real eigenvalues

    Eigen::Matrix2d U;
    schur_decomposition_2d(T, U);

    Eigen::Matrix2d Tsol;
    Tsol << -5, -6, 1, -5;

    Eigen::Matrix2d Usol;
    Usol << 1, 0, 0, 1;

    EXPECT_EQ_MATRIX(T, Tsol, 1e-4);
    EXPECT_EQ_MATRIX(U, Usol, 1e-4);
}

TEST_F(TestSchur, schur_decomposition_2d_a_minus_d_zero2)
{
    // special case: T(0,0) - T(1,1) = 0, but sign(1,b) == sign(1,c)
    Eigen::Matrix2d T;
    T << -5, -6, -1, -5;  // only real eigenvalues

    Eigen::Matrix2d Rsol = T;

    Eigen::Matrix2d U;
    schur_decomposition_2d(T, U);

    // in this example, our solution exhibits swapped eigenvalues compared to matlab.
    // so we check validity by transforming the problem back

    Eigen::MatrixXd R = U * T * U.adjoint();
    EXPECT_EQ_MATRIX(R, Rsol, 1e-4);

    //    Eigen::Matrix2d Tsol;
    //    Tsol << -2.5505, -5.0, 0, -7.4495;

    //    Eigen::Matrix2d Usol;
    //    Usol << 0.9258, 0.3780, -0.3780, 0.9258;

    //    EXPECT_EQ_MATRIX(T, Tsol, 1e-4);
    //    EXPECT_EQ_MATRIX(U, Usol, 1e-4);
}

TEST_F(TestSchur, swap_blocks_schur_form1)
{
    Eigen::Matrix2d T;
    T << -1.0000, 5.0000, 0, -2.0000;

    Eigen::Matrix2d Q;
    Q.setIdentity();
    swap_schur_blocks(T, 0, 1, 1, Q);
    EXPECT_EQ_MATRIX((Q.transpose() * Q).eval(), Eigen::Matrix2d::Identity(), 1e-5);

    EXPECT_NEAR(T(0, 0), -2, 1e-10);
    EXPECT_NEAR(T(1, 1), -1, 1e-10);
}

TEST_F(TestSchur, swap_blocks_schur_form2)
{
    // Eigen::Matrix3d A;
    // A << -2, -1.8125, -0.8125, 1, 0, 0, 0, 1, 0;

    Eigen::Matrix3d T;
    T << -1, -0.8378, 1.3184, 0, -0.5, 2.4397, 0, -0.2306, -0.5;

    Eigen::MatrixXd Q(3, 3);
    Q.setIdentity();
    swap_schur_blocks(T, 0, 1, 2, Q);
    EXPECT_EQ_MATRIX((Q.transpose() * Q).eval(), Eigen::Matrix3d::Identity(), 1e-5);

    // check eigenvalues A11
    Eigen::VectorXcd eig1 = T.block(0, 0, 2, 2).eigenvalues();
    EXPECT_TRUE(essentially_equal(eig1[0], std::complex<double>(-0.5, 0.75), 1e-3) ||
                essentially_equal(eig1[0], std::complex<double>(-0.5, -0.75), 1e-3))
        << "Eigenvalue: " << eig1[0];
    EXPECT_TRUE(essentially_equal(eig1[1], std::complex<double>(-0.5, 0.75), 1e-3) ||
                essentially_equal(eig1[1], std::complex<double>(-0.5, -0.75), 1e-3))
        << "Eigenvalue: " << eig1[1];

    EXPECT_NEAR(T(2, 2), -1, 1e-10);
}

TEST_F(TestSchur, swap_blocks_schur_form3_standardize)
{
    // Eigen::MatrixXd A(6, 6);
    // A << -17, -105.25, -332, -590, -536, -320, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0;

    // -> Schur-Decomposition
    //
    // T =
    //   -8.0000   30.7409  -57.5046 -247.0842 -169.7203  867.0578
    //         0   -4.0000    7.0149   30.1406   20.7053 -105.7696
    //         0         0   -2.0000   -7.5820   -5.1700   26.3801
    //         0         0    0.5276   -2.0000   -1.3243    6.2260
    //         0         0         0         0   -0.5000    3.4406
    //         0         0         0         0   -0.2906   -0.5000

    Eigen::MatrixXd T(6, 6);
    T << -8.0000, 30.7409, -57.5046, -247.0842, -169.7203, 867.0578, 0, -4.0000, 7.0149, 30.1406, 20.7053, -105.7696, 0, 0, -2.0000, -7.5820, -5.1700,
        26.3801, 0, 0, 0.5276, -2.0000, -1.3243, 6.2260, 0, 0, 0, 0, -0.5000, 3.4406, 0, 0, 0, 0, -0.2906, -0.5000;

    Eigen::MatrixXd T_original = T;  // backup for final test

    Eigen::MatrixXd Q1 = Eigen::MatrixXd::Identity(6, 6);
    swap_schur_blocks(T, 2, 2, 2, Q1, true);
    EXPECT_EQ_MATRIX((Q1.transpose() * Q1).eval(), Eigen::MatrixXd::Identity(6, 6), 1e-5);
    // check eigenvalues A11
    Eigen::VectorXcd eig1 = T.block(2, 2, 2, 2).eigenvalues();
    EXPECT_TRUE(essentially_equal(eig1[0], std::complex<double>(-0.5, 1), 1e-4) || essentially_equal(eig1[0], std::complex<double>(-0.5, -1), 1e-4))
        << "Eigenvalue: " << eig1[0];
    EXPECT_TRUE(essentially_equal(eig1[1], std::complex<double>(-0.5, 1), 1e-4) || essentially_equal(eig1[1], std::complex<double>(-0.5, -1), 1e-4))
        << "Eigenvalue: " << eig1[1];
    // check eigenvalues A22
    Eigen::VectorXcd eig2 = T.block(4, 4, 2, 2).eigenvalues();
    EXPECT_TRUE(essentially_equal(eig2[0], std::complex<double>(-2, 2), 1e-4) || essentially_equal(eig2[0], std::complex<double>(-2, -2), 1e-4))
        << "Eigenvalue: " << eig2[0];
    EXPECT_TRUE(essentially_equal(eig2[1], std::complex<double>(-2, 2), 1e-4) || essentially_equal(eig2[1], std::complex<double>(-2, -2), 1e-4))
        << "Eigenvalue: " << eig2[1];

    // check standard form
    auto is_standard_form = [](const Eigen::Ref<const Eigen::Matrix2d>& block) {
        if (essentially_equal(block(0, 0), block(1, 1), 1e-15) && block(0, 1) * block(1, 0) < 0) return true;
        return false;
    };
    EXPECT_TRUE(is_standard_form(T.block(2, 2, 2, 2))) << "T:\n" << T.block(2, 2, 2, 2);
    EXPECT_TRUE(is_standard_form(T.block(4, 4, 2, 2))) << "T:\n" << T.block(4, 4, 2, 2);

    // now swap first two eigenvalues
    Eigen::MatrixXd Q2 = Eigen::MatrixXd::Identity(6, 6);
    swap_schur_blocks(T, 0, 1, 1, Q2, true);
    EXPECT_EQ_MATRIX((Q2.transpose() * Q2).eval(), Eigen::MatrixXd::Identity(6, 6), 1e-5);
    EXPECT_NEAR(T(0, 0), -4, 1e-10);
    EXPECT_NEAR(T(1, 1), -8, 1e-10);

    // now swap 2nd and 3rd eigenvalues
    Eigen::MatrixXd Q3 = Eigen::MatrixXd::Identity(6, 6);
    swap_schur_blocks(T, 1, 1, 2, Q3, true);

    EXPECT_EQ_MATRIX((Q3.transpose() * Q3).eval(), Eigen::MatrixXd::Identity(6, 6), 1e-5);
    // check eigenvalues A11
    eig1 = T.block(1, 1, 2, 2).eigenvalues();
    EXPECT_TRUE(essentially_equal(eig1[0], std::complex<double>(-0.5, 1), 1e-4) || essentially_equal(eig1[0], std::complex<double>(-0.5, -1), 1e-4))
        << "Eigenvalue: " << eig1[0];
    EXPECT_TRUE(essentially_equal(eig1[1], std::complex<double>(-0.5, 1), 1e-4) || essentially_equal(eig1[1], std::complex<double>(-0.5, -1), 1e-4))
        << "Eigenvalue: " << eig1[1];
    EXPECT_NEAR(T(3, 3), -8, 1e-10);
    EXPECT_TRUE(is_standard_form(T.block(1, 1, 2, 2))) << "T:\n" << T.block(1, 1, 2, 2);

    // now swap again 2nd and 3rd to test p=2, q=1 case
    Eigen::MatrixXd Q4 = Eigen::MatrixXd::Identity(6, 6);
    swap_schur_blocks(T, 1, 2, 1, Q4, true);

    EXPECT_EQ_MATRIX((Q4.transpose() * Q4).eval(), Eigen::MatrixXd::Identity(6, 6), 1e-5);
    // check eigenvalues A11
    eig1 = T.block(2, 2, 2, 2).eigenvalues();
    EXPECT_TRUE(essentially_equal(eig1[0], std::complex<double>(-0.5, 1), 1e-4) || essentially_equal(eig1[0], std::complex<double>(-0.5, -1), 1e-4))
        << "Eigenvalue: " << eig1[0];
    EXPECT_TRUE(essentially_equal(eig1[1], std::complex<double>(-0.5, 1), 1e-4) || essentially_equal(eig1[1], std::complex<double>(-0.5, -1), 1e-4))
        << "Eigenvalue: " << eig1[1];
    EXPECT_NEAR(T(1, 1), -8, 1e-10);
    EXPECT_TRUE(is_standard_form(T.block(2, 2, 2, 2))) << "T:\n" << T.block(1, 1, 2, 2);

    // final test, transform T back to the original matrix
    Eigen::MatrixXd Tf = Q1 * Q2 * Q3 * Q4 * T * Q4.transpose() * Q3.transpose() * Q2.transpose() * Q1.transpose();
    EXPECT_EQ_MATRIX(Tf, T_original, 1e-4);
}

TEST_F(TestSchur, swap_blocks_schur_form3_non_standardize)
{
    Eigen::MatrixXd T(6, 6);
    T << -8.0000, 30.7409, -57.5046, -247.0842, -169.7203, 867.0578, 0, -4.0000, 7.0149, 30.1406, 20.7053, -105.7696, 0, 0, -2.0000, -7.5820, -5.1700,
        26.3801, 0, 0, 0.5276, -2.0000, -1.3243, 6.2260, 0, 0, 0, 0, -0.5000, 3.4406, 0, 0, 0, 0, -0.2906, -0.5000;

    Eigen::MatrixXd T_original = T;  // backup for final test

    Eigen::MatrixXd Q1 = Eigen::MatrixXd::Identity(6, 6);
    swap_schur_blocks(T, 2, 2, 2, Q1, false);
    EXPECT_EQ_MATRIX((Q1.transpose() * Q1).eval(), Eigen::MatrixXd::Identity(6, 6), 1e-5);
    // check eigenvalues A11
    Eigen::VectorXcd eig1 = T.block(2, 2, 2, 2).eigenvalues();
    EXPECT_TRUE(essentially_equal(eig1[0], std::complex<double>(-0.5, 1), 1e-4) || essentially_equal(eig1[0], std::complex<double>(-0.5, -1), 1e-4))
        << "Eigenvalue: " << eig1[0];
    EXPECT_TRUE(essentially_equal(eig1[1], std::complex<double>(-0.5, 1), 1e-4) || essentially_equal(eig1[1], std::complex<double>(-0.5, -1), 1e-4))
        << "Eigenvalue: " << eig1[1];
    // check eigenvalues A22
    Eigen::VectorXcd eig2 = T.block(4, 4, 2, 2).eigenvalues();
    EXPECT_TRUE(essentially_equal(eig2[0], std::complex<double>(-2, 2), 1e-4) || essentially_equal(eig2[0], std::complex<double>(-2, -2), 1e-4))
        << "Eigenvalue: " << eig2[0];
    EXPECT_TRUE(essentially_equal(eig2[1], std::complex<double>(-2, 2), 1e-4) || essentially_equal(eig2[1], std::complex<double>(-2, -2), 1e-4))
        << "Eigenvalue: " << eig2[1];

    // now swap first two eigenvalues
    Eigen::MatrixXd Q2 = Eigen::MatrixXd::Identity(6, 6);
    swap_schur_blocks(T, 0, 1, 1, Q2, false);
    EXPECT_EQ_MATRIX((Q2.transpose() * Q2).eval(), Eigen::MatrixXd::Identity(6, 6), 1e-5);
    EXPECT_NEAR(T(0, 0), -4, 1e-10);
    EXPECT_NEAR(T(1, 1), -8, 1e-10);

    // now swap 2nd and 3rd eigenvalues
    Eigen::MatrixXd Q3 = Eigen::MatrixXd::Identity(6, 6);
    swap_schur_blocks(T, 1, 1, 2, Q3, false);

    EXPECT_EQ_MATRIX((Q3.transpose() * Q3).eval(), Eigen::MatrixXd::Identity(6, 6), 1e-5);
    // check eigenvalues A11
    eig1 = T.block(1, 1, 2, 2).eigenvalues();
    EXPECT_TRUE(essentially_equal(eig1[0], std::complex<double>(-0.5, 1), 1e-4) || essentially_equal(eig1[0], std::complex<double>(-0.5, -1), 1e-4))
        << "Eigenvalue: " << eig1[0];
    EXPECT_TRUE(essentially_equal(eig1[1], std::complex<double>(-0.5, 1), 1e-4) || essentially_equal(eig1[1], std::complex<double>(-0.5, -1), 1e-4))
        << "Eigenvalue: " << eig1[1];
    EXPECT_NEAR(T(3, 3), -8, 1e-10);

    // now swap again 2nd and 3rd to test p=2, q=1 case
    Eigen::MatrixXd Q4 = Eigen::MatrixXd::Identity(6, 6);
    swap_schur_blocks(T, 1, 2, 1, Q4, false);

    EXPECT_EQ_MATRIX((Q4.transpose() * Q4).eval(), Eigen::MatrixXd::Identity(6, 6), 1e-5);
    // check eigenvalues A11
    eig1 = T.block(2, 2, 2, 2).eigenvalues();
    EXPECT_TRUE(essentially_equal(eig1[0], std::complex<double>(-0.5, 1), 1e-4) || essentially_equal(eig1[0], std::complex<double>(-0.5, -1), 1e-4))
        << "Eigenvalue: " << eig1[0];
    EXPECT_TRUE(essentially_equal(eig1[1], std::complex<double>(-0.5, 1), 1e-4) || essentially_equal(eig1[1], std::complex<double>(-0.5, -1), 1e-4))
        << "Eigenvalue: " << eig1[1];
    EXPECT_NEAR(T(1, 1), -8, 1e-10);

    // final test, transform T back to the original matrix
    Eigen::MatrixXd Tf = Q1 * Q2 * Q3 * Q4 * T * Q4.transpose() * Q3.transpose() * Q2.transpose() * Q1.transpose();
    EXPECT_EQ_MATRIX(Tf, T_original, 1e-4);
}

TEST_F(TestSchur, reorder_schur_blocks_real_2d)
{
    Eigen::Matrix2d T;
    T << 8, 31, 0, -4;

    Eigen::Matrix2d Q = Eigen::Matrix2d::Identity();

    int subspace_dim = 0;

    auto is_neg_real = [](const Eigen::Ref<const Eigen::MatrixXd>& block) { return block.eigenvalues()[0].real() < 0; };
    reorder_schur_blocks(T, Q, is_neg_real, &subspace_dim);

    EXPECT_NEAR(T(0, 0), -4, 1e-10);
    EXPECT_NEAR(T(1, 1), 8, 1e-10);
    EXPECT_NEQ_MATRIX(Q, Eigen::Matrix2d::Identity(), 1e-5);
    EXPECT_EQ_MATRIX((Q.transpose() * Q).eval(), Eigen::Matrix2d::Identity(), 1e-5);
    EXPECT_EQ(subspace_dim, 1);

    subspace_dim     = -10;
    auto is_pos_real = [](const Eigen::Ref<const Eigen::MatrixXd>& block) { return block.eigenvalues()[0].real() > 0; };
    reorder_schur_blocks(T, Q, is_pos_real);

    EXPECT_NEAR(T(0, 0), 8, 1e-10);
    EXPECT_NEAR(T(1, 1), -4, 1e-10);
    EXPECT_EQ(subspace_dim, -10);

    auto always_true = [](const Eigen::Ref<const Eigen::MatrixXd>& block) { return true; };
    reorder_schur_blocks(T, Q, always_true);

    EXPECT_NEAR(T(0, 0), 8, 1e-10);
    EXPECT_NEAR(T(1, 1), -4, 1e-10);

    auto always_false = [](const Eigen::Ref<const Eigen::MatrixXd>& block) { return false; };
    reorder_schur_blocks(T, Q, always_false);

    EXPECT_NEAR(T(0, 0), 8, 1e-10);
    EXPECT_NEAR(T(1, 1), -4, 1e-10);
}

TEST_F(TestSchur, reorder_schur_blocks_real_6d)
{
    Eigen::MatrixXd T(6, 6);
    T << 6, -28.9828, 114.6002, 334.8702, 642.8458, 858.6334, 0, 5, -18.9797, -55.4600, -106.4651, -142.2044, 0, 0, 4, 10.9599, 21.0520, 28.0873, 0,
        0, 0, -3, -5.1069, -6.8522, 0, 0, 0, 0, -2, -2.3436, 0, 0, 0, 0, 0, -1;

    Eigen::MatrixXd Q = Eigen::MatrixXd::Identity(6, 6);

    int subspace_dim = 0;

    auto select = [](const Eigen::Ref<const Eigen::MatrixXd>& block) { return block.eigenvalues()[0].real() < -1; };
    reorder_schur_blocks(T, Q, select, &subspace_dim);

    EXPECT_LT(T(0, 0), -1);
    EXPECT_LT(T(1, 1), -1);
    EXPECT_GE(T(2, 2), -1);
    EXPECT_GE(T(3, 3), -1);
    EXPECT_GE(T(5, 5), -1);
    EXPECT_GE(T(5, 5), -1);
    EXPECT_EQ(subspace_dim, 2);

    EXPECT_NEQ_MATRIX(Q, Eigen::MatrixXd::Identity(6, 6), 1e-5);
    EXPECT_EQ_MATRIX((Q.transpose() * Q).eval(), Eigen::MatrixXd::Identity(6, 6), 1e-5);
}

TEST_F(TestSchur, reorder_schur_blocks_complex_4d)
{
    /*
     T =
     2.0000   -5.6559    4.9142  -29.6122
     0.3978    2.0000   -1.4659    8.9097
     0         0        -1.0000    5.0024
     0         0        -0.7996   -1.0000
    */
    Eigen::Matrix4d T;
    T << 2, -5.6559, 4.9142, -29.6122, 0.3978, 2, -1.4659, 8.9097, 0, 0, -1, 5.0024, 0, 0, -0.7996, -1;

    Eigen::Matrix4d Q = Eigen::Matrix4d::Identity();

    int subspace_dim = 0;

    auto select = [](const Eigen::Ref<const Eigen::MatrixXd>& block) { return block.eigenvalues()[0].real() < 0; };
    reorder_schur_blocks(T, Q, select, &subspace_dim);

    double real1 = T.topLeftCorner(2, 2).eigenvalues()[0].real();
    double real2 = T.bottomRightCorner(2, 2).eigenvalues()[0].real();

    EXPECT_TRUE(essentially_equal(real1, -1, 1e-6)) << real1;
    EXPECT_TRUE(essentially_equal(real2, 2, 1e-6)) << real2;

    EXPECT_EQ(subspace_dim, 2);

    EXPECT_NEQ_MATRIX(Q, Eigen::Matrix4d::Identity(), 1e-5);
    EXPECT_EQ_MATRIX((Q.transpose() * Q).eval(), Eigen::Matrix4d::Identity(), 1e-5);
}

TEST_F(TestSchur, reorder_schur_blocks_mixed_3d)
{
    /*
     T =
     2.0000   -5.7972    3.4035
     0.3881    2.0000   -0.8847
     0         0        -1.0000
     */
    Eigen::Matrix3d T;
    T << 2, -5.7972, 3.4035, 0.3881, 2, -0.8847, 0, 0, -1;

    Eigen::Matrix3d Q = Eigen::Matrix3d::Identity();

    int subspace_dim = 0;

    auto select = [](const Eigen::Ref<const Eigen::MatrixXd>& block) { return block.eigenvalues()[0].real() < 0; };
    reorder_schur_blocks(T, Q, select, &subspace_dim);

    double real1 = T.bottomRightCorner(2, 2).eigenvalues()[0].real();
    EXPECT_TRUE(essentially_equal(real1, 2, 1e-6)) << real1;
    EXPECT_NEAR(T(0, 0), -1, 1e-10);
    EXPECT_EQ(subspace_dim, 1);

    EXPECT_NEQ_MATRIX(Q, Eigen::Matrix3d::Identity(), 1e-5);
    EXPECT_EQ_MATRIX((Q.transpose() * Q).eval(), Eigen::Matrix3d::Identity(), 1e-5);

    // swap back
    auto select2 = [](const Eigen::Ref<const Eigen::MatrixXd>& block) { return block.eigenvalues()[0].real() > 0; };
    reorder_schur_blocks(T, Q, select2, &subspace_dim);

    real1 = T.topLeftCorner(2, 2).eigenvalues()[0].real();
    EXPECT_TRUE(essentially_equal(real1, 2, 1e-6)) << real1;

    EXPECT_NEAR(T(2, 2), -1, 1e-10);
    EXPECT_EQ(subspace_dim, 2);

    EXPECT_NEQ_MATRIX(Q, Eigen::Matrix3d::Identity(), 1e-5);
    EXPECT_EQ_MATRIX((Q.transpose() * Q).eval(), Eigen::Matrix3d::Identity(), 1e-5);
}

TEST_F(TestSchur, reorder_schur_blocks_mixed_8d)
{
    /*
     T =
     4    6.8834  -22.8817  -68.4663  -60.1433 -350.9628  643.1324  989.0186
     0    2.0000   -5.6484  -16.7104  -14.6790  -85.6590  156.9690  241.3865
     0    0.3983    2.0000    5.0269    4.4159   25.7680  -47.2178  -72.6187
     0         0         0   -3.0000   -2.3424  -13.6697   25.0501   38.5195
     0         0         0         0   -1.0000   -4.8323    8.5495   13.1573
     0         0         0         0    0.8278   -1.0000    1.4519    2.2920
     0         0         0         0         0         0   -2.0000   -2.5705
     0         0         0         0         0         0         0   -1.0000
     */

    Eigen::MatrixXd T(8, 8);
    T << 4, 6.8834, -22.8817, -68.4663, -60.1433, -350.9628, 643.1324, 989.0186, 0, 2.0000, -5.6484, -16.7104, -14.6790, -85.6590, 156.9690, 241.3865,
        0, 0.3983, 2.0000, 5.0269, 4.4159, 25.7680, -47.2178, -72.6187, 0, 0, 0, -3.0000, -2.3424, -13.6697, 25.0501, 38.5195, 0, 0, 0, 0, -1.0000,
        -4.8323, 8.5495, 13.1573, 0, 0, 0, 0, 0.8278, -1.0000, 1.4519, 2.2920, 0, 0, 0, 0, 0, 0, -2.0000, -2.5705, 0, 0, 0, 0, 0, 0, 0, -1.0000;

    Eigen::MatrixXd Q = Eigen::MatrixXd::Identity(8, 8);

    int subspace_dim = 0;

    auto select = [](const Eigen::Ref<const Eigen::MatrixXd>& block) { return block.eigenvalues()[0].real() < 0; };
    reorder_schur_blocks(T, Q, select, &subspace_dim);

    // search for blocks
    int row = 0;
    while (row < T.rows())
    {
        // check, if we have a 1 x 1 or 2 x 2 block
        int p;
        if (row + 1 >= T.rows())
            p = 1;
        else if (approx_zero(T(row + 1, row), 1e-10))
            p = 1;
        else
            p = 2;

        if (p == 1)
        {
            if (row < 5)
                EXPECT_LT(T(row, row), 0);
            else
                EXPECT_GE(T(row, row), 0);
        }
        else
        {
            double real = T.block(row, row, 2, 2).eigenvalues()[0].real();
            if (row < 5)
                EXPECT_LT(real, 0);
            else
                EXPECT_GE(real, 0);
        }

        row += p;
    }
    EXPECT_EQ(subspace_dim, 5);

    EXPECT_NEQ_MATRIX(Q, Eigen::MatrixXd::Identity(8, 8), 1e-5);
    EXPECT_EQ_MATRIX((Q.transpose() * Q).eval(), Eigen::MatrixXd::Identity(8, 8), 1e-5);
}
