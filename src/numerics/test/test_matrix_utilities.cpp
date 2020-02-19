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

#include <corbo-numerics/matrix_utilities.h>

#include <corbo-core/console.h>
#include <corbo-core/macros.h>
#include <corbo-core/value_comparison.h>

#include <Eigen/Eigenvalues>

#include "gtest/gtest.h"

using corbo::is_square;
using corbo::have_equal_size;
using corbo::is_positive_definite;

class TestMatrixUtilities : public testing::Test
{
 protected:
    // You can do set-up work for each test here.
    TestMatrixUtilities() {}
    // You can do clean-up work that doesn't throw exceptions here.
    virtual ~TestMatrixUtilities() {}
    // If the constructor and destructor are not enough for setting up
    // and cleaning up each test, you can define the following methods:

    // Code here will be called immediately after the constructor (right
    // before each test).
    // virtual void SetUp() {}
    // Code here will be called immediately after each test (right
    // before the destructor).
    // virtual void TearDown();
};

TEST_F(TestMatrixUtilities, is_square)
{
    Eigen::MatrixXd square_mat = Eigen::MatrixXd::Random(15, 15);
    EXPECT_TRUE(is_square(square_mat));

    Eigen::MatrixXd non_square_mat = Eigen::MatrixXd::Random(16, 15);
    EXPECT_FALSE(is_square(non_square_mat));

    Eigen::VectorXd vector = Eigen::VectorXd::Random(10);
    EXPECT_FALSE(is_square(vector));
}

TEST_F(TestMatrixUtilities, have_equal_size)
{
    Eigen::MatrixXd mat1 = Eigen::MatrixXd::Random(15, 15);
    Eigen::MatrixXd mat2 = Eigen::MatrixXd::Random(14, 15);
    Eigen::VectorXd vec1 = Eigen::VectorXd::Random(15);
    Eigen::VectorXd vec2 = Eigen::VectorXd::Random(14);

    EXPECT_TRUE(have_equal_size(mat1, mat1));
    EXPECT_FALSE(have_equal_size(mat1, mat2));
    EXPECT_FALSE(have_equal_size(mat1, vec1));
    EXPECT_FALSE(have_equal_size(mat1, vec2));

    EXPECT_TRUE(have_equal_size(mat2, mat2));
    EXPECT_FALSE(have_equal_size(mat2, vec1));
    EXPECT_FALSE(have_equal_size(mat2, vec2));

    EXPECT_TRUE(have_equal_size(vec1, vec1));
    EXPECT_FALSE(have_equal_size(vec1, vec2));

    EXPECT_TRUE(have_equal_size(vec2, vec2));
}

TEST_F(TestMatrixUtilities, have_equal_size_multiple)
{
    Eigen::MatrixXd mat1 = Eigen::MatrixXd::Random(15, 15);
    Eigen::MatrixXd mat2 = Eigen::MatrixXd::Random(14, 15);
    Eigen::VectorXd vec1 = Eigen::VectorXd::Random(15);
    Eigen::VectorXd vec2 = Eigen::VectorXd::Random(14);

    EXPECT_TRUE(have_equal_size(mat1, mat1, mat1, mat1, mat1, mat1, mat1));
    EXPECT_FALSE(have_equal_size(mat1, mat1, mat1, mat2, mat1, mat1, mat1));
    EXPECT_FALSE(have_equal_size(mat1, mat1, mat1, mat1, vec1, mat1, mat1));
    EXPECT_FALSE(have_equal_size(mat1, mat1, mat1, mat1, mat1, vec2, mat1));

    EXPECT_TRUE(have_equal_size(mat2, mat2, mat2));
    EXPECT_FALSE(have_equal_size(mat2, mat2, mat1));
    EXPECT_FALSE(have_equal_size(mat2, mat1, mat2));
    EXPECT_FALSE(have_equal_size(mat1, mat2, mat2));

    EXPECT_FALSE(have_equal_size(vec1, mat1, mat1));
    EXPECT_FALSE(have_equal_size(vec1, mat2, mat2));

    EXPECT_TRUE(have_equal_size(vec1, vec1, vec1));
    EXPECT_TRUE(have_equal_size(vec2, vec2, vec2));
    EXPECT_FALSE(have_equal_size(vec1, vec2, vec1));

    EXPECT_FALSE(have_equal_size(mat1, mat2, vec1, vec2));
}

TEST_F(TestMatrixUtilities, is_positive_definite)
{
    Eigen::MatrixXd identity = Eigen::MatrixXd::Identity(15, 15);
    EXPECT_TRUE(is_positive_definite(identity));

    Eigen::MatrixXd zeros = Eigen::MatrixXd::Zero(30, 30);
    EXPECT_FALSE(is_positive_definite(zeros));

    Eigen::MatrixXd pos_semidef = Eigen::MatrixXd::Ones(15, 15);
    EXPECT_FALSE(is_positive_definite(pos_semidef));

    Eigen::MatrixXd non_square_mat = Eigen::MatrixXd::Random(12, 14);
    EXPECT_FALSE(is_positive_definite(non_square_mat));

    Eigen::VectorXd vector = Eigen::VectorXd::Random(10);
    EXPECT_FALSE(is_positive_definite(vector));

    Eigen::Matrix3d pos_def_1;
    pos_def_1 << 2, -1, 0, -1, 2, -1, 0, -1, 2;
    EXPECT_TRUE(is_positive_definite(pos_def_1));

    Eigen::Matrix2d pos_def_2;
    pos_def_2 << 1, 2, 2, 100;
    EXPECT_TRUE(is_positive_definite(pos_def_2));

    Eigen::Matrix2d pos_def_3;
    pos_def_3 << 1, -1, -1, 4;
    EXPECT_TRUE(is_positive_definite(pos_def_3));

    Eigen::Matrix4d pos_def_4;
    pos_def_4 << 9, 3, -6, 12, 3, 26, -7, -11, -6, -7, 9, 7, 12, -11, 7, 65;
    EXPECT_TRUE(is_positive_definite(pos_def_4));

    Eigen::Matrix4d neg_def_1;
    neg_def_1 << 9, 3, -6, 12, 3, 26, -7, -11, -6, -7, 9, 7, 12, -11, 7, 60;
    EXPECT_FALSE(is_positive_definite(neg_def_1));

    Eigen::Matrix2d neg_def_2;
    neg_def_2 << 1, 4, 4, 1;
    EXPECT_FALSE(is_positive_definite(neg_def_2));

    Eigen::Matrix3d pos_def_diag = Eigen::Vector3d(1, 3, 2).asDiagonal();
    EXPECT_TRUE(is_positive_definite(pos_def_diag));

    Eigen::Matrix3d neg_def_diag = Eigen::Vector3d(1, 3, -2).asDiagonal();
    EXPECT_FALSE(is_positive_definite(neg_def_diag));
}
