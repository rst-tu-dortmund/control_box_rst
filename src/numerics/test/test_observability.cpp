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

#include <corbo-numerics/observability.h>

#include <corbo-core/macros.h>

#include "gtest/gtest.h"
using corbo::Observability;

class TestObservability : public testing::Test
{
 protected:
    // You can do set-up work for each test here.
    TestObservability() {}
    // You can do clean-up work that doesn't throw exceptions here.
    virtual ~TestObservability() {}
    // If the constructor and destructor are not enough for setting up
    // and cleaning up each test, you can define the following methods:

    // Code here will be called immediately after the constructor (right
    // before each test).
    // virtual void SetUp() {}
    // Code here will be called immediately after each test (right
    // before the destructor).
    // virtual void TearDown();
};

TEST_F(TestObservability, observable_lti_system)
{
    Eigen::Matrix2d A;
    A << 1, 1, 4, -2;

    Eigen::Matrix2d C;
    C << 1, 0, 0, 1;

    EXPECT_TRUE(Observability::checkLinearTimeInvariantSystem(A, C));
}

TEST_F(TestObservability, observable_lti_system2)
{
    Eigen::Matrix2d A;
    A << 1, 2, 3, 4;

    Eigen::RowVector2d C;
    C << 1, 2;

    EXPECT_TRUE(Observability::checkLinearTimeInvariantSystem(A, C));
}

TEST_F(TestObservability, unobservable_lti_system)
{
    Eigen::Matrix2d A;
    A << 1, -2, -3, -4;

    Eigen::RowVector2d C;
    C << 1, 2;

    int obs_rank = 100;
    EXPECT_FALSE(Observability::checkLinearTimeInvariantSystem(A, C, &obs_rank));
    EXPECT_EQ(obs_rank, 1);
}

TEST_F(TestObservability, unobservable_lti_system2)
{
    Eigen::Matrix4d A;
    A.setZero();
    A(0, 1) = 1;
    A(2, 3) = 1;

    Eigen::RowVector4d C;
    C << 0, 1, 0, -1;

    int obs_rank = 100;
    EXPECT_FALSE(Observability::checkLinearTimeInvariantSystem(A, C, &obs_rank));
    EXPECT_EQ(obs_rank, 1);
}
