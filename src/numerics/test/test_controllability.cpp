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

#include <corbo-numerics/controllability.h>

#include <corbo-core/macros.h>

#include "gtest/gtest.h"

using corbo::Controllability;

class TestControllability : public testing::Test
{
 protected:
    // You can do set-up work for each test here.
    TestControllability() {}
    // You can do clean-up work that doesn't throw exceptions here.
    virtual ~TestControllability() {}
    // If the constructor and destructor are not enough for setting up
    // and cleaning up each test, you can define the following methods:

    // Code here will be called immediately after the constructor (right
    // before each test).
    // virtual void SetUp() {}
    // Code here will be called immediately after each test (right
    // before the destructor).
    // virtual void TearDown();
};

TEST_F(TestControllability, controllable_lti_system)
{
    Eigen::Matrix2d A;
    A << 0, 1, 0, 0;

    Eigen::Vector2d B;
    B << 0, 1;

    EXPECT_TRUE(Controllability::checkLinearTimeInvariantSystem(A, B));
}

TEST_F(TestControllability, uncontrollable_lti_system)
{
    Eigen::Matrix4d A;
    A.setZero();
    A(0, 1) = 1;
    A(2, 3) = 1;

    Eigen::Vector4d B;
    B << 0, 1, 0, -1;

    int ctrl_rank = 1000;
    EXPECT_FALSE(Controllability::checkLinearTimeInvariantSystem(A, B, &ctrl_rank));
    EXPECT_EQ(ctrl_rank, 2);
}

TEST_F(TestControllability, uncontrollable_lti_system2)
{
    Eigen::Matrix2d A;
    A << 1, 1, 4, -2;

    Eigen::Matrix2d B;
    B << 1, -1, 1, -1;

    EXPECT_FALSE(Controllability::checkLinearTimeInvariantSystem(A, B));
}
