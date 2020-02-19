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

#include <iostream>

#include "gtest/gtest.h"

// The fixture for testing class Foo.
class FooTest : public ::testing::Test
{
 protected:
    // You can do set-up work for each test here.
    FooTest() {}
    // You can do clean-up work that doesn't throw exceptions here.
    virtual ~FooTest() {}
    // If the constructor and destructor are not enough for setting up
    // and cleaning up each test, you can define the following methods:

    // Code here will be called immediately after the constructor (right
    // before each test).
    virtual void SetUp() { std::cout << "nice" << std::endl; }
    // Code here will be called immediately after each test (right
    // before the destructor).
    // virtual void TearDown();

    // The mock bar library shaed by all tests
    // int m_bar;
};

TEST_F(FooTest, blabla)
{
    int testval = 5;
    EXPECT_EQ(testval, 5);
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    int ret = RUN_ALL_TESTS();

    return ret;
    // return 1;
}

// int main()
//{
//  std::cout << "Hello World!" << std::endl;
//  return 0;
//}
