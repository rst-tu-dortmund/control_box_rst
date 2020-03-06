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

#include <corbo-controllers/pid_controller.h>
#include <corbo-systems/output_function_interface.h>

#include <corbo-core/console.h>
#include <thread>

#include "gtest/gtest.h"

// The fixture for testing class Foo.
class PidControllerTest : public testing::Test
{
 protected:
    // You can do set-up work for each test here.
    PidControllerTest() {}
    // You can do clean-up work that doesn't throw exceptions here.
    virtual ~PidControllerTest() {}
    // If the constructor and destructor are not enough for setting up
    // and cleaning up each test, you can define the following methods:

    // Code here will be called immediately after the constructor (right
    // before each test).
    virtual void SetUp() { _controller.reset(); }
    // Code here will be called immediately after each test (right
    // before the destructor).
    // virtual void TearDown();

    corbo::PidController _controller;
};

// class MockSystemOutput : public corbo::SystemOutputInterface
//{
//  public:
//    MOCK_METHOD2(output, void(const Eigen::Ref<const StateVector>& x, Eigen::Ref<OutputVector> y));
//};
// we expect the input values to be x and an output vector of size 1
// EXPECT_CALL(*sys_out, output(Eq(x), Property(&Eigen::Ref<OutputVector>::size, Eq(1)))).WillRepeatedly(SetArgReferee<1>(x));

TEST_F(PidControllerTest, zero_error)
{
    //    StateVector x(1);
    //    x[0] = 1;
    //    OutputVector yref(1);
    //    yref[0] = 1;
    //    corbo::StaticReference ref(yref);

    //    // setup system output
    //    std::shared_ptr<corbo::FullStateSystemOutput> sys_out = std::make_shared<corbo::FullStateSystemOutput>();
    //    _controller.initialize(1, 1, sys_out);

    //    ControlVector u(1);
    //    u[0] = std::numeric_limits<double>::max();
    //    _controller.step(x, ref, corbo::Duration(0), corbo::TimePoint(), u);

    //    EXPECT_EQ(u.size(), 1);
    //    EXPECT_NEAR(u[0], 0, 1e-3);
}

int main(int argc, char** argv)
{
    testing::InitGoogleTest(&argc, argv);
    int ret = RUN_ALL_TESTS();

    // std::this_thread::sleep_for(std::chrono::seconds(1));

    //
    // double test = time.toSec();

    // PRINT_INFO("current time:" << test);

    corbo::PidController::Ptr pid = corbo::ControllerFactory::instance().create("PidController");
    if (pid)
    {
        PRINT_INFO("yeees");
    }
    else
    {
        PRINT_ERROR("ooooh noooo");
    }

    return ret;
    // return 1;
}

// int main()
//{
//  std::cout << "Hello World!" << std::endl;
//  return 0;
//}
