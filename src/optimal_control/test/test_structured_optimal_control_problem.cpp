/*********************************************************************
 *
 * Software License Agreement
 *
 *  Copyright (c) 2018,
 *  TU Dortmund - Institute of Control Theory and Systems Engineering.
 *  All rights reserved.
 *
 *  This software is currently not released.
 *  Redistribution and use in source and binary forms,
 *  with or without modification, are prohibited.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 *  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 *  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 *  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 *  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 *  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 *  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 *  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 *
 * Authors: Christoph Rösmann
 *********************************************************************/

#include <corbo-optimization/optimal_control/structured_optimal_control_problem.h>

#include <corbo-core/reference_trajectory.h>

#include <corbo-core/console.h>
#include <corbo-core/macros.h>
#include <corbo-core/value_comparison.h>

#include "gtest/gtest.h"

// using corbo::FullDiscretizationGrid;
// using corbo::StaticReference;
// using corbo::ZeroReference;
// using corbo::VectorVertex;

class StructuredOptimalControlProblem : public testing::Test
{
 protected:
    // You can do set-up work for each test here.
    StructuredOptimalControlProblem() {}
    // You can do clean-up work that doesn't throw exceptions here.
    virtual ~StructuredOptimalControlProblem() {}
    // If the constructor and destructor are not enough for setting up
    // and cleaning up each test, you can define the following methods:

    // Code here will be called immediately after the constructor (right
    // before each test).
    void SetUp() override {}
    // Code here will be called immediately after each test (right
    // before the destructor).
    // virtual void TearDown();

    // FullDiscretizationGrid grid;
};

// TEST_F(TestFullDiscretizationGrid, initialize_grid)
//{
//    Eigen::Vector2d x0(0, 1);
//    Eigen::VectorXd u_ref(1);
//    u_ref[0] = -1;

//    StaticReference xref(Eigen::Vector2d(6, 4));
//    StaticReference uref(u_ref);

//    corbo::IntegratorSystem::Ptr system = std::make_shared<corbo::IntegratorSystem>(2);

//    grid.setDt(0.1);
//    grid.setNumStates(4);

//    grid.initialize(x0, xref, uref, system);

//    // check state inputs
//    EXPECT_EQ(grid.getStateSequence().size(), 3);
//    EXPECT_DOUBLE_EQ(grid.getStateSequence()[0].values()[0], 0.0);
//    EXPECT_DOUBLE_EQ(grid.getStateSequence()[0].values()[1], 1.0);
//    EXPECT_DOUBLE_EQ(grid.getStateSequence()[1].values()[0], 2.0);
//    EXPECT_DOUBLE_EQ(grid.getStateSequence()[1].values()[1], 2.0);
//    EXPECT_DOUBLE_EQ(grid.getStateSequence()[2].values()[0], 4.0);
//    EXPECT_DOUBLE_EQ(grid.getStateSequence()[2].values()[1], 3.0);
//    EXPECT_DOUBLE_EQ(grid.getFinalState().values()[0], 6.0);
//    EXPECT_DOUBLE_EQ(grid.getFinalState().values()[1], 4.0);

//    EXPECT_EQ(grid.getStateSequence().front().getDimensionUnfixed(), 0);
//    EXPECT_EQ(grid.getFinalState().getDimensionUnfixed(), 2);

//    // check control inputs
//    EXPECT_EQ(grid.getControlSequence().size(), 3);
//    for (const VectorVertex& vertex : grid.getControlSequence())
//    {
//        EXPECT_DOUBLE_EQ(vertex.values()[0], -1);
//    }
//}
