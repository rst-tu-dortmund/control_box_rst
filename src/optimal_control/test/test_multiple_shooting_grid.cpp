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

#include <corbo-optimal-control/structured_ocp/discretization_grids/multiple_shooting_grid.h>

#include <corbo-core/reference_trajectory.h>
#include <corbo-systems/benchmark/linear_benchmark_systems.h>

#include <corbo-core/console.h>
#include <corbo-core/macros.h>
#include <corbo-core/value_comparison.h>

#include "gtest/gtest.h"

using corbo::MultipleShootingGrid;
using corbo::StaticReference;
using corbo::ZeroReference;
using corbo::VectorVertex;

class TestMultipleShootingGrid : public testing::Test
{
 protected:
    // You can do set-up work for each test here.
    TestMultipleShootingGrid() {}
    // You can do clean-up work that doesn't throw exceptions here.
    virtual ~TestMultipleShootingGrid() {}
    // If the constructor and destructor are not enough for setting up
    // and cleaning up each test, you can define the following methods:

    // Code here will be called immediately after the constructor (right
    // before each test).
    void SetUp() override {}
    // Code here will be called immediately after each test (right
    // before the destructor).
    // virtual void TearDown();

    MultipleShootingGrid grid;
};

// #ifndef NDEBUG
//        // check if we have only unique vertices
//        std::vector<VertexInterface*> vertices_copy = _active_vertices;
//        std::sort(vertices_copy.begin(), vertices_copy.end());
//        auto it = std::adjacent_find(vertices_copy.begin(), vertices_copy.end());
//        if (it != vertices_copy.end())
//        {
//            PRINT_WARNING("MultipleShootingGrid::computeActiveVertices(): returns non");
//        }

// #endif

// TEST_F(TestMultipleShootingGrid, initialize_grid_start_state_only)
//{
//    Eigen::Vector2d x0(0, 5);

//    Eigen::Vector2d u0(0, 0);

//    int num_states = 5;
//    double dt      = 0.1;

//    grid.setHorizon(num_states, dt);
//    grid.initialize(x0, x0, u0);

//    corbo::TimeSeries::Ptr sn = grid.getShootingNodeTimeSeries();
//    EXPECT_DOUBLE_EQ(sn->getValueDimension(), 2);
//    EXPECT_DOUBLE_EQ(sn->getTimeDimension(), 5);

//    Eigen::MatrixXd s_seq = sn->getValuesMatrixView();

//    EXPECT_DOUBLE_EQ(s_seq(0, 0), 0.0);
//    EXPECT_DOUBLE_EQ(s_seq(1, 0), 5.0);
//    EXPECT_DOUBLE_EQ(s_seq(0, 1), 0.0);
//    EXPECT_DOUBLE_EQ(s_seq(1, 1), 5.0);
//    EXPECT_DOUBLE_EQ(s_seq(0, 2), 0.0);
//    EXPECT_DOUBLE_EQ(s_seq(1, 2), 5.0);
//    EXPECT_DOUBLE_EQ(s_seq(0, 3), 0.0);
//    EXPECT_DOUBLE_EQ(s_seq(1, 3), 5.0);
//    EXPECT_DOUBLE_EQ(s_seq(0, 4), 0.0);
//    EXPECT_DOUBLE_EQ(s_seq(1, 4), 5.0);

//    corbo::TimeSeries::Ptr un = grid.getControlInputTimeSeries();
//    EXPECT_DOUBLE_EQ(un->getValueDimension(), 2);
//    EXPECT_DOUBLE_EQ(un->getTimeDimension(), 5);

//    Eigen::MatrixXd u_seq = un->getValuesMatrixView();
//    EXPECT_TRUE(u_seq.isZero()) << u_seq;

//    // check fixed flags (first shooting node must be fixed, others not)
//    EXPECT_TRUE(grid.getShootingIntervals()[0].shooting_node->isFixed());
//    for (int i = 1; i < (int)grid.getShootingIntervals().size(); ++i)
//    {
//        EXPECT_FALSE(grid.getShootingIntervals()[i].shooting_node->isFixed());
//    }
//}

// TEST_F(TestMultipleShootingGrid, initialize_grid_start_state_end_state)
//{
//    Eigen::Vector2d x0(0, 5);
//    Eigen::Vector2d xf(5, 0);

//    Eigen::Vector2d u0(0, 0);

//    int num_states = 5;
//    double dt      = 0.1;

//    grid.setHorizon(num_states, dt);
//    grid.initialize(x0, xf, u0);

//    corbo::TimeSeries::Ptr sn = grid.getShootingNodeTimeSeries();
//    EXPECT_DOUBLE_EQ(sn->getValueDimension(), 2);
//    EXPECT_DOUBLE_EQ(sn->getTimeDimension(), 5);

//    Eigen::MatrixXd s_seq = sn->getValuesMatrixView();

//    EXPECT_DOUBLE_EQ(s_seq(0, 0), 0.0);
//    EXPECT_DOUBLE_EQ(s_seq(1, 0), 5.0);
//    EXPECT_DOUBLE_EQ(s_seq(0, 1), 1.25);
//    EXPECT_DOUBLE_EQ(s_seq(1, 1), 3.75);
//    EXPECT_DOUBLE_EQ(s_seq(0, 2), 2.5);
//    EXPECT_DOUBLE_EQ(s_seq(1, 2), 2.5);
//    EXPECT_DOUBLE_EQ(s_seq(0, 3), 3.75);
//    EXPECT_DOUBLE_EQ(s_seq(1, 3), 1.25);
//    EXPECT_DOUBLE_EQ(s_seq(0, 4), 5.0);
//    EXPECT_DOUBLE_EQ(s_seq(1, 4), 0.0);

//    corbo::TimeSeries::Ptr un = grid.getControlInputTimeSeries();
//    EXPECT_DOUBLE_EQ(un->getValueDimension(), 2);
//    EXPECT_DOUBLE_EQ(un->getTimeDimension(), 5);

//    Eigen::MatrixXd u_seq = un->getValuesMatrixView();
//    EXPECT_TRUE(u_seq.isZero()) << u_seq;

//    // check fixed flags (first shooting node must be fixed, others not)
//    EXPECT_TRUE(grid.getShootingIntervals()[0].shooting_node->isFixed());
//    for (int i = 1; i < (int)grid.getShootingIntervals().size(); ++i)
//    {
//        EXPECT_FALSE(grid.getShootingIntervals()[i].shooting_node->isFixed());
//    }
//}

// TEST_F(TestMultipleShootingGrid, grid_prune_trajectory)
//{
//    Eigen::Vector2d x0(0, 5);
//    Eigen::Vector2d xf(5, 0);

//    Eigen::Vector2d u0(0, 0);

//    int num_states = 5;
//    double dt      = 0.1;

//    grid.setHorizon(num_states, dt);
//    grid.initialize(x0, xf, u0);

//    grid.pruneTrajectory(Eigen::Vector2d(2.4, 2.4));  // erase the first two intervals

//    corbo::TimeSeries::Ptr sn = grid.getShootingNodeTimeSeries();
//    EXPECT_DOUBLE_EQ(sn->getValueDimension(), 2);
//    EXPECT_DOUBLE_EQ(sn->getTimeDimension(), 3);

//    Eigen::MatrixXd s_seq = sn->getValuesMatrixView();

//    // first states must be updated to (2.4, 2.4) during pruning
//    EXPECT_DOUBLE_EQ(s_seq(0, 0), 2.4);
//    EXPECT_DOUBLE_EQ(s_seq(1, 0), 2.4);
//    EXPECT_DOUBLE_EQ(s_seq(0, 1), 3.75);
//    EXPECT_DOUBLE_EQ(s_seq(1, 1), 1.25);
//    EXPECT_DOUBLE_EQ(s_seq(0, 2), 5.0);
//    EXPECT_DOUBLE_EQ(s_seq(1, 2), 0.0);

//    corbo::TimeSeries::Ptr un = grid.getControlInputTimeSeries();
//    EXPECT_DOUBLE_EQ(un->getValueDimension(), 2);
//    EXPECT_DOUBLE_EQ(un->getTimeDimension(), 3);

//    Eigen::MatrixXd u_seq = un->getValuesMatrixView();
//    EXPECT_TRUE(u_seq.isZero()) << u_seq;

//    // check fixed flags (first shooting node must be fixed, others not)
//    EXPECT_TRUE(grid.getShootingIntervals()[0].shooting_node->isFixed());
//    for (int i = 1; i < (int)grid.getShootingIntervals().size(); ++i)
//    {
//        EXPECT_FALSE(grid.getShootingIntervals()[i].shooting_node->isFixed());
//    }
//}

// TEST_F(TestMultipleShootingGrid, grid_prune_trajectory_keep_min)
//{
//    Eigen::Vector2d x0(0, 5);
//    Eigen::Vector2d xf(5, 0);

//    Eigen::Vector2d u0(0, 0);

//    int num_states = 5;
//    double dt      = 0.1;

//    grid.setHorizon(num_states, dt);
//    grid.initialize(x0, xf, u0);

//    // now specify final state as new initial state
//    grid.pruneTrajectory(xf);

//    corbo::TimeSeries::Ptr sn = grid.getShootingNodeTimeSeries();
//    EXPECT_DOUBLE_EQ(sn->getValueDimension(), 2);
//    EXPECT_DOUBLE_EQ(sn->getTimeDimension(), 2);  // we always keep 2 samples (1 interval) in our OCP

//    Eigen::MatrixXd s_seq = sn->getValuesMatrixView();

//    // first states must be updated to (2.4, 2.4) during pruning
//    EXPECT_DOUBLE_EQ(s_seq(0, 0), 5.0);
//    EXPECT_DOUBLE_EQ(s_seq(1, 0), 0.0);
//    EXPECT_DOUBLE_EQ(s_seq(0, 1), 5.0);
//    EXPECT_DOUBLE_EQ(s_seq(1, 1), 0.0);

//    corbo::TimeSeries::Ptr un = grid.getControlInputTimeSeries();
//    EXPECT_DOUBLE_EQ(un->getValueDimension(), 2);
//    EXPECT_DOUBLE_EQ(un->getTimeDimension(), 2);  // we always keep 2 samples (1 interval) in our OCP

//    Eigen::MatrixXd u_seq = un->getValuesMatrixView();
//    EXPECT_TRUE(u_seq.isZero()) << u_seq;

//    // check fixed flags (first shooting node must be fixed, others not)
//    EXPECT_TRUE(grid.getShootingIntervals()[0].shooting_node->isFixed());
//    for (int i = 1; i < (int)grid.getShootingIntervals().size(); ++i)
//    {
//        EXPECT_FALSE(grid.getShootingIntervals()[i].shooting_node->isFixed());
//    }
//}

// TEST_F(TestMultipleShootingGrid, grid_extrapolate_trajectory)
//{
//    Eigen::Vector2d x0(0, 5);
//    Eigen::Vector2d xf(5, 0);

//    Eigen::Vector2d u0(0, 0);

//    int num_states = 5;
//    double dt      = 0.1;

//    grid.setHorizon(num_states, dt);
//    grid.initialize(x0, xf, u0);

//    grid.extrapolateTrajectory(2);

//    corbo::TimeSeries::Ptr sn = grid.getShootingNodeTimeSeries();
//    EXPECT_DOUBLE_EQ(sn->getValueDimension(), 2);
//    EXPECT_DOUBLE_EQ(sn->getTimeDimension(), 7);

//    Eigen::MatrixXd s_seq = sn->getValuesMatrixView();

//    EXPECT_DOUBLE_EQ(s_seq(0, 0), 0.0);
//    EXPECT_DOUBLE_EQ(s_seq(1, 0), 5.0);
//    EXPECT_DOUBLE_EQ(s_seq(0, 1), 1.25);
//    EXPECT_DOUBLE_EQ(s_seq(1, 1), 3.75);
//    EXPECT_DOUBLE_EQ(s_seq(0, 2), 2.5);
//    EXPECT_DOUBLE_EQ(s_seq(1, 2), 2.5);
//    EXPECT_DOUBLE_EQ(s_seq(0, 3), 3.75);
//    EXPECT_DOUBLE_EQ(s_seq(1, 3), 1.25);
//    EXPECT_DOUBLE_EQ(s_seq(0, 4), 5.0);
//    EXPECT_DOUBLE_EQ(s_seq(1, 4), 0.0);
//    EXPECT_DOUBLE_EQ(s_seq(0, 5), 6.25);
//    EXPECT_DOUBLE_EQ(s_seq(1, 5), -1.25);
//    EXPECT_DOUBLE_EQ(s_seq(0, 6), 7.5);
//    EXPECT_DOUBLE_EQ(s_seq(1, 6), -2.5);

//    corbo::TimeSeries::Ptr un = grid.getControlInputTimeSeries();
//    EXPECT_DOUBLE_EQ(un->getValueDimension(), 2);
//    EXPECT_DOUBLE_EQ(un->getTimeDimension(), 7);

//    Eigen::MatrixXd u_seq = un->getValuesMatrixView();
//    EXPECT_TRUE(u_seq.isZero()) << u_seq;

//    // check fixed flags (first shooting node must be fixed, others not)
//    EXPECT_TRUE(grid.getShootingIntervals()[0].shooting_node->isFixed());
//    for (int i = 1; i < (int)grid.getShootingIntervals().size(); ++i)
//    {
//        EXPECT_FALSE(grid.getShootingIntervals()[i].shooting_node->isFixed());
//    }
//}
