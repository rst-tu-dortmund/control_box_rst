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

#include <corbo-optimal-control/structured_ocp/discretization_grids/full_discretization_grid.h>

#include <corbo-core/reference_trajectory.h>

#include <corbo-core/console.h>
#include <corbo-core/macros.h>
#include <corbo-core/value_comparison.h>

#include "gtest/gtest.h"

using corbo::FullDiscretizationGrid;
using corbo::StaticReference;
using corbo::VectorVertex;
using corbo::ZeroReference;

class TestFullDiscretizationGrid : public testing::Test
{
 protected:
    // You can do set-up work for each test here.
    TestFullDiscretizationGrid() {}
    // You can do clean-up work that doesn't throw exceptions here.
    virtual ~TestFullDiscretizationGrid() {}
    // If the constructor and destructor are not enough for setting up
    // and cleaning up each test, you can define the following methods:

    // Code here will be called immediately after the constructor (right
    // before each test).
    void SetUp() override {}
    // Code here will be called immediately after each test (right
    // before the destructor).
    // virtual void TearDown();

    FullDiscretizationGrid grid;
};

TEST_F(TestFullDiscretizationGrid, initialize_grid_start_state_only)
{
    Eigen::Vector2d x0(0, 5);

    Eigen::Vector2d u0(0, 0);

    int num_states = 5;
    double dt      = 0.1;

    grid.setHorizon(num_states, dt);
    grid.initialize(x0, x0, u0);

    corbo::TimeSeries::Ptr sn = std::make_shared<corbo::TimeSeries>();
    grid.getShootingNodeTimeSeries(sn);
    EXPECT_DOUBLE_EQ(sn->getValueDimension(), 2);
    EXPECT_DOUBLE_EQ(sn->getTimeDimension(), 5);

    Eigen::MatrixXd s_seq = sn->getValuesMatrixView();

    EXPECT_DOUBLE_EQ(s_seq(0, 0), 0.0);
    EXPECT_DOUBLE_EQ(s_seq(1, 0), 5.0);
    EXPECT_DOUBLE_EQ(s_seq(0, 1), 0.0);
    EXPECT_DOUBLE_EQ(s_seq(1, 1), 5.0);
    EXPECT_DOUBLE_EQ(s_seq(0, 2), 0.0);
    EXPECT_DOUBLE_EQ(s_seq(1, 2), 5.0);
    EXPECT_DOUBLE_EQ(s_seq(0, 3), 0.0);
    EXPECT_DOUBLE_EQ(s_seq(1, 3), 5.0);
    EXPECT_DOUBLE_EQ(s_seq(0, 4), 0.0);
    EXPECT_DOUBLE_EQ(s_seq(1, 4), 5.0);

    corbo::TimeSeries::Ptr un = std::make_shared<corbo::TimeSeries>();
    grid.getControlInputTimeSeries(un);
    EXPECT_DOUBLE_EQ(un->getValueDimension(), 2);
    EXPECT_DOUBLE_EQ(un->getTimeDimension(), 5);

    Eigen::MatrixXd u_seq = un->getValuesMatrixView();
    EXPECT_TRUE(u_seq.isZero()) << u_seq;

    // check fixed flags (first shooting node must be fixed, others not)
    EXPECT_TRUE(grid.getShootingIntervals()[0].shooting_node->isFixed());
    for (int i = 1; i < (int)grid.getShootingIntervals().size(); ++i)
    {
        EXPECT_FALSE(grid.getShootingIntervals()[i].shooting_node->isFixed());
    }
}

TEST_F(TestFullDiscretizationGrid, initialize_grid_start_state_end_state)
{
    Eigen::Vector2d x0(0, 5);
    Eigen::Vector2d xf(5, 0);

    Eigen::Vector2d u0(0, 0);

    int num_states = 5;
    double dt      = 0.1;

    grid.setHorizon(num_states, dt);
    grid.initialize(x0, xf, u0);

    corbo::TimeSeries::Ptr sn = std::make_shared<corbo::TimeSeries>();
    grid.getShootingNodeTimeSeries(sn);
    EXPECT_DOUBLE_EQ(sn->getValueDimension(), 2);
    EXPECT_DOUBLE_EQ(sn->getTimeDimension(), 5);

    Eigen::MatrixXd s_seq = sn->getValuesMatrixView();

    EXPECT_DOUBLE_EQ(s_seq(0, 0), 0.0);
    EXPECT_DOUBLE_EQ(s_seq(1, 0), 5.0);
    EXPECT_DOUBLE_EQ(s_seq(0, 1), 1.25);
    EXPECT_DOUBLE_EQ(s_seq(1, 1), 3.75);
    EXPECT_DOUBLE_EQ(s_seq(0, 2), 2.5);
    EXPECT_DOUBLE_EQ(s_seq(1, 2), 2.5);
    EXPECT_DOUBLE_EQ(s_seq(0, 3), 3.75);
    EXPECT_DOUBLE_EQ(s_seq(1, 3), 1.25);
    EXPECT_DOUBLE_EQ(s_seq(0, 4), 5.0);
    EXPECT_DOUBLE_EQ(s_seq(1, 4), 0.0);

    corbo::TimeSeries::Ptr un = std::make_shared<corbo::TimeSeries>();
    grid.getControlInputTimeSeries(un);
    EXPECT_DOUBLE_EQ(un->getValueDimension(), 2);
    EXPECT_DOUBLE_EQ(un->getTimeDimension(), 5);

    Eigen::MatrixXd u_seq = un->getValuesMatrixView();
    EXPECT_TRUE(u_seq.isZero()) << u_seq;

    // check fixed flags (first shooting node must be fixed, others not)
    EXPECT_TRUE(grid.getShootingIntervals()[0].shooting_node->isFixed());
    for (int i = 1; i < (int)grid.getShootingIntervals().size(); ++i)
    {
        EXPECT_FALSE(grid.getShootingIntervals()[i].shooting_node->isFixed());
    }
}

TEST_F(TestFullDiscretizationGrid, grid_prune_trajectory)
{
    Eigen::Vector2d x0(0, 5);
    Eigen::Vector2d xf(5, 0);

    Eigen::Vector2d u0(0, 0);

    int num_states = 5;
    double dt      = 0.1;

    grid.setHorizon(num_states, dt);
    grid.initialize(x0, xf, u0);

    grid.pruneTrajectory(Eigen::Vector2d(2.4, 2.4));  // erase the first two intervals

    corbo::TimeSeries::Ptr sn = std::make_shared<corbo::TimeSeries>();
    grid.getShootingNodeTimeSeries(sn);
    EXPECT_DOUBLE_EQ(sn->getValueDimension(), 2);
    EXPECT_DOUBLE_EQ(sn->getTimeDimension(), 3);

    Eigen::MatrixXd s_seq = sn->getValuesMatrixView();

    // first states must be updated to (2.4, 2.4) during pruning
    EXPECT_DOUBLE_EQ(s_seq(0, 0), 2.4);
    EXPECT_DOUBLE_EQ(s_seq(1, 0), 2.4);
    EXPECT_DOUBLE_EQ(s_seq(0, 1), 3.75);
    EXPECT_DOUBLE_EQ(s_seq(1, 1), 1.25);
    EXPECT_DOUBLE_EQ(s_seq(0, 2), 5.0);
    EXPECT_DOUBLE_EQ(s_seq(1, 2), 0.0);

    corbo::TimeSeries::Ptr un = std::make_shared<corbo::TimeSeries>();
    grid.getControlInputTimeSeries(un);
    EXPECT_DOUBLE_EQ(un->getValueDimension(), 2);
    EXPECT_DOUBLE_EQ(un->getTimeDimension(), 3);

    Eigen::MatrixXd u_seq = un->getValuesMatrixView();
    EXPECT_TRUE(u_seq.isZero()) << u_seq;

    // check fixed flags (first shooting node must be fixed, others not)
    EXPECT_TRUE(grid.getShootingIntervals()[0].shooting_node->isFixed());
    for (int i = 1; i < (int)grid.getShootingIntervals().size(); ++i)
    {
        EXPECT_FALSE(grid.getShootingIntervals()[i].shooting_node->isFixed());
    }
}

TEST_F(TestFullDiscretizationGrid, grid_prune_trajectory_keep_min)
{
    Eigen::Vector2d x0(0, 5);
    Eigen::Vector2d xf(5, 0);

    Eigen::Vector2d u0(0, 0);

    int num_states = 5;
    double dt      = 0.1;

    grid.setHorizon(num_states, dt);
    grid.initialize(x0, xf, u0);

    // now specify final state as new initial state
    grid.pruneTrajectory(xf);

    corbo::TimeSeries::Ptr sn = std::make_shared<corbo::TimeSeries>();
    grid.getShootingNodeTimeSeries(sn);
    EXPECT_DOUBLE_EQ(sn->getValueDimension(), 2);
    EXPECT_DOUBLE_EQ(sn->getTimeDimension(), 2);  // we always keep 2 samples (1 interval) in our OCP

    Eigen::MatrixXd s_seq = sn->getValuesMatrixView();

    // first states must be updated to (2.4, 2.4) during pruning
    EXPECT_DOUBLE_EQ(s_seq(0, 0), 5.0);
    EXPECT_DOUBLE_EQ(s_seq(1, 0), 0.0);
    EXPECT_DOUBLE_EQ(s_seq(0, 1), 5.0);
    EXPECT_DOUBLE_EQ(s_seq(1, 1), 0.0);

    corbo::TimeSeries::Ptr un = std::make_shared<corbo::TimeSeries>();
    grid.getControlInputTimeSeries(un);
    EXPECT_DOUBLE_EQ(un->getValueDimension(), 2);
    EXPECT_DOUBLE_EQ(un->getTimeDimension(), 2);  // we always keep 2 samples (1 interval) in our OCP

    Eigen::MatrixXd u_seq = un->getValuesMatrixView();
    EXPECT_TRUE(u_seq.isZero()) << u_seq;

    // check fixed flags (first shooting node must be fixed, others not)
    EXPECT_TRUE(grid.getShootingIntervals()[0].shooting_node->isFixed());
    for (int i = 1; i < (int)grid.getShootingIntervals().size(); ++i)
    {
        EXPECT_FALSE(grid.getShootingIntervals()[i].shooting_node->isFixed());
    }
}

TEST_F(TestFullDiscretizationGrid, grid_extrapolate_trajectory)
{
    Eigen::Vector2d x0(0, 5);
    Eigen::Vector2d xf(5, 0);

    Eigen::Vector2d u0(0, 0);

    int num_states = 5;
    double dt      = 0.1;

    grid.setHorizon(num_states, dt);
    grid.initialize(x0, xf, u0);

    grid.extrapolateTrajectory(2);

    corbo::TimeSeries::Ptr sn = std::make_shared<corbo::TimeSeries>();
    grid.getShootingNodeTimeSeries(sn);
    EXPECT_DOUBLE_EQ(sn->getValueDimension(), 2);
    EXPECT_DOUBLE_EQ(sn->getTimeDimension(), 7);

    Eigen::MatrixXd s_seq = sn->getValuesMatrixView();

    EXPECT_DOUBLE_EQ(s_seq(0, 0), 0.0);
    EXPECT_DOUBLE_EQ(s_seq(1, 0), 5.0);
    EXPECT_DOUBLE_EQ(s_seq(0, 1), 1.25);
    EXPECT_DOUBLE_EQ(s_seq(1, 1), 3.75);
    EXPECT_DOUBLE_EQ(s_seq(0, 2), 2.5);
    EXPECT_DOUBLE_EQ(s_seq(1, 2), 2.5);
    EXPECT_DOUBLE_EQ(s_seq(0, 3), 3.75);
    EXPECT_DOUBLE_EQ(s_seq(1, 3), 1.25);
    EXPECT_DOUBLE_EQ(s_seq(0, 4), 5.0);
    EXPECT_DOUBLE_EQ(s_seq(1, 4), 0.0);
    EXPECT_DOUBLE_EQ(s_seq(0, 5), 6.25);
    EXPECT_DOUBLE_EQ(s_seq(1, 5), -1.25);
    EXPECT_DOUBLE_EQ(s_seq(0, 6), 7.5);
    EXPECT_DOUBLE_EQ(s_seq(1, 6), -2.5);

    corbo::TimeSeries::Ptr un = std::make_shared<corbo::TimeSeries>();
    grid.getControlInputTimeSeries(un);
    EXPECT_DOUBLE_EQ(un->getValueDimension(), 2);
    EXPECT_DOUBLE_EQ(un->getTimeDimension(), 7);

    Eigen::MatrixXd u_seq = un->getValuesMatrixView();
    EXPECT_TRUE(u_seq.isZero()) << u_seq;

    // check fixed flags (first shooting node must be fixed, others not)
    EXPECT_TRUE(grid.getShootingIntervals()[0].shooting_node->isFixed());
    for (int i = 1; i < (int)grid.getShootingIntervals().size(); ++i)
    {
        EXPECT_FALSE(grid.getShootingIntervals()[i].shooting_node->isFixed());
    }
}

TEST_F(TestFullDiscretizationGrid, resample_trajectory_single_dt)
{
    Eigen::Vector2d x0(0, 0);
    Eigen::Vector2d xf(1, -1);
    Eigen::Vector2d u0(0, 0);

    int num_states = 10;
    double dt      = 0.1;

    grid.setSingleDt(true);

    grid.setHorizon(num_states, dt);
    grid.initialize(x0, xf, u0);

    // Now resample to 5 states
    grid.resampleTrajectory(5);

    // states
    corbo::TimeSeries::Ptr sn = std::make_shared<corbo::TimeSeries>();
    grid.getShootingNodeTimeSeries(sn);
    EXPECT_DOUBLE_EQ(sn->getValueDimension(), 2);
    EXPECT_DOUBLE_EQ(sn->getTimeDimension(), 5);
    Eigen::MatrixXd s_seq = sn->getValuesMatrixView();

    // Check values
    for (int i = 0; i < s_seq.cols(); ++i)
    {
        EXPECT_DOUBLE_EQ(x0[0] + (double)i * (xf[0] - x0[0]) / 4.0, s_seq(0, i));
        EXPECT_DOUBLE_EQ(x0[1] + (double)i * (xf[1] - x0[1]) / 4.0, s_seq(1, i));
    }

    // controls
    corbo::TimeSeries::Ptr un = std::make_shared<corbo::TimeSeries>();
    grid.getControlInputTimeSeries(un);
    EXPECT_DOUBLE_EQ(un->getValueDimension(), 2);
    EXPECT_DOUBLE_EQ(un->getTimeDimension(), 5);

    Eigen::MatrixXd u_seq = un->getValuesMatrixView();
    EXPECT_TRUE(u_seq.isZero()) << u_seq;

    // dt
    EXPECT_DOUBLE_EQ(4.0 * grid.getFirstDt(), 0.9) << "dt: " << grid.getFirstDt();

    // check fixed flags (first shooting node must be fixed, others not)
    EXPECT_TRUE(grid.getShootingIntervals()[0].shooting_node->isFixed());
    for (int i = 1; i < (int)grid.getShootingIntervals().size(); ++i)
    {
        EXPECT_FALSE(grid.getShootingIntervals()[i].shooting_node->isFixed());
    }

    // Now resample to 15 states
    grid.resampleTrajectory(15);

    grid.getShootingNodeTimeSeries(sn);
    EXPECT_DOUBLE_EQ(sn->getValueDimension(), 2);
    EXPECT_DOUBLE_EQ(sn->getTimeDimension(), 15);
    s_seq = sn->getValuesMatrixView();

    // Check values
    for (int i = 0; i < s_seq.cols(); ++i)
    {
        EXPECT_DOUBLE_EQ(x0[0] + (double)i * (xf[0] - x0[0]) / 14.0, s_seq(0, i));
        EXPECT_DOUBLE_EQ(x0[1] + (double)i * (xf[1] - x0[1]) / 14.0, s_seq(1, i));
    }

    // controls
    grid.getControlInputTimeSeries(un);
    EXPECT_DOUBLE_EQ(un->getValueDimension(), 2);
    EXPECT_DOUBLE_EQ(un->getTimeDimension(), 15);

    u_seq = un->getValuesMatrixView();
    EXPECT_TRUE(u_seq.isZero()) << u_seq;

    // dt
    EXPECT_DOUBLE_EQ(14.0 * grid.getFirstDt(), 0.9) << "dt: " << grid.getFirstDt();

    // check fixed flags (first shooting node must be fixed, others not)
    EXPECT_TRUE(grid.getShootingIntervals()[0].shooting_node->isFixed());
    for (int i = 1; i < (int)grid.getShootingIntervals().size(); ++i)
    {
        EXPECT_FALSE(grid.getShootingIntervals()[i].shooting_node->isFixed());
    }
}

TEST_F(TestFullDiscretizationGrid, resample_trajectory_single_dt2)
{
    Eigen::Vector2d x0(0, 0);
    Eigen::Vector2d xf(1, -1);
    Eigen::Vector2d u0(0, 0);

    int num_states = 10;
    double dt      = 0.1;

    grid.setSingleDt(true);

    grid.setHorizon(num_states, dt);
    grid.initialize(x0, xf, u0);

    // Now resample to 11 states
    grid.resampleTrajectory(11);

    // states
    corbo::TimeSeries::Ptr sn = std::make_shared<corbo::TimeSeries>();
    grid.getShootingNodeTimeSeries(sn);
    EXPECT_DOUBLE_EQ(sn->getValueDimension(), 2);
    EXPECT_DOUBLE_EQ(sn->getTimeDimension(), 11);
    Eigen::MatrixXd s_seq = sn->getValuesMatrixView();

    // Check values
    for (int i = 0; i < s_seq.cols(); ++i)
    {
        EXPECT_DOUBLE_EQ(x0[0] + (double)i * (xf[0] - x0[0]) / 10.0, s_seq(0, i));
        EXPECT_DOUBLE_EQ(x0[1] + (double)i * (xf[1] - x0[1]) / 10.0, s_seq(1, i));
    }

    // controls
    corbo::TimeSeries::Ptr un = std::make_shared<corbo::TimeSeries>();
    grid.getControlInputTimeSeries(un);
    EXPECT_DOUBLE_EQ(un->getValueDimension(), 2);
    EXPECT_DOUBLE_EQ(un->getTimeDimension(), 11);

    Eigen::MatrixXd u_seq = un->getValuesMatrixView();
    EXPECT_TRUE(u_seq.isZero()) << u_seq;

    // dt
    EXPECT_DOUBLE_EQ(10.0 * grid.getFirstDt(), 0.9) << "dt: " << grid.getFirstDt();

    // check fixed flags (first shooting node must be fixed, others not)
    EXPECT_TRUE(grid.getShootingIntervals()[0].shooting_node->isFixed());
    for (int i = 1; i < (int)grid.getShootingIntervals().size(); ++i)
    {
        EXPECT_FALSE(grid.getShootingIntervals()[i].shooting_node->isFixed());
    }

    // Now resample to 10 states
    grid.resampleTrajectory(10);

    grid.getShootingNodeTimeSeries(sn);
    EXPECT_DOUBLE_EQ(sn->getValueDimension(), 2);
    EXPECT_DOUBLE_EQ(sn->getTimeDimension(), 10);
    s_seq = sn->getValuesMatrixView();

    // Check values
    for (int i = 0; i < s_seq.cols(); ++i)
    {
        EXPECT_DOUBLE_EQ(x0[0] + (double)i * (xf[0] - x0[0]) / 9.0, s_seq(0, i));
        EXPECT_DOUBLE_EQ(x0[1] + (double)i * (xf[1] - x0[1]) / 9.0, s_seq(1, i));
    }

    // controls
    grid.getControlInputTimeSeries(un);
    EXPECT_DOUBLE_EQ(un->getValueDimension(), 2);
    EXPECT_DOUBLE_EQ(un->getTimeDimension(), 10);

    u_seq = un->getValuesMatrixView();
    EXPECT_TRUE(u_seq.isZero()) << u_seq;

    // dt
    EXPECT_DOUBLE_EQ(9.0 * grid.getFirstDt(), 0.9) << "dt: " << grid.getFirstDt();

    // check fixed flags (first shooting node must be fixed, others not)
    EXPECT_TRUE(grid.getShootingIntervals()[0].shooting_node->isFixed());
    for (int i = 1; i < (int)grid.getShootingIntervals().size(); ++i)
    {
        EXPECT_FALSE(grid.getShootingIntervals()[i].shooting_node->isFixed());
    }
}

TEST_F(TestFullDiscretizationGrid, grid_shift_trajectory1)
{
    Eigen::Vector2d x0(3, 0);
    Eigen::Vector2d xf(0, 3);

    Eigen::Vector2d u0(0, 0);

    int num_states = 4;
    double dt      = 0.1;

    grid.setHorizon(num_states, dt);
    grid.initialize(x0, xf, u0);

    // set some values for u
    grid.getShootingIntervalsRaw()[0].controls.front()->values()[0] = 1.5;
    grid.getShootingIntervalsRaw()[1].controls.front()->values()[0] = 2.5;
    grid.getShootingIntervalsRaw()[2].controls.front()->values()[0] = 3.5;

    // now shift grid by 1
    EXPECT_TRUE(grid.shiftGridValues(1));

    // get grid values
    corbo::TimeSeries::Ptr sn = std::make_shared<corbo::TimeSeries>();
    grid.getShootingNodeTimeSeries(sn);
    EXPECT_DOUBLE_EQ(sn->getValueDimension(), 2);
    EXPECT_DOUBLE_EQ(sn->getTimeDimension(), 4);

    Eigen::MatrixXd s_seq = sn->getValuesMatrixView();

    EXPECT_DOUBLE_EQ(s_seq(0, 0), 2);
    EXPECT_DOUBLE_EQ(s_seq(1, 0), 1);
    EXPECT_DOUBLE_EQ(s_seq(0, 1), 1);
    EXPECT_DOUBLE_EQ(s_seq(1, 1), 2);
    EXPECT_DOUBLE_EQ(s_seq(0, 2), 0);
    EXPECT_DOUBLE_EQ(s_seq(1, 2), 3);
    EXPECT_DOUBLE_EQ(s_seq(0, 3), -1);
    EXPECT_DOUBLE_EQ(s_seq(1, 3), 4);

    // controls
    corbo::TimeSeries::Ptr un = std::make_shared<corbo::TimeSeries>();
    grid.getControlInputTimeSeries(un);
    EXPECT_DOUBLE_EQ(un->getValueDimension(), 2);
    EXPECT_DOUBLE_EQ(un->getTimeDimension(), 4);

    Eigen::MatrixXd u_seq = un->getValuesMatrixView();
    EXPECT_DOUBLE_EQ(u_seq(0, 0), 2.5);
    EXPECT_DOUBLE_EQ(u_seq(1, 0), 0);
    EXPECT_DOUBLE_EQ(u_seq(0, 1), 3.5);
    EXPECT_DOUBLE_EQ(u_seq(1, 1), 0);
    EXPECT_DOUBLE_EQ(u_seq(0, 2), 3.5);  // zero-order hold (extrapolation)
    EXPECT_DOUBLE_EQ(u_seq(1, 2), 0);
    // last control is duplicated in getControlInputTimeSeries()...
}

TEST_F(TestFullDiscretizationGrid, grid_shift_trajectory2)
{
    Eigen::Vector2d x0(3, 0);
    Eigen::Vector2d xf(0, 3);

    Eigen::Vector2d u0(0, 0);

    int num_states = 4;
    double dt      = 0.1;

    grid.setHorizon(num_states, dt);
    grid.initialize(x0, xf, u0);

    // set some values for u
    grid.getShootingIntervalsRaw()[0].controls.front()->values()[0] = 1.5;
    grid.getShootingIntervalsRaw()[1].controls.front()->values()[0] = 2.5;
    grid.getShootingIntervalsRaw()[2].controls.front()->values()[0] = 3.5;

    // now shift grid by 1
    EXPECT_TRUE(grid.shiftGridValues(2));

    // get grid values
    corbo::TimeSeries::Ptr sn = std::make_shared<corbo::TimeSeries>();
    grid.getShootingNodeTimeSeries(sn);
    EXPECT_DOUBLE_EQ(sn->getValueDimension(), 2);
    EXPECT_DOUBLE_EQ(sn->getTimeDimension(), 4);

    Eigen::MatrixXd s_seq = sn->getValuesMatrixView();

    EXPECT_DOUBLE_EQ(s_seq(0, 0), 1);
    EXPECT_DOUBLE_EQ(s_seq(1, 0), 2);
    EXPECT_DOUBLE_EQ(s_seq(0, 1), 0);
    EXPECT_DOUBLE_EQ(s_seq(1, 1), 3);
    EXPECT_DOUBLE_EQ(s_seq(0, 2), -1);
    EXPECT_DOUBLE_EQ(s_seq(1, 2), 4);
    EXPECT_DOUBLE_EQ(s_seq(0, 3), -2);
    EXPECT_DOUBLE_EQ(s_seq(1, 3), 5);

    // controls
    corbo::TimeSeries::Ptr un = std::make_shared<corbo::TimeSeries>();
    grid.getControlInputTimeSeries(un);
    EXPECT_DOUBLE_EQ(un->getValueDimension(), 2);
    EXPECT_DOUBLE_EQ(un->getTimeDimension(), 4);

    Eigen::MatrixXd u_seq = un->getValuesMatrixView();
    EXPECT_DOUBLE_EQ(u_seq(0, 0), 3.5);
    EXPECT_DOUBLE_EQ(u_seq(1, 0), 0);
    EXPECT_DOUBLE_EQ(u_seq(0, 1), 3.5);
    EXPECT_DOUBLE_EQ(u_seq(1, 1), 0);
    EXPECT_DOUBLE_EQ(u_seq(0, 2), 3.5);  // zero-order hold (extrapolation)
    EXPECT_DOUBLE_EQ(u_seq(1, 2), 0);
    // last control is duplicated in getControlInputTimeSeries()...
}

TEST_F(TestFullDiscretizationGrid, grid_shift_trajectory_x0)
{
    Eigen::Vector2d x0(3, 0);
    Eigen::Vector2d xf(0, 3);

    Eigen::Vector2d u0(0, 0);

    int num_states = 4;
    double dt      = 0.1;

    grid.setHorizon(num_states, dt);
    grid.initialize(x0, xf, u0);

    // set some values for u
    grid.getShootingIntervalsRaw()[0].controls.front()->values()[0] = 1.5;
    grid.getShootingIntervalsRaw()[1].controls.front()->values()[0] = 2.5;
    grid.getShootingIntervalsRaw()[2].controls.front()->values()[0] = 3.5;

    // now shift grid by 1
    EXPECT_TRUE(grid.shiftGridValues(Eigen::Vector2d(2.1, 0.9)));

    // get grid values
    corbo::TimeSeries::Ptr sn = std::make_shared<corbo::TimeSeries>();
    grid.getShootingNodeTimeSeries(sn);
    EXPECT_DOUBLE_EQ(sn->getValueDimension(), 2);
    EXPECT_DOUBLE_EQ(sn->getTimeDimension(), 4);

    Eigen::MatrixXd s_seq = sn->getValuesMatrixView();

    EXPECT_DOUBLE_EQ(s_seq(0, 0), 2);
    EXPECT_DOUBLE_EQ(s_seq(1, 0), 1);
    EXPECT_DOUBLE_EQ(s_seq(0, 1), 1);
    EXPECT_DOUBLE_EQ(s_seq(1, 1), 2);
    EXPECT_DOUBLE_EQ(s_seq(0, 2), 0);
    EXPECT_DOUBLE_EQ(s_seq(1, 2), 3);
    EXPECT_DOUBLE_EQ(s_seq(0, 3), -1);
    EXPECT_DOUBLE_EQ(s_seq(1, 3), 4);

    // controls
    corbo::TimeSeries::Ptr un = std::make_shared<corbo::TimeSeries>();
    grid.getControlInputTimeSeries(un);
    EXPECT_DOUBLE_EQ(un->getValueDimension(), 2);
    EXPECT_DOUBLE_EQ(un->getTimeDimension(), 4);

    Eigen::MatrixXd u_seq = un->getValuesMatrixView();
    EXPECT_DOUBLE_EQ(u_seq(0, 0), 2.5);
    EXPECT_DOUBLE_EQ(u_seq(1, 0), 0);
    EXPECT_DOUBLE_EQ(u_seq(0, 1), 3.5);
    EXPECT_DOUBLE_EQ(u_seq(1, 1), 0);
    EXPECT_DOUBLE_EQ(u_seq(0, 2), 3.5);  // zero-order hold (extrapolation)
    EXPECT_DOUBLE_EQ(u_seq(1, 2), 0);
    // last control is duplicated in getControlInputTimeSeries()...
}

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
