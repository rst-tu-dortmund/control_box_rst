/*********************************************************************
 *
 * Software License Agreement
 *
 *  Copyright (c) 2017,
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

#include <corbo-optimization/solver/levenberg_marquardt_sparse.h>

#include <corbo-optimization/standard_optimization_problem.h>

#include <corbo-core/console.h>
#include <corbo-core/macros.h>
#include <corbo-core/utilities.h>
#include <corbo-core/value_comparison.h>

#include <array>
#include <functional>

#include "gtest/gtest.h"

using corbo::LevenbergMarquardtSparse;
using corbo::StandardOptimizationProblemWithCallbacks;

class TestLevenbergMarquardtSparse : public testing::Test
{
 protected:
    // You can do set-up work for each test here.
    TestLevenbergMarquardtSparse() {}
    // You can do clean-up work that doesn't throw exceptions here.
    virtual ~TestLevenbergMarquardtSparse() {}
    // If the constructor and destructor are not enough for setting up
    // and cleaning up each test, you can define the following methods:

    // Code here will be called immediately after the constructor (right
    // before each test).
    void SetUp() override
    {
        // configure solver
        solver.setIterations(100);
    }
    // Code here will be called immediately after each test (right
    // before the destructor).
    // virtual void TearDown()

    StandardOptimizationProblemWithCallbacks optim;

    LevenbergMarquardtSparse solver;
};

TEST_F(TestLevenbergMarquardtSparse, solve_unconstr_1)
{
    // parameters
    optim.resizeParameterVector(1);
    Eigen::VectorXd x(1);
    x.setOnes();
    optim.setX(x);

    // create objective function
    // (x-2)^2
    auto objective = [](const Eigen::VectorXd& x, Eigen::Ref<Eigen::VectorXd> values) { values[0] = x[0] - 2; };
    optim.setObjectiveFunction(objective, 1, true);  // least-squares type

    EXPECT_TRUE(solver.initialize(optim));

    bool success = solver.solve(optim, true);

    EXPECT_TRUE(success);
    EXPECT_NEAR(optim.getX()[0], 2.0, 1e-6);
}

TEST_F(TestLevenbergMarquardtSparse, solve_unconstr_2)
{
    // parameters
    optim.resizeParameterVector(3);
    Eigen::VectorXd x(3);
    x.setOnes();
    optim.setX(x);

    // create objective function
    auto objective = [](const Eigen::VectorXd& x, Eigen::Ref<Eigen::VectorXd> values) {
        values[0] = x[0] - 5;
        values[1] = x[1] + 3;
        values[2] = x[2];
    };
    optim.setObjectiveFunction(objective, 3, true);  // least-squares type

    EXPECT_TRUE(solver.initialize(optim));

    bool success = solver.solve(optim, true);

    EXPECT_TRUE(success);
    EXPECT_NEAR(optim.getX()[0], 5.0, 1e-6);
    EXPECT_NEAR(optim.getX()[1], -3.0, 1e-6);
    EXPECT_NEAR(optim.getX()[2], 0.0, 1e-6);
}

TEST_F(TestLevenbergMarquardtSparse, solve_rosenbrock_unconstr)
{
    // parameters
    optim.resizeParameterVector(2);
    Eigen::VectorXd x(2);
    x.setOnes();
    optim.setX(x);

    // create objective function
    // Create edges for minimizing 0.5 * (100*(x2-x1^2)^2 + (1-x1)^2 ) => f1^2 + f2^2
    auto objective = [](const Eigen::VectorXd& x, Eigen::Ref<Eigen::VectorXd> values) {
        values[0] = std::sqrt(100) * (x[1] - x[0] * x[0]);
        values[1] = 1 - x[0];
    };
    optim.setObjectiveFunction(objective, 2, true);  // least-squares type

    EXPECT_TRUE(solver.initialize(optim));

    // now solve
    bool success = solver.solve(optim, true);

    EXPECT_TRUE(success);
    EXPECT_NEAR(optim.getX()[0], 1.0, 1e-3);
    EXPECT_NEAR(optim.getX()[1], 1.0, 1e-3);
}

TEST_F(TestLevenbergMarquardtSparse, solve_eqconstr_1)
{
    // parameters
    optim.resizeParameterVector(1);
    Eigen::VectorXd x(1);
    x.setOnes();
    optim.setX(x);

    // create objective function
    auto objective = [](const Eigen::VectorXd& x, Eigen::Ref<Eigen::VectorXd> values) { values[0] = x[0] - 2; };
    optim.setObjectiveFunction(objective, 1, true);  // least-squares type

    // create equality constraint function
    auto equality = [](const Eigen::VectorXd& x, Eigen::Ref<Eigen::VectorXd> values) { values[0] = x[0] - 3; };
    optim.setEqualityConstraint(equality, 1);

    solver.setPenaltyWeights(100, 100, 100);

    EXPECT_TRUE(solver.initialize(optim));

    bool success = solver.solve(optim, true);

    EXPECT_TRUE(success);
    EXPECT_NEAR(optim.getX()[0], 3.0, 1e-4);
}

TEST_F(TestLevenbergMarquardtSparse, solve_ineqconstr_1)
{
    // parameters
    optim.resizeParameterVector(1);
    Eigen::VectorXd x(1);
    x.setOnes();
    optim.setX(x);

    // create objective function
    auto objective = [](const Eigen::VectorXd& x, Eigen::Ref<Eigen::VectorXd> values) { values[0] = x[0] - 2; };  // min -> x = 2
    optim.setObjectiveFunction(objective, 1, true);                                                               // least-squares type

    // create inequality constraint function
    auto inequality = [](const Eigen::VectorXd& x, Eigen::Ref<Eigen::VectorXd> values) { values[0] = -x[0] + 3; };  // x > 3
    optim.setInequalityConstraint(inequality, 1);

    solver.setPenaltyWeights(100, 100, 100);

    EXPECT_TRUE(solver.initialize(optim));

    bool success = solver.solve(optim, true);

    EXPECT_TRUE(success);
    EXPECT_NEAR(optim.getX()[0], 3.0, 1e-4);
}

TEST_F(TestLevenbergMarquardtSparse, solve_lower_bounds)
{
    // parameters
    optim.resizeParameterVector(1);
    Eigen::VectorXd x(1);
    x.setOnes();
    optim.setX(x);

    // set bounds
    Eigen::VectorXd lb(1);
    lb[0] = 5;
    optim.setLowerBounds(lb);

    // create objective function
    auto objective = [](const Eigen::VectorXd& x, Eigen::Ref<Eigen::VectorXd> values) { values[0] = x[0] - 2; };  // min -> x = 2
    optim.setObjectiveFunction(objective, 1, true);                                                               // least-squares type

    solver.setPenaltyWeights(100, 100, 100);

    EXPECT_TRUE(solver.initialize(optim));

    bool success = solver.solve(optim, true);

    EXPECT_TRUE(success);
    EXPECT_NEAR(optim.getX()[0], 5.0, 1e-3);
}

TEST_F(TestLevenbergMarquardtSparse, solve_upper_bounds)
{
    // parameters
    optim.resizeParameterVector(1);
    Eigen::VectorXd x(1);
    x.setOnes();
    optim.setX(x);

    // set bounds
    Eigen::VectorXd ub(1);
    ub[0] = -1;
    optim.setUpperBounds(ub);

    // create objective function
    auto objective = [](const Eigen::VectorXd& x, Eigen::Ref<Eigen::VectorXd> values) { values[0] = x[0] - 2; };  // min -> x = 2
    optim.setObjectiveFunction(objective, 1, true);                                                               // least-squares type

    solver.setPenaltyWeights(100, 100, 100);

    EXPECT_TRUE(solver.initialize(optim));

    bool success = solver.solve(optim, true);

    EXPECT_TRUE(success);
    EXPECT_NEAR(optim.getX()[0], -1.0, 1e-3);
}

TEST_F(TestLevenbergMarquardtSparse, solve_betts_fun_constr)
{
    // parameters
    optim.resizeParameterVector(2);

    // create bounds on x1 and x2
    optim.setLowerBound(0, 2);
    optim.setUpperBound(0, 50);
    optim.setLowerBound(1, -50);
    optim.setUpperBound(1, 50);

    // Problem definition
    // min 0.01 * x1^2 + x2^2 - 100
    // s.t. 2<=x1<=50,
    //      -50<=x2<=50,
    //      10*x1-x2>=10

    // Create objective for minimizing
    // 0.01 * x1^2 + x2^2 - 100 = f1^2 + f^2
    auto objective = [](const Eigen::VectorXd& x, Eigen::Ref<Eigen::VectorXd> values) {
        values[0] = std::sqrt(0.01) * x[0];
        values[1] = x[1];
    };
    optim.setObjectiveFunction(objective, 2, true);  // least-squares type

    // Create inequality for satisfying
    // 10*x1-x2>=10 -> 10*x1-x2-10>=0 -> x2-10*x1+10<=0
    auto inequality = [](const Eigen::VectorXd& x, Eigen::Ref<Eigen::VectorXd> values) {
        values[0] = x[1] - 10.0 * x[0] + 10.0;  // c(x)<=0 convention
    };
    optim.setInequalityConstraint(inequality, 1);

    // feasible start
    optim.setParameterValue(0, 5);
    optim.setParameterValue(0, -5);

    // now solve
    bool success = solver.solve(optim, true);

    EXPECT_TRUE(success);
    EXPECT_NEAR(optim.getX()[0], 2.0, 1e-2);
    EXPECT_NEAR(optim.getX()[1], 0.0, 1e-2);

    // infeasible start
    optim.setParameterValue(0, -1);
    optim.setParameterValue(0, -1);

    solver.setPenaltyWeights(1, 10, 10);
    solver.setIterations(5000);
    EXPECT_TRUE(solver.initialize(optim));

    // now solve
    success = solver.solve(optim, true);

    EXPECT_TRUE(success);
    EXPECT_NEAR(optim.getX()[0], 2.0, 1e-2);
    EXPECT_NEAR(optim.getX()[1], 0.0, 1e-2);
}

TEST_F(TestLevenbergMarquardtSparse, solve_betts_fun_constr_weight_adapt)
{
    // parameters
    optim.resizeParameterVector(2);

    // create bounds on x1 and x2
    optim.setLowerBound(0, 2);
    optim.setUpperBound(0, 50);
    optim.setLowerBound(1, -50);
    optim.setUpperBound(1, 50);

    // Problem definition
    // min 0.01 * x1^2 + x2^2 - 100
    // s.t. 2<=x1<=50,
    //      -50<=x2<=50,
    //      10*x1-x2>=10

    // Create objective for minimizing
    // 0.01 * x1^2 + x2^2 - 100 = f1^2 + f^2
    auto objective = [](const Eigen::VectorXd& x, Eigen::Ref<Eigen::VectorXd> values) {
        values[0] = std::sqrt(0.01) * x[0];
        values[1] = x[1];
    };
    optim.setObjectiveFunction(objective, 2, true);  // least-squares type

    // Create inequality for satisfying
    // 10*x1-x2>=10 -> 10*x1-x2-10>=0 -> x2-10*x1+10<=0
    auto inequality = [](const Eigen::VectorXd& x, Eigen::Ref<Eigen::VectorXd> values) {
        values[0] = x[1] - 10.0 * x[0] + 10.0;  // c(x)<=0 convention
    };
    optim.setInequalityConstraint(inequality, 1);

    // feasible start
    optim.setParameterValue(0, 5);
    optim.setParameterValue(0, -5);

    // set weight adaptation
    solver.setPenaltyWeights(2, 2, 2);
    solver.setWeightAdapation(5, 5, 5, 500, 500, 500);
    solver.setIterations(5);

    EXPECT_TRUE(solver.initialize(optim));

    // now solve
    bool success;
    for (int i = 0; i < 5; ++i) success = solver.solve(optim, true, i == 0 ? true : false);

    EXPECT_TRUE(success);
    EXPECT_NEAR(optim.getX()[0], 2.0, 1e-2);
    EXPECT_NEAR(optim.getX()[1], 0.0, 1e-2);

    // infeasible start
    optim.setParameterValue(0, -1);
    optim.setParameterValue(0, -1);

    // now solve
    for (int i = 0; i < 5; ++i) success = solver.solve(optim, true, i == 0 ? true : false);

    EXPECT_TRUE(success);
    EXPECT_NEAR(optim.getX()[0], 2.0, 1e-2);
    EXPECT_NEAR(optim.getX()[1], 0.0, 1e-2);
}
