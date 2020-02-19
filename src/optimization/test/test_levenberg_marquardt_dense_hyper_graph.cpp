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

#include <corbo-optimization/solver/levenberg_marquardt_dense.h>

#include <corbo-optimization/hyper_graph/base_edge.h>
#include <corbo-optimization/hyper_graph/generic_edge.h>
#include <corbo-optimization/hyper_graph/hyper_graph.h>

#include <corbo-core/console.h>
#include <corbo-core/macros.h>
#include <corbo-core/utilities.h>
#include <corbo-core/value_comparison.h>

#include <array>
#include <functional>

#include "gtest/gtest.h"

using corbo::LevenbergMarquardtDense;
using corbo::HyperGraph2;
using corbo::VectorVertex;
using pfVectorVertex = corbo::PartiallyFixedVectorVertex;
using corbo::EdgeInterface;
using corbo::EdgeGenericScalarFun;
using corbo::EdgeGenericVectorFun;

class TestLevenbergMarquardtDenseHyperGraph : public testing::Test
{
 protected:
    // You can do set-up work for each test here.
    TestLevenbergMarquardtDenseHyperGraph() {}
    // You can do clean-up work that doesn't throw exceptions here.
    virtual ~TestLevenbergMarquardtDenseHyperGraph() {}
    // If the constructor and destructor are not enough for setting up
    // and cleaning up each test, you can define the following methods:

    // Code here will be called immediately after the constructor (right
    // before each test).
    void SetUp() override
    {
        // configure hyper graph
        graph.setEdgeIterationStrategy(false);
        // set retain lsq_form to true (required by LM)
        graph.setRetainLsqForm(true);

        // configure solver
        solver.setIterations(100);
    }
    // Code here will be called immediately after each test (right
    // before the destructor).
    // virtual void TearDown()

    HyperGraph2 graph;

    LevenbergMarquardtDense solver;
};

TEST_F(TestLevenbergMarquardtDenseHyperGraph, solve_unconstr_1)
{
    // Vertices
    VectorVertex v1(1);
    v1.values().setOnes();
    graph.addVertex(&v1);

    // Edges
    using EdgeT = EdgeGenericScalarFun<VectorVertex>;

    // create some edges types
    auto edge_fun = [](const EdgeT::VertexContainer& vertices) { return vertices[0]->getData()[0] - 2; };

    EdgeInterface::UPtr edge(new EdgeT(edge_fun, true, v1));
    graph.addObjectiveEdge(std::move(edge));

    EXPECT_TRUE(solver.initialize(graph));

    bool success = solver.solve(graph, true);

    EXPECT_TRUE(success);
    EXPECT_NEAR(v1.values()[0], 2.0, 1e-6);
}

TEST_F(TestLevenbergMarquardtDenseHyperGraph, solve_unconstr_2)
{
    // Vertices
    VectorVertex v1(1);
    v1.values().setOnes();
    VectorVertex v2(Eigen::Vector2d(1, 1));
    graph.addVertex(&v1);
    graph.addVertex(&v2);

    // Edges
    using EdgeT = EdgeGenericVectorFun<3, VectorVertex, VectorVertex>;

    auto edge_fun = [](const EdgeT::VertexContainer& vertices, Eigen::Ref<EdgeT::ErrorVector> values) {
        values[0] = vertices[0]->getData()[0] - 5;
        values[1] = vertices[1]->getData()[0] + 3;
        values[2] = vertices[1]->getData()[1];
    };

    EdgeInterface::UPtr edge(new EdgeT(edge_fun, true, v1, v2));
    graph.addObjectiveEdge(std::move(edge));

    EXPECT_TRUE(solver.initialize(graph));

    bool success = solver.solve(graph, true);

    EXPECT_TRUE(success);
    EXPECT_NEAR(v1.values()[0], 5.0, 1e-6);
    EXPECT_NEAR(v2.values()[0], -3.0, 1e-6);
    EXPECT_NEAR(v2.values()[1], 0.0, 1e-6);
}

TEST_F(TestLevenbergMarquardtDenseHyperGraph, solve_unconstr_2_fixed)
{
    // Vertices
    VectorVertex v1(1);
    v1.values().setOnes();
    pfVectorVertex v2(Eigen::Vector2d(1, 1));
    graph.addVertex(&v1);
    graph.addVertex(&v2);

    v1.setFixed(true);
    v2.setFixed(0, true);

    // Edges
    using EdgeT = EdgeGenericVectorFun<3, VectorVertex, pfVectorVertex>;

    auto edge_fun = [](const EdgeT::VertexContainer& vertices, Eigen::Ref<EdgeT::ErrorVector> values) {
        values[0] = vertices[0]->getData()[0] - 5;
        values[1] = vertices[1]->getData()[0] + 3;
        values[2] = vertices[1]->getData()[1];
    };

    EdgeInterface::UPtr edge(new EdgeT(edge_fun, true, v1, v2));
    graph.addObjectiveEdge(std::move(edge));

    EXPECT_TRUE(solver.initialize(graph));

    bool success = solver.solve(graph, true);

    EXPECT_TRUE(success);
    EXPECT_NEAR(v1.values()[0], 1.0, 1e-6);
    EXPECT_NEAR(v2.values()[0], 1.0, 1e-6);
    EXPECT_NEAR(v2.values()[1], 0.0, 1e-6);
}

TEST_F(TestLevenbergMarquardtDenseHyperGraph, solve_rosenbrock_unconstr)
{
    // Add variables
    VectorVertex x1(1);
    VectorVertex x2(1);

    x1.values()[0] = 5;
    x2.values()[0] = -5;

    graph.addVertex(&x1);
    graph.addVertex(&x2);

    // Create edges for minimizing 0.5 * (100*(x2-x1^2)^2 + (1-x1)^2 ) => f1^2 + f2^2

    // start with f1
    using Edge1T   = EdgeGenericScalarFun<VectorVertex, VectorVertex>;
    auto edge_fun1 = [](const Edge1T::VertexContainer& vertices) {
        double x1 = vertices[0]->getData()[0];
        double x2 = vertices[1]->getData()[0];
        return std::sqrt(100) * (x2 - x1 * x1);
    };

    EdgeInterface::UPtr edge1(new Edge1T(edge_fun1, true, x1, x2));
    graph.addObjectiveEdge(std::move(edge1));

    // now create f2
    using Edge2T   = EdgeGenericScalarFun<VectorVertex>;
    auto edge_fun2 = [](const Edge2T::VertexContainer& vertices) {
        double x1 = vertices[0]->getData()[0];
        return 1 - x1;
    };

    EdgeInterface::UPtr edge2(new Edge2T(edge_fun2, true, x1));
    graph.addObjectiveEdge(std::move(edge2));

    EXPECT_TRUE(solver.initialize(graph));

    // now solve
    bool success = solver.solve(graph, true);

    EXPECT_TRUE(success);
    EXPECT_NEAR(x1.values()[0], 1.0, 1e-3);
    EXPECT_NEAR(x2.values()[0], 1.0, 1e-3);
}

TEST_F(TestLevenbergMarquardtDenseHyperGraph, solve_eqconstr_1)
{
    // Vertices
    VectorVertex v1(1);
    v1.values().setOnes();
    graph.addVertex(&v1);

    // Edges
    using EdgeT = EdgeGenericScalarFun<VectorVertex>;

    // create some edges types
    auto edge_fun = [](const EdgeT::VertexContainer& vertices) { return vertices[0]->getData()[0] - 2; };
    EdgeInterface::UPtr edge1(new EdgeT(edge_fun, true, v1));
    graph.addObjectiveEdge(std::move(edge1));

    auto eq_constr = [](const EdgeT::VertexContainer& vertices) { return vertices[0]->getData()[0] - 3; };
    EdgeInterface::UPtr edge2(new EdgeT(eq_constr, false, v1));
    graph.addEqualityConstraintEdge(std::move(edge2));

    solver.setPenaltyWeights(100, 100, 100);

    EXPECT_TRUE(solver.initialize(graph));

    bool success = solver.solve(graph, true);

    EXPECT_TRUE(success);
    EXPECT_NEAR(v1.values()[0], 3.0, 1e-4);
}

TEST_F(TestLevenbergMarquardtDenseHyperGraph, solve_ineqconstr_1)
{
    // Vertices
    VectorVertex v1(1);
    v1.values().setOnes();
    graph.addVertex(&v1);

    // Edges
    using EdgeT = EdgeGenericScalarFun<VectorVertex>;

    // create some edges types
    auto edge_fun = [](const EdgeT::VertexContainer& vertices) { return vertices[0]->getData()[0] - 2; };  // min -> x = 2
    EdgeInterface::UPtr edge1(new EdgeT(edge_fun, true, v1));
    graph.addObjectiveEdge(std::move(edge1));

    auto eq_constr = [](const EdgeT::VertexContainer& vertices) { return -vertices[0]->getData()[0] + 3; };  // x > 3
    EdgeInterface::UPtr edge2(new EdgeT(eq_constr, false, v1));
    graph.addInequalityConstraintEdge(std::move(edge2));

    solver.setPenaltyWeights(100, 100, 100);

    EXPECT_TRUE(solver.initialize(graph));

    bool success = solver.solve(graph, true);

    EXPECT_TRUE(success);
    EXPECT_NEAR(v1.values()[0], 3.0, 1e-4);
}

TEST_F(TestLevenbergMarquardtDenseHyperGraph, solve_lower_bounds)
{
    // Vertices
    VectorVertex v1(1);
    v1.values().setOnes();

    Eigen::VectorXd lb(1);
    lb[0] = 5;
    v1.setLowerBounds(lb);

    graph.addVertex(&v1);

    // Edges
    using EdgeT = EdgeGenericScalarFun<VectorVertex>;

    // create some edges types
    auto edge_fun = [](const EdgeT::VertexContainer& vertices) { return vertices[0]->getData()[0] - 2; };

    EdgeInterface::UPtr edge(new EdgeT(edge_fun, true, v1));
    graph.addObjectiveEdge(std::move(edge));

    solver.setPenaltyWeights(100, 100, 100);

    EXPECT_TRUE(solver.initialize(graph));

    bool success = solver.solve(graph, true);

    EXPECT_TRUE(success);
    EXPECT_NEAR(v1.values()[0], 5.0, 1e-3);
}

TEST_F(TestLevenbergMarquardtDenseHyperGraph, solve_upper_bounds)
{
    // Vertices
    VectorVertex v1(1);
    v1.values().setOnes();

    Eigen::VectorXd ub(1);
    ub[0] = -1;
    v1.setUpperBounds(ub);

    graph.addVertex(&v1);

    // Edges
    using EdgeT = EdgeGenericScalarFun<VectorVertex>;

    // create some edges types
    auto edge_fun = [](const EdgeT::VertexContainer& vertices) { return vertices[0]->getData()[0] - 2; };

    EdgeInterface::UPtr edge(new EdgeT(edge_fun, true, v1));
    graph.addObjectiveEdge(std::move(edge));

    solver.setPenaltyWeights(100, 100, 100);

    EXPECT_TRUE(solver.initialize(graph));

    bool success = solver.solve(graph, true);

    EXPECT_TRUE(success);
    EXPECT_NEAR(v1.values()[0], -1.0, 1e-3);
}

TEST_F(TestLevenbergMarquardtDenseHyperGraph, solve_betts_fun_constr)
{
    // Add variables
    VectorVertex x1(1);
    VectorVertex x2(1);

    graph.addVertex(&x1);
    graph.addVertex(&x2);

    // Create edges for minimizing 0.01 * x1^2 + x2^2 - 100 = f1^2 + f^2  s.t. 2<=x1<=50, -50<=x2<=50, 10*x1-x2>=10

    // start with f1
    using Edge1T   = EdgeGenericScalarFun<VectorVertex>;
    auto edge_fun1 = [](const Edge1T::VertexContainer& vertices) {
        double x1 = vertices[0]->getData()[0];
        return std::sqrt(0.01) * x1;
    };
    EdgeInterface::UPtr edge1(new Edge1T(edge_fun1, true, x1));
    graph.addObjectiveEdge(std::move(edge1));

    // now create f2
    using Edge2T   = EdgeGenericScalarFun<VectorVertex>;
    auto edge_fun2 = [](const Edge2T::VertexContainer& vertices) {
        double x2 = vertices[0]->getData()[0];
        return x2;
    };
    EdgeInterface::UPtr edge2(new Edge2T(edge_fun2, true, x2));
    graph.addObjectiveEdge(std::move(edge2));

    // create bounds on x1 and x2
    x1.setLowerBound(0, 2);
    x1.setUpperBound(0, 50);
    x2.setLowerBound(0, -50);
    x2.setUpperBound(0, 50);

    // create linear inequality constraint 10*x1-x2>=10 -> 10*x1-x2-10>=0 -> x2-10*x1+10<=0
    using Edge3T   = EdgeGenericScalarFun<VectorVertex, VectorVertex>;
    auto edge_fun3 = [](const Edge3T::VertexContainer& vertices) {
        double x1 = vertices.at(0)->getData()[0];
        double x2 = vertices.at(1)->getData()[0];
        return x2 - 10.0 * x1 + 10.0;  // c(x)<=0 convention
    };
    EdgeInterface::UPtr edge3(new Edge3T(edge_fun3, false, x1, x2));
    graph.addInequalityConstraintEdge(std::move(edge3));

    // feasible start
    x1.values()[0] = 5;
    x2.values()[0] = -5;

    // now solve
    bool success = solver.solve(graph, true);

    EXPECT_TRUE(success);
    EXPECT_NEAR(x1.values()[0], 2.0, 1e-2);
    EXPECT_NEAR(x2.values()[0], 0.0, 1e-2);

    // infeasible start
    x1.values()[0] = -1;
    x2.values()[0] = -1;

    solver.setPenaltyWeights(1, 10, 10);
    solver.setIterations(5000);
    EXPECT_TRUE(solver.initialize(graph));

    // now solve
    success = solver.solve(graph, true);

    EXPECT_TRUE(success);
    EXPECT_NEAR(x1.values()[0], 2.0, 1e-2);
    EXPECT_NEAR(x2.values()[0], 0.0, 1e-2);
}

TEST_F(TestLevenbergMarquardtDenseHyperGraph, solve_betts_fun_constr_weight_adapt)
{
    // Add variables
    VectorVertex x1(1);
    VectorVertex x2(1);

    graph.addVertex(&x1);
    graph.addVertex(&x2);

    // Create edges for minimizing 0.01 * x1^2 + x2^2 - 100 = f1^2 + f^2  s.t. 2<=x1<=50, -50<=x2<=50, 10*x1-x2>=10

    // start with f1
    using Edge1T   = EdgeGenericScalarFun<VectorVertex>;
    auto edge_fun1 = [](const Edge1T::VertexContainer& vertices) {
        double x1 = vertices[0]->getData()[0];
        return std::sqrt(0.01) * x1;
    };
    EdgeInterface::UPtr edge1(new Edge1T(edge_fun1, true, x1));
    graph.addObjectiveEdge(std::move(edge1));

    // now create f2
    using Edge2T   = EdgeGenericScalarFun<VectorVertex>;
    auto edge_fun2 = [](const Edge2T::VertexContainer& vertices) {
        double x2 = vertices[0]->getData()[0];
        return x2;
    };
    EdgeInterface::UPtr edge2(new Edge2T(edge_fun2, true, x2));
    graph.addObjectiveEdge(std::move(edge2));

    // create bounds on x1 and x2
    x1.setLowerBound(0, 2);
    x1.setUpperBound(0, 50);
    x2.setLowerBound(0, -50);
    x2.setUpperBound(0, 50);

    // create linear inequality constraint 10*x1-x2>=10 -> 10*x1-x2-10>=0 -> x2-10*x1+10<=0
    using Edge3T   = EdgeGenericScalarFun<VectorVertex, VectorVertex>;
    auto edge_fun3 = [](const Edge3T::VertexContainer& vertices) {
        double x1 = vertices.at(0)->getData()[0];
        double x2 = vertices.at(1)->getData()[0];
        return x2 - 10.0 * x1 + 10.0;  // c(x)<=0 convention
    };
    EdgeInterface::UPtr edge3(new Edge3T(edge_fun3, false, x1, x2));
    graph.addInequalityConstraintEdge(std::move(edge3));

    // feasible start
    x1.values()[0] = 5;
    x2.values()[0] = -5;

    // set weight adaptation
    solver.setPenaltyWeights(2, 2, 2);
    solver.setWeightAdapation(5, 5, 5, 500, 500, 500);
    solver.setIterations(5);

    EXPECT_TRUE(solver.initialize(graph));

    // now solve
    bool success;
    for (int i = 0; i < 5; ++i) success = solver.solve(graph, true, i == 0 ? true : false);

    EXPECT_TRUE(success);
    EXPECT_NEAR(x1.values()[0], 2.0, 1e-2);
    EXPECT_NEAR(x2.values()[0], 0.0, 1e-2);

    // infeasible start
    x1.values()[0] = -1;
    x2.values()[0] = -1;

    // now solve
    for (int i = 0; i < 5; ++i) success = solver.solve(graph, true, i == 0 ? true : false);

    EXPECT_TRUE(success);
    EXPECT_NEAR(x1.values()[0], 2.0, 1e-2);
    EXPECT_NEAR(x2.values()[0], 0.0, 1e-2);
}

TEST_F(TestLevenbergMarquardtDenseHyperGraph, solve_betts_fun_constr_fixed)
{
    // Add variables
    pfVectorVertex x1(1);
    VectorVertex x2(1);

    x1.setFixed(0, true);

    graph.addVertex(&x1);
    graph.addVertex(&x2);

    // Create edges for minimizing 0.01 * x1^2 + x2^2 - 100 = f1^2 + f^2  s.t. 2<=x1<=50, -50<=x2<=50, 10*x1-x2>=10

    // start with f1
    using Edge1T   = EdgeGenericScalarFun<VectorVertex>;
    auto edge_fun1 = [](const Edge1T::VertexContainer& vertices) {
        double x1 = vertices[0]->getData()[0];
        return std::sqrt(0.01) * x1;
    };
    EdgeInterface::UPtr edge1(new Edge1T(edge_fun1, true, x1));
    graph.addObjectiveEdge(std::move(edge1));

    // now create f2
    using Edge2T   = EdgeGenericScalarFun<VectorVertex>;
    auto edge_fun2 = [](const Edge2T::VertexContainer& vertices) {
        double x2 = vertices[0]->getData()[0];
        return x2;
    };
    EdgeInterface::UPtr edge2(new Edge2T(edge_fun2, true, x2));
    graph.addObjectiveEdge(std::move(edge2));

    // create bounds on x1 and x2
    x1.setLowerBound(0, 2);
    x1.setUpperBound(0, 50);
    x2.setLowerBound(0, -50);
    x2.setUpperBound(0, 50);

    // create linear inequality constraint 10*x1-x2>=10 -> 10*x1-x2-10>=0 -> x2-10*x1+10<=0
    using Edge3T   = EdgeGenericScalarFun<VectorVertex, VectorVertex>;
    auto edge_fun3 = [](const Edge3T::VertexContainer& vertices) {
        double x1 = vertices.at(0)->getData()[0];
        double x2 = vertices.at(1)->getData()[0];
        return x2 - 10.0 * x1 + 10.0;  // c(x)<=0 convention
    };
    EdgeInterface::UPtr edge3(new Edge3T(edge_fun3, false, x1, x2));
    graph.addInequalityConstraintEdge(std::move(edge3));

    // feasible start
    x1.values()[0] = 2;
    x2.values()[0] = -5;

    // now solve
    bool success = solver.solve(graph, true);

    EXPECT_TRUE(success);
    EXPECT_NEAR(x1.values()[0], 2.0, 1e-2);
    EXPECT_NEAR(x2.values()[0], 0.0, 1e-2);

    // infeasible start
    x1.values()[0] = 2;
    x2.values()[0] = -4;

    solver.setPenaltyWeights(1, 10, 10);
    solver.setIterations(100);
    EXPECT_TRUE(solver.initialize(graph));

    // now solve
    success = solver.solve(graph, true);

    EXPECT_TRUE(success);
    EXPECT_NEAR(x1.values()[0], 2.0, 1e-2);
    EXPECT_NEAR(x2.values()[0], 0.0, 1e-2);
}
