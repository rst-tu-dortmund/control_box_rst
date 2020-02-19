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

#include <corbo-optimization/hyper_graph/edge.h>
#include <corbo-optimization/hyper_graph/generic_edge.h>
#include <corbo-optimization/hyper_graph/hyper_graph.h>

#include <corbo-core/console.h>
#include <corbo-core/macros.h>
#include <corbo-core/types.h>
#include <corbo-core/utilities.h>
#include <corbo-core/value_comparison.h>

#include <array>
#include <functional>

#include "gtest/gtest.h"

using corbo::HyperGraph;
using corbo::VertexSet;
using corbo::VertexInterface;
using corbo::OptimizationEdgeSet;
using pFVectorVertex = corbo::PartiallyFixedVectorVertex;
using corbo::VectorVertex;
using corbo::EdgeInterface;
using corbo::EdgeGenericScalarFun;
using corbo::EdgeGenericVectorFun;
using corbo::BaseEdge;
using corbo::BaseMixedEdge;
using corbo::CORBO_INF_DBL;

class TestHyperGraph : public testing::Test
{
 protected:
    // You can do set-up work for each test here.
    TestHyperGraph()
    {
        vertices = std::make_shared<VertexSet>();
        edges    = std::make_shared<OptimizationEdgeSet>();
        graph.setVertexSet(vertices);
        graph.setEdgeSet(edges);

        v1 = std::make_shared<pFVectorVertex>(1);
        v1->values().setOnes();
        v2 = std::make_shared<pFVectorVertex>(Eigen::Vector2d(1, 1));
    }

    // You can do clean-up work that doesn't throw exceptions here.
    virtual ~TestHyperGraph() {}
    // If the constructor and destructor are not enough for setting up
    // and cleaning up each test, you can define the following methods:

    // Code here will be called immediately after the constructor (right
    // before each test).
    void SetUp() override {}

    // Code here will be called immediately after each test (right
    // before the destructor).
    // virtual void TearDown()

    using Edge1T = EdgeGenericScalarFun<pFVectorVertex, pFVectorVertex>;

    // create some edges types
    static double edge1_fun(const Edge1T::VertexContainer& vertices)
    {
        return vertices[0]->getData()[0] + 2 * vertices[1]->getData()[0] + 3 * vertices[1]->getData()[1];
    }

    using Edge2T = EdgeGenericVectorFun<2, pFVectorVertex>;

    static void edge2_fun(const Edge2T::VertexContainer& vertices, Eigen::Ref<Edge2T::ErrorVector> values)
    {
        values[0] = 4 * vertices[0]->getData()[0];
        values[1] = 5 * vertices[0]->getData()[0] + 6 * vertices[0]->getData()[1];
    }

    using Edge3T = EdgeGenericVectorFun<2, pFVectorVertex, pFVectorVertex>;

    static void edge3_fun(const Edge1T::VertexContainer& vertices, Eigen::Ref<Edge2T::ErrorVector> values)
    {
        // x^2 + 2*y^2 + 3*z^2
        values[0] = vertices[0]->getData()[0] * vertices[0]->getData()[0] + 2 * vertices[1]->getData()[0] * vertices[1]->getData()[0] +
                    3 * vertices[1]->getData()[1] * vertices[1]->getData()[1];
        // x*y + 2*x*z
        values[1] = vertices[0]->getData()[0] * vertices[1]->getData()[0] + 2 * vertices[0]->getData()[0] * vertices[1]->getData()[1];
    }

    using Edge4T = EdgeGenericVectorFun<3, pFVectorVertex, pFVectorVertex>;

    static void edge4_fun(const Edge4T::VertexContainer& vertices, Eigen::Ref<Edge4T::ErrorVector> values)
    {
        values[0] = vertices[0]->getData()[0] - 5;
        values[1] = vertices[1]->getData()[0] + 3;
        values[2] = vertices[1]->getData()[1];
    }

    HyperGraph graph;
    VertexSet::Ptr vertices;
    OptimizationEdgeSet::Ptr edges;

    pFVectorVertex::Ptr v1;
    pFVectorVertex::Ptr v2;
};

TEST_F(TestHyperGraph, vertex_set_active_vertices)
{
    // Vertices
    pFVectorVertex::Ptr p1 = std::make_shared<pFVectorVertex>(Eigen::Vector2d(2, 3));
    p1->setFixed(1, true);
    VectorVertex::Ptr p2 = std::make_shared<VectorVertex>(Eigen::Vector3d(4, 5, 3));
    p2->setFixed(true);
    pFVectorVertex::Ptr p3 = std::make_shared<pFVectorVertex>(1);
    p3->setData(0, -1);

    // Create Vertex Set
    VertexSet::Ptr vertices = VertexSet::Ptr(new VertexSet({p1, p2, p3}));

    std::vector<VertexInterface*>& active_vertices = vertices->getActiveVertices();

    EXPECT_EQ(active_vertices.size(), 2);
    EXPECT_EQ(active_vertices[0], p1.get());
    EXPECT_EQ(active_vertices[1], p3.get());
}

TEST_F(TestHyperGraph, vertex_set_increment)
{
    // Vertices
    pFVectorVertex::Ptr p1 = std::make_shared<pFVectorVertex>(1);
    p1->setData(0, 1);
    pFVectorVertex::Ptr p2 = std::make_shared<pFVectorVertex>(Eigen::Vector2d(1, 1));

    // Create Vertex Set
    VertexSet::Ptr vertices = VertexSet::Ptr(new VertexSet({p1, p2}));

    EXPECT_EQ(vertices->getParameterDimension(), 3);

    Eigen::Vector3d inc(5, -1, 7);
    vertices->applyIncrementNonFixed(inc);

    EXPECT_DOUBLE_EQ(p1->values()[0], 6);
    EXPECT_DOUBLE_EQ(p2->values()[0], 0);
    EXPECT_DOUBLE_EQ(p2->values()[1], 8);
}

TEST_F(TestHyperGraph, vertex_set_increment_fixed)
{
    // Vertices
    pFVectorVertex::Ptr p1 = std::make_shared<pFVectorVertex>(1);
    p1->setData(0, 1);
    p1->setFixed(0, true);
    pFVectorVertex::Ptr p2 = std::make_shared<pFVectorVertex>(Eigen::Vector2d(1, 1));
    p2->setFixed(1, true);

    // Create Vertex Set
    VertexSet::Ptr vertices = VertexSet::Ptr(new VertexSet({p1, p2}));

    EXPECT_EQ(vertices->getParameterDimension(), 1);

    Eigen::Matrix<double, 1, 1> inc;
    inc[0] = -1;
    vertices->applyIncrementNonFixed(inc);

    EXPECT_DOUBLE_EQ(p1->values()[0], 1);
    EXPECT_DOUBLE_EQ(p2->values()[0], 0);
    EXPECT_DOUBLE_EQ(p2->values()[1], 1);
}

TEST_F(TestHyperGraph, vertex_set_increment_fixed2)
{
    // Vertices
    pFVectorVertex::Ptr p1 = std::make_shared<pFVectorVertex>(1);
    p1->setData(0, 1);
    p1->setFixed(0, true);
    pFVectorVertex::Ptr p2 = std::make_shared<pFVectorVertex>(Eigen::Vector2d(1, 1));
    p2->setFixed(0, true);

    // Create Vertex Set
    VertexSet::Ptr vertices = VertexSet::Ptr(new VertexSet({p1, p2}));

    EXPECT_EQ(vertices->getParameterDimension(), 1);

    Eigen::Matrix<double, 1, 1> inc;
    inc[0] = -1;
    vertices->applyIncrementNonFixed(inc);

    EXPECT_DOUBLE_EQ(p1->values()[0], 1);
    EXPECT_DOUBLE_EQ(p2->values()[0], 1);
    EXPECT_DOUBLE_EQ(p2->values()[1], 0);
}

TEST_F(TestHyperGraph, graph_consistency)
{
    vertices->addVertex(v1);

    // Edges
    BaseEdge::Ptr edge1(new Edge1T(edge1_fun, false, *v1, *v2));
    BaseEdge::Ptr edge2(new Edge2T(edge2_fun, false, *v2));
    edges->addEdges({edge1, edge2}, {}, {}, {}, {});

    EXPECT_FALSE(graph.checkGraphConsistency());

    vertices->addVertex(v2);

    EXPECT_TRUE(graph.checkGraphConsistency());
}

TEST_F(TestHyperGraph, numerical_jacobian1)
{
    vertices->addVertex(v1);
    vertices->addVertex(v2);

    // Edge1
    BaseEdge::Ptr edge1(new Edge1T(edge1_fun, false, *v1, *v2));

    Eigen::MatrixXd jac1_v1(edge1->getDimension(), v1->getDimensionUnfixed());
    edge1->computeJacobian(0, jac1_v1);
    EXPECT_NEAR(jac1_v1(0, 0), 1, 1e-6);

    Eigen::MatrixXd jac1_v2(edge1->getDimension(), v2->getDimensionUnfixed());
    edge1->computeJacobian(1, jac1_v2);
    EXPECT_NEAR(jac1_v2(0, 0), 2, 1e-6);
    EXPECT_NEAR(jac1_v2(0, 1), 3, 1e-6);
}

TEST_F(TestHyperGraph, numerical_jacobian2)
{
    vertices->addVertex(v1);
    vertices->addVertex(v2);

    // Edge3
    BaseEdge::Ptr edge3(new Edge3T(edge3_fun, false, *v1, *v2));

    Eigen::MatrixXd jac3_v1(edge3->getDimension(), v1->getDimensionUnfixed());
    edge3->computeJacobian(0, jac3_v1);

    Eigen::MatrixXd jac3_v1_sol(2, 1);
    jac3_v1_sol << 2, 3;
    EXPECT_EQ_MATRIX(jac3_v1, jac3_v1_sol, 1e-6);

    Eigen::MatrixXd jac3_v2(edge3->getDimension(), v2->getDimensionUnfixed());
    edge3->computeJacobian(1, jac3_v2);

    Eigen::MatrixXd jac3_v2_sol(2, 2);
    jac3_v2_sol << 4, 6, 1, 2;
    EXPECT_EQ_MATRIX(jac3_v2, jac3_v2_sol, 1e-6);
}

TEST_F(TestHyperGraph, numerical_hessian)
{
    vertices->addVertex(v1);
    vertices->addVertex(v2);

    // Edge1
    BaseEdge::Ptr edge1(new Edge1T(edge1_fun, false, *v1, *v2));

    Eigen::MatrixXd hes_v1_v1(v1->getDimensionUnfixed(), v1->getDimensionUnfixed());
    hes_v1_v1.setZero();
    edge1->computeHessianInc(0, 0, hes_v1_v1);
    EXPECT_TRUE(hes_v1_v1.isZero(1e-3)) << "hes_v1_v1: " << hes_v1_v1;

    Eigen::MatrixXd hes_v1_v2(v1->getDimensionUnfixed(), v2->getDimensionUnfixed());
    hes_v1_v2.setZero();
    edge1->computeHessianInc(0, 1, hes_v1_v2);
    EXPECT_TRUE(hes_v1_v2.isZero(1e-3)) << "hes_v1_v2: " << hes_v1_v2;

    Eigen::MatrixXd hes_v2_v1(v2->getDimensionUnfixed(), v1->getDimensionUnfixed());
    hes_v2_v1.setZero();
    edge1->computeHessianInc(1, 0, hes_v2_v1);
    EXPECT_TRUE(hes_v2_v1.isZero(1e-3)) << "hes_v2_v1: " << hes_v2_v1;

    Eigen::MatrixXd hes_v2_v2(v2->getDimensionUnfixed(), v2->getDimensionUnfixed());
    hes_v2_v2.setZero();
    edge1->computeHessianInc(1, 1, hes_v2_v2);
    EXPECT_TRUE(hes_v2_v2.isZero(1e-3)) << "hes_v2_v2: " << hes_v2_v2;
}

TEST_F(TestHyperGraph, numerical_hessian2)
{
    vertices->addVertex(v1);
    vertices->addVertex(v2);

    // Edge1
    BaseEdge::Ptr edge3(new Edge3T(edge3_fun, false, *v1, *v2));

    Eigen::MatrixXd hes_v1_v1(v1->getDimensionUnfixed(), v1->getDimensionUnfixed());
    hes_v1_v1.setZero();
    edge3->computeHessianInc(0, 0, hes_v1_v1);
    Eigen::MatrixXd hes_v1_v1_sol(v1->getDimensionUnfixed(), v1->getDimensionUnfixed());
    hes_v1_v1_sol << 2;
    EXPECT_EQ_MATRIX(hes_v1_v1, hes_v1_v1_sol, 1e-3);

    Eigen::MatrixXd hes_v1_v2(v1->getDimensionUnfixed(), v2->getDimensionUnfixed());
    hes_v1_v2.setZero();
    edge3->computeHessianInc(0, 1, hes_v1_v2);
    Eigen::MatrixXd hes_v1_v2_sol(v1->getDimensionUnfixed(), v2->getDimensionUnfixed());
    hes_v1_v2_sol << 1, 2;
    EXPECT_EQ_MATRIX(hes_v1_v2, hes_v1_v2_sol, 1e-3);

    Eigen::MatrixXd hes_v2_v1(v2->getDimensionUnfixed(), v1->getDimensionUnfixed());
    hes_v2_v1.setZero();
    edge3->computeHessianInc(1, 0, hes_v2_v1);
    Eigen::MatrixXd hes_v2_v1_sol(v2->getDimensionUnfixed(), v1->getDimensionUnfixed());
    hes_v2_v1_sol << 1, 2;
    EXPECT_EQ_MATRIX(hes_v2_v1, hes_v2_v1_sol, 1e-3);

    Eigen::MatrixXd hes_v2_v2(v2->getDimensionUnfixed(), v2->getDimensionUnfixed());
    hes_v2_v2.setZero();
    edge3->computeHessianInc(1, 1, hes_v2_v2);
    Eigen::MatrixXd hes_v2_v2_sol(v2->getDimensionUnfixed(), v2->getDimensionUnfixed());
    hes_v2_v2_sol << 4, 0, 0, 6;
    EXPECT_EQ_MATRIX(hes_v2_v2, hes_v2_v2_sol, 1e-3);
}
// TEST_F(TestHyperGraph, optim_edge_set)
//{
// Vertices
//   pFVectorVertex::Ptr p1 = std::make_shared<pFVectorVertex>(Eigen::Vector2d(2, 3));
//  p1->setFixed(1, true);
//   VectorVertex::Ptr p2 = std::make_shared<VectorVertex>(Eigen::Vector3d(4, 5, 3));
//   p2->setFixed(true);
//    pFVectorVertex::Ptr p3 = std::make_shared<pFVectorVertex>(1);
//     p3->setData(0, -1);

//    // Edge
//    using Edge1T   = EdgeGenericScalarFun<pFVectorVertex, pFVectorVertex>;
//    auto edge1_fun = [](const Edge1T::VertexContainer& vertices) {
//        return vertices[0]->getData()[0] + 2 * vertices[1]->getData()[0] + 3 * vertices[1]->getData()[1];
//    };

//    // Create Edge Set
//    OptimizationEdgeSet::Ptr vertices = VertexSet::Ptr(new VertexSet({p1, p2, p3}));
//}
