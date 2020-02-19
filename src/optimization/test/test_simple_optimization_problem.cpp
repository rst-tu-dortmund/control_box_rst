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

#include <corbo-optimization/simple_optimization_problem.h>

#include <corbo-core/console.h>
#include <corbo-core/macros.h>
#include <corbo-core/types.h>
#include <corbo-core/utilities.h>
#include <corbo-core/value_comparison.h>
#include <corbo-optimization/misc.h>

#include <array>
#include <functional>

#include "gtest/gtest.h"

using corbo::SimpleOptimizationProblemWithCallbacks;
using corbo::CORBO_INF_DBL;

class TestStandardOptimizationProblem : public testing::Test
{
 protected:
    // You can do set-up work for each test here.
    TestStandardOptimizationProblem() {}
    // You can do clean-up work that doesn't throw exceptions here.
    virtual ~TestStandardOptimizationProblem() {}
    // If the constructor and destructor are not enough for setting up
    // and cleaning up each test, you can define the following methods:

    // Code here will be called immediately after the constructor (right
    // before each test).
    void SetUp() override {}
    // Code here will be called immediately after each test (right
    // before the destructor).
    // virtual void TearDown()

    SimpleOptimizationProblemWithCallbacks optim;
};

TEST_F(TestStandardOptimizationProblem, compute_objective)
{
    Eigen::VectorXd x(1);
    x[0] = 2;
    optim.setX(x);

    // define objective function
    // f = x^2 +1
    auto objective = [](const Eigen::VectorXd& x, Eigen::Ref<Eigen::VectorXd> values) { values[0] = x[0] * x[0] + 1; };
    optim.setObjectiveFunction(objective, 1, false);

    EXPECT_EQ(optim.getObjectiveDimension(), 1);
    EXPECT_EQ(optim.getNonLsqObjectiveDimension(), 1);
    EXPECT_EQ(optim.getLsqObjectiveDimension(), 0);

    // compute values
    double value = optim.computeValueObjective();
    EXPECT_EQ(value, 5) << "x*x+1 with x=2";

    value = optim.computeValueNonLsqObjective();
    EXPECT_EQ(value, 5) << "x*x+1 with x=2";

    // compute gradient
    Eigen::VectorXd gradient(optim.getParameterDimension());
    optim.computeGradientObjective(gradient);

    EXPECT_NEAR(gradient(0), 4, 1e-5) << "2*x = 2*2 = 4";

    // compute hessian
    Eigen::MatrixXd hessian(optim.getParameterDimension(), optim.getParameterDimension());
    optim.computeDenseHessianObjective(hessian);

    EXPECT_NEAR(hessian(0, 0), 2, 1e-5);
}

TEST_F(TestStandardOptimizationProblem, compute_objective_lsq)
{
    Eigen::VectorXd x(1);
    x[0] = 2;
    optim.setX(x);

    // define lsq objective function
    // f = (x-1)^2 + (6-x)^2
    auto objective = [](const Eigen::VectorXd& x, Eigen::Ref<Eigen::VectorXd> values) {
        values[0] = x[0] - 1;
        values[1] = 6 - x[0];
    };
    optim.setObjectiveFunction(objective, 2, true);

    EXPECT_EQ(optim.getObjectiveDimension(), 1);
    EXPECT_EQ(optim.getNonLsqObjectiveDimension(), 0);
    EXPECT_EQ(optim.getLsqObjectiveDimension(), 2);

    // compute values
    Eigen::Vector2d lsq_values;
    optim.computeValuesLsqObjective(lsq_values);
    EXPECT_EQ(lsq_values(0), 1) << "(2-1)";
    EXPECT_EQ(lsq_values(1), 4) << "(6-2)";

    double value = optim.computeValueObjective();

    EXPECT_EQ(value, lsq_values.squaredNorm()) << "(2-1)^2 + (6-4)^2 = 1 + 16 = 17";

    // compute lsq jacobian
    Eigen::MatrixXd lsq_jacob(optim.getLsqObjectiveDimension(), optim.getParameterDimension());
    optim.computeDenseJacobianLsqObjective(lsq_jacob);
    EXPECT_NEAR(lsq_jacob(0, 0), 1, 1e-5);
    EXPECT_NEAR(lsq_jacob(1, 0), -1, 1e-5);

    // compute overall objective gradient
    // df/dx = 2*(x-1) - 2*(6-x) = 2*1 - 2*4 = -6;
    Eigen::VectorXd gradient(optim.getParameterDimension());
    optim.computeGradientObjective(gradient);
    EXPECT_NEAR(gradient(0), -6, 1e-5);

    // compute overall objective hessian
    // df/dx^2 = 2 + 2 = 4
    Eigen::MatrixXd hessian(optim.getParameterDimension(), optim.getParameterDimension());
    optim.computeDenseHessianObjective(hessian);
    EXPECT_NEAR(hessian(0, 0), 4, 1e-4);
}

TEST_F(TestStandardOptimizationProblem, compute_objective_sparse_lsq)
{
    Eigen::VectorXd x(1);
    x[0] = 2;
    optim.setX(x);

    // define lsq objective function
    // f = (x-1)^2 + (6-x)^2
    auto objective = [](const Eigen::VectorXd& x, Eigen::Ref<Eigen::VectorXd> values) {
        values[0] = x[0] - 1;
        values[1] = 6 - x[0];
    };
    optim.setObjectiveFunction(objective, 2, true);

    // compute sparse jacobian
    int nnz = optim.computeSparseJacobianLsqObjectiveNNZ();
    EXPECT_GE(nnz, 2);

    Eigen::VectorXi irow(nnz);
    Eigen::VectorXi icol(nnz);
    Eigen::VectorXd values(nnz);
    optim.computeSparseJacobianLsqObjectiveStructure(irow, icol);
    optim.computeSparseJacobianLsqObjectiveValues(values);

    // convert to Eigen type for convenient comparison
    Eigen::SparseMatrix<double> jacob_sparse(optim.getLsqObjectiveDimension(), optim.getParameterDimension());
    std::vector<Eigen::Triplet<double>> triplets;
    corbo::convert_triplet(irow, icol, values, triplets);
    if (!triplets.empty()) jacob_sparse.setFromTriplets(triplets.begin(), triplets.end());

    Eigen::MatrixXd jacob_sol(2, 1);
    jacob_sol << 1, -1;

    EXPECT_EQ_MATRIX(jacob_sparse, jacob_sol, 1e-5);

    // compute sparse eigen jacobian directly
    Eigen::SparseMatrix<double> jacob_sparse_eigen(optim.getLsqObjectiveDimension(), optim.getParameterDimension());
    optim.computeSparseJacobianLsqObjective(jacob_sparse_eigen);
    EXPECT_EQ_MATRIX(jacob_sparse_eigen, jacob_sol, 1e-5);

    // compute sparse hessian
    int hess_nnz = optim.computeSparseHessianObjectiveNNZ();

    Eigen::VectorXi hess_irow(hess_nnz);
    Eigen::VectorXi hess_icol(hess_nnz);
    Eigen::VectorXd hess_values(hess_nnz);
    optim.computeSparseHessianObjectiveStructure(hess_irow, hess_icol);
    optim.computeSparseHessianObjectiveValues(hess_values);

    // convert to Eigen type for convenient comparison
    triplets.clear();
    Eigen::SparseMatrix<double> sparse_hessian(optim.getParameterDimension(), optim.getParameterDimension());
    corbo::convert_triplet(hess_irow, hess_icol, hess_values, triplets);
    if (!triplets.empty()) sparse_hessian.setFromTriplets(triplets.begin(), triplets.end());

    Eigen::MatrixXd hessian_sol(1, 1);
    hessian_sol << 4;

    EXPECT_EQ_MATRIX(sparse_hessian, hessian_sol, 1e-4);

    // compute sparse eigen hessian directly
    Eigen::SparseMatrix<double> sparse_hessian_eigen(optim.getParameterDimension(), optim.getParameterDimension());
    optim.computeSparseHessianObjective(sparse_hessian_eigen);
    EXPECT_EQ_MATRIX(sparse_hessian_eigen, hessian_sol, 1e-4);
}

TEST_F(TestStandardOptimizationProblem, compute_equality_constr)
{
    Eigen::VectorXd x(3);
    x[0] = 2;
    x[1] = 3;
    x[2] = 4;
    optim.setX(x);

    // define equality constraint function
    // ceq = [x1^2+1; x2+x3-10]
    auto ceq = [](const Eigen::VectorXd& x, Eigen::Ref<Eigen::VectorXd> values) {
        values[0] = x[0] * x[0] + 1;
        values[1] = x[1] + x[2] - 10;
    };
    optim.setEqualityConstraint(ceq, 2);

    // compute values
    Eigen::VectorXd values(optim.getEqualityDimension());
    optim.computeValuesEquality(values);

    EXPECT_EQ(values[0], 5) << "x*x+1 with x=2";
    EXPECT_EQ(values[1], -3) << "3+4-10";

    // compute jacobian
    Eigen::MatrixXd jacobian(optim.getEqualityDimension(), optim.getParameterDimension());
    optim.computeDenseJacobianEqualities(jacobian, nullptr);

    Eigen::MatrixXd jacob_sol(2, 3);
    jacob_sol << 4, 0, 0, 0, 1, 1;

    EXPECT_EQ_MATRIX(jacobian, jacob_sol, 1e-5);

    // compute hessian
    Eigen::MatrixXd hessian(optim.getParameterDimension(), optim.getParameterDimension());
    optim.computeDenseHessianEqualities(hessian, nullptr);

    Eigen::MatrixXd hessian_sol(3, 3);
    hessian_sol << 2, 0, 0, 0, 0, 0, 0, 0, 0;

    EXPECT_EQ_MATRIX(hessian, hessian_sol, 1e-5);
}

TEST_F(TestStandardOptimizationProblem, compute_equality_constr_sparse)
{
    Eigen::VectorXd x(3);
    x[0] = 2;
    x[1] = 3;
    x[2] = 4;
    optim.setX(x);

    // define equality constraint
    // ceq = [x1^2+1; x2+x3-10]
    auto ceq = [](const Eigen::VectorXd& x, Eigen::Ref<Eigen::VectorXd> values) {
        values[0] = x[0] * x[0] + 1;
        values[1] = x[1] + x[2] - 10;
    };
    optim.setEqualityConstraint(ceq, 2);

    // compute sparse jacobian
    int nnz = optim.computeSparseJacobianEqualitiesNNZ();
    // This formulation does not support sparse matrices actually, so we check if default conversions are working
    EXPECT_EQ(nnz, optim.getParameterDimension() * optim.getEqualityDimension());

    Eigen::VectorXi irow(nnz);
    Eigen::VectorXi icol(nnz);
    Eigen::VectorXd values(nnz);
    optim.computeSparseJacobianEqualitiesStructure(irow, icol);
    optim.computeSparseJacobianEqualitiesValues(values);

    // convert to Eigen type for convenient comparison
    std::vector<Eigen::Triplet<double>> triplet_list;
    for (int i = 0; i < nnz; ++i) triplet_list.emplace_back(irow[i], icol[i], values[i]);
    Eigen::SparseMatrix<double> sparse_jacobian(optim.getEqualityDimension(), optim.getParameterDimension());
    sparse_jacobian.setFromTriplets(triplet_list.begin(), triplet_list.end());

    Eigen::MatrixXd jacob_sol(2, 3);
    jacob_sol << 4, 0, 0, 0, 1, 1;

    EXPECT_EQ_MATRIX(sparse_jacobian, jacob_sol, 1e-5);

    // compute sparse eigen jacobian directly
    Eigen::SparseMatrix<double> sparse_jacob_eigen(optim.getEqualityDimension(), optim.getParameterDimension());
    optim.computeSparseJacobianEqualities(sparse_jacob_eigen);
    EXPECT_EQ_MATRIX(sparse_jacob_eigen, jacob_sol, 1e-5);

    // compute sparse hessian
    int hess_nnz = optim.computeSparseHessianEqualitiesNNZ();

    Eigen::VectorXi hess_irow(hess_nnz);
    Eigen::VectorXi hess_icol(hess_nnz);
    Eigen::VectorXd hess_values(hess_nnz);
    optim.computeSparseHessianEqualitiesStructure(hess_irow, hess_icol);
    optim.computeSparseHessianEqualitiesValues(hess_values);

    // convert to Eigen type for convenient comparison
    triplet_list.clear();
    for (int i = 0; i < hess_nnz; ++i) triplet_list.emplace_back(hess_irow[i], hess_icol[i], hess_values[i]);
    Eigen::SparseMatrix<double> sparse_hessian(optim.getParameterDimension(), optim.getParameterDimension());
    sparse_hessian.setFromTriplets(triplet_list.begin(), triplet_list.end());

    Eigen::MatrixXd hessian_sol(3, 3);
    hessian_sol << 2, 0, 0, 0, 0, 0, 0, 0, 0;

    EXPECT_EQ_MATRIX(sparse_hessian, hessian_sol, 1e-5);

    // compute sparse eigen hessian directly
    Eigen::SparseMatrix<double> sparse_hessian_eigen(optim.getParameterDimension(), optim.getParameterDimension());
    optim.computeSparseHessianEqualities(sparse_hessian_eigen);
    EXPECT_EQ_MATRIX(sparse_hessian_eigen, hessian_sol, 1e-5);
}

TEST_F(TestStandardOptimizationProblem, compute_inequality_constr)
{
    Eigen::VectorXd x(3);
    x[0] = 2;
    x[1] = 3;
    x[2] = -4;
    optim.setX(x);

    // define inequality constraint function
    // c = [x1^2+1; x2+x3-9]
    auto c = [](const Eigen::VectorXd& x, Eigen::Ref<Eigen::VectorXd> values) {
        values[0] = x[0] * x[0] + 1;
        values[1] = x[1] + x[2] - 9;
    };
    optim.setInequalityConstraint(c, 2);

    // compute values
    Eigen::VectorXd values(optim.getInequalityDimension());
    optim.computeValuesInequality(values);

    EXPECT_EQ(values[0], 5) << "x*x+1 with x=2";
    EXPECT_EQ(values[1], -10) << "3-4-9";

    // compute only active constraints max(0,c(x))
    optim.computeValuesActiveInequality(values);
    EXPECT_EQ(values[0], 5);
    EXPECT_EQ(values[1], 0);

    // compute jacobian
    Eigen::MatrixXd jacobian(optim.getInequalityDimension(), optim.getParameterDimension());
    optim.computeDenseJacobianInequalities(jacobian, nullptr);

    Eigen::MatrixXd jacob_sol(2, 3);
    jacob_sol << 4, 0, 0, 0, 1, 1;
    EXPECT_EQ_MATRIX(jacobian, jacob_sol, 1e-5);

    // compute jacobian with zero-rows for inactive constraints
    optim.computeDenseJacobianActiveInequalities(jacobian);

    Eigen::MatrixXd active_ineq_jacob_sol(2, 3);
    active_ineq_jacob_sol << 4, 0, 0, 0, 0, 0;
    EXPECT_EQ_MATRIX(jacobian, active_ineq_jacob_sol, 1e-5);

    // compute hessian
    Eigen::MatrixXd hessian(optim.getParameterDimension(), optim.getParameterDimension());
    optim.computeDenseHessianInequalities(hessian, nullptr);

    Eigen::MatrixXd hessian_sol(3, 3);
    hessian_sol << 2, 0, 0, 0, 0, 0, 0, 0, 0;

    EXPECT_EQ_MATRIX(hessian, hessian_sol, 1e-5);
}

TEST_F(TestStandardOptimizationProblem, compute_inequality_constr_sparse)
{
    Eigen::VectorXd x(3);
    x[0] = 2;
    x[1] = 3;
    x[2] = -4;
    optim.setX(x);

    // define inequality constraint function
    // c = [x1^2+1; x2+x3-9]
    auto c = [](const Eigen::VectorXd& x, Eigen::Ref<Eigen::VectorXd> values) {
        values[0] = x[0] * x[0] + 1;
        values[1] = x[1] + x[2] - 9;
    };
    optim.setInequalityConstraint(c, 2);

    // compute sparse jacobian
    int nnz = optim.computeSparseJacobianInequalitiesNNZ();
    // This formulation does not support sparse matrices actually, so we check if default conversions are working
    EXPECT_EQ(nnz, optim.getParameterDimension() * optim.getInequalityDimension());

    Eigen::VectorXi irow(nnz);
    Eigen::VectorXi icol(nnz);
    Eigen::VectorXd values(nnz);
    optim.computeSparseJacobianInequalitiesStructure(irow, icol);
    optim.computeSparseJacobianInequalitiesValues(values);

    // convert to Eigen type for convenient comparison
    std::vector<Eigen::Triplet<double>> triplet_list;
    for (int i = 0; i < nnz; ++i) triplet_list.emplace_back(irow[i], icol[i], values[i]);
    Eigen::SparseMatrix<double> sparse_jacobian(optim.getInequalityDimension(), optim.getParameterDimension());
    sparse_jacobian.setFromTriplets(triplet_list.begin(), triplet_list.end());

    Eigen::MatrixXd jacob_sol(2, 3);
    jacob_sol << 4, 0, 0, 0, 1, 1;

    EXPECT_EQ_MATRIX(sparse_jacobian, jacob_sol, 1e-5);

    // compute sparse eigen matrix directly
    Eigen::SparseMatrix<double> sparse_jacob_eigen(optim.getInequalityDimension(), optim.getParameterDimension());
    optim.computeSparseJacobianInequalities(sparse_jacob_eigen);
    EXPECT_EQ_MATRIX(sparse_jacob_eigen, jacob_sol, 1e-5);

    // compute jacobian with zero-rows for inactive constraints
    optim.computeSparseJacobianActiveInequalitiesValues(values);
    triplet_list.clear();
    for (int i = 0; i < nnz; ++i) triplet_list.emplace_back(irow[i], icol[i], values[i]);
    sparse_jacobian.setFromTriplets(triplet_list.begin(), triplet_list.end());

    Eigen::MatrixXd active_ineq_jacob_sol(2, 3);
    active_ineq_jacob_sol << 4, 0, 0, 0, 0, 0;
    EXPECT_EQ_MATRIX(sparse_jacobian, active_ineq_jacob_sol, 1e-5);

    // compute sparse eigen jacobian directly
    optim.computeSparseJacobianActiveInequalities(sparse_jacob_eigen);
    EXPECT_EQ_MATRIX(sparse_jacob_eigen, active_ineq_jacob_sol, 1e-5);

    // compute sparse hessian
    int hess_nnz = optim.computeSparseHessianInequalitiesNNZ();

    Eigen::VectorXi hess_irow(hess_nnz);
    Eigen::VectorXi hess_icol(hess_nnz);
    Eigen::VectorXd hess_values(hess_nnz);
    optim.computeSparseHessianInequalitiesStructure(hess_irow, hess_icol);
    optim.computeSparseHessianInequalitiesValues(hess_values);

    // convert to Eigen type for convenient comparison
    triplet_list.clear();
    for (int i = 0; i < hess_nnz; ++i) triplet_list.emplace_back(hess_irow[i], hess_icol[i], hess_values[i]);
    Eigen::SparseMatrix<double> sparse_hessian(optim.getParameterDimension(), optim.getParameterDimension());
    sparse_hessian.setFromTriplets(triplet_list.begin(), triplet_list.end());

    Eigen::MatrixXd hessian_sol(3, 3);
    hessian_sol << 2, 0, 0, 0, 0, 0, 0, 0, 0;

    EXPECT_EQ_MATRIX(sparse_hessian, hessian_sol, 1e-5);

    // compute sparse eigen hessian directly
    Eigen::SparseMatrix<double> sparse_hessian_eigen(optim.getParameterDimension(), optim.getParameterDimension());
    optim.computeSparseHessianInequalities(sparse_hessian_eigen);
    EXPECT_EQ_MATRIX(sparse_hessian_eigen, hessian_sol, 1e-5);
}

TEST_F(TestStandardOptimizationProblem, test_bounds)
{
    optim.resizeParameterVector(3);

    Eigen::VectorXd x(3);
    x[0] = 2;
    x[1] = 3;
    x[2] = 4;
    optim.setParameterVector(x);

    // set bounds
    Eigen::VectorXd lb(3);
    lb[0] = -CORBO_INF_DBL;
    lb[1] = -2;
    lb[2] = -CORBO_INF_DBL;
    optim.setLowerBounds(lb);

    // optim.setUpperBound(0, CORBO_INF_DBL); //<- Should be initialized with CORBO_INF_DBL
    optim.setUpperBound(1, 5);
    optim.setUpperBound(2, 6);

    // check bounds
    Eigen::VectorXd lb_return(optim.getParameterDimension());
    Eigen::VectorXd ub_return(optim.getParameterDimension());
    optim.getBounds(lb_return, ub_return);

    EXPECT_EQ_MATRIX(lb_return, lb, 0.01);
    EXPECT_GE(ub_return[0], CORBO_INF_DBL);
    EXPECT_EQ(ub_return[1], 5);
    EXPECT_EQ(ub_return[2], 6);

    // check dimension of active bounds
    EXPECT_EQ(optim.finiteBoundsDimension(), 3);
    EXPECT_EQ(optim.finiteCombinedBoundsDimension(), 2);

    // Jacobian for bounds

    // Reduced Jacobian for active bounds only
    Eigen::MatrixXd active_bound_jacobian(optim.finiteCombinedBoundsDimension(), optim.getParameterDimension());
    optim.computeDenseJacobianFiniteCombinedBounds(active_bound_jacobian);

    EXPECT_ZERO_MATRIX(active_bound_jacobian, 1e-5);

    // now change parameters such that bounds are violated
    optim.setParameterValue(1, -3);    // violating bounds now
    optim.setParameterValue(2, 6.01);  // here as well
    Eigen::MatrixXd active_bound_jacobian_sol(2, 3);
    active_bound_jacobian_sol << 0, -1, 0, 0, 0, 1;

    optim.computeDenseJacobianFiniteCombinedBounds(active_bound_jacobian);
    EXPECT_EQ_MATRIX(active_bound_jacobian, active_bound_jacobian_sol, 1e-5);
}

TEST_F(TestStandardOptimizationProblem, test_bounds_sparse)
{
    optim.resizeParameterVector(3);

    Eigen::VectorXd x(3);
    x[0] = 2;
    x[1] = 3;
    x[2] = 4;
    optim.setParameterVector(x);

    // set bounds
    Eigen::VectorXd lb(3);
    lb[0] = -CORBO_INF_DBL;
    lb[1] = -2;
    lb[2] = -CORBO_INF_DBL;
    optim.setLowerBounds(lb);

    // optim.setUpperBound(0, CORBO_INF_DBL); //<- Should be initialized with CORBO_INF_DBL
    optim.setUpperBound(1, 5);
    optim.setUpperBound(2, 6);

    // Reduced Jacobian for active bounds only
    // compute sparse jacobian
    int nnz = optim.computeSparseJacobianFiniteCombinedBoundsNNZ();

    Eigen::VectorXi irow(nnz);
    Eigen::VectorXi icol(nnz);
    Eigen::VectorXd values(nnz);
    optim.computeSparseJacobianFiniteCombinedBoundsStructure(irow, icol);
    optim.computeSparseJacobianFiniteCombinedBoundsValues(values);

    // convert to Eigen type for convenient comparison
    std::vector<Eigen::Triplet<double>> triplet_list;
    for (int i = 0; i < nnz; ++i) triplet_list.emplace_back(irow[i], icol[i], values[i]);
    Eigen::SparseMatrix<double> sparse_jacobian(optim.finiteCombinedBoundsDimension(), optim.getParameterDimension());
    sparse_jacobian.setFromTriplets(triplet_list.begin(), triplet_list.end());
    EXPECT_ZERO_MATRIX(sparse_jacobian.toDense(), 1e-5);

    // compute sparse eigen matrix directly
    Eigen::SparseMatrix<double> sparse_jacob_eigen(optim.finiteCombinedBoundsDimension(), optim.getParameterDimension());
    optim.computeSparseJacobianFiniteCombinedBounds(sparse_jacob_eigen);
    EXPECT_ZERO_MATRIX(sparse_jacobian.toDense(), 1e-5);

    // now change parameters such that bounds are violated
    optim.setParameterValue(1, -3);    // violating bounds now
    optim.setParameterValue(2, 6.01);  // here as well

    optim.computeSparseJacobianFiniteCombinedBoundsValues(values);

    triplet_list.clear();
    for (int i = 0; i < nnz; ++i) triplet_list.emplace_back(irow[i], icol[i], values[i]);
    sparse_jacobian.setFromTriplets(triplet_list.begin(), triplet_list.end());

    Eigen::MatrixXd active_bound_jacobian_sol(2, 3);
    active_bound_jacobian_sol << 0, -1, 0, 0, 0, 1;

    EXPECT_EQ_MATRIX(sparse_jacobian, active_bound_jacobian_sol, 1e-5);

    // compute sparse eigen matrix directly
    optim.computeSparseJacobianFiniteCombinedBounds(sparse_jacob_eigen);
    EXPECT_EQ_MATRIX(sparse_jacob_eigen, active_bound_jacobian_sol, 1e-5);
}

TEST_F(TestStandardOptimizationProblem, test_combined_jacobian_sparse)
{
    optim.resizeParameterVector(3);

    Eigen::VectorXd x(3);
    x[0] = 2;
    x[1] = -3;
    x[2] = 7;
    optim.setParameterVector(x);

    // set bounds
    Eigen::VectorXd lb(3);
    lb[0] = -CORBO_INF_DBL;
    lb[1] = -2;
    lb[2] = -CORBO_INF_DBL;
    optim.setLowerBounds(lb);

    // optim.setUpperBound(0, CORBO_INF_DBL); //<- Should be initialized with CORBO_INF_DBL
    optim.setUpperBound(1, 5);
    optim.setUpperBound(2, 6);

    // define objective function
    // f = (x^2 +1)^2
    auto objective = [](const Eigen::VectorXd& x, Eigen::Ref<Eigen::VectorXd> values) { values[0] = x[0] * x[0] + 1; };
    optim.setObjectiveFunction(objective, 1, true);

    // define equality constraint function
    // ceq = [x1^2+1; x2+x3-10]
    auto ceq = [](const Eigen::VectorXd& x, Eigen::Ref<Eigen::VectorXd> values) {
        values[0] = x[0] * x[0] + 1;
        values[1] = x[1] + x[2] - 10;
    };
    optim.setEqualityConstraint(ceq, 2);

    // define inequality constraint function
    // c = [x1^2+1; x2+x3-9]
    auto c = [](const Eigen::VectorXd& x, Eigen::Ref<Eigen::VectorXd> values) {
        values[0] = x[0] * x[0] + 1;
        values[1] = x[1] + x[2] - 9;
    };
    optim.setInequalityConstraint(c, 2);

    int dim_obj_lsq = optim.getLsqObjectiveDimension();
    int dim_eq      = optim.getEqualityDimension();
    int dim_ineq    = optim.getInequalityDimension();
    int dim_bounds  = optim.finiteCombinedBoundsDimension();

    int dim_total = dim_obj_lsq + dim_eq + dim_ineq + dim_bounds;

    Eigen::SparseMatrix<double> jacobian(dim_total, optim.getParameterDimension());
    optim.computeCombinedSparseJacobian(jacobian, true, true, true, true, true, 2, 3, 4);

    Eigen::MatrixXd sol1(7, 3);
    sol1 << 4, 0, 0, 4, 0, 0, 0, 1, 1, 4, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 1;

    sol1.middleRows(dim_obj_lsq, dim_eq) *= 2;
    sol1.middleRows(dim_obj_lsq + dim_eq, dim_ineq) *= 3;
    sol1.middleRows(dim_obj_lsq + dim_eq + dim_ineq, dim_bounds) *= 4;
    EXPECT_EQ_MATRIX(jacobian, sol1, 1e-5);

    optim.computeCombinedSparseJacobian(jacobian, true, true, true, true, false, 2, 3, 4);

    Eigen::MatrixXd sol2 = sol1;
    sol2(4, 1) = 1 * 3;
    sol2(4, 2) = 1 * 3;
    EXPECT_EQ_MATRIX(jacobian, sol2, 1e-5);

    // test triplet version
    Eigen::MatrixXd sol3(5, 3);
    sol3 << 4, 0, 0, 4, 0, 0, 0, 1, 1, 4, 0, 0, 0, 1, 1;

    int jac_nnz = optim.computeCombinedSparseJacobiansNNZ();

    Eigen::VectorXi jac_irow(jac_nnz);
    Eigen::VectorXi jac_icol(jac_nnz);
    Eigen::VectorXd jac_values(jac_nnz);
    optim.computeCombinedSparseJacobiansStructure(jac_irow, jac_icol);
    optim.computeCombinedSparseJacobiansValues(jac_values);

    // convert to Eigen type for convenient comparison
    std::vector<Eigen::Triplet<double>> triplet_list;
    corbo::convert_triplet(jac_irow, jac_icol, jac_values, triplet_list);
    Eigen::SparseMatrix<double> sparse_hessian(dim_obj_lsq + dim_eq + dim_ineq, optim.getParameterDimension());
    sparse_hessian.setFromTriplets(triplet_list.begin(), triplet_list.end());
    EXPECT_EQ_MATRIX(sparse_hessian, sol3, 1e-5);

    // now with multipliers
    double multiplier_obj_lsq = 5;
    Eigen::Vector2d multiplier_eq(2, 3);
    Eigen::Vector2d multiplier_ineq(3, 4);

    sol3.row(0) *= multiplier_obj_lsq;
    sol3.row(1) *= multiplier_eq(0);
    sol3.row(2) *= multiplier_eq(1);
    sol3.row(3) *= multiplier_ineq(0);
    sol3.row(4) *= multiplier_ineq(1);

    jac_values.setZero();
    optim.computeCombinedSparseJacobiansValues(jac_values, true, true, true, &multiplier_obj_lsq, multiplier_eq.data(), multiplier_ineq.data());
    corbo::convert_triplet(jac_irow, jac_icol, jac_values, triplet_list);
    sparse_hessian.setZero();
    sparse_hessian.setFromTriplets(triplet_list.begin(), triplet_list.end());
    EXPECT_EQ_MATRIX(sparse_hessian, sol3, 1e-5);
}
