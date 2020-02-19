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

#include <corbo-core/console.h>

#include "gtest/gtest.h"

#ifdef OSQP

#include <corbo-optimization/simple_optimization_problem.h>
#include <corbo-optimization/solver/qp_solver_osqp.h>

#include <corbo-core/macros.h>
#include <corbo-core/utilities.h>
#include <corbo-core/value_comparison.h>

#include <array>
#include <functional>
#include <memory>

// using corbo::SimpleOptimizationProblemWithCallbacks;
// using corbo::SolverOsqp;

class TestSolverOsqp : public testing::Test
{
 protected:
    // You can do set-up work for each test here.
    TestSolverOsqp() {}
    // You can do clean-up work that doesn't throw exceptions here.
    virtual ~TestSolverOsqp() {}
    // If the constructor and destructor are not enough for setting up
    // and cleaning up each test, you can define the following methods:

    // Code here will be called immediately after the constructor (right
    // before each test).
    void SetUp() override
    {
        // configure solver
        // EXPECT_TRUE(solver.initialize(&optim));  // we need to initialize the sover before setting parameters
    }
    // Code here will be called immediately after each test (right
    // before the destructor).
    // virtual void TearDown()

    // SimpleOptimizationProblemWithCallbacks optim;

    // SolverOsqp solver;
};
TEST_F(TestSolverOsqp, osqp_demo)
{
    using VectorXl = Eigen::Matrix<c_int, Eigen::Dynamic, 1>;

    // Load problem data (sparse format)
    // 0.5*x^T * P * x + q*x  s.t. l <= A <= u  with P: n x n and A: n x m
    c_int n = 2;
    c_int m = 3;

    c_int P_nnz = 3;
    Eigen::Vector3d P_x;  // value vector
    P_x << 4, 1, 2;
    VectorXl P_i(3);  // row indices
    P_i << 0, 0, 1;
    VectorXl P_p(3);  // column pointers
    P_p << 0, 1, 3;
    csc* P = csc_matrix(n, n, P_nnz, P_x.data(), P_i.data(), P_p.data());

    Eigen::Vector2d q(1, 1);

    c_int A_nnz = 4;
    Eigen::Vector4d A_x;
    A_x << 1, 1, 1, 1;
    VectorXl A_i(4);
    A_i << 0, 1, 0, 2;
    VectorXl A_p(3);
    A_p << 0, 2, 4;
    csc* A = csc_matrix(m, n, A_nnz, A_x.data(), A_i.data(), A_p.data());

    Eigen::Vector3d l(1, 0, 0);
    Eigen::Vector3d u(1, 0.7, 0.7);

    // Problem settings
    std::unique_ptr<OSQPSettings> settings = std::unique_ptr<OSQPSettings>(new OSQPSettings);

    // Populate data
    std::unique_ptr<OSQPData> data =  std::unique_ptr<OSQPData>(new OSQPData);
    data->n        = n;
    data->m        = m;
    data->P        = P;
    data->q        = q.data();
    data->A        = A;
    data->l        = l.data();
    data->u        = u.data();

    // Define Solver settings as default
    osqp_set_default_settings(settings.get());
    settings->verbose = 0;
    settings->alpha   = 1.0;  // Change alpha parameter

    // Setup workspace
    OSQPWorkspace* work;
    c_int exitflag = osqp_setup(&work, data.get(), settings.get());
    EXPECT_EQ(exitflag, 0);

    // Solve Problem
    osqp_solve(work);

    // Access solution
    const Eigen::Map<const Eigen::VectorXd> x_sol(work->solution->x, n);

    EXPECT_NEAR(x_sol[0], 0.3, 1e-2);
    EXPECT_NEAR(x_sol[1], 0.7, 1e-2);

    // Cleanup
    osqp_cleanup(work);
    c_free(A);
    c_free(P);
}

TEST_F(TestSolverOsqp, osqp_demo_eigen_sparse_interface)
{
    using VectorXl = Eigen::Matrix<c_int, Eigen::Dynamic, 1>;

    // Load problem data (sparse format)
    // 0.5*x^T * P * x + q*x  s.t. l <= A <= u  with P: n x n and A: n x m

    c_int n = 2;
    c_int m = 3;

    Eigen::MatrixXd P(n, n);
    P << 4, 1, 1, 2;

    Eigen::VectorXd q(n);
    q << 1, 1;

    Eigen::MatrixXd A(m, n);
    A << 1, 1, 1, 0, 0, 1;

    Eigen::VectorXd l(m);
    l << 1, 0, 0;
    Eigen::VectorXd u(m);
    u << 1, 0.7, 0.7;

    Eigen::SparseMatrix<double, Eigen::ColMajor, c_int> P_sp;
    P_sp = P.sparseView();
    P_sp = P_sp.triangularView<Eigen::Upper>();
    P_sp.makeCompressed();

    Eigen::SparseMatrix<double, Eigen::ColMajor, c_int> A_sp;
    A_sp = A.sparseView();
    A_sp.makeCompressed();

    // just for testing:
    Eigen::Map<Eigen::VectorXd> P_v(P_sp.valuePtr(), P_sp.nonZeros());
    EXPECT_EQ(P_v.size(), 3);
    EXPECT_EQ(P_v[0], 4);
    EXPECT_EQ(P_v[1], 1);
    EXPECT_EQ(P_v[2], 2);
    Eigen::Map<VectorXl> P_i(P_sp.innerIndexPtr(), P_sp.nonZeros());
    EXPECT_EQ(P_i.size(), 3);
    EXPECT_EQ(P_i[0], 0);
    EXPECT_EQ(P_i[1], 0);
    EXPECT_EQ(P_i[2], 1);
    Eigen::Map<VectorXl> P_p(P_sp.outerIndexPtr(), P_sp.outerSize() + 1);
    EXPECT_EQ(P_p.size(), 3);
    EXPECT_EQ(P_p[0], 0);
    EXPECT_EQ(P_p[1], 1);
    EXPECT_EQ(P_p[2], 3);

    Eigen::Map<Eigen::VectorXd> A_v(A_sp.valuePtr(), A_sp.nonZeros());
    EXPECT_EQ(A_v.size(), 4);
    EXPECT_EQ(A_v[0], 1);
    EXPECT_EQ(A_v[1], 1);
    EXPECT_EQ(A_v[2], 1);
    EXPECT_EQ(A_v[3], 1);
    Eigen::Map<VectorXl> A_i(A_sp.innerIndexPtr(), A_sp.nonZeros());
    EXPECT_EQ(A_i.size(), 4);
    EXPECT_EQ(A_i[0], 0);
    EXPECT_EQ(A_i[1], 1);
    EXPECT_EQ(A_i[2], 0);
    EXPECT_EQ(A_i[3], 2);
    Eigen::Map<VectorXl> A_p(A_sp.outerIndexPtr(), A_sp.outerSize() + 1);
    EXPECT_EQ(A_p.size(), 3);
    EXPECT_EQ(A_p[0], 0);
    EXPECT_EQ(A_p[1], 2);
    EXPECT_EQ(A_p[2], 4);

    csc* P_csc = csc_matrix(n, n, P_sp.nonZeros(), P_sp.valuePtr(), P_sp.innerIndexPtr(), P_sp.outerIndexPtr());

    csc* A_csc = csc_matrix(m, n, A_sp.nonZeros(), A_sp.valuePtr(), A_sp.innerIndexPtr(), A_sp.outerIndexPtr());

    // Problem settings
    std::unique_ptr<OSQPSettings> settings = std::unique_ptr<OSQPSettings>(new OSQPSettings);

    // Populate data
    std::unique_ptr<OSQPData> data =  std::unique_ptr<OSQPData>(new OSQPData);
    data->n        = n;
    data->m        = m;
    data->P        = P_csc;
    data->q        = q.data();
    data->A        = A_csc;
    data->l        = l.data();
    data->u        = u.data();

    // Define Solver settings as default
    osqp_set_default_settings(settings.get());
    settings->verbose = 0;
    settings->alpha   = 1.0;  // Change alpha parameter

    // Setup workspace
    OSQPWorkspace* work;
    c_int exitflag = osqp_setup(&work, data.get(), settings.get());
    EXPECT_EQ(exitflag, 0);

    // Solve Problem
    osqp_solve(work);

    // Access solution
    const Eigen::Map<const Eigen::VectorXd> x_sol(work->solution->x, n);

    EXPECT_NEAR(x_sol[0], 0.3, 1e-4);
    EXPECT_NEAR(x_sol[1], 0.7, 1e-4);

    // Cleanup
    osqp_cleanup(work);
    c_free(P_csc);
    c_free(A_csc);
}

#else  // OSQP

class TestSolverOsqp : public testing::Test
{
 protected:
    // You can do set-up work for each test here.
    TestSolverOsqp() {}
    // You can do clean-up work that doesn't throw exceptions here.
    virtual ~TestSolverOsqp() {}
};

TEST_F(TestSolverOsqp, osqp_not_found) { PRINT_WARNING("Skipping OSQP tests, since OSQP is not found."); }

#endif  // OSQP
