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

#include <corbo-optimization/hyper_graph/hyper_graph_optimization_problem_vertex_based.h>

#include <corbo-optimization/hyper_graph/generic_edge.h>
#include <corbo-optimization/hyper_graph/vector_vertex.h>
#include <corbo-optimization/misc.h>

#include <corbo-core/console.h>
#include <corbo-core/macros.h>
#include <corbo-core/types.h>
#include <corbo-core/utilities.h>
#include <corbo-core/value_comparison.h>

#include <array>
#include <functional>

#include "gtest/gtest.h"

using HyperGraphOptimizationProblem = corbo::HyperGraphOptimizationProblemVertexBased;
using corbo::VertexSet;
using corbo::OptimizationEdgeSet;
using pFVectorVertex = corbo::PartiallyFixedVectorVertex;
using corbo::EdgeGenericScalarFun;
using corbo::EdgeGenericVectorFun;
using corbo::BaseEdge;
using corbo::BaseMixedEdge;
using corbo::MixedEdgeGenericVectorFun;
using corbo::CORBO_INF_DBL;

class TestHyperGraphOptimizationProblemVertexBased : public testing::Test
{
 protected:
    // You can do set-up work for each test here.
    TestHyperGraphOptimizationProblemVertexBased()
    {
        vertices = std::make_shared<VertexSet>();
        edges    = std::make_shared<OptimizationEdgeSet>();
    }
    // You can do clean-up work that doesn't throw exceptions here.
    virtual ~TestHyperGraphOptimizationProblemVertexBased() {}
    // If the constructor and destructor are not enough for setting up
    // and cleaning up each test, you can define the following methods:

    // Code here will be called immediately after the constructor (right
    // before each test).
    void SetUp() override
    {
        // setup vertices
        v1 = std::make_shared<pFVectorVertex>(1);
        v1->values().setOnes();
        v2 = std::make_shared<pFVectorVertex>(Eigen::Vector2d(1, 1));

        optim.getGraph().setVertexSet(vertices);
        optim.getGraph().setEdgeSet(edges);
    }
    // Code here will be called immediately after each test (right
    // before the destructor).
    // virtual void TearDown()

    void computeSparseJacobianFiniteCombinedBoundsViaTriplets(int nnz, int n, int m, Eigen::SparseMatrix<double>& jacobian, double weight = 1.0)
    {
        jacobian.resize(n, m);
        jacobian.setZero();
        Eigen::VectorXi i_row(nnz), j_col(nnz);
        Eigen::VectorXd jac_values(nnz);
        optim.computeSparseJacobianFiniteCombinedBoundsStructure(i_row, j_col);
        optim.computeSparseJacobianFiniteCombinedBoundsValues(jac_values, weight);
        std::vector<Eigen::Triplet<double>> triplets;
        corbo::convert_triplet(i_row, j_col, jac_values, triplets);
        jacobian.setFromTriplets(triplets.begin(), triplets.end());
    }

    void computeSparseJacobianObjectiveViaTriplets(int nnz, int n, int m, Eigen::SparseMatrix<double>& jacobian, const double* multipliers = nullptr)
    {
        jacobian.resize(n, m);
        jacobian.setZero();
        Eigen::VectorXi i_row(nnz), j_col(nnz);
        Eigen::VectorXd jac_values(nnz);
        optim.computeSparseJacobianLsqObjectiveStructure(i_row, j_col);
        optim.computeSparseJacobianLsqObjectiveValues(jac_values, multipliers);
        std::vector<Eigen::Triplet<double>> triplets;
        corbo::convert_triplet(i_row, j_col, jac_values, triplets);
        jacobian.setFromTriplets(triplets.begin(), triplets.end());
    }

    void computeSparseJacobianEqualitiesViaTriplets(int nnz, int n, int m, Eigen::SparseMatrix<double>& jacobian, const double* multipliers = nullptr)
    {
        jacobian.resize(n, m);
        jacobian.setZero();
        Eigen::VectorXi i_row(nnz), j_col(nnz);
        Eigen::VectorXd jac_values(nnz);
        optim.computeSparseJacobianEqualitiesStructure(i_row, j_col);
        optim.computeSparseJacobianEqualitiesValues(jac_values, multipliers);
        std::vector<Eigen::Triplet<double>> triplets;
        corbo::convert_triplet(i_row, j_col, jac_values, triplets);
        jacobian.setFromTriplets(triplets.begin(), triplets.end());
    }

    void computeSparseJacobianInequalitiesViaTriplets(int nnz, int n, int m, Eigen::SparseMatrix<double>& jacobian,
                                                      const double* multipliers = nullptr)
    {
        jacobian.resize(n, m);
        jacobian.setZero();
        Eigen::VectorXi i_row(nnz), j_col(nnz);
        Eigen::VectorXd jac_values(nnz);
        optim.computeSparseJacobianInequalitiesStructure(i_row, j_col);
        optim.computeSparseJacobianInequalitiesValues(jac_values, multipliers);
        std::vector<Eigen::Triplet<double>> triplets;
        corbo::convert_triplet(i_row, j_col, jac_values, triplets);
        jacobian.setFromTriplets(triplets.begin(), triplets.end());
    }

    void computeSparseJacobianActiveInequalitiesViaTriplets(int nnz, int n, int m, Eigen::SparseMatrix<double>& jacobian, double weight = 1.0)
    {
        jacobian.resize(n, m);
        jacobian.setZero();
        Eigen::VectorXi i_row(nnz), j_col(nnz);
        Eigen::VectorXd jac_values(nnz);
        optim.computeSparseJacobianInequalitiesStructure(i_row, j_col);
        optim.computeSparseJacobianActiveInequalitiesValues(jac_values, weight);
        std::vector<Eigen::Triplet<double>> triplets;
        corbo::convert_triplet(i_row, j_col, jac_values, triplets);
        jacobian.setFromTriplets(triplets.begin(), triplets.end());
    }

    void computeSparseJacobiansViaTriplets(int nnz_obj, int n_obj, int m_obj, Eigen::SparseMatrix<double>& jacobian_obj, int nnz_eq, int n_eq,
                                           int m_eq, Eigen::SparseMatrix<double>& jacobian_eq, int nnz_ineq, int n_ineq, int m_ineq,
                                           Eigen::SparseMatrix<double>& jacobian_ineq, const double* multipliers_obj = nullptr,
                                           const double* multipliers_eq = nullptr, const double* multipliers_ineq = nullptr, bool active_ineq = false,
                                           double active_ineq_weight = 1.0)
    {
        jacobian_obj.resize(n_obj, m_obj);
        jacobian_eq.resize(n_eq, m_eq);
        jacobian_ineq.resize(n_ineq, m_ineq);
        jacobian_obj.setZero();
        jacobian_eq.setZero();
        jacobian_ineq.setZero();
        Eigen::VectorXi i_row_obj(nnz_obj), j_col_obj(nnz_obj);
        Eigen::VectorXi i_row_eq(nnz_eq), j_col_eq(nnz_eq);
        Eigen::VectorXi i_row_ineq(nnz_ineq), j_col_ineq(nnz_ineq);
        Eigen::VectorXd jac_values_obj(nnz_obj);
        Eigen::VectorXd jac_values_eq(nnz_eq);
        Eigen::VectorXd jac_values_ineq(nnz_ineq);
        optim.computeSparseJacobiansStructure(i_row_obj, j_col_obj, i_row_eq, j_col_eq, i_row_ineq, j_col_ineq);
        optim.computeSparseJacobiansValues(jac_values_obj, jac_values_eq, jac_values_ineq, multipliers_obj, multipliers_eq, multipliers_ineq,
                                           active_ineq, active_ineq_weight);
        std::vector<Eigen::Triplet<double>> triplets_obj, triplets_eq, triplets_ineq;
        corbo::convert_triplet(i_row_obj, j_col_obj, jac_values_obj, triplets_obj);
        corbo::convert_triplet(i_row_eq, j_col_eq, jac_values_eq, triplets_eq);
        corbo::convert_triplet(i_row_ineq, j_col_ineq, jac_values_ineq, triplets_ineq);
        if (!triplets_obj.empty()) jacobian_obj.setFromTriplets(triplets_obj.begin(), triplets_obj.end());
        if (!triplets_eq.empty()) jacobian_eq.setFromTriplets(triplets_eq.begin(), triplets_eq.end());
        if (!triplets_ineq.empty()) jacobian_ineq.setFromTriplets(triplets_ineq.begin(), triplets_ineq.end());
    }

    void computeSparseHessianObjectiveViaTriplets(int nnz, int n, Eigen::SparseMatrix<double>& hessian, double multiplier = 1.0,
                                                  bool lower_part_only = false)
    {
        hessian.resize(n, n);
        hessian.setZero();
        Eigen::VectorXi i_row(nnz), j_col(nnz);
        Eigen::VectorXd hessian_values(nnz);
        optim.computeSparseHessianObjectiveStructure(i_row, j_col, lower_part_only);
        optim.computeSparseHessianObjectiveValues(hessian_values, multiplier, lower_part_only);
        std::vector<Eigen::Triplet<double>> triplets;
        corbo::convert_triplet(i_row, j_col, hessian_values, triplets);
        hessian.setFromTriplets(triplets.begin(), triplets.end());
    }

    void computeSparseHessianEqualitiesViaTriplets(int nnz, int n, Eigen::SparseMatrix<double>& hessian, const double* multipliers = nullptr,
                                                   bool lower_part_only = false)
    {
        hessian.resize(n, n);
        hessian.setZero();
        Eigen::VectorXi i_row(nnz), j_col(nnz);
        Eigen::VectorXd hessian_values(nnz);
        optim.computeSparseHessianEqualitiesStructure(i_row, j_col, lower_part_only);
        optim.computeSparseHessianEqualitiesValues(hessian_values, multipliers, lower_part_only);
        std::vector<Eigen::Triplet<double>> triplets;
        corbo::convert_triplet(i_row, j_col, hessian_values, triplets);
        hessian.setFromTriplets(triplets.begin(), triplets.end());
    }

    void computeSparseHessianInequalitiesViaTriplets(int nnz, int n, Eigen::SparseMatrix<double>& hessian, const double* multipliers = nullptr,
                                                     bool lower_part_only = false)
    {
        hessian.resize(n, n);
        hessian.setZero();
        Eigen::VectorXi i_row(nnz), j_col(nnz);
        Eigen::VectorXd hessian_values(nnz);
        optim.computeSparseHessianInequalitiesStructure(i_row, j_col, lower_part_only);
        optim.computeSparseHessianInequalitiesValues(hessian_values, multipliers, lower_part_only);
        std::vector<Eigen::Triplet<double>> triplets;
        corbo::convert_triplet(i_row, j_col, hessian_values, triplets);
        hessian.setFromTriplets(triplets.begin(), triplets.end());
    }

    void computeSparseHessiansViaTriplets(int nnz_obj, int n_obj, Eigen::SparseMatrix<double>& hessian_obj, int nnz_eq, int n_eq,
                                          Eigen::SparseMatrix<double>& hessian_eq, int nnz_ineq, int n_ineq,
                                          Eigen::SparseMatrix<double>& hessian_ineq, double multiplier_obj = 1.0,
                                          const double* multipliers_eq = nullptr, const double* multipliers_ineq = nullptr,
                                          bool lower_part_only = false)
    {
        hessian_obj.resize(n_obj, n_obj);
        hessian_eq.resize(n_eq, n_eq);
        hessian_ineq.resize(n_ineq, n_ineq);
        hessian_obj.setZero();
        hessian_eq.setZero();
        hessian_ineq.setZero();
        Eigen::VectorXi i_row_obj(nnz_obj), j_col_obj(nnz_obj);
        Eigen::VectorXi i_row_eq(nnz_eq), j_col_eq(nnz_eq);
        Eigen::VectorXi i_row_ineq(nnz_ineq), j_col_ineq(nnz_ineq);
        Eigen::VectorXd hessian_values_obj(nnz_obj);
        Eigen::VectorXd hessian_values_eq(nnz_eq);
        Eigen::VectorXd hessian_values_ineq(nnz_ineq);
        optim.computeSparseHessiansStructure(i_row_obj, j_col_obj, i_row_eq, j_col_eq, i_row_ineq, j_col_ineq, lower_part_only);
        optim.computeSparseHessiansValues(hessian_values_obj, hessian_values_eq, hessian_values_ineq, multiplier_obj, multipliers_eq,
                                          multipliers_ineq, lower_part_only);
        std::vector<Eigen::Triplet<double>> triplets_obj, triplets_eq, triplets_ineq;
        corbo::convert_triplet(i_row_obj, j_col_obj, hessian_values_obj, triplets_obj);
        corbo::convert_triplet(i_row_eq, j_col_eq, hessian_values_eq, triplets_eq);
        corbo::convert_triplet(i_row_ineq, j_col_ineq, hessian_values_ineq, triplets_ineq);
        if (!triplets_obj.empty()) hessian_obj.setFromTriplets(triplets_obj.begin(), triplets_obj.end());
        if (!triplets_eq.empty()) hessian_eq.setFromTriplets(triplets_eq.begin(), triplets_eq.end());
        if (!triplets_ineq.empty()) hessian_ineq.setFromTriplets(triplets_ineq.begin(), triplets_ineq.end());
    }

    //    void testObjectiveAndEqualityAndInequalityJacobians(const Eigen::MatrixXd& jacobian_sol, int nnz_min, double tol,
    //                                                        const double* multipliers = nullptr)
    //    {
    //        testObjectiveAndEqualityAndInequalityJacobians(jacobian_sol, jacobian_sol, jacobian_sol, nnz_min, nnz_min, nnz_min, tol, multipliers,
    //                                                       multipliers, multipliers);
    //    }

    void testObjectiveAndEqualityAndInequalityJacobians(const Eigen::VectorXd& gradient_obj_sol, const Eigen::VectorXd& gradient_non_lsq_obj_sol,
                                                        const Eigen::MatrixXd& jacobian_lsq_obj_sol, const Eigen::MatrixXd& jacobian_eq_sol,
                                                        const Eigen::MatrixXd& jacobian_ineq_sol, int nnz_min_obj, int nnz_min_eq, int nnz_min_ineq,
                                                        double tol, const double* multipliers_lsq_obj = nullptr,
                                                        const double* multipliers_eq = nullptr, const double* multipliers_ineq = nullptr)
    {
        bool has_gradient_obj     = gradient_obj_sol.size() > 0;
        bool has_gradient_non_lsq = gradient_non_lsq_obj_sol.size() > 0;
        bool has_jacobian_lsq_obj = jacobian_lsq_obj_sol.rows() > 0;
        bool has_jacobian_eq      = jacobian_eq_sol.rows() > 0;
        bool has_jacobian_ineq    = jacobian_ineq_sol.rows() > 0;

        // compute objective gradient
        Eigen::VectorXd gradient_non_lsq_obj = Eigen::VectorXd::Zero(optim.getParameterDimension());
        if (has_gradient_non_lsq)
        {
            optim.computeGradientNonLsqObjective(gradient_non_lsq_obj);
            EXPECT_EQ_MATRIX(gradient_non_lsq_obj, gradient_non_lsq_obj_sol, tol);
        }
        Eigen::VectorXd gradient_obj = Eigen::VectorXd::Zero(optim.getParameterDimension());
        if (has_gradient_obj)
        {
            optim.computeGradientObjective(gradient_obj);
            EXPECT_EQ_MATRIX(gradient_obj, gradient_obj_sol, tol);
        }

        // compute dense jacobian
        Eigen::MatrixXd jacobian_lsq_obj = Eigen::MatrixXd::Zero(optim.getLsqObjectiveDimension(), optim.getParameterDimension());
        if (has_jacobian_lsq_obj)
        {
            optim.computeDenseJacobianLsqObjective(jacobian_lsq_obj, multipliers_lsq_obj);
            EXPECT_EQ_MATRIX(jacobian_lsq_obj, jacobian_lsq_obj_sol, tol);
        }

        Eigen::MatrixXd jacobian_eq = Eigen::MatrixXd::Zero(optim.getEqualityDimension(), optim.getParameterDimension());
        if (has_jacobian_eq)
        {
            optim.computeDenseJacobianEqualities(jacobian_eq, multipliers_eq);
            EXPECT_EQ_MATRIX(jacobian_eq, jacobian_eq_sol, tol);
        }

        Eigen::MatrixXd jacobian_ineq = Eigen::MatrixXd::Zero(optim.getInequalityDimension(), optim.getParameterDimension());
        if (has_jacobian_ineq)
        {
            optim.computeDenseJacobianInequalities(jacobian_ineq, multipliers_ineq);
            EXPECT_EQ_MATRIX(jacobian_ineq, jacobian_ineq_sol, tol);
        }
        // combined dense computation
        jacobian_lsq_obj.setZero();
        gradient_non_lsq_obj.setZero();
        jacobian_eq.setZero();
        jacobian_ineq.setZero();

        optim.computeDenseJacobians(gradient_non_lsq_obj, jacobian_lsq_obj, jacobian_eq, jacobian_ineq, multipliers_lsq_obj, multipliers_eq,
                                    multipliers_ineq);
        if (has_gradient_non_lsq) EXPECT_EQ_MATRIX(gradient_non_lsq_obj, gradient_non_lsq_obj_sol, tol);
        if (has_jacobian_lsq_obj) EXPECT_EQ_MATRIX(jacobian_lsq_obj, jacobian_lsq_obj_sol, tol);
        if (has_jacobian_eq) EXPECT_EQ_MATRIX(jacobian_eq, jacobian_eq_sol, tol);
        if (has_jacobian_ineq) EXPECT_EQ_MATRIX(jacobian_ineq, jacobian_ineq_sol, tol);

        // compute sparse jacobian
        Eigen::SparseMatrix<double> sparse_jacobian_lsq_obj(optim.getLsqObjectiveDimension(), optim.getParameterDimension());
        Eigen::SparseMatrix<double> sparse_jacobian_eq(optim.getEqualityDimension(), optim.getParameterDimension());
        Eigen::SparseMatrix<double> sparse_jacobian_ineq(optim.getInequalityDimension(), optim.getParameterDimension());

        if (has_jacobian_lsq_obj)
        {
            optim.computeSparseJacobianLsqObjective(sparse_jacobian_lsq_obj, multipliers_lsq_obj);
            EXPECT_EQ_MATRIX(sparse_jacobian_lsq_obj, jacobian_lsq_obj_sol, tol);
        }
        if (has_jacobian_eq)
        {
            optim.computeSparseJacobianEqualities(sparse_jacobian_eq, multipliers_eq);
            EXPECT_EQ_MATRIX(sparse_jacobian_eq, jacobian_eq_sol, tol);
        }
        if (has_jacobian_ineq)
        {
            optim.computeSparseJacobianInequalities(sparse_jacobian_ineq, multipliers_ineq);
            EXPECT_EQ_MATRIX(sparse_jacobian_ineq, jacobian_ineq_sol, tol);
        }

        // combined sparse computation
        sparse_jacobian_lsq_obj.setZero();
        sparse_jacobian_eq.setZero();
        sparse_jacobian_ineq.setZero();
        optim.computeSparseJacobians(sparse_jacobian_lsq_obj, sparse_jacobian_eq, sparse_jacobian_ineq, multipliers_lsq_obj, multipliers_eq,
                                     multipliers_ineq);
        if (has_jacobian_lsq_obj) EXPECT_EQ_MATRIX(sparse_jacobian_lsq_obj, jacobian_lsq_obj_sol, tol);
        if (has_jacobian_eq) EXPECT_EQ_MATRIX(sparse_jacobian_eq, jacobian_eq_sol, tol);
        if (has_jacobian_ineq) EXPECT_EQ_MATRIX(sparse_jacobian_ineq, jacobian_ineq_sol, tol);

        // compute other sparse jacobian
        int nnz_lsq_obj = optim.computeSparseJacobianLsqObjectiveNNZ();
        if (has_jacobian_lsq_obj)
        {
            EXPECT_LE(nnz_lsq_obj, jacobian_lsq_obj_sol.rows() * jacobian_lsq_obj_sol.cols());
            EXPECT_GE(nnz_lsq_obj, nnz_min_obj);

            computeSparseJacobianObjectiveViaTriplets(nnz_lsq_obj, jacobian_lsq_obj.rows(), jacobian_lsq_obj.cols(), sparse_jacobian_lsq_obj,
                                                      multipliers_lsq_obj);
            EXPECT_EQ_MATRIX(sparse_jacobian_lsq_obj, jacobian_lsq_obj_sol, tol);
        }
        int nnz_eq = optim.computeSparseJacobianEqualitiesNNZ();
        if (has_jacobian_eq)
        {
            EXPECT_LE(nnz_eq, jacobian_eq_sol.rows() * jacobian_eq_sol.cols());
            EXPECT_GE(nnz_eq, nnz_min_eq);

            computeSparseJacobianEqualitiesViaTriplets(nnz_eq, jacobian_eq.rows(), jacobian_eq.cols(), sparse_jacobian_eq, multipliers_eq);
            EXPECT_EQ_MATRIX(sparse_jacobian_eq, jacobian_eq_sol, tol);
        }
        int nnz_ineq = optim.computeSparseJacobianInequalitiesNNZ();
        if (has_jacobian_ineq)
        {
            EXPECT_LE(nnz_ineq, jacobian_ineq_sol.rows() * jacobian_ineq_sol.cols());
            EXPECT_GE(nnz_ineq, nnz_min_ineq);

            computeSparseJacobianInequalitiesViaTriplets(nnz_ineq, jacobian_ineq.rows(), jacobian_ineq.cols(), sparse_jacobian_ineq,
                                                         multipliers_ineq);
            EXPECT_EQ_MATRIX(sparse_jacobian_ineq, jacobian_ineq_sol, tol);
        }

        // combined sparse computation
        optim.computeSparseJacobiansNNZ(nnz_lsq_obj, nnz_eq, nnz_ineq);
        computeSparseJacobiansViaTriplets(nnz_lsq_obj, jacobian_lsq_obj.rows(), jacobian_lsq_obj.cols(), sparse_jacobian_lsq_obj, nnz_eq,
                                          jacobian_eq.rows(), jacobian_eq.cols(), sparse_jacobian_eq, nnz_ineq, jacobian_ineq.rows(),
                                          jacobian_ineq.cols(), sparse_jacobian_ineq, multipliers_lsq_obj, multipliers_eq, multipliers_ineq);
        if (has_jacobian_lsq_obj) EXPECT_EQ_MATRIX(sparse_jacobian_lsq_obj, jacobian_lsq_obj_sol, tol);
        if (has_jacobian_eq) EXPECT_EQ_MATRIX(sparse_jacobian_eq, jacobian_eq_sol, tol);
        if (has_jacobian_ineq) EXPECT_EQ_MATRIX(sparse_jacobian_ineq, jacobian_ineq_sol, tol);
    }

    void testCombinedSparseJacobian(const Eigen::MatrixXd* jacobian_lsq_obj_sol, const Eigen::MatrixXd* jacobian_eq_sol,
                                    const Eigen::MatrixXd* jacobian_ineq_sol, double tol, const Eigen::MatrixXd* finite_combined_bounds = nullptr,
                                    bool active_ineq = false, double weight_eq = 1.0, double weight_ineq = 1.0, double weight_bounds = 1.0)
    {
        // create composed solution
        int dim_sol = 0;
        if (jacobian_lsq_obj_sol) dim_sol += jacobian_lsq_obj_sol->rows();
        if (jacobian_eq_sol) dim_sol += jacobian_eq_sol->rows();
        if (jacobian_ineq_sol) dim_sol += jacobian_ineq_sol->rows();
        if (finite_combined_bounds) dim_sol += finite_combined_bounds->rows();

        Eigen::MatrixXd combined_sol(dim_sol, optim.getParameterDimension());
        int row_idx = 0;
        if (jacobian_lsq_obj_sol)
        {
            combined_sol.topRows(jacobian_lsq_obj_sol->rows()) = *jacobian_lsq_obj_sol;
            row_idx += jacobian_lsq_obj_sol->rows();
        }
        if (jacobian_eq_sol)
        {
            combined_sol.middleRows(row_idx, jacobian_eq_sol->rows()) = *jacobian_eq_sol;
            row_idx += jacobian_eq_sol->rows();
        }
        if (jacobian_ineq_sol)
        {
            combined_sol.middleRows(row_idx, jacobian_ineq_sol->rows()) = *jacobian_ineq_sol;
            row_idx += jacobian_ineq_sol->rows();
        }
        if (finite_combined_bounds)
        {
            combined_sol.bottomRows(finite_combined_bounds->rows()) = *finite_combined_bounds;
        }

        // single matrix sparse jacobian
        int combined_dim = (jacobian_lsq_obj_sol != 0 ? optim.getLsqObjectiveDimension() : 0) +
                           (jacobian_eq_sol != 0 ? optim.getEqualityDimension() : 0) + (jacobian_ineq_sol != 0 ? optim.getInequalityDimension() : 0);
        if (finite_combined_bounds) combined_dim += optim.finiteCombinedBoundsDimension();
        Eigen::SparseMatrix<double> combined_sparse_jacobian(combined_dim, optim.getParameterDimension());
        optim.computeCombinedSparseJacobian(combined_sparse_jacobian, (bool)jacobian_lsq_obj_sol, (bool)jacobian_eq_sol, (bool)jacobian_ineq_sol,
                                            (bool)finite_combined_bounds, active_ineq, weight_eq, weight_ineq, weight_bounds);
        EXPECT_EQ_MATRIX(combined_sparse_jacobian, combined_sol, tol);

        // triplet based computation
        if (!active_ineq)  // TODO(roesmann) currently not supported by sparse computation
        {
            int jac_nnz = optim.computeCombinedSparseJacobiansNNZ(jacobian_lsq_obj_sol, jacobian_eq_sol, jacobian_ineq_sol);

            Eigen::VectorXd multplier_eq   = Eigen::VectorXd::Constant(optim.getEqualityDimension(), weight_eq);
            Eigen::VectorXd multplier_ineq = Eigen::VectorXd::Constant(optim.getInequalityDimension(), weight_ineq);

            Eigen::VectorXi jac_irow(jac_nnz);
            Eigen::VectorXi jac_icol(jac_nnz);
            Eigen::VectorXd jac_values(jac_nnz);
            optim.computeCombinedSparseJacobiansStructure(jac_irow, jac_icol, jacobian_lsq_obj_sol, jacobian_eq_sol, jacobian_ineq_sol);
            optim.computeCombinedSparseJacobiansValues(jac_values, jacobian_lsq_obj_sol, jacobian_eq_sol, jacobian_ineq_sol, nullptr,
                                                       weight_eq == 1.0 ? nullptr : multplier_eq.data(),
                                                       weight_ineq == 1.0 ? nullptr : multplier_ineq.data());

            dim_sol = 0;
            if (jacobian_lsq_obj_sol) dim_sol += jacobian_lsq_obj_sol->rows();
            if (jacobian_eq_sol) dim_sol += jacobian_eq_sol->rows();
            if (jacobian_ineq_sol) dim_sol += jacobian_ineq_sol->rows();

            std::vector<Eigen::Triplet<double>> triplet_list;
            corbo::convert_triplet(jac_irow, jac_icol, jac_values, triplet_list);
            Eigen::SparseMatrix<double> sparse_hessian(dim_sol, optim.getParameterDimension());
            sparse_hessian.setFromTriplets(triplet_list.begin(), triplet_list.end());
            EXPECT_EQ_MATRIX(sparse_hessian, combined_sol.topRows(dim_sol), 1e-5);
        }
    }

    //    void testObjectiveAndEqualityAndInequalityHessians(const Eigen::MatrixXd& hessian_sol, int dim_values, int nnz_min, double tol,
    //                                                       const double* multipliers = nullptr)
    //    {
    //        testObjectiveAndEqualityAndInequalityHessians(hessian_sol, hessian_sol, hessian_sol, dim_values, dim_values, dim_values, nnz_min,
    //        nnz_min,
    //                                                      nnz_min, tol, multipliers, multipliers, multipliers);
    //    }

    void testObjectiveAndEqualityAndInequalityHessians(const Eigen::MatrixXd& hessian_obj_sol, const Eigen::MatrixXd& hessian_eq_sol,
                                                       const Eigen::MatrixXd& hessian_ineq_sol, int dim_eq, int dim_ineq, int nnz_min_obj,
                                                       int nnz_min_eq, int nnz_min_ineq, double tol, double multiplier_obj = 1.0,
                                                       const double* multipliers_eq = nullptr, const double* multipliers_ineq = nullptr)
    {
        // compute dense jacobian
        Eigen::MatrixXd hessian_obj = Eigen::MatrixXd::Zero(optim.getParameterDimension(), optim.getParameterDimension());
        optim.computeDenseHessianObjective(hessian_obj, multiplier_obj);
        EXPECT_EQ_MATRIX(hessian_obj, hessian_obj_sol, tol);
        Eigen::MatrixXd hessian_eq = Eigen::MatrixXd::Zero(optim.getParameterDimension(), optim.getParameterDimension());
        optim.computeDenseHessianEqualities(hessian_eq, multipliers_eq);
        EXPECT_EQ_MATRIX(hessian_eq, hessian_eq_sol, tol);
        Eigen::MatrixXd hessian_ineq = Eigen::MatrixXd::Zero(optim.getParameterDimension(), optim.getParameterDimension());
        optim.computeDenseHessianInequalities(hessian_ineq, multipliers_ineq);
        EXPECT_EQ_MATRIX(hessian_ineq, hessian_ineq_sol, tol);

        // combined dense computation
        hessian_obj.setZero();
        hessian_eq.setZero();
        hessian_ineq.setZero();
        optim.computeDenseHessians(hessian_obj, hessian_eq, hessian_ineq, multiplier_obj, multipliers_eq, multipliers_ineq);
        EXPECT_EQ_MATRIX(hessian_obj, hessian_obj_sol, tol);
        EXPECT_EQ_MATRIX(hessian_eq, hessian_eq_sol, tol);
        EXPECT_EQ_MATRIX(hessian_ineq, hessian_ineq_sol, tol);

        // compute sparse hessian
        Eigen::SparseMatrix<double> sparse_hessian_obj(optim.getParameterDimension(), optim.getParameterDimension());
        optim.computeSparseHessianObjective(sparse_hessian_obj, multiplier_obj);
        EXPECT_EQ_MATRIX(sparse_hessian_obj, hessian_obj_sol, tol);
        Eigen::SparseMatrix<double> sparse_hessian_eq(optim.getParameterDimension(), optim.getParameterDimension());
        optim.computeSparseHessianEqualities(sparse_hessian_eq, multipliers_eq);
        EXPECT_EQ_MATRIX(sparse_hessian_eq, hessian_eq_sol, tol);
        Eigen::SparseMatrix<double> sparse_hessian_ineq(optim.getParameterDimension(), optim.getParameterDimension());
        optim.computeSparseHessianInequalities(sparse_hessian_ineq, multipliers_ineq);
        EXPECT_EQ_MATRIX(sparse_hessian_ineq, hessian_ineq_sol, tol);

        // combined sparse computation
        sparse_hessian_obj.setZero();
        sparse_hessian_eq.setZero();
        sparse_hessian_ineq.setZero();
        optim.computeSparseHessians(sparse_hessian_obj, sparse_hessian_eq, sparse_hessian_ineq, multiplier_obj, multipliers_eq, multipliers_ineq);
        EXPECT_EQ_MATRIX(sparse_hessian_obj, hessian_obj_sol, tol);
        EXPECT_EQ_MATRIX(sparse_hessian_eq, hessian_eq_sol, tol);
        EXPECT_EQ_MATRIX(sparse_hessian_ineq, hessian_ineq_sol, tol);

        // compute other sparse jacobian
        int nnz_obj = optim.computeSparseHessianObjectiveNNZ();
        EXPECT_GE(nnz_obj, nnz_min_obj);
        int nnz_eq = optim.computeSparseHessianEqualitiesNNZ();
        EXPECT_GE(nnz_eq, nnz_min_eq);
        int nnz_ineq = optim.computeSparseHessianInequalitiesNNZ();
        EXPECT_GE(nnz_ineq, nnz_min_ineq);
        computeSparseHessianObjectiveViaTriplets(nnz_obj, hessian_obj.rows(), sparse_hessian_obj, multiplier_obj);
        EXPECT_EQ_MATRIX(sparse_hessian_obj, hessian_obj_sol, tol);
        computeSparseHessianEqualitiesViaTriplets(nnz_eq, hessian_obj.rows(), sparse_hessian_eq, multipliers_eq);
        EXPECT_EQ_MATRIX(sparse_hessian_eq, hessian_eq_sol, tol);
        computeSparseHessianInequalitiesViaTriplets(nnz_ineq, hessian_obj.rows(), sparse_hessian_ineq, multipliers_ineq);
        EXPECT_EQ_MATRIX(sparse_hessian_ineq, hessian_ineq_sol, tol);

        // lower part only
        computeSparseHessianObjectiveViaTriplets(optim.computeSparseHessianObjectiveNNZ(true), hessian_obj.rows(), sparse_hessian_obj, multiplier_obj,
                                                 true);
        Eigen::SparseMatrix<double> hessian_obj_lower = sparse_hessian_obj.selfadjointView<Eigen::Lower>();
        EXPECT_EQ_MATRIX(hessian_obj_lower, hessian_obj_sol, tol);
        computeSparseHessianEqualitiesViaTriplets(optim.computeSparseHessianEqualitiesNNZ(true), hessian_obj.rows(), sparse_hessian_eq,
                                                  multipliers_eq, true);
        Eigen::SparseMatrix<double> hessian_eq_lower = sparse_hessian_eq.selfadjointView<Eigen::Lower>();
        EXPECT_EQ_MATRIX(hessian_eq_lower, hessian_eq_sol, tol);
        computeSparseHessianInequalitiesViaTriplets(optim.computeSparseHessianInequalitiesNNZ(true), hessian_obj.rows(), sparse_hessian_ineq,
                                                    multipliers_ineq, true);
        Eigen::SparseMatrix<double> hessian_ineq_lower = sparse_hessian_ineq.selfadjointView<Eigen::Lower>();
        EXPECT_EQ_MATRIX(hessian_ineq_lower, hessian_ineq_sol, tol);

        // combined sparse computation
        optim.computeSparseHessiansNNZ(nnz_obj, nnz_eq, nnz_ineq);
        computeSparseHessiansViaTriplets(nnz_obj, hessian_obj.rows(), sparse_hessian_obj, nnz_eq, hessian_eq.rows(), sparse_hessian_eq, nnz_ineq,
                                         hessian_ineq.rows(), sparse_hessian_ineq, multiplier_obj, multipliers_eq, multipliers_ineq);
        EXPECT_EQ_MATRIX(sparse_hessian_obj, hessian_obj_sol, tol);
        EXPECT_EQ_MATRIX(sparse_hessian_eq, hessian_eq_sol, tol);
        EXPECT_EQ_MATRIX(sparse_hessian_ineq, hessian_ineq_sol, tol);

        // lower part only
        computeSparseHessiansViaTriplets(optim.computeSparseHessianObjectiveNNZ(true), hessian_obj.rows(), sparse_hessian_obj,
                                         optim.computeSparseHessianEqualitiesNNZ(true), hessian_eq.rows(), sparse_hessian_eq,
                                         optim.computeSparseHessianInequalitiesNNZ(true), hessian_ineq.rows(), sparse_hessian_ineq, multiplier_obj,
                                         multipliers_eq, multipliers_ineq, true);
        hessian_obj_lower = sparse_hessian_obj.selfadjointView<Eigen::Lower>();
        EXPECT_EQ_MATRIX(hessian_obj_lower, hessian_obj_sol, tol);
        hessian_eq_lower = sparse_hessian_eq.selfadjointView<Eigen::Lower>();
        EXPECT_EQ_MATRIX(hessian_eq_lower, hessian_eq_sol, tol);
        hessian_ineq_lower = sparse_hessian_ineq.selfadjointView<Eigen::Lower>();
        EXPECT_EQ_MATRIX(hessian_ineq_lower, hessian_ineq_sol, tol);
    }

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

    using MixedEdge1T = MixedEdgeGenericVectorFun<pFVectorVertex, pFVectorVertex>;

    void mixed_edge1_precompute(const Edge4T::VertexContainer& vertices) { mixed_edge_aux = 2.0; }
    void mixed_edge1_obj(const MixedEdge1T::VertexContainer& vertices, Eigen::Ref<Eigen::VectorXd> values_obj)
    {
        values_obj[0] = vertices[0]->getData()[0] + 2 * vertices[1]->getData()[0] + 3 * vertices[1]->getData()[1];
        values_obj *= mixed_edge_aux;
    }
    void mixed_edge1_eq(const MixedEdge1T::VertexContainer& vertices, Eigen::Ref<Eigen::VectorXd> values_eq)
    {
        values_eq[0] = 4 * vertices[0]->getData()[0];
        values_eq[1] = 5 * vertices[1]->getData()[0] + 6 * vertices[1]->getData()[1];
        values_eq *= mixed_edge_aux;
    }
    void mixed_edge1_ineq(const MixedEdge1T::VertexContainer& vertices, Eigen::Ref<Eigen::VectorXd> values_ineq)
    {
        values_ineq[0] = vertices[0]->getData()[0] - 5;
        values_ineq[1] = vertices[1]->getData()[0] + 3;
        values_ineq[2] = vertices[1]->getData()[1];
        values_ineq *= mixed_edge_aux;
    }

    using MixedEdge2T = MixedEdgeGenericVectorFun<pFVectorVertex, pFVectorVertex>;

    void mixed_edge2_precompute(const Edge4T::VertexContainer& vertices) { mixed_edge_aux = 0.5; }
    void mixed_edge2_obj(const MixedEdge1T::VertexContainer& vertices, Eigen::Ref<Eigen::VectorXd> values_obj)
    {
        // x^2 + 2*y^2 + 3*z^2
        values_obj[0] = vertices[0]->getData()[0] * vertices[0]->getData()[0] + 2 * vertices[1]->getData()[0] * vertices[1]->getData()[0] +
                        3 * vertices[1]->getData()[1] * vertices[1]->getData()[1];
        // x*y + 2*x*z
        values_obj[1] = vertices[0]->getData()[0] * vertices[1]->getData()[0] + 2 * vertices[0]->getData()[0] * vertices[1]->getData()[1];
        values_obj *= mixed_edge_aux;
    }
    void mixed_edge2_eq(const MixedEdge1T::VertexContainer& vertices, Eigen::Ref<Eigen::VectorXd> values_eq)
    {
        // x^2 + 2*y^2 + 3*z^2
        values_eq[0] = vertices[0]->getData()[0] * vertices[0]->getData()[0] + 2 * vertices[1]->getData()[0] * vertices[1]->getData()[0] +
                       3 * vertices[1]->getData()[1] * vertices[1]->getData()[1];
        // x*y + 2*x*z
        values_eq[1] = vertices[0]->getData()[0] * vertices[1]->getData()[0] + 2 * vertices[0]->getData()[0] * vertices[1]->getData()[1];
        values_eq *= mixed_edge_aux;
    }
    void mixed_edge2_ineq(const MixedEdge1T::VertexContainer& vertices, Eigen::Ref<Eigen::VectorXd> values_ineq)
    {
        // x^2 + 2*y^2 + 3*z^2
        values_ineq[0] = vertices[0]->getData()[0] * vertices[0]->getData()[0] + 2 * vertices[1]->getData()[0] * vertices[1]->getData()[0] +
                         3 * vertices[1]->getData()[1] * vertices[1]->getData()[1];
        // x*y + 2*x*z
        values_ineq[1] = vertices[0]->getData()[0] * vertices[1]->getData()[0] + 2 * vertices[0]->getData()[0] * vertices[1]->getData()[1];
        values_ineq *= mixed_edge_aux;
    }

    HyperGraphOptimizationProblem optim;
    VertexSet::Ptr vertices;
    OptimizationEdgeSet::Ptr edges;

    pFVectorVertex::Ptr v1;
    pFVectorVertex::Ptr v2;

    double mixed_edge_aux = 1.0;
};

TEST_F(TestHyperGraphOptimizationProblemVertexBased, values)
{
    // Vertices
    vertices->addVertices({v1, v2});

    // Edges
    BaseEdge::Ptr edge1(new Edge1T(edge1_fun, false, *v1, *v2));
    edges->addObjectiveEdge(edge1);
    edges->addEqualityEdge(edge1);
    edges->addInequalityEdge(edge1);

    BaseEdge::Ptr edge2(new Edge2T(edge2_fun, false, *v2));
    edges->addObjectiveEdge(edge2);
    edges->addEqualityEdge(edge2);
    edges->addInequalityEdge(edge2);

    // specify solution
    Eigen::VectorXd values_sol(3);
    values_sol << 1 + 2 + 3, 4, 11;

    // compute values
    double value_obj = optim.computeValueObjective();
    EXPECT_DOUBLE_EQ(value_obj, values_sol.sum());

    Eigen::VectorXd values_eq = Eigen::VectorXd::Zero(optim.getEqualityDimension());
    optim.computeValuesEquality(values_eq);
    EXPECT_EQ_MATRIX(values_eq, values_sol, 1e-7);
    Eigen::VectorXd values_ineq = Eigen::VectorXd::Zero(optim.getInequalityDimension());
    optim.computeValuesInequality(values_ineq);
    EXPECT_EQ_MATRIX(values_ineq, values_sol, 1e-7);
    // combined computation
    double non_lsq_value = 0;
    Eigen::VectorXd values_lsq_obj_dummy;
    values_eq.setZero();
    values_ineq.setZero();
    optim.computeValues(non_lsq_value, values_lsq_obj_dummy, values_eq, values_ineq);
    EXPECT_DOUBLE_EQ(non_lsq_value, values_sol.sum());
    EXPECT_EQ_MATRIX(values_eq, values_sol, 1e-7);
    EXPECT_EQ_MATRIX(values_ineq, values_sol, 1e-7);
}

TEST_F(TestHyperGraphOptimizationProblemVertexBased, values_mixed)
{
    // Vertices
    vertices->addVertices({v1, v2});

    // Edges
    BaseEdge::Ptr edge1(new Edge1T(edge1_fun, false, *v1, *v2));
    edges->addObjectiveEdge(edge1);
    edges->addEqualityEdge(edge1);
    edges->addInequalityEdge(edge1);

    auto e2_precompute = [this](const MixedEdge1T::VertexContainer& vertices) { mixed_edge1_precompute(vertices); };
    auto e2_obj        = [this](const MixedEdge1T::VertexContainer& vertices, Eigen::Ref<Eigen::VectorXd> values_obj) {
        mixed_edge1_obj(vertices, values_obj);
    };
    auto e2_eq = [this](const MixedEdge1T::VertexContainer& vertices, Eigen::Ref<Eigen::VectorXd> values_eq) { mixed_edge1_eq(vertices, values_eq); };
    auto e2_ineq = [this](const MixedEdge1T::VertexContainer& vertices, Eigen::Ref<Eigen::VectorXd> values_ineq) {
        mixed_edge1_ineq(vertices, values_ineq);
    };
    BaseMixedEdge::Ptr edge2(new MixedEdge1T(1, 2, 3, e2_precompute, e2_obj, e2_eq, e2_ineq, *v1, *v2));
    edges->addMixedEdge(edge2);

    // specify solution
    Eigen::VectorXd values_obj_sol(2);
    values_obj_sol << 1 + 2 + 3, (1 + 2 + 3) * 2;
    Eigen::VectorXd values_eq_sol(3);
    values_eq_sol << 1 + 2 + 3, 4 * 2, (5 + 6) * 2;
    Eigen::VectorXd values_ineq_sol(4);
    values_ineq_sol << 1 + 2 + 3, (1 - 5) * 2, (1 + 3) * 2, 1 * 2;

    // compute values
    double value_obj = optim.computeValueObjective();
    EXPECT_DOUBLE_EQ(value_obj, values_obj_sol.sum());
    Eigen::VectorXd values_eq = Eigen::VectorXd::Zero(optim.getEqualityDimension());
    optim.computeValuesEquality(values_eq);
    EXPECT_EQ_MATRIX(values_eq, values_eq_sol, 1e-7);
    Eigen::VectorXd values_ineq = Eigen::VectorXd::Zero(optim.getInequalityDimension());
    optim.computeValuesInequality(values_ineq);
    EXPECT_EQ_MATRIX(values_ineq, values_ineq_sol, 1e-7);
    // combined computation
    double non_lsq_value = 0;
    Eigen::VectorXd values_lsq_obj_dummy;
    values_eq.setZero();
    values_ineq.setZero();
    optim.computeValues(non_lsq_value, values_lsq_obj_dummy, values_eq, values_ineq);
    EXPECT_DOUBLE_EQ(value_obj, values_obj_sol.sum());
    EXPECT_EQ_MATRIX(values_eq, values_eq_sol, 1e-7);
    EXPECT_EQ_MATRIX(values_ineq, values_ineq_sol, 1e-7);
}

TEST_F(TestHyperGraphOptimizationProblemVertexBased, values_lsq)
{
    // Vertices
    v2->setFixed(0, true);

    vertices->addVertices({v1, v2});

    // Edges
    edges->addObjectiveEdge(std::make_shared<Edge1T>(edge1_fun, false, *v1, *v2));
    edges->addEqualityEdge(std::make_shared<Edge1T>(edge1_fun, false, *v1, *v2));
    edges->addInequalityEdge(std::make_shared<Edge1T>(edge1_fun, false, *v1, *v2));

    edges->addLsqObjectiveEdge(std::make_shared<Edge2T>(edge2_fun, true, *v2));
    edges->addEqualityEdge(std::make_shared<Edge2T>(edge2_fun, true, *v2));
    edges->addInequalityEdge(std::make_shared<Edge2T>(edge2_fun, true, *v2));

    // specify solution
    Eigen::VectorXd values_lsq_sol(2);
    values_lsq_sol << 4, 11;

    Eigen::VectorXd values_constraint_sol(3);
    values_constraint_sol << 1 + 2 + 3, 4, 11;

    // compute values
    Eigen::VectorXd values_lsq_obj = Eigen::VectorXd::Zero(optim.getLsqObjectiveDimension());
    optim.computeValuesLsqObjective(values_lsq_obj);
    EXPECT_EQ_MATRIX(values_lsq_obj, values_lsq_sol, 1e-7);

    double values_non_lsq_obj = optim.computeValueNonLsqObjective();
    EXPECT_DOUBLE_EQ(values_non_lsq_obj, 1 + 2 + 3);

    double value_obj = optim.computeValueObjective();
    EXPECT_DOUBLE_EQ(value_obj, values_non_lsq_obj + values_lsq_obj.squaredNorm());

    Eigen::VectorXd values_eq = Eigen::VectorXd::Zero(optim.getEqualityDimension());
    optim.computeValuesEquality(values_eq);
    EXPECT_EQ_MATRIX(values_eq, values_constraint_sol, 1e-7);
    Eigen::VectorXd values_ineq = Eigen::VectorXd::Zero(optim.getInequalityDimension());
    optim.computeValuesInequality(values_ineq);
    EXPECT_EQ_MATRIX(values_ineq, values_constraint_sol, 1e-7);
    // combined computation
    values_lsq_obj.setZero();
    values_eq.setZero();
    values_ineq.setZero();
    optim.computeValues(values_non_lsq_obj, values_lsq_obj, values_eq, values_ineq);
    EXPECT_DOUBLE_EQ(values_non_lsq_obj, 1 + 2 + 3);
    EXPECT_EQ_MATRIX(values_lsq_obj, values_lsq_sol, 1e-7);
    EXPECT_EQ_MATRIX(values_eq, values_constraint_sol, 1e-7);
    EXPECT_EQ_MATRIX(values_ineq, values_constraint_sol, 1e-7);
}

TEST_F(TestHyperGraphOptimizationProblemVertexBased, values_lsq_mixed)
{
    // Vertices
    v2->setFixed(0, true);

    vertices->addVertices({v1, v2});

    // Edges
    BaseEdge::Ptr edge1o(new Edge1T(edge1_fun, false, *v1, *v2));
    BaseEdge::Ptr edge1e(new Edge1T(edge1_fun, false, *v1, *v2));
    BaseEdge::Ptr edge1i(new Edge1T(edge1_fun, false, *v1, *v2));
    edges->addObjectiveEdge(edge1o);
    edges->addEqualityEdge(edge1e);
    edges->addInequalityEdge(edge1i);

    BaseEdge::Ptr edge2o(new Edge2T(edge2_fun, true, *v2));
    BaseEdge::Ptr edge2e(new Edge2T(edge2_fun, false, *v2));  // for equalities, the LSQ flag is ignored
    BaseEdge::Ptr edge2i(new Edge2T(edge2_fun, false, *v2));  // for inequalities, the LSQ flag is ignored
    edges->addLsqObjectiveEdge(edge2o);
    edges->addEqualityEdge(edge2e);
    edges->addInequalityEdge(edge2i);

    auto e3_precompute = [this](const MixedEdge1T::VertexContainer& vertices) { mixed_edge1_precompute(vertices); };
    auto e3_obj        = [this](const MixedEdge1T::VertexContainer& vertices, Eigen::Ref<Eigen::VectorXd> values_obj) {
        mixed_edge1_obj(vertices, values_obj);
    };
    auto e3_eq = [this](const MixedEdge1T::VertexContainer& vertices, Eigen::Ref<Eigen::VectorXd> values_eq) { mixed_edge1_eq(vertices, values_eq); };
    auto e3_ineq = [this](const MixedEdge1T::VertexContainer& vertices, Eigen::Ref<Eigen::VectorXd> values_ineq) {
        mixed_edge1_ineq(vertices, values_ineq);
    };
    MixedEdge1T::Ptr edge3(new MixedEdge1T(1, 2, 3, e3_precompute, e3_obj, e3_eq, e3_ineq, *v1, *v2));
    edge3->setObjectiveLsqForm(true);  // active LSQ form for underlying objective function
    edges->addMixedEdge(edge3);

    // specify solution
    Eigen::VectorXd values_lsq_obj_sol(3);
    values_lsq_obj_sol << 4, 11, 2 * (1 + 2 + 3);

    Eigen::VectorXd values_obj_sol(3);
    values_obj_sol << 1 + 2 + 3, 4 * 4 + 11 * 11, std::pow((1 + 2 + 3) * 2, 2);
    //    Eigen::VectorXd values_eq_sol(3);
    //    values_eq_sol << 1 + 2 + 3, 4 * 4 + 11 * 11, std::pow(4 * 2, 2) + std::pow((5 + 6) * 2, 2);
    //    Eigen::VectorXd values_ineq_sol(3);
    //    values_ineq_sol << 1 + 2 + 3, 4 * 4 + 11 * 11, std::pow((1 - 5) * 2, 2) + std::pow((1 + 3) * 2, 2) + std::pow(1 * 2, 2);
    Eigen::VectorXd values_eq_sol(5);
    values_eq_sol << 1 + 2 + 3, 4, 11, 4 * 2, 2 * (5 + 6);
    Eigen::VectorXd values_ineq_sol(6);
    values_ineq_sol << 1 + 2 + 3, 4, 11, (1 - 5) * 2, (1 + 3) * 2, 1 * 2;

    // compute values
    Eigen::VectorXd values_lsq_obj = Eigen::VectorXd::Zero(optim.getLsqObjectiveDimension());
    optim.computeValuesLsqObjective(values_lsq_obj);
    EXPECT_EQ_MATRIX(values_lsq_obj, values_lsq_obj_sol, 1e-7);

    double values_non_lsq_obj = optim.computeValueNonLsqObjective();
    EXPECT_DOUBLE_EQ(values_non_lsq_obj, 1 + 2 + 3);

    double value_obj = optim.computeValueObjective();
    EXPECT_DOUBLE_EQ(value_obj, values_non_lsq_obj + values_lsq_obj.squaredNorm());

    Eigen::VectorXd values_eq = Eigen::VectorXd::Zero(optim.getEqualityDimension());
    optim.computeValuesEquality(values_eq);
    EXPECT_EQ_MATRIX(values_eq, values_eq_sol, 1e-7);
    Eigen::VectorXd values_ineq = Eigen::VectorXd::Zero(optim.getInequalityDimension());
    optim.computeValuesInequality(values_ineq);
    EXPECT_EQ_MATRIX(values_ineq, values_ineq_sol, 1e-7);
    // combined computation
    values_lsq_obj.setZero();
    values_eq.setZero();
    values_ineq.setZero();
    optim.computeValues(values_non_lsq_obj, values_lsq_obj, values_eq, values_ineq);
    EXPECT_DOUBLE_EQ(values_non_lsq_obj, 1 + 2 + 3);
    EXPECT_EQ_MATRIX(values_lsq_obj, values_lsq_obj_sol, 1e-7);
    EXPECT_EQ_MATRIX(values_eq, values_eq_sol, 1e-7);
    EXPECT_EQ_MATRIX(values_ineq, values_ineq_sol, 1e-7);
}

TEST_F(TestHyperGraphOptimizationProblemVertexBased, values_active_ineq)
{
    v1->values()[0] = -2;
    v2->values()[0] = -1;
    v2->values()[1] = 2;

    // Vertices
    vertices->addVertices({v1, v2});

    // Edges
    BaseEdge::Ptr edge1(new Edge1T(edge1_fun, false, *v1, *v2));

    Eigen::VectorXd values1(edge1->getDimension());
    edge1->computeValues(values1);
    bool active1 = values1[0] >= 0;
    edges->addInequalityEdge(edge1);

    BaseEdge::Ptr edge2(new Edge2T(edge2_fun, false, *v2));
    Eigen::VectorXd values2(edge2->getDimension());
    edge2->computeValues(values2);
    bool active2 = values2[0] >= 0;
    bool active3 = values2[1] >= 0;
    edges->addInequalityEdge(edge2);

    // make sure to have both active and inactive constraints
    EXPECT_TRUE(active1 || active2 || active3);
    EXPECT_FALSE(active1 && active2 && active3);

    Eigen::VectorXd values(3);
    optim.computeValuesActiveInequality(values);

    if (!active1)
        EXPECT_DOUBLE_EQ(values[0], 0);
    else
        EXPECT_GT(values[0], 0);
    if (!active2)
        EXPECT_DOUBLE_EQ(values[1], 0);
    else
        EXPECT_GT(values[1], 0);
    if (!active3)
        EXPECT_DOUBLE_EQ(values[2], 0);
    else
        EXPECT_GT(values[2], 0);

    Eigen::VectorXd values_scaled(3);
    optim.computeValuesActiveInequality(values_scaled, 5);

    if (!active1)
        EXPECT_DOUBLE_EQ(values_scaled[0], 0);
    else
        EXPECT_GT(values_scaled[0], 0);
    if (!active2)
        EXPECT_DOUBLE_EQ(values_scaled[1], 0);
    else
        EXPECT_GT(values_scaled[1], 0);
    if (!active3)
        EXPECT_DOUBLE_EQ(values_scaled[2], 0);
    else
        EXPECT_GT(values_scaled[2], 0);

    values *= 5.0;

    EXPECT_EQ_MATRIX(values, values_scaled, 1e-3);
}

TEST_F(TestHyperGraphOptimizationProblemVertexBased, values_active_ineq_mixed)
{
    v1->values()[0] = -2;
    v2->values()[0] = -1;
    v2->values()[1] = 2;

    // Vertices
    vertices->addVertices({v1, v2});

    // Edges
    BaseEdge::Ptr edge1(new Edge1T(edge1_fun, false, *v1, *v2));

    Eigen::VectorXd values1(edge1->getDimension());
    edge1->computeValues(values1);
    bool active1 = values1[0] >= 0;
    edges->addInequalityEdge(edge1);

    BaseEdge::Ptr edge2(new Edge2T(edge2_fun, false, *v2));
    Eigen::VectorXd values2(edge2->getDimension());
    edge2->computeValues(values2);
    bool active2 = values2[0] >= 0;
    bool active3 = values2[1] >= 0;
    edges->addInequalityEdge(edge2);

    auto e3_precompute = [this](const MixedEdge1T::VertexContainer& vertices) { mixed_edge1_precompute(vertices); };
    auto e3_obj        = [this](const MixedEdge1T::VertexContainer& vertices, Eigen::Ref<Eigen::VectorXd> values_obj) {
        mixed_edge1_obj(vertices, values_obj);
    };
    auto e3_eq = [this](const MixedEdge1T::VertexContainer& vertices, Eigen::Ref<Eigen::VectorXd> values_eq) { mixed_edge1_eq(vertices, values_eq); };
    auto e3_ineq = [this](const MixedEdge1T::VertexContainer& vertices, Eigen::Ref<Eigen::VectorXd> values_ineq) {
        mixed_edge1_ineq(vertices, values_ineq);
    };
    MixedEdge1T::Ptr edge3(new MixedEdge1T(1, 2, 3, e3_precompute, e3_obj, e3_eq, e3_ineq, *v1, *v2));
    edge3->precompute();
    Eigen::VectorXd values3(edge3->getInequalityDimension());
    edge3->computeInequalityValues(values3);
    bool active4 = values3[0] >= 0;
    bool active5 = values3[1] >= 0;
    bool active6 = values3[2] >= 0;
    edges->addMixedEdge(edge3);

    // make sure to have both active and inactive constraints in the mixed constraint set
    EXPECT_TRUE(active4 || active5 || active6);
    EXPECT_FALSE(active4 && active5 && active6);

    Eigen::VectorXd values(6);
    optim.computeValuesActiveInequality(values);

    if (!active1)
        EXPECT_DOUBLE_EQ(values[0], 0);
    else
        EXPECT_GT(values[0], 0);
    if (!active2)
        EXPECT_DOUBLE_EQ(values[1], 0);
    else
        EXPECT_GT(values[1], 0);
    if (!active3)
        EXPECT_DOUBLE_EQ(values[2], 0);
    else
        EXPECT_GT(values[2], 0);
    if (!active4)
        EXPECT_DOUBLE_EQ(values[3], 0);
    else
        EXPECT_GT(values[3], 0);
    if (!active5)
        EXPECT_DOUBLE_EQ(values[4], 0);
    else
        EXPECT_GT(values[4], 0);
    if (!active6)
        EXPECT_DOUBLE_EQ(values[5], 0);
    else
        EXPECT_GT(values[5], 0);

    Eigen::VectorXd values_scaled(6);
    optim.computeValuesActiveInequality(values_scaled, 5);

    if (!active1)
        EXPECT_DOUBLE_EQ(values_scaled[0], 0);
    else
        EXPECT_GT(values_scaled[0], 0);
    if (!active2)
        EXPECT_DOUBLE_EQ(values_scaled[1], 0);
    else
        EXPECT_GT(values_scaled[1], 0);
    if (!active3)
        EXPECT_DOUBLE_EQ(values_scaled[2], 0);
    else
        EXPECT_GT(values_scaled[2], 0);
    if (!active4)
        EXPECT_DOUBLE_EQ(values_scaled[3], 0);
    else
        EXPECT_GT(values_scaled[3], 0);
    if (!active5)
        EXPECT_DOUBLE_EQ(values_scaled[4], 0);
    else
        EXPECT_GT(values_scaled[4], 0);
    if (!active6)
        EXPECT_DOUBLE_EQ(values_scaled[5], 0);
    else
        EXPECT_GT(values_scaled[5], 0);

    values *= 5.0;

    EXPECT_EQ_MATRIX(values, values_scaled, 1e-3);
}

TEST_F(TestHyperGraphOptimizationProblemVertexBased, bound_values)
{
    v1->values()[0] = -2;
    v2->values()[0] = -1;
    v2->values()[1] = 5;

    Eigen::Vector2d lb(-5, std::numeric_limits<double>::lowest());
    Eigen::Vector2d ub(0, 4);  // the component violets bound (5 > 4)

    // -inf < -2 < inf
    // -5 < -1 < 0
    // -inf < 5 < 4

    v2->setLowerBounds(lb);
    v2->setUpperBounds(ub);

    // Vertices
    vertices->addVertices({v1, v2});

    EXPECT_EQ(optim.finiteBoundsDimension(), 3);
    EXPECT_EQ(optim.finiteCombinedBoundsDimension(), 2);

    Eigen::Vector3d lb_values;
    Eigen::Vector3d ub_values;
    optim.getBounds(lb_values, ub_values);

    EXPECT_LT(lb_values[0], -1e20);
    EXPECT_DOUBLE_EQ(lb_values[1], -5);
    EXPECT_LT(lb_values[2], -1e20);

    EXPECT_GT(ub_values[0], 1e20);
    EXPECT_DOUBLE_EQ(ub_values[1], 0);
    EXPECT_DOUBLE_EQ(ub_values[2], 4);

    Eigen::Vector2d distances;
    optim.computeDistanceFiniteCombinedBounds(distances);
    EXPECT_DOUBLE_EQ(distances[0], 0);
    EXPECT_DOUBLE_EQ(distances[1], 1);

    // specify jacobian solution
    Eigen::MatrixXd jacobian_sol(2, 3);
    jacobian_sol << 0, 0, 0, 0, 0, 1;

    // compute dense jacobian
    Eigen::MatrixXd jacobian(optim.finiteCombinedBoundsDimension(), optim.getParameterDimension());
    jacobian.setZero();
    optim.computeDenseJacobianFiniteCombinedBounds(jacobian);

    EXPECT_EQ_MATRIX(jacobian, jacobian_sol, 1e-7);

    // compute sparse jacobian
    Eigen::SparseMatrix<double> sparse_jacobian(optim.finiteCombinedBoundsDimension(), optim.getParameterDimension());
    optim.computeSparseJacobianFiniteCombinedBounds(sparse_jacobian);

    EXPECT_EQ_MATRIX(sparse_jacobian, jacobian_sol, 1e-7);

    // compute other sparse jacobian
    int nnz = optim.computeSparseJacobianFiniteCombinedBoundsNNZ();
    EXPECT_LE(nnz, jacobian_sol.rows() * jacobian_sol.cols());
    EXPECT_GE(nnz, 1);
    computeSparseJacobianFiniteCombinedBoundsViaTriplets(nnz, jacobian.rows(), jacobian.cols(), sparse_jacobian);

    EXPECT_EQ_MATRIX(sparse_jacobian, jacobian_sol, 1e-7);

    // combined jacobian
    testCombinedSparseJacobian(nullptr, nullptr, nullptr, 1e-6, &jacobian_sol);

    // now violate the lower bound as well
    v2->values()[0] = -10;

    optim.computeDistanceFiniteCombinedBounds(distances);
    EXPECT_DOUBLE_EQ(distances[0], 5);
    EXPECT_DOUBLE_EQ(distances[1], 1);

    // specify jacobian solution
    jacobian_sol << 0, -1, 0, 0, 0, 1;

    // compute dense jacobian
    jacobian.setZero();
    optim.computeDenseJacobianFiniteCombinedBounds(jacobian);

    EXPECT_EQ_MATRIX(jacobian, jacobian_sol, 1e-7);

    // compute sparse jacobian
    sparse_jacobian.setZero();
    optim.computeSparseJacobianFiniteCombinedBounds(sparse_jacobian);

    EXPECT_EQ_MATRIX(sparse_jacobian, jacobian_sol, 1e-7);

    // compute other sparse jacobian
    nnz = optim.computeSparseJacobianFiniteCombinedBoundsNNZ();
    EXPECT_LE(nnz, jacobian_sol.rows() * jacobian_sol.cols());
    EXPECT_GE(nnz, 1);
    computeSparseJacobianFiniteCombinedBoundsViaTriplets(nnz, jacobian.rows(), jacobian.cols(), sparse_jacobian);

    EXPECT_EQ_MATRIX(sparse_jacobian, jacobian_sol, 1e-7);

    // combined jacobian
    testCombinedSparseJacobian(nullptr, nullptr, nullptr, 1e-6, &jacobian_sol);

    // check universal bound setter and getter methods
    EXPECT_LE(optim.getLowerBound(0), -CORBO_INF_DBL);
    EXPECT_DOUBLE_EQ(optim.getLowerBound(1), -5);
    EXPECT_LE(optim.getLowerBound(2), -CORBO_INF_DBL);
    EXPECT_DOUBLE_EQ(optim.getUpperBound(1), 0);
    EXPECT_DOUBLE_EQ(optim.getUpperBound(2), 4);

    optim.setLowerBound(0, -10.5);
    optim.setUpperBound(2, 10.5);
    EXPECT_DOUBLE_EQ(optim.getLowerBound(0), -10.5);
    EXPECT_DOUBLE_EQ(optim.getUpperBound(2), 10.5);
}

TEST_F(TestHyperGraphOptimizationProblemVertexBased, bound_values_jacobian_fixed)
{
    v1->values()[0] = -2;
    v2->values()[0] = -1;
    v2->values()[1] = 5;

    v1->setFixed(0, true);
    v2->setFixed(1, true);

    Eigen::Vector2d lb(-5, std::numeric_limits<double>::lowest());
    Eigen::Vector2d ub(0, 4);  // the component violets bound (5 > 4)

    // -inf < -2 < inf
    // -5 < -1 < 0
    // -inf < 5 < 4

    v2->setLowerBounds(lb);
    v2->setUpperBounds(ub);

    // Vertices
    vertices->addVertices({v1, v2});

    EXPECT_EQ(optim.getParameterDimension(), 1);

    // specify jacobian solution
    Eigen::MatrixXd jacobian_sol(1, 1);
    jacobian_sol << 0;

    // compute dense jacobian
    Eigen::MatrixXd jacobian(optim.finiteCombinedBoundsDimension(), optim.getParameterDimension());
    jacobian.setZero();
    optim.computeDenseJacobianFiniteCombinedBounds(jacobian);

    EXPECT_EQ_MATRIX(jacobian, jacobian_sol, 1e-7);

    // compute sparse jacobian
    Eigen::SparseMatrix<double> sparse_jacobian(optim.finiteCombinedBoundsDimension(), optim.getParameterDimension());
    optim.computeSparseJacobianFiniteCombinedBounds(sparse_jacobian);

    EXPECT_EQ_MATRIX(sparse_jacobian, jacobian_sol, 1e-7);

    // compute other sparse jacobian
    int nnz = optim.computeSparseJacobianFiniteCombinedBoundsNNZ();
    EXPECT_LE(nnz, jacobian_sol.rows() * jacobian_sol.cols());
    EXPECT_GE(nnz, 1);
    computeSparseJacobianFiniteCombinedBoundsViaTriplets(nnz, jacobian.rows(), jacobian.cols(), sparse_jacobian);

    EXPECT_EQ_MATRIX(sparse_jacobian, jacobian_sol, 1e-7);

    // combined jacobian
    testCombinedSparseJacobian(nullptr, nullptr, nullptr, 1e-6, &jacobian_sol);

    // now violate the lower bound as well
    v2->values()[0] = -10;

    // specify jacobian solution
    jacobian_sol << -1;

    // compute dense jacobian
    jacobian.setZero();
    optim.computeDenseJacobianFiniteCombinedBounds(jacobian);

    EXPECT_EQ_MATRIX(jacobian, jacobian_sol, 1e-7);

    // compute sparse jacobian
    sparse_jacobian.setZero();
    optim.computeSparseJacobianFiniteCombinedBounds(sparse_jacobian);

    EXPECT_EQ_MATRIX(sparse_jacobian, jacobian_sol, 1e-7);

    // compute other sparse jacobian
    nnz = optim.computeSparseJacobianFiniteCombinedBoundsNNZ();
    EXPECT_LE(nnz, jacobian_sol.rows() * jacobian_sol.cols());
    EXPECT_GE(nnz, 1);
    computeSparseJacobianFiniteCombinedBoundsViaTriplets(nnz, jacobian.rows(), jacobian.cols(), sparse_jacobian);

    EXPECT_EQ_MATRIX(sparse_jacobian, jacobian_sol, 1e-7);

    // combined jacobian
    testCombinedSparseJacobian(nullptr, nullptr, nullptr, 1e-6, &jacobian_sol);

    // check universal bound setter and getter methods
    EXPECT_DOUBLE_EQ(optim.getLowerBound(0), -5);
    optim.setLowerBound(0, 5);
    optim.setUpperBound(0, 10);
    EXPECT_DOUBLE_EQ(optim.getLowerBound(0), 5);
    EXPECT_DOUBLE_EQ(optim.getUpperBound(0), 10);
}

/****************
 *   JACOBIANS
 * **************/

TEST_F(TestHyperGraphOptimizationProblemVertexBased, jacobian)
{
    // Vertices
    vertices->addVertices({v1, v2});

    // Edges
    BaseEdge::Ptr edge1o(new Edge1T(edge1_fun, false, *v1, *v2));
    BaseEdge::Ptr edge1e(new Edge1T(edge1_fun, false, *v1, *v2));
    BaseEdge::Ptr edge1i(new Edge1T(edge1_fun, false, *v1, *v2));
    edges->addObjectiveEdge(edge1o);
    edges->addEqualityEdge(edge1e);
    edges->addInequalityEdge(edge1i);

    BaseEdge::Ptr edge2o(new Edge2T(edge2_fun, false, *v2));
    BaseEdge::Ptr edge2e(new Edge2T(edge2_fun, false, *v2));
    BaseEdge::Ptr edge2i(new Edge2T(edge2_fun, false, *v2));
    edges->addObjectiveEdge(edge2o);
    edges->addEqualityEdge(edge2e);
    edges->addInequalityEdge(edge2i);

    // specify gradient solution
    Eigen::VectorXd gradient_sol(3);
    gradient_sol << 1, 11, 9;

    Eigen::VectorXd gradient_non_lsq_sol(3);
    gradient_non_lsq_sol << 1, 11, 9;

    // specify jacobian solution
    Eigen::MatrixXd jacobian_sol(3, 3);
    jacobian_sol << 1, 2, 3, 0, 4, 0, 0, 5, 6;

    Eigen::MatrixXd jacobian_lsq_obj_sol;  // emtpy -> no lsq objectives

    testObjectiveAndEqualityAndInequalityJacobians(gradient_sol, gradient_non_lsq_sol, jacobian_lsq_obj_sol, jacobian_sol, jacobian_sol, -1, 6, 6,
                                                   1e-7);

    testCombinedSparseJacobian(nullptr, &jacobian_sol, &jacobian_sol, 1e-7);
}

TEST_F(TestHyperGraphOptimizationProblemVertexBased, jacobian_mixed_edge)
{
    // Vertices
    vertices->addVertices({v1, v2});

    // Edges
    BaseEdge::Ptr edge1o(new Edge1T(edge1_fun, false, *v1, *v2));
    BaseEdge::Ptr edge1e(new Edge1T(edge1_fun, false, *v1, *v2));
    BaseEdge::Ptr edge1i(new Edge1T(edge1_fun, false, *v1, *v2));
    edges->addObjectiveEdge(edge1o);
    edges->addEqualityEdge(edge1e);
    edges->addInequalityEdge(edge1i);

    BaseEdge::Ptr edge2o(new Edge2T(edge2_fun, false, *v2));
    BaseEdge::Ptr edge2e(new Edge2T(edge2_fun, false, *v2));
    BaseEdge::Ptr edge2i(new Edge2T(edge2_fun, false, *v2));
    edges->addObjectiveEdge(edge2o);
    edges->addEqualityEdge(edge2e);
    edges->addInequalityEdge(edge2i);

    auto e3_precompute = [this](const MixedEdge1T::VertexContainer& vertices) { mixed_edge1_precompute(vertices); };
    auto e3_obj        = [this](const MixedEdge1T::VertexContainer& vertices, Eigen::Ref<Eigen::VectorXd> values_obj) {
        mixed_edge1_obj(vertices, values_obj);
    };
    auto e3_eq = [this](const MixedEdge1T::VertexContainer& vertices, Eigen::Ref<Eigen::VectorXd> values_eq) { mixed_edge1_eq(vertices, values_eq); };
    auto e3_ineq = [this](const MixedEdge1T::VertexContainer& vertices, Eigen::Ref<Eigen::VectorXd> values_ineq) {
        mixed_edge1_ineq(vertices, values_ineq);
    };
    MixedEdge1T::Ptr edge3(new MixedEdge1T(1, 2, 3, e3_precompute, e3_obj, e3_eq, e3_ineq, *v1, *v2));
    edges->addMixedEdge(edge3);

    // specify gradient solution
    Eigen::VectorXd gradient_sol(3);
    gradient_sol << 3, 15, 15;

    Eigen::VectorXd gradient_non_lsq_sol = gradient_sol;

    // specify jacobian solution
    // Eigen::MatrixXd jacobian_obj_sol(4, 3);
    // jacobian_obj_sol << 1, 2, 3, 0, 4, 0, 0, 5, 6, 1 * 2, 2 * 2, 3 * 2;

    Eigen::MatrixXd jacobian_lsq_obj_sol;  // emtpy -> no lsq objectives

    Eigen::MatrixXd jacobian_eq_sol(5, 3);
    jacobian_eq_sol << 1, 2, 3, 0, 4, 0, 0, 5, 6, 4 * 2, 0, 0, 0, 5 * 2, 6 * 2;

    Eigen::MatrixXd jacobian_ineq_sol(6, 3);
    jacobian_ineq_sol << 1, 2, 3, 0, 4, 0, 0, 5, 6, 2, 0, 0, 0, 2, 0, 0, 0, 2;

    testObjectiveAndEqualityAndInequalityJacobians(gradient_sol, gradient_non_lsq_sol, jacobian_lsq_obj_sol, jacobian_eq_sol, jacobian_ineq_sol, -1,
                                                   9, 9, 1e-6);
    testCombinedSparseJacobian(nullptr, &jacobian_eq_sol, &jacobian_ineq_sol, 1e-7);
}

TEST_F(TestHyperGraphOptimizationProblemVertexBased, jacobian_fixed)
{
    // Vertices
    vertices->addVertices({v1, v2});

    v1->setFixed(0, true);
    v2->setFixed(1, true);

    // Edges
    BaseEdge::Ptr edge1o(new Edge1T(edge1_fun, false, *v1, *v2));
    BaseEdge::Ptr edge1e(new Edge1T(edge1_fun, false, *v1, *v2));
    BaseEdge::Ptr edge1i(new Edge1T(edge1_fun, false, *v1, *v2));
    edges->addObjectiveEdge(edge1o);
    edges->addEqualityEdge(edge1e);
    edges->addInequalityEdge(edge1i);

    BaseEdge::Ptr edge2o(new Edge2T(edge2_fun, false, *v2));
    BaseEdge::Ptr edge2e(new Edge2T(edge2_fun, false, *v2));
    BaseEdge::Ptr edge2i(new Edge2T(edge2_fun, false, *v2));
    edges->addObjectiveEdge(edge2o);
    edges->addEqualityEdge(edge2e);
    edges->addInequalityEdge(edge2i);

    // specify gradient solution
    Eigen::VectorXd gradient_sol(1);
    gradient_sol << 11;

    Eigen::VectorXd gradient_non_lsq_sol = gradient_sol;

    // specify jacobian solution
    Eigen::MatrixXd jacobian_lsq_obj_sol;  // emtpy -> no lsq objectives

    Eigen::MatrixXd jacobian_sol(3, 1);
    jacobian_sol << 2, 4, 5;

    testObjectiveAndEqualityAndInequalityJacobians(gradient_sol, gradient_non_lsq_sol, jacobian_lsq_obj_sol, jacobian_sol, jacobian_sol, -1, 3, 3,
                                                   1e-6);
    testCombinedSparseJacobian(nullptr, &jacobian_sol, &jacobian_sol, 1e-7);
}

TEST_F(TestHyperGraphOptimizationProblemVertexBased, jacobian_fixed_mixed)
{
    // Vertices
    vertices->addVertices({v1, v2});

    v1->setFixed(0, true);
    v2->setFixed(1, true);

    // Edges
    BaseEdge::Ptr edge1o(new Edge1T(edge1_fun, false, *v1, *v2));
    BaseEdge::Ptr edge1e(new Edge1T(edge1_fun, false, *v1, *v2));
    BaseEdge::Ptr edge1i(new Edge1T(edge1_fun, false, *v1, *v2));
    edges->addObjectiveEdge(edge1o);
    edges->addEqualityEdge(edge1e);
    edges->addInequalityEdge(edge1i);

    BaseEdge::Ptr edge2o(new Edge2T(edge2_fun, false, *v2));
    BaseEdge::Ptr edge2e(new Edge2T(edge2_fun, false, *v2));
    BaseEdge::Ptr edge2i(new Edge2T(edge2_fun, false, *v2));
    edges->addObjectiveEdge(edge2o);
    edges->addEqualityEdge(edge2e);
    edges->addInequalityEdge(edge2i);

    auto e3_precompute = [this](const MixedEdge1T::VertexContainer& vertices) { mixed_edge1_precompute(vertices); };
    auto e3_obj        = [this](const MixedEdge1T::VertexContainer& vertices, Eigen::Ref<Eigen::VectorXd> values_obj) {
        mixed_edge1_obj(vertices, values_obj);
    };
    auto e3_eq = [this](const MixedEdge1T::VertexContainer& vertices, Eigen::Ref<Eigen::VectorXd> values_eq) { mixed_edge1_eq(vertices, values_eq); };
    auto e3_ineq = [this](const MixedEdge1T::VertexContainer& vertices, Eigen::Ref<Eigen::VectorXd> values_ineq) {
        mixed_edge1_ineq(vertices, values_ineq);
    };
    MixedEdge1T::Ptr edge3(new MixedEdge1T(1, 2, 3, e3_precompute, e3_obj, e3_eq, e3_ineq, *v1, *v2));
    edges->addMixedEdge(edge3);

    // specify gradient solution
    Eigen::VectorXd gradient_sol(1);
    gradient_sol << 15;

    Eigen::VectorXd gradient_non_lsq_sol = gradient_sol;

    // specify jacobian solution
    // Eigen::MatrixXd jacobian_obj_sol(4, 1);
    // jacobian_obj_sol << 2, 4, 5, 2 * 2;
    Eigen::MatrixXd jacobian_lsq_obj_sol;  // emtpy -> no lsq objectives

    Eigen::MatrixXd jacobian_eq_sol(5, 1);
    jacobian_eq_sol << 2, 4, 5, 0, 5 * 2;

    Eigen::MatrixXd jacobian_ineq_sol(6, 1);
    jacobian_ineq_sol << 2, 4, 5, 0, 2, 0;

    testObjectiveAndEqualityAndInequalityJacobians(gradient_sol, gradient_non_lsq_sol, jacobian_lsq_obj_sol, jacobian_eq_sol, jacobian_ineq_sol, -1,
                                                   4, 4, 1e-6);
    testCombinedSparseJacobian(nullptr, &jacobian_eq_sol, &jacobian_ineq_sol, 1e-6);
}

TEST_F(TestHyperGraphOptimizationProblemVertexBased, jacobian_fixed2)
{
    // Vertices
    vertices->addVertex(v1);

    Eigen::Vector4d val(1, 2, 3, 4);
    pFVectorVertex::Ptr v3 = std::make_shared<pFVectorVertex>(val);

    vertices->addVertex(v3);

    v3->setFixed(1, true);
    v3->setFixed(2, true);

    // Edges
    BaseEdge::Ptr edge1o(new Edge1T(edge1_fun, false, *v1, *v3));
    BaseEdge::Ptr edge1e(new Edge1T(edge1_fun, false, *v1, *v3));
    BaseEdge::Ptr edge1i(new Edge1T(edge1_fun, false, *v1, *v3));
    edges->addObjectiveEdge(edge1o);
    edges->addEqualityEdge(edge1e);
    edges->addInequalityEdge(edge1i);

    BaseEdge::Ptr edge2o(new Edge2T(edge2_fun, false, *v3));
    BaseEdge::Ptr edge2e(new Edge2T(edge2_fun, false, *v3));
    BaseEdge::Ptr edge2i(new Edge2T(edge2_fun, false, *v3));
    edges->addObjectiveEdge(edge2o);
    edges->addEqualityEdge(edge2e);
    edges->addInequalityEdge(edge2i);

    // specify gradient solution
    Eigen::VectorXd gradient_sol(3);
    gradient_sol << 1, 11, 0;

    Eigen::VectorXd gradient_non_lsq_sol = gradient_sol;

    // specify jacobian solution
    Eigen::MatrixXd jacobian_lsq_obj_sol;  // emtpy -> no lsq objectives

    Eigen::MatrixXd jacobian_sol(3, 3);
    jacobian_sol << 1, 2, 0, 0, 4, 0, 0, 5, 0;

    testObjectiveAndEqualityAndInequalityJacobians(gradient_sol, gradient_non_lsq_sol, jacobian_lsq_obj_sol, jacobian_sol, jacobian_sol, -1, 4, 4,
                                                   1e-6);
    testCombinedSparseJacobian(nullptr, &jacobian_sol, &jacobian_sol, 1e-7);
}

TEST_F(TestHyperGraphOptimizationProblemVertexBased, jacobian_fixed2_mixed)
{
    // Vertices
    vertices->addVertex(v1);

    Eigen::Vector4d val(1, 2, 3, 4);
    pFVectorVertex::Ptr v3 = std::make_shared<pFVectorVertex>(val);

    vertices->addVertex(v3);

    v3->setFixed(1, true);
    v3->setFixed(2, true);

    // Edges
    BaseEdge::Ptr edge1o(new Edge1T(edge1_fun, false, *v1, *v3));
    BaseEdge::Ptr edge1e(new Edge1T(edge1_fun, false, *v1, *v3));
    BaseEdge::Ptr edge1i(new Edge1T(edge1_fun, false, *v1, *v3));
    edges->addObjectiveEdge(edge1o);
    edges->addEqualityEdge(edge1e);
    edges->addInequalityEdge(edge1i);

    BaseEdge::Ptr edge2o(new Edge2T(edge2_fun, false, *v3));
    BaseEdge::Ptr edge2e(new Edge2T(edge2_fun, false, *v3));
    BaseEdge::Ptr edge2i(new Edge2T(edge2_fun, false, *v3));
    edges->addObjectiveEdge(edge2e);
    edges->addEqualityEdge(edge2e);
    edges->addInequalityEdge(edge2i);

    auto e3_precompute = [this](const MixedEdge1T::VertexContainer& vertices) { mixed_edge1_precompute(vertices); };
    auto e3_obj        = [this](const MixedEdge1T::VertexContainer& vertices, Eigen::Ref<Eigen::VectorXd> values_obj) {
        mixed_edge1_obj(vertices, values_obj);
    };
    auto e3_eq = [this](const MixedEdge1T::VertexContainer& vertices, Eigen::Ref<Eigen::VectorXd> values_eq) { mixed_edge1_eq(vertices, values_eq); };
    auto e3_ineq = [this](const MixedEdge1T::VertexContainer& vertices, Eigen::Ref<Eigen::VectorXd> values_ineq) {
        mixed_edge1_ineq(vertices, values_ineq);
    };
    MixedEdge1T::Ptr edge3(new MixedEdge1T(1, 2, 3, e3_precompute, e3_obj, e3_eq, e3_ineq, *v1, *v3));
    edges->addMixedEdge(edge3);

    // specify gradient solution
    Eigen::VectorXd gradient_sol(3);
    gradient_sol << 1 + 0 + 0 + 2, 2 + 4 + 5 + 4, 0 + 0 + 0 + 0;

    Eigen::VectorXd gradient_non_lsq_sol = gradient_sol;

    // specify jacobian solution
    Eigen::MatrixXd jacobian_lsq_obj_sol;  // emtpy -> no lsq objectives

    // Eigen::MatrixXd jacobian_obj_sol(4, 3);
    // jacobian_obj_sol << 1, 2, 0, 0, 4, 0, 0, 5, 0, 1 * 2, 2 * 2, 0;

    Eigen::MatrixXd jacobian_eq_sol(5, 3);
    jacobian_eq_sol << 1, 2, 0, 0, 4, 0, 0, 5, 0, 4 * 2, 0, 0, 0, 5 * 2, 0;

    Eigen::MatrixXd jacobian_ineq_sol(6, 3);
    jacobian_ineq_sol << 1, 2, 0, 0, 4, 0, 0, 5, 0, 2, 0, 0, 0, 2, 0, 0, 0, 0;

    testObjectiveAndEqualityAndInequalityJacobians(gradient_sol, gradient_non_lsq_sol, jacobian_lsq_obj_sol, jacobian_eq_sol, jacobian_ineq_sol, -1,
                                                   9, 9, 1e-6);
    testCombinedSparseJacobian(nullptr, &jacobian_eq_sol, &jacobian_ineq_sol, 1e-6);
}

TEST_F(TestHyperGraphOptimizationProblemVertexBased, jacobian_modified_values)
{
    // Vertices
    vertices->addVertices({v1, v2});

    // Edges
    BaseEdge::Ptr edge1o(new Edge1T(edge1_fun, false, *v1, *v2));
    BaseEdge::Ptr edge1e(new Edge1T(edge1_fun, false, *v1, *v2));
    BaseEdge::Ptr edge1i(new Edge1T(edge1_fun, false, *v1, *v2));
    edges->addObjectiveEdge(edge1o);
    edges->addEqualityEdge(edge1e);
    edges->addInequalityEdge(edge1i);

    BaseEdge::Ptr edge2o(new Edge2T(edge2_fun, false, *v2));
    BaseEdge::Ptr edge2e(new Edge2T(edge2_fun, false, *v2));
    BaseEdge::Ptr edge2i(new Edge2T(edge2_fun, false, *v2));
    edges->addObjectiveEdge(edge2o);
    edges->addEqualityEdge(edge2e);
    edges->addInequalityEdge(edge2i);

    std::array<double, 3> multipliers = {{2.0, 3.0, 4.0}};

    // specify gradient solution
    Eigen::VectorXd gradient_sol(3);
    gradient_sol << 1, 2 + 4 + 5, 3 + 6;

    Eigen::VectorXd gradient_non_lsq_sol = gradient_sol;

    // Specify jacobian solution
    Eigen::MatrixXd jacobian_lsq_obj_sol;  // emtpy -> no lsq objectives

    Eigen::MatrixXd jacobian_sol(3, 3);
    jacobian_sol << 1, 2, 3, 0, 4, 0, 0, 5, 6;
    jacobian_sol.row(0) *= multipliers[0];
    jacobian_sol.row(1) *= multipliers[1];
    jacobian_sol.row(2) *= multipliers[2];

    testObjectiveAndEqualityAndInequalityJacobians(gradient_sol, gradient_non_lsq_sol, jacobian_lsq_obj_sol, jacobian_sol, jacobian_sol, -1, 6, 6,
                                                   1e-6, nullptr, multipliers.data(), multipliers.data());
}

TEST_F(TestHyperGraphOptimizationProblemVertexBased, jacobian_modified_values_mixed)
{
    // Vertices
    vertices->addVertices({v1, v2});

    // Edges
    BaseEdge::Ptr edge1o(new Edge1T(edge1_fun, false, *v1, *v2));
    BaseEdge::Ptr edge1e(new Edge1T(edge1_fun, false, *v1, *v2));
    BaseEdge::Ptr edge1i(new Edge1T(edge1_fun, false, *v1, *v2));
    edges->addObjectiveEdge(edge1o);
    edges->addEqualityEdge(edge1e);
    edges->addInequalityEdge(edge1i);

    BaseEdge::Ptr edge2o(new Edge2T(edge2_fun, false, *v2));
    BaseEdge::Ptr edge2e(new Edge2T(edge2_fun, false, *v2));
    BaseEdge::Ptr edge2i(new Edge2T(edge2_fun, false, *v2));
    edges->addObjectiveEdge(edge2o);
    edges->addEqualityEdge(edge2e);
    edges->addInequalityEdge(edge2i);

    auto e3_precompute = [this](const MixedEdge1T::VertexContainer& vertices) { mixed_edge1_precompute(vertices); };
    auto e3_obj        = [this](const MixedEdge1T::VertexContainer& vertices, Eigen::Ref<Eigen::VectorXd> values_obj) {
        mixed_edge1_obj(vertices, values_obj);
    };
    auto e3_eq = [this](const MixedEdge1T::VertexContainer& vertices, Eigen::Ref<Eigen::VectorXd> values_eq) { mixed_edge1_eq(vertices, values_eq); };
    auto e3_ineq = [this](const MixedEdge1T::VertexContainer& vertices, Eigen::Ref<Eigen::VectorXd> values_ineq) {
        mixed_edge1_ineq(vertices, values_ineq);
    };
    MixedEdge1T::Ptr edge3(new MixedEdge1T(1, 2, 3, e3_precompute, e3_obj, e3_eq, e3_ineq, *v1, *v2));
    edges->addMixedEdge(edge3);

    // std::array<double, 4> multipliers_obj  = {{2.0, 3.0, 4.0, 5.0}};
    std::array<double, 5> multipliers_eq   = {{2.0, 3.0, 4.0, 3.0, 2.0}};
    std::array<double, 6> multipliers_ineq = {{2.0, 3.0, 4.0, 5.0, 3.0, 2.0}};

    // specify jacobian solution
    Eigen::MatrixXd jacobian_lsq_obj_sol;  // emtpy -> no lsq objectives

    Eigen::MatrixXd jacobian_obj_sol(4, 3);
    jacobian_obj_sol << 1, 2, 3, 0, 4, 0, 0, 5, 6, 1 * 2, 2 * 2, 3 * 2;
    //    jacobian_obj_sol.row(0) *= multipliers_obj[0];
    //    jacobian_obj_sol.row(1) *= multipliers_obj[1];
    //    jacobian_obj_sol.row(2) *= multipliers_obj[2];
    //    jacobian_obj_sol.row(3) *= multipliers_obj[3];

    Eigen::MatrixXd jacobian_eq_sol(5, 3);
    jacobian_eq_sol << 1, 2, 3, 0, 4, 0, 0, 5, 6, 4 * 2, 0, 0, 0, 5 * 2, 6 * 2;
    jacobian_eq_sol.row(0) *= multipliers_eq[0];
    jacobian_eq_sol.row(1) *= multipliers_eq[1];
    jacobian_eq_sol.row(2) *= multipliers_eq[2];
    jacobian_eq_sol.row(3) *= multipliers_eq[3];
    jacobian_eq_sol.row(4) *= multipliers_eq[4];

    Eigen::MatrixXd jacobian_ineq_sol(6, 3);
    jacobian_ineq_sol << 1, 2, 3, 0, 4, 0, 0, 5, 6, 2, 0, 0, 0, 2, 0, 0, 0, 2;
    jacobian_ineq_sol.row(0) *= multipliers_ineq[0];
    jacobian_ineq_sol.row(1) *= multipliers_ineq[1];
    jacobian_ineq_sol.row(2) *= multipliers_ineq[2];
    jacobian_ineq_sol.row(3) *= multipliers_ineq[3];
    jacobian_ineq_sol.row(4) *= multipliers_ineq[4];
    jacobian_ineq_sol.row(5) *= multipliers_ineq[5];

    // specify gradient solution
    Eigen::VectorXd gradient_sol(3);
    gradient_sol(0) = jacobian_obj_sol.col(0).sum();
    gradient_sol(1) = jacobian_obj_sol.col(1).sum();
    gradient_sol(2) = jacobian_obj_sol.col(2).sum();

    Eigen::VectorXd gradient_non_lsq_sol = gradient_sol;

    testObjectiveAndEqualityAndInequalityJacobians(gradient_sol, gradient_non_lsq_sol, jacobian_lsq_obj_sol, jacobian_eq_sol, jacobian_ineq_sol, -1,
                                                   9, 9, 1e-6, nullptr, multipliers_eq.data(), multipliers_ineq.data());
}

TEST_F(TestHyperGraphOptimizationProblemVertexBased, jacobian_lsq)
{
    // second edge is a least-squares edge:
    // v_initial = [1 1]'
    // f = [4*1   5*1+6*1]' = [4 11]'
    // J(f) = [4 0; 5 6]
    // J(f^2) = 2 * f' * J(f) = 2 * [16+55 11*6] = [142 132]

    // Vertices
    vertices->addVertices({v1, v2});

    // Edges
    BaseEdge::Ptr edge1o(new Edge1T(edge1_fun, false, *v1, *v2));
    BaseEdge::Ptr edge1e(new Edge1T(edge1_fun, false, *v1, *v2));
    BaseEdge::Ptr edge1i(new Edge1T(edge1_fun, false, *v1, *v2));
    BaseEdge::Ptr edge2o(new Edge2T(edge2_fun, true, *v2));  // least squares edge -> error is squared
    edges->addEdges({edge1o}, {edge2o}, {edge1e}, {edge1i}, {});

    // specify gradient solution
    Eigen::VectorXd gradient_sol(3);
    gradient_sol << 1, 144, 135;

    Eigen::VectorXd gradient_non_lsq_sol(3);
    gradient_non_lsq_sol << 1, 2, 3;

    // specify jacobian solution
    Eigen::MatrixXd jacobian_sol(1, 3);
    jacobian_sol << 1, 2, 3;  //, 0, 142, 132;

    Eigen::MatrixXd jacobian_lsq_obj_sol(2, 3);
    jacobian_lsq_obj_sol << 0, 4, 0, 0, 5, 6;

    testObjectiveAndEqualityAndInequalityJacobians(gradient_sol, gradient_non_lsq_sol, jacobian_lsq_obj_sol, jacobian_sol, jacobian_sol, -1, 3, 3,
                                                   1e-6);

    testCombinedSparseJacobian(&jacobian_lsq_obj_sol, &jacobian_sol, &jacobian_sol, 1e-6);
}

TEST_F(TestHyperGraphOptimizationProblemVertexBased, jacobian_lsq_mixed)
{
    // Vertices
    vertices->addVertices({v1, v2});

    // Edges
    BaseEdge::Ptr edge1o(new Edge1T(edge1_fun, false, *v1, *v2));
    BaseEdge::Ptr edge1e(new Edge1T(edge1_fun, false, *v1, *v2));
    BaseEdge::Ptr edge1i(new Edge1T(edge1_fun, false, *v1, *v2));
    BaseEdge::Ptr edge2o(new Edge2T(edge2_fun, true, *v2));  // least squares edge -> error is squared

    auto e3_precompute = [this](const MixedEdge1T::VertexContainer& vertices) { mixed_edge1_precompute(vertices); };
    auto e3_obj        = [this](const MixedEdge1T::VertexContainer& vertices, Eigen::Ref<Eigen::VectorXd> values_obj) {
        mixed_edge1_obj(vertices, values_obj);
    };
    auto e3_eq = [this](const MixedEdge1T::VertexContainer& vertices, Eigen::Ref<Eigen::VectorXd> values_eq) { mixed_edge1_eq(vertices, values_eq); };
    auto e3_ineq = [this](const MixedEdge1T::VertexContainer& vertices, Eigen::Ref<Eigen::VectorXd> values_ineq) {
        mixed_edge1_ineq(vertices, values_ineq);
    };
    MixedEdge1T::Ptr edge3(new MixedEdge1T(1, 2, 3, e3_precompute, e3_obj, e3_eq, e3_ineq, *v1, *v2));
    edge3->setObjectiveLsqForm(true);  // set lsq form for objective

    edges->addEdges({edge1o, edge2o}, {edge1e}, {edge1i}, {edge3});

    // specify jacobian solution
    Eigen::MatrixXd jacobian_obj_sol(3, 3);
    jacobian_obj_sol << 1, 2, 3, 0, 142, 132, 2 * 12 * 2, 2 * 12 * 4, 2 * 12 * 6;

    Eigen::MatrixXd jacobian_lsq_obj_sol(3, 3);
    jacobian_lsq_obj_sol << 0, 4, 0, 0, 5, 6, 2, 4, 6;

    Eigen::MatrixXd jacobian_eq_sol(3, 3);
    // jacobian_eq_sol << 1, 2, 3, 0, 142, 132, 8 * 8 * 2, 2 * 22 * 10, 2 * 22 * 12; // lsq edge case
    jacobian_eq_sol << 1, 2, 3, 8, 0, 0, 0, 10, 12;

    Eigen::MatrixXd jacobian_ineq_sol(4, 3);
    // jacobian_ineq_sol << 1, 2, 3, 0, 142, 132, 2 * 2 * (1 - 5) * 2, 2 * 2 * (1 + 3) * 2, 2 * 2 * 2; // lsq edge case
    jacobian_ineq_sol << 1, 2, 3, 2, 0, 0, 0, 2, 0, 0, 0, 2;

    // specify gradient solution
    Eigen::VectorXd gradient_sol(3);
    gradient_sol(0) = jacobian_obj_sol.col(0).sum();
    gradient_sol(1) = jacobian_obj_sol.col(1).sum();
    gradient_sol(2) = jacobian_obj_sol.col(2).sum();

    Eigen::VectorXd gradient_non_lsq_sol(3);
    gradient_non_lsq_sol << 1, 2, 3;

    testObjectiveAndEqualityAndInequalityJacobians(gradient_sol, gradient_non_lsq_sol, jacobian_lsq_obj_sol, jacobian_eq_sol, jacobian_ineq_sol, 6, 6,
                                                   6, 1e-6);
    testCombinedSparseJacobian(&jacobian_lsq_obj_sol, &jacobian_eq_sol, &jacobian_ineq_sol, 1e-6);
}

TEST_F(TestHyperGraphOptimizationProblemVertexBased, jacobian_lsq_fixed)
{
    // second edge is a least-squares edge:
    // v_initial = [1 1]'
    // f = [4*1   5*1+6*1]' = [4 11]'
    // J(f) = [4 0; 5 6]
    // J(f^2) = 2 * f' * J(f) = 2 * [16+55 11*6] = [142 132]

    // Vertices
    vertices->addVertices({v1, v2});

    v1->setFixed(0, true);
    v2->setFixed(1, true);

    // Edges
    BaseEdge::Ptr edge1o(new Edge1T(edge1_fun, false, *v1, *v2));
    BaseEdge::Ptr edge1e(new Edge1T(edge1_fun, false, *v1, *v2));
    BaseEdge::Ptr edge1i(new Edge1T(edge1_fun, false, *v1, *v2));
    BaseEdge::Ptr edge2o(new Edge2T(edge2_fun, true, *v2));  // least squares edge -> error is squared
    edges->addEdges({edge1o, edge2o}, {edge1e}, {edge1i}, {});

    // specify jacobian solution
    Eigen::MatrixXd jacobian_obj_sol(2, 1);
    jacobian_obj_sol << 2, 142;

    Eigen::MatrixXd jacobian_sol(1, 1);
    jacobian_sol << 2;

    Eigen::MatrixXd jacobian_lsq_obj_sol(2, 1);
    jacobian_lsq_obj_sol << 4, 5;

    // specify gradient solution
    Eigen::VectorXd gradient_sol(1);
    gradient_sol(0) = jacobian_obj_sol.col(0).sum();

    Eigen::VectorXd gradient_non_lsq_sol(1);
    gradient_non_lsq_sol << 2;

    testObjectiveAndEqualityAndInequalityJacobians(gradient_sol, gradient_non_lsq_sol, jacobian_lsq_obj_sol, jacobian_sol, jacobian_sol, 2, 1, 1,
                                                   1e-6);
    testCombinedSparseJacobian(&jacobian_lsq_obj_sol, &jacobian_sol, &jacobian_sol, 1e-6);
}

TEST_F(TestHyperGraphOptimizationProblemVertexBased, jacobian_lsq_fixed_mixed)
{
    // Vertices
    vertices->addVertices({v1, v2});

    v1->setFixed(0, true);
    v2->setFixed(1, true);

    // Edges
    BaseEdge::Ptr edge1o(new Edge1T(edge1_fun, false, *v1, *v2));
    BaseEdge::Ptr edge1e(new Edge1T(edge1_fun, false, *v1, *v2));
    BaseEdge::Ptr edge1i(new Edge1T(edge1_fun, false, *v1, *v2));
    BaseEdge::Ptr edge2o(new Edge2T(edge2_fun, true, *v2));  // least squares edge -> error is squared

    auto e3_precompute = [this](const MixedEdge1T::VertexContainer& vertices) { mixed_edge1_precompute(vertices); };
    auto e3_obj        = [this](const MixedEdge1T::VertexContainer& vertices, Eigen::Ref<Eigen::VectorXd> values_obj) {
        mixed_edge1_obj(vertices, values_obj);
    };
    auto e3_eq = [this](const MixedEdge1T::VertexContainer& vertices, Eigen::Ref<Eigen::VectorXd> values_eq) { mixed_edge1_eq(vertices, values_eq); };
    auto e3_ineq = [this](const MixedEdge1T::VertexContainer& vertices, Eigen::Ref<Eigen::VectorXd> values_ineq) {
        mixed_edge1_ineq(vertices, values_ineq);
    };
    MixedEdge1T::Ptr edge3(new MixedEdge1T(1, 2, 3, e3_precompute, e3_obj, e3_eq, e3_ineq, *v1, *v2));
    edge3->setObjectiveLsqForm(true);  // set lsq form for objective

    edges->addEdges({edge1o, edge2o}, {edge1e}, {edge1i}, {edge3});

    // specify jacobian solution
    Eigen::MatrixXd jacobian_obj_sol(3, 1);
    jacobian_obj_sol << 2, 142, 2 * 12 * 4;

    Eigen::MatrixXd jacobian_lsq_obj_sol(3, 1);
    jacobian_lsq_obj_sol << 4, 5, 4;

    Eigen::MatrixXd jacobian_eq_sol(3, 1);
    jacobian_eq_sol << 2, 0, 10;

    Eigen::MatrixXd jacobian_ineq_sol(4, 1);
    jacobian_ineq_sol << 2, 0, 2, 0;

    // specify gradient solution
    Eigen::VectorXd gradient_sol(1);
    gradient_sol(0) = jacobian_obj_sol.col(0).sum();

    Eigen::VectorXd gradient_non_lsq_sol(1);
    gradient_non_lsq_sol << 2;

    testObjectiveAndEqualityAndInequalityJacobians(gradient_sol, gradient_non_lsq_sol, jacobian_lsq_obj_sol, jacobian_eq_sol, jacobian_ineq_sol, 3, 2,
                                                   2, 1e-6);
    testCombinedSparseJacobian(&jacobian_lsq_obj_sol, &jacobian_eq_sol, &jacobian_ineq_sol, 1e-6);
}

TEST_F(TestHyperGraphOptimizationProblemVertexBased, jacobian_lsq_modified_values)
{
    // second edge is a least-squares edge:
    // v_initial = [1 1]'
    // f = [4*1   5*1+6*1]' = [4 11]'
    // J(f) = [4 0; 5 6]
    // J(f^2) = 2 * f' * J(f) = 2 * [16+55 11*6] = [142 132]

    // Vertices
    vertices->addVertices({v1, v2});

    // Edges
    BaseEdge::Ptr edge1o(new Edge1T(edge1_fun, false, *v1, *v2));
    BaseEdge::Ptr edge1e(new Edge1T(edge1_fun, false, *v1, *v2));
    BaseEdge::Ptr edge1i(new Edge1T(edge1_fun, false, *v1, *v2));
    BaseEdge::Ptr edge2o(new Edge2T(edge2_fun, true, *v2));  // least squares edge -> error is squared
    edges->addEdges({edge1o, edge2o}, {edge1e}, {edge1i}, {});

    std::array<double, 2> multipliers_lsq_obj = {{2.0, 3.0}};
    std::array<double, 1> multipliers_constr  = {{4.0}};

    // specify jacobian solution
    // Eigen::MatrixXd jacobian_sol(2, 3);
    // jacobian_sol << 1, 2, 3, 0, 142, 132;
    // jacobian_sol.row(0) *= multipliers[0];
    // jacobian_sol.row(1) *= multipliers[1];

    // specify jacobian solution
    Eigen::MatrixXd jacobian_lsq_obj_sol(2, 3);
    jacobian_lsq_obj_sol << 0, 8, 0, 0, 15, 18;

    Eigen::MatrixXd jacobian_eq_sol(1, 3);
    jacobian_eq_sol << 4, 8, 12;

    Eigen::MatrixXd jacobian_ineq_sol(1, 3);
    jacobian_ineq_sol << 4, 8, 12;

    // specify gradient solution
    Eigen::VectorXd gradient_sol(3);
    gradient_sol << 1, 144, 135;

    Eigen::VectorXd gradient_non_lsq_sol(3);
    gradient_non_lsq_sol << 1, 2, 3;

    testObjectiveAndEqualityAndInequalityJacobians(gradient_sol, gradient_non_lsq_sol, jacobian_lsq_obj_sol, jacobian_eq_sol, jacobian_ineq_sol, 3, 3,
                                                   3, 1e-6, multipliers_lsq_obj.data(), multipliers_constr.data(), multipliers_constr.data());
}

TEST_F(TestHyperGraphOptimizationProblemVertexBased, jacobian_lsq_modified_values_mixed)
{
    // second edge is a least-squares edge:
    // v_initial = [1 1]'
    // f = [4*1   5*1+6*1]' = [4 11]'
    // J(f) = [4 0; 5 6]
    // J(f^2) = 2 * f' * J(f) = 2 * [16+55 11*6] = [142 132]

    // Vertices
    vertices->addVertices({v1, v2});

    // Edges
    BaseEdge::Ptr edge1o(new Edge1T(edge1_fun, false, *v1, *v2));
    BaseEdge::Ptr edge1e(new Edge1T(edge1_fun, false, *v1, *v2));
    BaseEdge::Ptr edge1i(new Edge1T(edge1_fun, false, *v1, *v2));
    BaseEdge::Ptr edge2o(new Edge2T(edge2_fun, true, *v2));  // least squares edge -> error is squared

    auto e3_precompute = [this](const MixedEdge1T::VertexContainer& vertices) { mixed_edge1_precompute(vertices); };
    auto e3_obj        = [this](const MixedEdge1T::VertexContainer& vertices, Eigen::Ref<Eigen::VectorXd> values_obj) {
        mixed_edge1_obj(vertices, values_obj);
    };
    auto e3_eq = [this](const MixedEdge1T::VertexContainer& vertices, Eigen::Ref<Eigen::VectorXd> values_eq) { mixed_edge1_eq(vertices, values_eq); };
    auto e3_ineq = [this](const MixedEdge1T::VertexContainer& vertices, Eigen::Ref<Eigen::VectorXd> values_ineq) {
        mixed_edge1_ineq(vertices, values_ineq);
    };
    MixedEdge1T::Ptr edge3(new MixedEdge1T(1, 2, 3, e3_precompute, e3_obj, e3_eq, e3_ineq, *v1, *v2));
    edge3->setObjectiveLsqForm(true);  // set lsq form for objective

    edges->addEdges({edge1o, edge2o}, {edge1e}, {edge1i}, {edge3});

    std::array<double, 3> multipliers_lsq_obj = {{2.0, 3.0, 4.0}};
    std::array<double, 3> multipliers_eq      = {{2.0, 3.0, 5.0}};
    std::array<double, 4> multipliers_ineq    = {{2.0, 3.0, 6.0, 7.0}};

    // specify jacobian solution
    //    Eigen::MatrixXd jacobian_obj_sol(3, 3);
    //    jacobian_obj_sol << 1, 2, 3, 0, 142, 132, 2 * 12 * 2, 2 * 12 * 4, 2 * 12 * 6;
    //    jacobian_obj_sol.row(0) *= multipliers_obj[0];
    //    jacobian_obj_sol.row(1) *= multipliers_obj[1];
    //    jacobian_obj_sol.row(2) *= multipliers_obj[2];

    Eigen::MatrixXd jacobian_lsq_obj_sol(3, 3);
    jacobian_lsq_obj_sol << 0, 8, 0, 0, 15, 18, 8, 16, 24;

    Eigen::MatrixXd jacobian_eq_sol(3, 3);
    // jacobian_eq_sol << 1, 2, 3, 0, 142, 132, 8 * 8 * 2, 2 * 22 * 10, 2 * 22 * 12;
    jacobian_eq_sol << 2, 4, 6, 24, 0, 0, 0, 50, 60;
    // jacobian_eq_sol.row(0) *= multipliers_eq[0];
    // jacobian_eq_sol.row(1) *= multipliers_eq[1];
    // jacobian_eq_sol.row(2) *= multipliers_eq[2];

    Eigen::MatrixXd jacobian_ineq_sol(4, 3);
    // jacobian_ineq_sol << 1, 2, 3, 0, 142, 132, 2 * 2 * (1 - 5) * 2, 2 * 2 * (1 + 3) * 2, 2 * 2 * 2;
    jacobian_ineq_sol << 2, 4, 6, 6, 0, 0, 0, 12, 0, 0, 0, 14;
    // jacobian_ineq_sol.row(0) *= multipliers_ineq[0];
    // jacobian_ineq_sol.row(1) *= multipliers_ineq[1];
    // jacobian_ineq_sol.row(2) *= multipliers_ineq[2];

    // specify gradient solution
    Eigen::VectorXd gradient_sol(3);
    gradient_sol << 49, 240, 279;

    Eigen::VectorXd gradient_non_lsq_sol(3);
    gradient_non_lsq_sol << 1, 2, 3;

    testObjectiveAndEqualityAndInequalityJacobians(gradient_sol, gradient_non_lsq_sol, jacobian_lsq_obj_sol, jacobian_eq_sol, jacobian_ineq_sol, 6, 6,
                                                   6, 1e-6, multipliers_lsq_obj.data(), multipliers_eq.data(), multipliers_ineq.data());
}

TEST_F(TestHyperGraphOptimizationProblemVertexBased, jacobian_ineq_active_only)
{
    v1->values()[0] = -2;
    v2->values()[0] = -1;
    v2->values()[1] = 2;

    // Vertices
    vertices->addVertices({v1, v2});

    // Edges
    BaseEdge::Ptr edge1(new Edge1T(edge1_fun, false, *v1, *v2));
    Eigen::VectorXd values1(edge1->getDimension());
    edge1->computeValues(values1);
    bool active1 = values1[0] >= 0;
    edges->addInequalityEdge(edge1);

    BaseEdge::Ptr edge2(new Edge2T(edge2_fun, false, *v2));
    Eigen::VectorXd values2(edge2->getDimension());
    edge2->computeValues(values2);
    bool active2 = values2[0] >= 0;
    bool active3 = values2[1] >= 0;
    edges->addInequalityEdge(edge2);

    // make sure to have both active and inactive constraints
    EXPECT_TRUE(active1 || active2 || active3);
    EXPECT_FALSE(active1 && active2 && active3);

    // specify jacobian solution
    Eigen::MatrixXd jacobian_sol(3, 3);
    jacobian_sol << 1, 2, 3, 0, 4, 0, 0, 5, 6;
    if (!active1) jacobian_sol.row(0).setZero();
    if (!active2) jacobian_sol.row(1).setZero();
    if (!active3) jacobian_sol.row(2).setZero();

    // compute dense jacobian
    Eigen::MatrixXd jacobian(optim.getInequalityDimension(), optim.getParameterDimension());
    jacobian.setZero();
    optim.computeDenseJacobianActiveInequalities(jacobian);
    EXPECT_EQ_MATRIX(jacobian, jacobian_sol, 1e-7);

    // combined dense jacobian
    Eigen::MatrixXd jacob_lsq_obj(optim.getLsqObjectiveDimension(), optim.getParameterDimension());
    Eigen::MatrixXd jacob_eq(optim.getEqualityDimension(), optim.getParameterDimension());
    Eigen::VectorXd gradient_non_lsq_obj(optim.getParameterDimension());
    jacobian.setZero();
    optim.computeDenseJacobians(gradient_non_lsq_obj, jacob_lsq_obj, jacob_eq, jacobian, nullptr, nullptr, nullptr, true);
    EXPECT_EQ_MATRIX(jacobian, jacobian_sol, 1e-7);

    // compute sparse jacobian
    Eigen::SparseMatrix<double> sparse_jacobian(optim.getInequalityDimension(), optim.getParameterDimension());
    optim.computeSparseJacobianActiveInequalities(sparse_jacobian);
    EXPECT_EQ_MATRIX(sparse_jacobian, jacobian_sol, 1e-6);

    // combined sparse jacobian
    Eigen::SparseMatrix<double> sparse_jacob_obj(optim.getLsqObjectiveDimension(), optim.getParameterDimension());
    Eigen::SparseMatrix<double> sparse_jacob_eq(optim.getEqualityDimension(), optim.getParameterDimension());
    sparse_jacobian.setZero();
    optim.computeSparseJacobians(sparse_jacob_obj, sparse_jacob_eq, sparse_jacobian, nullptr, nullptr, nullptr, true);
    EXPECT_EQ_MATRIX(sparse_jacobian, jacobian_sol, 1e-7);

    // compute other sparse jacobian
    int nnz = optim.computeSparseJacobianInequalitiesNNZ();
    EXPECT_LE(nnz, jacobian_sol.rows() * jacobian_sol.cols());
    EXPECT_GE(nnz, 6);
    computeSparseJacobianActiveInequalitiesViaTriplets(nnz, jacobian.rows(), jacobian.cols(), sparse_jacobian);
    EXPECT_EQ_MATRIX(sparse_jacobian, jacobian_sol, 1e-6);

    // combined sparse jacobian
    computeSparseJacobiansViaTriplets(0, jacob_lsq_obj.rows(), jacob_lsq_obj.cols(), sparse_jacob_obj, 0, jacob_eq.rows(), jacob_eq.cols(),
                                      sparse_jacob_eq, nnz, jacobian.rows(), jacobian.cols(), sparse_jacobian, nullptr, nullptr, nullptr, true);
    EXPECT_EQ_MATRIX(sparse_jacobian, jacobian_sol, 1e-6);

    testCombinedSparseJacobian(nullptr, nullptr, &jacobian_sol, 1e-6, nullptr, true);
}

TEST_F(TestHyperGraphOptimizationProblemVertexBased, jacobian_ineq_active_only_mixed)
{
    v1->values()[0] = -2;
    v2->values()[0] = -1;
    v2->values()[1] = 2;

    // Vertices
    vertices->addVertices({v1, v2});

    // Edges
    BaseEdge::Ptr edge1(new Edge1T(edge1_fun, false, *v1, *v2));

    Eigen::VectorXd values1(edge1->getDimension());
    edge1->computeValues(values1);
    bool active1 = values1[0] >= 0;
    edges->addInequalityEdge(edge1);

    BaseEdge::Ptr edge2(new Edge2T(edge2_fun, false, *v2));
    Eigen::VectorXd values2(edge2->getDimension());
    edge2->computeValues(values2);
    bool active2 = values2[0] >= 0;
    bool active3 = values2[1] >= 0;
    edges->addInequalityEdge(edge2);

    auto e3_precompute = [this](const MixedEdge1T::VertexContainer& vertices) { mixed_edge1_precompute(vertices); };
    auto e3_obj        = [this](const MixedEdge1T::VertexContainer& vertices, Eigen::Ref<Eigen::VectorXd> values_obj) {
        mixed_edge1_obj(vertices, values_obj);
    };
    auto e3_eq = [this](const MixedEdge1T::VertexContainer& vertices, Eigen::Ref<Eigen::VectorXd> values_eq) { mixed_edge1_eq(vertices, values_eq); };
    auto e3_ineq = [this](const MixedEdge1T::VertexContainer& vertices, Eigen::Ref<Eigen::VectorXd> values_ineq) {
        mixed_edge1_ineq(vertices, values_ineq);
    };
    MixedEdge1T::Ptr edge3(new MixedEdge1T(1, 2, 3, e3_precompute, e3_obj, e3_eq, e3_ineq, *v1, *v2));
    edge3->precompute();
    Eigen::VectorXd values3(edge3->getDimension());
    edge3->computeInequalityValues(values3);
    bool active4 = values3[0] >= 0;
    bool active5 = values3[1] >= 0;
    bool active6 = values3[2] >= 0;
    edges->addMixedEdge(edge3);

    // make sure to have both active and inactive constraints in the mixed constraint set
    EXPECT_TRUE(active4 || active5 || active6);
    EXPECT_FALSE(active4 && active5 && active6);

    // specify jacobian solution
    Eigen::MatrixXd jacobian_sol(6, 3);
    jacobian_sol << 1, 2, 3, 0, 4, 0, 0, 5, 6, 2, 0, 0, 0, 2, 0, 0, 0, 2;
    if (!active1) jacobian_sol.row(0).setZero();
    if (!active2) jacobian_sol.row(1).setZero();
    if (!active3) jacobian_sol.row(2).setZero();
    if (!active4) jacobian_sol.row(3).setZero();
    if (!active5) jacobian_sol.row(4).setZero();
    if (!active6) jacobian_sol.row(5).setZero();

    // compute dense jacobian
    Eigen::MatrixXd jacobian(optim.getInequalityDimension(), optim.getParameterDimension());
    jacobian.setZero();
    optim.computeDenseJacobianActiveInequalities(jacobian);
    EXPECT_EQ_MATRIX(jacobian, jacobian_sol, 1e-7);

    // combined dense jacobian
    Eigen::MatrixXd jacob_lsq_obj(optim.getLsqObjectiveDimension(), optim.getParameterDimension());
    Eigen::MatrixXd jacob_eq(optim.getEqualityDimension(), optim.getParameterDimension());
    Eigen::VectorXd gradient_non_lsq_obj(optim.getParameterDimension());
    jacobian.setZero();
    optim.computeDenseJacobians(gradient_non_lsq_obj, jacob_lsq_obj, jacob_eq, jacobian, nullptr, nullptr, nullptr, true);
    EXPECT_EQ_MATRIX(jacobian, jacobian_sol, 1e-7);

    // compute sparse jacobian
    Eigen::SparseMatrix<double> sparse_jacobian(optim.getInequalityDimension(), optim.getParameterDimension());
    optim.computeSparseJacobianActiveInequalities(sparse_jacobian);
    EXPECT_EQ_MATRIX(sparse_jacobian, jacobian_sol, 1e-6);

    // combined sparse jacobian
    Eigen::SparseMatrix<double> sparse_jacob_obj(1, 3), sparse_jacob_eq(2, 3);
    sparse_jacobian.setZero();
    optim.computeSparseJacobians(sparse_jacob_obj, sparse_jacob_eq, sparse_jacobian, nullptr, nullptr, nullptr, true);
    EXPECT_EQ_MATRIX(sparse_jacobian, jacobian_sol, 1e-7);

    // compute other sparse jacobian
    int nnz = optim.computeSparseJacobianInequalitiesNNZ();
    EXPECT_LE(nnz, jacobian_sol.rows() * jacobian_sol.cols());
    EXPECT_GE(nnz, 6);
    computeSparseJacobianActiveInequalitiesViaTriplets(nnz, jacobian.rows(), jacobian.cols(), sparse_jacobian);
    EXPECT_EQ_MATRIX(sparse_jacobian, jacobian_sol, 1e-6);

    // combined sparse jacobian
    computeSparseJacobiansViaTriplets(optim.computeSparseJacobianLsqObjectiveNNZ(), jacob_lsq_obj.rows(), jacob_lsq_obj.cols(), sparse_jacob_obj,
                                      optim.computeSparseJacobianEqualitiesNNZ(), jacob_eq.rows(), jacob_eq.cols(), sparse_jacob_eq, nnz,
                                      jacobian.rows(), jacobian.cols(), sparse_jacobian, nullptr, nullptr, nullptr, true);
    EXPECT_EQ_MATRIX(sparse_jacobian, jacobian_sol, 1e-6);

    // combined single sparse jacobian
    testCombinedSparseJacobian(nullptr, nullptr, &jacobian_sol, 1e-6, nullptr, true);
}

TEST_F(TestHyperGraphOptimizationProblemVertexBased, jacobian_ineq_active_only2)
{
    // Add variables
    pFVectorVertex::Ptr x1 = std::make_shared<pFVectorVertex>(1);
    pFVectorVertex::Ptr x2 = std::make_shared<pFVectorVertex>(1);

    vertices->addVertices({x1, x2});

    using EdgeT   = EdgeGenericScalarFun<pFVectorVertex, pFVectorVertex>;
    auto edge_fun = [](const EdgeT::VertexContainer& vertices) {
        double x1 = vertices.at(0)->getData()[0];
        double x2 = vertices.at(1)->getData()[0];
        return x2 - x1 + 10;
    };
    BaseEdge::Ptr edge(new EdgeT(edge_fun, false, *x1, *x2));
    edges->addInequalityEdge(edge);

    x1->values()[0] = 2;
    x2->values()[0] = -10;

    Eigen::Matrix<double, 1, 1> values;
    optim.computeValuesActiveInequality(values);

    EXPECT_DOUBLE_EQ(values[0], 0.0);

    // specify jacobian solution
    Eigen::MatrixXd jacobian_sol(1, 2);
    jacobian_sol << 0, 0;

    // compute dense jacobian
    Eigen::MatrixXd jacobian(optim.getInequalityDimension(), optim.getParameterDimension());
    jacobian.setZero();
    optim.computeDenseJacobianActiveInequalities(jacobian);
    EXPECT_EQ_MATRIX(jacobian, jacobian_sol, 1e-7);

    // combined dense jacobian
    Eigen::MatrixXd jacob_lsq_obj(optim.getLsqObjectiveDimension(), optim.getParameterDimension());
    Eigen::MatrixXd jacob_eq(optim.getEqualityDimension(), optim.getParameterDimension());
    Eigen::VectorXd gradient_non_lsq_obj(optim.getParameterDimension());
    jacobian.setZero();
    optim.computeDenseJacobians(gradient_non_lsq_obj, jacob_lsq_obj, jacob_eq, jacobian, nullptr, nullptr, nullptr, true);
    EXPECT_EQ_MATRIX(jacobian, jacobian_sol, 1e-7);

    // compute sparse jacobian
    Eigen::SparseMatrix<double> sparse_jacobian(optim.getInequalityDimension(), optim.getParameterDimension());
    optim.computeSparseJacobianActiveInequalities(sparse_jacobian);
    EXPECT_EQ_MATRIX(sparse_jacobian, jacobian_sol, 1e-6);

    // combined sparse jacobian
    Eigen::SparseMatrix<double> sparse_jacob_obj(optim.getLsqObjectiveDimension(), optim.getParameterDimension());
    Eigen::SparseMatrix<double> sparse_jacob_eq(optim.getEqualityDimension(), optim.getParameterDimension());
    sparse_jacobian.setZero();
    optim.computeSparseJacobians(sparse_jacob_obj, sparse_jacob_eq, sparse_jacobian, nullptr, nullptr, nullptr, true);
    EXPECT_EQ_MATRIX(sparse_jacobian, jacobian_sol, 1e-7);

    // compute other sparse jacobian
    int nnz = optim.computeSparseJacobianInequalitiesNNZ();
    EXPECT_LE(nnz, jacobian_sol.rows() * jacobian_sol.cols());
    EXPECT_GE(nnz, 0);  // 2
    computeSparseJacobianActiveInequalitiesViaTriplets(nnz, jacobian.rows(), jacobian.cols(), sparse_jacobian);
    EXPECT_EQ_MATRIX(sparse_jacobian, jacobian_sol, 1e-6);

    // combined sparse jacobian
    computeSparseJacobiansViaTriplets(0, jacob_lsq_obj.rows(), jacob_lsq_obj.cols(), sparse_jacob_obj, 0, jacob_eq.rows(), jacob_eq.cols(),
                                      sparse_jacob_eq, nnz, jacobian.rows(), jacobian.cols(), sparse_jacobian, nullptr, nullptr, nullptr, true);
    EXPECT_EQ_MATRIX(sparse_jacobian, jacobian_sol, 1e-6);

    // combined single matrix jajacobian
    testCombinedSparseJacobian(nullptr, nullptr, &jacobian_sol, 1e-6, nullptr, true);

    // Change values in order to activate inequalities
    x1->values()[0] = 2;
    x2->values()[0] = 20;

    optim.computeValuesActiveInequality(values);

    EXPECT_DOUBLE_EQ(values[0], 28.0);

    // specify jacobian solution
    jacobian_sol(0, 0) = -1.0;
    jacobian_sol(0, 1) = 1.0;

    // compute dense jacobian
    jacobian.setZero();
    optim.computeDenseJacobianActiveInequalities(jacobian);
    EXPECT_EQ_MATRIX(jacobian, jacobian_sol, 1e-7);

    // combined dense jacobian
    jacobian.setZero();
    optim.computeDenseJacobians(gradient_non_lsq_obj, jacob_lsq_obj, jacob_eq, jacobian, nullptr, nullptr, nullptr, true);
    EXPECT_EQ_MATRIX(jacobian, jacobian_sol, 1e-7);

    // compute sparse jacobian
    sparse_jacobian.setZero();
    optim.computeSparseJacobianActiveInequalities(sparse_jacobian);
    EXPECT_EQ_MATRIX(sparse_jacobian, jacobian_sol, 1e-6);

    // combined sparse jacobian
    sparse_jacobian.setZero();
    optim.computeSparseJacobians(sparse_jacob_obj, sparse_jacob_eq, sparse_jacobian, nullptr, nullptr, nullptr, true);
    EXPECT_EQ_MATRIX(sparse_jacobian, jacobian_sol, 1e-7);

    // compute other sparse jacobian
    nnz = optim.computeSparseJacobianInequalitiesNNZ();
    EXPECT_LE(nnz, jacobian_sol.rows() * jacobian_sol.cols());
    EXPECT_GE(nnz, 1);  // 2
    computeSparseJacobianActiveInequalitiesViaTriplets(nnz, jacobian.rows(), jacobian.cols(), sparse_jacobian);
    EXPECT_EQ_MATRIX(sparse_jacobian, jacobian_sol, 1e-6);

    computeSparseJacobiansViaTriplets(0, jacob_lsq_obj.rows(), jacob_lsq_obj.cols(), sparse_jacob_obj, 0, jacob_eq.rows(), jacob_eq.cols(),
                                      sparse_jacob_eq, nnz, jacobian.rows(), jacobian.cols(), sparse_jacobian, nullptr, nullptr, nullptr, true);
    EXPECT_EQ_MATRIX(sparse_jacobian, jacobian_sol, 1e-6);

    // combined single matrix jajacobian
    testCombinedSparseJacobian(nullptr, nullptr, &jacobian_sol, 1e-6, nullptr, true);
}

TEST_F(TestHyperGraphOptimizationProblemVertexBased, jacobian_ineq_active_only2_mixed)
{
    // Add variables
    pFVectorVertex::Ptr x1 = std::make_shared<pFVectorVertex>(1);
    pFVectorVertex::Ptr x2 = std::make_shared<pFVectorVertex>(1);

    vertices->addVertices({x1, x2});

    auto e_precompute = [this](const MixedEdge1T::VertexContainer& vertices) {};
    auto e_obj        = [this](const MixedEdge1T::VertexContainer& vertices, Eigen::Ref<Eigen::VectorXd> values_obj) {};
    auto e_eq         = [this](const MixedEdge1T::VertexContainer& vertices, Eigen::Ref<Eigen::VectorXd> values_eq) {};
    auto e_ineq       = [this](const MixedEdge1T::VertexContainer& vertices, Eigen::Ref<Eigen::VectorXd> values_ineq) {
        double x1      = vertices.at(0)->getData()[0];
        double x2      = vertices.at(1)->getData()[0];
        values_ineq[0] = x2 - x1 + 10;
    };
    MixedEdge1T::Ptr edge(new MixedEdge1T(0, 0, 1, e_precompute, e_obj, e_eq, e_ineq, *x1, *x2));

    edges->addMixedEdge(edge);

    x1->values()[0] = 2;
    x2->values()[0] = -10;

    Eigen::Matrix<double, 1, 1> values;
    optim.computeValuesActiveInequality(values);

    EXPECT_DOUBLE_EQ(values[0], 0.0);

    // specify jacobian solution
    Eigen::MatrixXd jacobian_sol(1, 2);
    jacobian_sol << 0, 0;

    // compute dense jacobian
    Eigen::MatrixXd jacobian(optim.getInequalityDimension(), optim.getParameterDimension());
    jacobian.setZero();
    optim.computeDenseJacobianActiveInequalities(jacobian);
    EXPECT_EQ_MATRIX(jacobian, jacobian_sol, 1e-7);

    // combined dense jacobian
    Eigen::MatrixXd jacob_obj(optim.getLsqObjectiveDimension(), optim.getParameterDimension());
    Eigen::MatrixXd jacob_eq(optim.getEqualityDimension(), optim.getParameterDimension());
    Eigen::VectorXd gradient_non_lsq_obj(optim.getParameterDimension());
    jacobian.setZero();
    optim.computeDenseJacobians(gradient_non_lsq_obj, jacob_obj, jacob_eq, jacobian, nullptr, nullptr, nullptr, true);
    EXPECT_EQ_MATRIX(jacobian, jacobian_sol, 1e-7);

    // compute sparse jacobian
    Eigen::SparseMatrix<double> sparse_jacobian(optim.getInequalityDimension(), optim.getParameterDimension());
    optim.computeSparseJacobianActiveInequalities(sparse_jacobian);
    EXPECT_EQ_MATRIX(sparse_jacobian, jacobian_sol, 1e-6);

    // combined sparse jacobian
    Eigen::SparseMatrix<double> sparse_jacob_obj(optim.getLsqObjectiveDimension(), optim.getParameterDimension());
    Eigen::SparseMatrix<double> sparse_jacob_eq(optim.getEqualityDimension(), optim.getParameterDimension());
    sparse_jacobian.setZero();
    optim.computeSparseJacobians(sparse_jacob_obj, sparse_jacob_eq, sparse_jacobian, nullptr, nullptr, nullptr, true);
    EXPECT_EQ_MATRIX(sparse_jacobian, jacobian_sol, 1e-7);

    // compute other sparse jacobian
    int nnz = optim.computeSparseJacobianInequalitiesNNZ();
    EXPECT_LE(nnz, jacobian_sol.rows() * jacobian_sol.cols());
    EXPECT_GE(nnz, 0);  // 2
    computeSparseJacobianActiveInequalitiesViaTriplets(nnz, jacobian.rows(), jacobian.cols(), sparse_jacobian);
    EXPECT_EQ_MATRIX(sparse_jacobian, jacobian_sol, 1e-6);

    // combined sparse jacobian
    computeSparseJacobiansViaTriplets(0, jacob_obj.rows(), jacob_obj.cols(), sparse_jacob_obj, 0, jacob_eq.rows(), jacob_eq.cols(), sparse_jacob_eq,
                                      nnz, jacobian.rows(), jacobian.cols(), sparse_jacobian, nullptr, nullptr, nullptr, true);
    EXPECT_EQ_MATRIX(sparse_jacobian, jacobian_sol, 1e-6);

    // combined single sparse jacobian
    testCombinedSparseJacobian(nullptr, nullptr, &jacobian_sol, 1e-6, nullptr, true);

    // Change values in order to activate inequalities
    x1->values()[0] = 2;
    x2->values()[0] = 20;

    optim.computeValuesActiveInequality(values);

    EXPECT_DOUBLE_EQ(values[0], 28.0);

    // specify jacobian solution
    jacobian_sol(0, 0) = -1.0;
    jacobian_sol(0, 1) = 1.0;

    // compute dense jacobian
    jacobian.setZero();
    optim.computeDenseJacobianActiveInequalities(jacobian);
    EXPECT_EQ_MATRIX(jacobian, jacobian_sol, 1e-7);

    // combined dense jacobian
    jacobian.setZero();
    optim.computeDenseJacobians(gradient_non_lsq_obj, jacob_obj, jacob_eq, jacobian, nullptr, nullptr, nullptr, true);
    EXPECT_EQ_MATRIX(jacobian, jacobian_sol, 1e-7);

    // compute sparse jacobian
    sparse_jacobian.setZero();
    optim.computeSparseJacobianActiveInequalities(sparse_jacobian);
    EXPECT_EQ_MATRIX(sparse_jacobian, jacobian_sol, 1e-6);

    // combined sparse jacobian
    sparse_jacobian.setZero();
    optim.computeSparseJacobians(sparse_jacob_obj, sparse_jacob_eq, sparse_jacobian, nullptr, nullptr, nullptr, true);
    EXPECT_EQ_MATRIX(sparse_jacobian, jacobian_sol, 1e-7);

    // compute other sparse jacobian
    nnz = optim.computeSparseJacobianInequalitiesNNZ();
    EXPECT_LE(nnz, jacobian_sol.rows() * jacobian_sol.cols());
    EXPECT_GE(nnz, 1);  // 2
    computeSparseJacobianActiveInequalitiesViaTriplets(nnz, jacobian.rows(), jacobian.cols(), sparse_jacobian);
    EXPECT_EQ_MATRIX(sparse_jacobian, jacobian_sol, 1e-6);

    computeSparseJacobiansViaTriplets(0, jacob_obj.rows(), jacob_obj.cols(), sparse_jacob_obj, 0, jacob_eq.rows(), jacob_eq.cols(), sparse_jacob_eq,
                                      nnz, jacobian.rows(), jacobian.cols(), sparse_jacobian, nullptr, nullptr, nullptr, true);
    EXPECT_EQ_MATRIX(sparse_jacobian, jacobian_sol, 1e-6);

    // combined single sparse jacobian
    testCombinedSparseJacobian(nullptr, nullptr, &jacobian_sol, 1e-6, nullptr, true);
}

TEST_F(TestHyperGraphOptimizationProblemVertexBased, jacobian_ineq_active_only2_fixed)
{
    // Add variables
    pFVectorVertex::Ptr x1 = std::make_shared<pFVectorVertex>(1);
    pFVectorVertex::Ptr x2 = std::make_shared<pFVectorVertex>(1);

    vertices->addVertices({x1, x2});

    x1->setFixed(0, true);

    using EdgeT   = EdgeGenericScalarFun<pFVectorVertex, pFVectorVertex>;
    auto edge_fun = [](const EdgeT::VertexContainer& vertices) {
        double x1 = vertices.at(0)->getData()[0];
        double x2 = vertices.at(1)->getData()[0];
        return x2 - x1 + 10;
    };
    BaseEdge::Ptr edge(new EdgeT(edge_fun, false, *x1, *x2));
    edges->addInequalityEdge(edge);

    x1->values()[0] = 2;
    x2->values()[0] = -10;

    Eigen::Matrix<double, 1, 1> values;
    optim.computeValuesActiveInequality(values);

    EXPECT_DOUBLE_EQ(values[0], 0.0);

    // specify jacobian solution
    Eigen::MatrixXd jacobian_sol(1, 1);
    jacobian_sol << 0;

    // compute dense jacobian
    Eigen::MatrixXd jacobian(optim.getInequalityDimension(), optim.getParameterDimension());
    jacobian.setZero();
    optim.computeDenseJacobianActiveInequalities(jacobian);
    EXPECT_EQ_MATRIX(jacobian, jacobian_sol, 1e-6);

    // combined dense jacobian
    Eigen::MatrixXd jacob_obj(optim.getLsqObjectiveDimension(), optim.getParameterDimension());
    Eigen::MatrixXd jacob_eq(optim.getEqualityDimension(), optim.getParameterDimension());
    Eigen::VectorXd gradient_non_lsq_obj(optim.getParameterDimension());
    jacobian.setZero();
    optim.computeDenseJacobians(gradient_non_lsq_obj, jacob_obj, jacob_eq, jacobian, nullptr, nullptr, nullptr, true);
    EXPECT_EQ_MATRIX(jacobian, jacobian_sol, 1e-7);

    // compute sparse jacobian
    Eigen::SparseMatrix<double> sparse_jacobian(optim.getInequalityDimension(), optim.getParameterDimension());
    optim.computeSparseJacobianActiveInequalities(sparse_jacobian);
    EXPECT_EQ_MATRIX(sparse_jacobian, jacobian_sol, 1e-6);

    // combined sparse jacobian
    Eigen::SparseMatrix<double> sparse_jacob_obj(optim.getLsqObjectiveDimension(), optim.getParameterDimension());
    Eigen::SparseMatrix<double> sparse_jacob_eq(optim.getEqualityDimension(), optim.getParameterDimension());
    sparse_jacobian.setZero();
    optim.computeSparseJacobians(sparse_jacob_obj, sparse_jacob_eq, sparse_jacobian, nullptr, nullptr, nullptr, true);
    EXPECT_EQ_MATRIX(sparse_jacobian, jacobian_sol, 1e-7);

    // compute other sparse jacobian
    int nnz = optim.computeSparseJacobianInequalitiesNNZ();
    EXPECT_LE(nnz, jacobian_sol.rows() * jacobian_sol.cols());
    EXPECT_GE(nnz, 0);  // 2
    computeSparseJacobianActiveInequalitiesViaTriplets(nnz, jacobian.rows(), jacobian.cols(), sparse_jacobian);
    EXPECT_EQ_MATRIX(sparse_jacobian, jacobian_sol, 1e-6);

    // combined sparse jacobian
    computeSparseJacobiansViaTriplets(0, jacob_obj.rows(), jacob_obj.cols(), sparse_jacob_obj, 0, jacob_eq.rows(), jacob_eq.cols(), sparse_jacob_eq,
                                      nnz, jacobian.rows(), jacobian.cols(), sparse_jacobian, nullptr, nullptr, nullptr, true);
    EXPECT_EQ_MATRIX(sparse_jacobian, jacobian_sol, 1e-6);

    // combined single matrix jajacobian
    testCombinedSparseJacobian(nullptr, nullptr, &jacobian_sol, 1e-6, nullptr, true);

    // Change values in order to activate inequalities
    x1->values()[0] = 2;
    x2->values()[0] = 20;

    optim.computeValuesActiveInequality(values);

    EXPECT_DOUBLE_EQ(values[0], 28.0);

    // specify jacobian solution
    jacobian_sol(0, 0) = 1.0;

    // compute dense jacobian
    jacobian.setZero();
    optim.computeDenseJacobianActiveInequalities(jacobian);
    EXPECT_EQ_MATRIX(jacobian, jacobian_sol, 1e-6);

    // combined dense jacobian
    jacobian.setZero();
    optim.computeDenseJacobians(gradient_non_lsq_obj, jacob_obj, jacob_eq, jacobian, nullptr, nullptr, nullptr, true);
    EXPECT_EQ_MATRIX(jacobian, jacobian_sol, 1e-7);

    // compute sparse jacobian
    sparse_jacobian.setZero();
    optim.computeSparseJacobianActiveInequalities(sparse_jacobian);
    EXPECT_EQ_MATRIX(sparse_jacobian, jacobian_sol, 1e-6);

    // combined sparse jacobian
    sparse_jacobian.setZero();
    optim.computeSparseJacobians(sparse_jacob_obj, sparse_jacob_eq, sparse_jacobian, nullptr, nullptr, nullptr, true);
    EXPECT_EQ_MATRIX(sparse_jacobian, jacobian_sol, 1e-7);

    // compute other sparse jacobian
    nnz = optim.computeSparseJacobianInequalitiesNNZ();
    EXPECT_LE(nnz, jacobian_sol.rows() * jacobian_sol.cols());
    EXPECT_GE(nnz, 1);  // 2
    computeSparseJacobianActiveInequalitiesViaTriplets(nnz, jacobian.rows(), jacobian.cols(), sparse_jacobian);
    EXPECT_EQ_MATRIX(sparse_jacobian, jacobian_sol, 1e-6);

    // combined sparse jacobian
    computeSparseJacobiansViaTriplets(0, jacob_obj.rows(), jacob_obj.cols(), sparse_jacob_obj, 0, jacob_eq.rows(), jacob_eq.cols(), sparse_jacob_eq,
                                      nnz, jacobian.rows(), jacobian.cols(), sparse_jacobian, nullptr, nullptr, nullptr, true);
    EXPECT_EQ_MATRIX(sparse_jacobian, jacobian_sol, 1e-6);

    // combined single matrix jajacobian
    testCombinedSparseJacobian(nullptr, nullptr, &jacobian_sol, 1e-6, nullptr, true);
}

TEST_F(TestHyperGraphOptimizationProblemVertexBased, jacobian_ineq_active_only2_fixed_mixed)
{
    // Add variables
    pFVectorVertex::Ptr x1 = std::make_shared<pFVectorVertex>(1);
    pFVectorVertex::Ptr x2 = std::make_shared<pFVectorVertex>(1);

    vertices->addVertices({x1, x2});

    x1->setFixed(0, true);

    auto e_precompute = [this](const MixedEdge1T::VertexContainer& vertices) {};
    auto e_obj        = [this](const MixedEdge1T::VertexContainer& vertices, Eigen::Ref<Eigen::VectorXd> values_obj) {};
    auto e_eq         = [this](const MixedEdge1T::VertexContainer& vertices, Eigen::Ref<Eigen::VectorXd> values_eq) {};
    auto e_ineq       = [this](const MixedEdge1T::VertexContainer& vertices, Eigen::Ref<Eigen::VectorXd> values_ineq) {
        double x1      = vertices.at(0)->getData()[0];
        double x2      = vertices.at(1)->getData()[0];
        values_ineq[0] = x2 - x1 + 10;
    };
    MixedEdge1T::Ptr edge(new MixedEdge1T(0, 0, 1, e_precompute, e_obj, e_eq, e_ineq, *x1, *x2));

    edges->addMixedEdge(edge);

    x1->values()[0] = 2;
    x2->values()[0] = -10;

    Eigen::Matrix<double, 1, 1> values;
    optim.computeValuesActiveInequality(values);

    EXPECT_DOUBLE_EQ(values[0], 0.0);

    // specify jacobian solution
    Eigen::MatrixXd jacobian_sol(1, 1);
    jacobian_sol << 0;

    // compute dense jacobian
    Eigen::MatrixXd jacobian(optim.getInequalityDimension(), optim.getParameterDimension());
    jacobian.setZero();
    optim.computeDenseJacobianActiveInequalities(jacobian);
    EXPECT_EQ_MATRIX(jacobian, jacobian_sol, 1e-6);

    // combined dense jacobian
    Eigen::MatrixXd jacob_obj(optim.getLsqObjectiveDimension(), optim.getParameterDimension());
    Eigen::MatrixXd jacob_eq(optim.getEqualityDimension(), optim.getParameterDimension());
    Eigen::VectorXd gradient_non_lsq_obj(optim.getParameterDimension());
    jacobian.setZero();
    optim.computeDenseJacobians(gradient_non_lsq_obj, jacob_obj, jacob_eq, jacobian, nullptr, nullptr, nullptr, true);
    EXPECT_EQ_MATRIX(jacobian, jacobian_sol, 1e-7);

    // compute sparse jacobian
    Eigen::SparseMatrix<double> sparse_jacobian(optim.getInequalityDimension(), optim.getParameterDimension());
    optim.computeSparseJacobianActiveInequalities(sparse_jacobian);
    EXPECT_EQ_MATRIX(sparse_jacobian, jacobian_sol, 1e-6);

    // combined sparse jacobian
    Eigen::SparseMatrix<double> sparse_jacob_obj(optim.getLsqObjectiveDimension(), optim.getParameterDimension());
    Eigen::SparseMatrix<double> sparse_jacob_eq(optim.getEqualityDimension(), optim.getParameterDimension());
    sparse_jacobian.setZero();
    optim.computeSparseJacobians(sparse_jacob_obj, sparse_jacob_eq, sparse_jacobian, nullptr, nullptr, nullptr, true);
    EXPECT_EQ_MATRIX(sparse_jacobian, jacobian_sol, 1e-7);

    // compute other sparse jacobian
    int nnz = optim.computeSparseJacobianInequalitiesNNZ();
    EXPECT_LE(nnz, jacobian_sol.rows() * jacobian_sol.cols());
    EXPECT_GE(nnz, 0);  // 2
    computeSparseJacobianActiveInequalitiesViaTriplets(nnz, jacobian.rows(), jacobian.cols(), sparse_jacobian);
    EXPECT_EQ_MATRIX(sparse_jacobian, jacobian_sol, 1e-6);

    // combined sparse jacobian
    computeSparseJacobiansViaTriplets(0, jacob_obj.rows(), jacob_obj.cols(), sparse_jacob_obj, 0, jacob_eq.rows(), jacob_eq.cols(), sparse_jacob_eq,
                                      nnz, jacobian.rows(), jacobian.cols(), sparse_jacobian, nullptr, nullptr, nullptr, true);
    EXPECT_EQ_MATRIX(sparse_jacobian, jacobian_sol, 1e-6);

    // combined single sparse jacobian
    testCombinedSparseJacobian(nullptr, nullptr, &jacobian_sol, 1e-6, nullptr, true);

    // Change values in order to activate inequalities
    x1->values()[0] = 2;
    x2->values()[0] = 20;

    optim.computeValuesActiveInequality(values);

    EXPECT_DOUBLE_EQ(values[0], 28.0);

    // specify jacobian solution
    jacobian_sol(0, 0) = 1.0;

    // compute dense jacobian
    jacobian.setZero();
    optim.computeDenseJacobianActiveInequalities(jacobian);
    EXPECT_EQ_MATRIX(jacobian, jacobian_sol, 1e-6);

    // combined dense jacobian
    jacobian.setZero();
    optim.computeDenseJacobians(gradient_non_lsq_obj, jacob_obj, jacob_eq, jacobian, nullptr, nullptr, nullptr, true);
    EXPECT_EQ_MATRIX(jacobian, jacobian_sol, 1e-7);

    // compute sparse jacobian
    sparse_jacobian.setZero();
    optim.computeSparseJacobianActiveInequalities(sparse_jacobian);
    EXPECT_EQ_MATRIX(sparse_jacobian, jacobian_sol, 1e-6);

    // combined sparse jacobian
    sparse_jacobian.setZero();
    optim.computeSparseJacobians(sparse_jacob_obj, sparse_jacob_eq, sparse_jacobian, nullptr, nullptr, nullptr, true);
    EXPECT_EQ_MATRIX(sparse_jacobian, jacobian_sol, 1e-7);

    // compute other sparse jacobian
    nnz = optim.computeSparseJacobianInequalitiesNNZ();
    EXPECT_LE(nnz, jacobian_sol.rows() * jacobian_sol.cols());
    EXPECT_GE(nnz, 1);  // 2
    computeSparseJacobianActiveInequalitiesViaTriplets(nnz, jacobian.rows(), jacobian.cols(), sparse_jacobian);
    EXPECT_EQ_MATRIX(sparse_jacobian, jacobian_sol, 1e-6);

    // combined sparse jacobian
    computeSparseJacobiansViaTriplets(0, jacob_obj.rows(), jacob_obj.cols(), sparse_jacob_obj, 0, jacob_eq.rows(), jacob_eq.cols(), sparse_jacob_eq,
                                      nnz, jacobian.rows(), jacobian.cols(), sparse_jacobian, nullptr, nullptr, nullptr, true);
    EXPECT_EQ_MATRIX(sparse_jacobian, jacobian_sol, 1e-6);

    // combined single sparse jacobian
    testCombinedSparseJacobian(nullptr, nullptr, &jacobian_sol, 1e-6, nullptr, true);
}

// TEST_F(TestHyperGraph, dense_user_defined_jacobian1)
//{
//    // first part of the jacobian user defined
//
//    // Vertices
//    graph.addVertex(&v1);
//    graph.addVertex(&v2);
//
//    v1.setFixed(0, true);
//    v2.setFixed(1, true);
//
//    // Edges
//    EdgeInterface::UPtr edge1(new Edge1T(edge1_fun, edge1_jacob, false, v1, v2));
//    graph.addObjectiveEdge(std::move(edge1));
//
//    EdgeInterface::UPtr edge2(new Edge2T(edge2_fun, false, v2));
//    graph.addObjectiveEdge(std::move(edge2));
//
//    Eigen::MatrixXd jacobian(graph.objectiveDimension(), graph.getParameterDimension());
//    jacobian.setZero();
//    graph.computeDenseJacobianObjective(jacobian);
//
//    Eigen::MatrixXd jacobian_sol(3, 1);
//    jacobian_sol << 2, 4, 5;
//
//    EXPECT_EQ_MATRIX(jacobian, jacobian_sol, 1e-7);
//}
//
// TEST_F(TestHyperGraph, dense_user_defined_jacobian1_fixed)
//{
//    // first part of the jacobian user defined
//
//    // Vertices
//    graph.addVertex(&v1);
//    graph.addVertex(&v2);
//
//    // Edges
//    EdgeInterface::UPtr edge1(new Edge1T(edge1_fun, edge1_jacob, false, v1, v2));
//    graph.addObjectiveEdge(std::move(edge1));
//
//    EdgeInterface::UPtr edge2(new Edge2T(edge2_fun, false, v2));
//    graph.addObjectiveEdge(std::move(edge2));
//
//    Eigen::MatrixXd jacobian(graph.objectiveDimension(), graph.getParameterDimension());
//    jacobian.setZero();
//    graph.computeDenseJacobianObjective(jacobian);
//
//    Eigen::MatrixXd jacobian_sol(3, 3);
//    jacobian_sol << 1, 2, 3, 0, 4, 0, 0, 5, 6;
//
//    EXPECT_EQ_MATRIX(jacobian, jacobian_sol, 1e-7);
//}
//
// TEST_F(TestHyperGraph, dense_user_defined_jacobian2)
//{
//    // second part of the jacobian user defined
//
//    // Vertices
//    graph.addVertex(&v1);
//    graph.addVertex(&v2);
//
//    // Edges
//    EdgeInterface::UPtr edge1(new Edge1T(edge1_fun, false, v1, v2));
//    graph.addObjectiveEdge(std::move(edge1));
//
//    EdgeInterface::UPtr edge2(new Edge2T(edge2_fun, edge2_jacob, false, v2));
//    graph.addObjectiveEdge(std::move(edge2));
//
//    Eigen::MatrixXd jacobian(graph.objectiveDimension(), graph.getParameterDimension());
//    jacobian.setZero();
//    graph.computeDenseJacobianObjective(jacobian);
//
//    Eigen::MatrixXd jacobian_sol(3, 3);
//    jacobian_sol << 1, 2, 3, 0, 4, 0, 0, 5, 6;
//
//    EXPECT_EQ_MATRIX(jacobian, jacobian_sol, 1e-7);
//}
//
// TEST_F(TestHyperGraph, dense_user_defined_jacobian2_fixed)
//{
//    // second part of the jacobian user defined
//
//    // Vertices
//    graph.addVertex(&v1);
//    graph.addVertex(&v2);
//
//    v1.setFixed(0, true);
//    v2.setFixed(1, true);
//
//    // Edges
//    EdgeInterface::UPtr edge1(new Edge1T(edge1_fun, false, v1, v2));
//    graph.addObjectiveEdge(std::move(edge1));
//
//    EdgeInterface::UPtr edge2(new Edge2T(edge2_fun, edge2_jacob, false, v2));
//    graph.addObjectiveEdge(std::move(edge2));
//
//    Eigen::MatrixXd jacobian(graph.objectiveDimension(), graph.getParameterDimension());
//    jacobian.setZero();
//    graph.computeDenseJacobianObjective(jacobian);
//
//    Eigen::MatrixXd jacobian_sol(3, 1);
//    jacobian_sol << 2, 4, 5;
//
//    EXPECT_EQ_MATRIX(jacobian, jacobian_sol, 1e-7);
//}
//
// TEST_F(TestHyperGraph, dense_user_defined_jacobian3)
//{
//    // complete jacobian user defined
//
//    // Vertices
//    graph.addVertex(&v1);
//    graph.addVertex(&v2);
//
//    // Edges
//    EdgeInterface::UPtr edge1(new Edge1T(edge1_fun, edge1_jacob, false, v1, v2));
//    graph.addObjectiveEdge(std::move(edge1));
//
//    EdgeInterface::UPtr edge2(new Edge2T(edge2_fun, edge2_jacob, false, v2));
//    graph.addObjectiveEdge(std::move(edge2));
//
//    Eigen::MatrixXd jacobian(graph.objectiveDimension(), graph.getParameterDimension());
//    jacobian.setZero();
//
//    graph.computeDenseJacobianObjective(jacobian);
//
//    Eigen::MatrixXd jacobian_sol(3, 3);
//    jacobian_sol << 1, 2, 3, 0, 4, 0, 0, 5, 6;
//
//    EXPECT_EQ_MATRIX(jacobian, jacobian_sol, 1e-7);
//}

//**********************************************************
// *
// *  HESSIAN COMPUTATION
// *
// *********************************************************/

TEST_F(TestHyperGraphOptimizationProblemVertexBased, hessian)
{
    // Vertices
    vertices->addVertices({v1, v2});

    // Edges
    BaseEdge::Ptr edge3o(new Edge3T(edge3_fun, false, *v1, *v2));
    BaseEdge::Ptr edge3e(new Edge3T(edge3_fun, false, *v1, *v2));
    BaseEdge::Ptr edge3i(new Edge3T(edge3_fun, false, *v1, *v2));
    edges->addEdges({edge3o}, {edge3e}, {edge3i}, {});

    // Eigen::MatrixXd jacobian(optim.objectiveDimension(), optim.getParameterDimension());
    // jacobian.setZero();
    // optim.computeDenseJacobianObjective(jacobian);

    // Eigen::MatrixXd hessian(optim.getParameterDimension(), optim.getParameterDimension());
    // hessian.setZero();
    // optim.computeDenseHessianObjective(jacobian, hessian); // NOT YET IMPLEMENTED
    // optim.computeDenseHessianObjective(hessian);

    Eigen::MatrixXd hessian_obj_sol(3, 3);
    hessian_obj_sol << 2, 1, 2, 1, 4, 0, 2, 0, 6;

    Eigen::MatrixXd hessian_constr_sol(3, 3);
    hessian_constr_sol << 2, 1, 2, 1, 4, 0, 2, 0, 6;

    testObjectiveAndEqualityAndInequalityHessians(hessian_obj_sol, hessian_constr_sol, hessian_constr_sol, optim.getEqualityDimension(),
                                                  optim.getInequalityDimension(), 7, 7, 7, 1e-4);
}

TEST_F(TestHyperGraphOptimizationProblemVertexBased, hessian_mixed)
{
    // Vertices
    vertices->addVertices({v1, v2});

    // Edges
    BaseEdge::Ptr edge3o(new Edge3T(edge3_fun, false, *v1, *v2));
    BaseEdge::Ptr edge3e(new Edge3T(edge3_fun, false, *v1, *v2));
    BaseEdge::Ptr edge3i(new Edge3T(edge3_fun, false, *v1, *v2));

    auto e2_precompute = [this](const MixedEdge2T::VertexContainer& vertices) { mixed_edge2_precompute(vertices); };
    auto e2_obj        = [this](const MixedEdge2T::VertexContainer& vertices, Eigen::Ref<Eigen::VectorXd> values_obj) {
        mixed_edge2_obj(vertices, values_obj);
    };
    auto e2_eq = [this](const MixedEdge2T::VertexContainer& vertices, Eigen::Ref<Eigen::VectorXd> values_eq) { mixed_edge2_eq(vertices, values_eq); };
    auto e2_ineq = [this](const MixedEdge2T::VertexContainer& vertices, Eigen::Ref<Eigen::VectorXd> values_ineq) {
        mixed_edge2_ineq(vertices, values_ineq);
    };
    BaseMixedEdge::Ptr edge2(new MixedEdge2T(2, 2, 2, e2_precompute, e2_obj, e2_eq, e2_ineq, *v1, *v2));

    edges->addEdges({edge3o}, {edge3e}, {edge3i}, {edge2});

    Eigen::MatrixXd hessian_sol(3, 3);
    hessian_sol << 2, 1, 2, 1, 4, 0, 2, 0, 6;
    hessian_sol *= 1.5;  // the mixed edge results in the same hessian scaled by 0.5, hence 1.5*hessian_sol

    testObjectiveAndEqualityAndInequalityHessians(hessian_sol, hessian_sol, hessian_sol, optim.getEqualityDimension(), optim.getInequalityDimension(),
                                                  7, 7, 7, 1e-4);
}

TEST_F(TestHyperGraphOptimizationProblemVertexBased, hessian_fixed)
{
    // Vertices
    vertices->addVertices({v1, v2});

    v1->setFixed(0, true);
    v2->setFixed(1, true);

    // Edges
    BaseEdge::Ptr edge3o(new Edge3T(edge3_fun, false, *v1, *v2));
    BaseEdge::Ptr edge3e(new Edge3T(edge3_fun, false, *v1, *v2));
    BaseEdge::Ptr edge3i(new Edge3T(edge3_fun, false, *v1, *v2));
    edges->addEdges({edge3o}, {edge3e}, {edge3i}, {});

    Eigen::MatrixXd hessian_sol(1, 1);
    hessian_sol << 4;

    testObjectiveAndEqualityAndInequalityHessians(hessian_sol, hessian_sol, hessian_sol, optim.getEqualityDimension(), optim.getInequalityDimension(),
                                                  1, 1, 1, 1e-4);
}

TEST_F(TestHyperGraphOptimizationProblemVertexBased, hessian_fixed_mixed)
{
    // Vertices
    vertices->addVertices({v1, v2});

    v1->setFixed(0, true);
    v2->setFixed(1, true);

    // Edges
    BaseEdge::Ptr edge3o(new Edge3T(edge3_fun, false, *v1, *v2));
    BaseEdge::Ptr edge3e(new Edge3T(edge3_fun, false, *v1, *v2));
    BaseEdge::Ptr edge3i(new Edge3T(edge3_fun, false, *v1, *v2));

    auto e2_precompute = [this](const MixedEdge2T::VertexContainer& vertices) { mixed_edge2_precompute(vertices); };
    auto e2_obj        = [this](const MixedEdge2T::VertexContainer& vertices, Eigen::Ref<Eigen::VectorXd> values_obj) {
        mixed_edge2_obj(vertices, values_obj);
    };
    auto e2_eq = [this](const MixedEdge2T::VertexContainer& vertices, Eigen::Ref<Eigen::VectorXd> values_eq) { mixed_edge2_eq(vertices, values_eq); };
    auto e2_ineq = [this](const MixedEdge2T::VertexContainer& vertices, Eigen::Ref<Eigen::VectorXd> values_ineq) {
        mixed_edge2_ineq(vertices, values_ineq);
    };
    BaseMixedEdge::Ptr edge2(new MixedEdge2T(2, 2, 2, e2_precompute, e2_obj, e2_eq, e2_ineq, *v1, *v2));

    edges->addEdges({edge3o}, {edge3e}, {edge3i}, {edge2});

    Eigen::MatrixXd hessian_sol(1, 1);
    hessian_sol << 4;
    hessian_sol *= 1.5;  // the mixed edge results in the same hessian scaled by 0.5, hence 1.5*hessian_sol

    testObjectiveAndEqualityAndInequalityHessians(hessian_sol, hessian_sol, hessian_sol, optim.getEqualityDimension(), optim.getInequalityDimension(),
                                                  1, 1, 1, 1e-4);
}

TEST_F(TestHyperGraphOptimizationProblemVertexBased, hessian_modified_values)
{
    // Vertices
    vertices->addVertices({v1, v2});

    // Edges
    BaseEdge::Ptr edge3o(new Edge3T(edge3_fun, false, *v1, *v2));
    BaseEdge::Ptr edge3e(new Edge3T(edge3_fun, false, *v1, *v2));
    BaseEdge::Ptr edge3i(new Edge3T(edge3_fun, false, *v1, *v2));
    edges->addEdges({edge3o}, {edge3e}, {edge3i}, {});

    double multiplier_obj = 3;
    Eigen::Vector2d multipliers_eq(10, 20);
    Eigen::Vector2d multipliers_ineq(20, 30);

    // x^2 + 2*y^2 + 3*z^2
    // Eigen::MatrixXd hessian_sol1(3, 3);
    // hessian_sol1 << 2, 0, 0, 0, 4, 0, 0, 0, 6;

    // x*y + 2*x*z
    // Eigen::MatrixXd hessian_sol2(3, 3);
    // hessian_sol2 << 0, 1, 2, 1, 0, 0, 2, 0, 0;

    Eigen::MatrixXd hessian_sol_obj(3, 3);
    hessian_sol_obj << 6, 3, 6, 3, 12, 0, 6, 0, 18;

    Eigen::MatrixXd hessian_sol_eq(3, 3);
    hessian_sol_eq << 20, 20, 40, 20, 40, 0, 40, 0, 60;

    Eigen::MatrixXd hessian_sol_ineq(3, 3);
    hessian_sol_ineq << 40, 30, 60, 30, 80, 0, 60, 0, 120;

    testObjectiveAndEqualityAndInequalityHessians(hessian_sol_obj, hessian_sol_eq, hessian_sol_ineq, optim.getEqualityDimension(),
                                                  optim.getInequalityDimension(), 7, 7, 7, 1e-4, multiplier_obj, multipliers_eq.data(),
                                                  multipliers_ineq.data());
}

TEST_F(TestHyperGraphOptimizationProblemVertexBased, hessian_modified_values_mixed)
{
    // Vertices
    vertices->addVertices({v1, v2});

    // Edges
    BaseEdge::Ptr edge3o(new Edge3T(edge3_fun, false, *v1, *v2));
    BaseEdge::Ptr edge3e(new Edge3T(edge3_fun, false, *v1, *v2));
    BaseEdge::Ptr edge3i(new Edge3T(edge3_fun, false, *v1, *v2));

    auto e2_precompute = [this](const MixedEdge2T::VertexContainer& vertices) { mixed_edge2_precompute(vertices); };
    auto e2_obj        = [this](const MixedEdge2T::VertexContainer& vertices, Eigen::Ref<Eigen::VectorXd> values_obj) {
        mixed_edge2_obj(vertices, values_obj);
    };
    auto e2_eq = [this](const MixedEdge2T::VertexContainer& vertices, Eigen::Ref<Eigen::VectorXd> values_eq) { mixed_edge2_eq(vertices, values_eq); };
    auto e2_ineq = [this](const MixedEdge2T::VertexContainer& vertices, Eigen::Ref<Eigen::VectorXd> values_ineq) {
        mixed_edge2_ineq(vertices, values_ineq);
    };
    BaseMixedEdge::Ptr edge2(new MixedEdge2T(2, 2, 2, e2_precompute, e2_obj, e2_eq, e2_ineq, *v1, *v2));

    edges->addEdges({edge3o}, {edge3e}, {edge3i}, {edge2});

    double multiplier_obj = 3;
    Eigen::Vector4d multipliers_eq(10, 20, 30, 40);
    Eigen::Vector4d multipliers_ineq(20, 30, 40, 50);

    // x^2 + 2*y^2 + 3*z^2
    // Eigen::MatrixXd hessian_sol1(3, 3);
    // hessian_sol1 << 2, 0, 0, 0, 4, 0, 0, 0, 6;

    // x*y + 2*x*z
    // Eigen::MatrixXd hessian_sol2(3, 3);
    // hessian_sol2 << 0, 1, 2, 1, 0, 0, 2, 0, 0;

    // Eigen::MatrixXd hessian_sol = multipliers[0] * hessian_sol1 + multipliers[1] * hessian_sol2;
    // hessian_sol *= 1.5;  // the mixed edge results in the same hessian scaled by 0.5, hence 1.5*hessian_sol

    Eigen::MatrixXd hessian_sol_obj(3, 3);
    hessian_sol_obj << 9, 4.5, 9, 4.5, 18, 0, 9, 0, 27;

    Eigen::MatrixXd hessian_sol_eq(3, 3);
    hessian_sol_eq << 50, 40, 80, 40, 100, 0, 80, 0, 150;

    Eigen::MatrixXd hessian_sol_ineq(3, 3);
    hessian_sol_ineq << 80, 55, 110, 55, 160, 0, 110, 0, 240;

    testObjectiveAndEqualityAndInequalityHessians(hessian_sol_obj, hessian_sol_eq, hessian_sol_ineq, optim.getEqualityDimension(),
                                                  optim.getInequalityDimension(), 7, 7, 7, 1e-4, multiplier_obj, multipliers_eq.data(),
                                                  multipliers_ineq.data());
}

TEST_F(TestHyperGraphOptimizationProblemVertexBased, hessian_lsq)
{
    // Vertices
    vertices->addVertices({v1, v2});

    // Edges
    BaseEdge::Ptr edge2o(new Edge2T(edge2_fun, true, *v2));  // least squares edge -> error is squared

    BaseEdge::Ptr edge3o(new Edge3T(edge3_fun, false, *v1, *v2));
    BaseEdge::Ptr edge3e(new Edge3T(edge3_fun, false, *v1, *v2));
    BaseEdge::Ptr edge3i(new Edge3T(edge3_fun, false, *v1, *v2));
    edges->addEdges({edge2o, edge3o}, {edge3e}, {edge3i}, {});

    // 16x^2 + (5x + 6y)^2
    // J_x = 32*x + 10*(5x + 6y)
    // J_y = 12*(5x + 6y)
    // Eigen::MatrixXd hessian_sol1(3, 3);
    // hessian_sol1 << 0, 0, 0, 0, 32 + 50, 60, 0, 60, 12 * 6;

    // Eigen::MatrixXd hessian_sol2(3, 3);
    // hessian_sol2 << 2, 1, 2, 1, 4, 0, 2, 0, 6;

    // Eigen::MatrixXd hessian_sol = hessian_sol1 + hessian_sol2;

    Eigen::MatrixXd hessian_sol_obj(3, 3);
    hessian_sol_obj << 2, 1, 2, 1, 86, 60, 2, 60, 78;

    Eigen::MatrixXd hessian_sol_eq(3, 3);
    hessian_sol_eq << 2, 1, 2, 1, 4, 0, 2, 0, 6;

    Eigen::MatrixXd hessian_sol_ineq(3, 3);
    hessian_sol_ineq << 2, 1, 2, 1, 4, 0, 2, 0, 6;

    testObjectiveAndEqualityAndInequalityHessians(hessian_sol_obj, hessian_sol_eq, hessian_sol_ineq, optim.getEqualityDimension(),
                                                  optim.getInequalityDimension(), 9, 7, 7, 1e-4);
}

TEST_F(TestHyperGraphOptimizationProblemVertexBased, hessian_lsq_mixed)
{
    // Vertices
    vertices->addVertices({v1, v2});

    // Edges
    BaseEdge::Ptr edge1o(new Edge2T(edge2_fun, true, *v2));  // least squares edge -> error is squared

    BaseEdge::Ptr edge2o(new Edge3T(edge3_fun, false, *v1, *v2));
    BaseEdge::Ptr edge2e(new Edge3T(edge3_fun, false, *v1, *v2));
    BaseEdge::Ptr edge2i(new Edge3T(edge3_fun, false, *v1, *v2));

    auto e3_precompute = [this](const MixedEdge1T::VertexContainer& vertices) { mixed_edge1_precompute(vertices); };
    auto e3_obj        = [this](const MixedEdge1T::VertexContainer& vertices, Eigen::Ref<Eigen::VectorXd> values_obj) {
        mixed_edge1_obj(vertices, values_obj);
    };
    auto e3_eq = [this](const MixedEdge1T::VertexContainer& vertices, Eigen::Ref<Eigen::VectorXd> values_eq) { mixed_edge1_eq(vertices, values_eq); };
    auto e3_ineq = [this](const MixedEdge1T::VertexContainer& vertices, Eigen::Ref<Eigen::VectorXd> values_ineq) {
        mixed_edge1_ineq(vertices, values_ineq);
    };
    MixedEdge1T::Ptr edge3(new MixedEdge1T(1, 2, 3, e3_precompute, e3_obj, e3_eq, e3_ineq, *v1, *v2));
    edge3->setObjectiveLsqForm(true);

    edges->addEdges({edge1o, edge2o}, {edge2e}, {edge2i}, {edge3});

    // 16x^2 + (5x + 6y)^2
    // J_x = 32*x + 10*(5x + 6y)
    // J_y = 12*(5x + 6y)
    // Eigen::MatrixXd hessian_sol1(3, 3);
    // hessian_sol1 << 0, 0, 0, 0, 32 + 50, 60, 0, 60, 12 * 6;

    // Eigen::MatrixXd hessian_sol2(3, 3);
    // hessian_sol2 << 2, 1, 2, 1, 4, 0, 2, 0, 6;

    // obj: 4*(x + 2y + 3z)^2
    // Jx = 8*(x + 2y + 3z)
    // Jy = 16*(x + 2y + 3z)
    // Jz = 24*(x + 2y + 3z)
    // Eigen::MatrixXd hessian_sol3_obj(3, 3);
    // hessian_sol3_obj << 8, 16, 24, 16, 32, 48, 24, 48, 72;

    // eq: 4*(16x^2 + (5y + 6z)^2)
    // Jx = 128x
    // Jy = 40*(5y + 6z)
    // Jz = 48*(5y + 6z)
    // Eigen::MatrixXd hessian_sol3_eq(3, 3);
    // hessian_sol3_eq << 128, 0, 0, 0, 200, 240, 0, 240, 288;

    // ineq: 4*( (x-5)^2 + (y+3)^2 + z^2 )
    // Jx = 8*(x-5)
    // Jy = 8*(y+3)
    // Jz = 8*z
    // Eigen::MatrixXd hessian_sol3_ineq(3, 3);
    // hessian_sol3_ineq << 8, 0, 0, 0, 8, 0, 0, 0, 8;

    // Eigen::MatrixXd hessian_sol_obj  = hessian_sol1 + hessian_sol2 + hessian_sol3_obj;
    // Eigen::MatrixXd hessian_sol_eq   = hessian_sol1 + hessian_sol2 + hessian_sol3_eq;
    // Eigen::MatrixXd hessian_sol_ineq = hessian_sol1 + hessian_sol2 + hessian_sol3_ineq;

    Eigen::MatrixXd hessian_sol_obj(3, 3);
    hessian_sol_obj << 10, 17, 26, 17, 118, 108, 26, 108, 150;

    Eigen::MatrixXd hessian_sol_eq(3, 3);
    hessian_sol_eq << 2, 1, 2, 1, 4, 0, 2, 0, 6;

    Eigen::MatrixXd hessian_sol_ineq(3, 3);
    hessian_sol_ineq << 2, 1, 2, 1, 4, 0, 2, 0, 6;

    testObjectiveAndEqualityAndInequalityHessians(hessian_sol_obj, hessian_sol_eq, hessian_sol_ineq, optim.getEqualityDimension(),
                                                  optim.getInequalityDimension(), 9, 7, 7, 1e-4);
}

TEST_F(TestHyperGraphOptimizationProblemVertexBased, hessian_lsq_fixed)
{
    // Vertices
    vertices->addVertices({v1, v2});

    v1->setFixed(0, true);
    v2->setFixed(1, true);

    // Edges
    BaseEdge::Ptr edge2o(new Edge2T(edge2_fun, true, *v2));  // least squares edge -> error is squared

    BaseEdge::Ptr edge3o(new Edge3T(edge3_fun, false, *v1, *v2));
    BaseEdge::Ptr edge3e(new Edge3T(edge3_fun, false, *v1, *v2));
    BaseEdge::Ptr edge3i(new Edge3T(edge3_fun, false, *v1, *v2));

    edges->addEdges({edge2o, edge3o}, {edge3e}, {edge3i}, {});

    // 16x^2 + (5x + 6y)^2
    // J_x = 32*x + 10*(5x + 6y)
    // J_y = 12*(5x + 6y)
    // Eigen::MatrixXd hessian_sol1(3, 3);
    // hessian_sol1 << 0, 0, 0, 0, 32 + 50, 60, 0, 60, 12 * 6;

    // Eigen::MatrixXd hessian_sol2(3, 3);
    // hessian_sol2 << 2, 1, 2, 1, 4, 0, 2, 0, 6;

    // Eigen::MatrixXd hessian_sol = (hessian_sol1 + hessian_sol2).block(1, 1, 1, 1);

    Eigen::MatrixXd hessian_sol_obj(1, 1);
    hessian_sol_obj << 86;

    Eigen::MatrixXd hessian_sol_eq(1, 1);
    hessian_sol_eq << 4;

    Eigen::MatrixXd hessian_sol_ineq(1, 1);
    hessian_sol_ineq << 4;

    testObjectiveAndEqualityAndInequalityHessians(hessian_sol_obj, hessian_sol_eq, hessian_sol_ineq, optim.getEqualityDimension(),
                                                  optim.getInequalityDimension(), 1, 1, 1, 1e-4);
}

TEST_F(TestHyperGraphOptimizationProblemVertexBased, hessian_lsq_fixed_mixed)
{
    // Vertices
    vertices->addVertices({v1, v2});

    v1->setFixed(0, true);
    v2->setFixed(1, true);

    // Edges
    BaseEdge::Ptr edge1o(new Edge2T(edge2_fun, true, *v2));  // least squares edge -> error is squared

    BaseEdge::Ptr edge2o(new Edge3T(edge3_fun, false, *v1, *v2));
    BaseEdge::Ptr edge2e(new Edge3T(edge3_fun, false, *v1, *v2));
    BaseEdge::Ptr edge2i(new Edge3T(edge3_fun, false, *v1, *v2));

    auto e3_precompute = [this](const MixedEdge1T::VertexContainer& vertices) { mixed_edge1_precompute(vertices); };
    auto e3_obj        = [this](const MixedEdge1T::VertexContainer& vertices, Eigen::Ref<Eigen::VectorXd> values_obj) {
        mixed_edge1_obj(vertices, values_obj);
    };
    auto e3_eq = [this](const MixedEdge1T::VertexContainer& vertices, Eigen::Ref<Eigen::VectorXd> values_eq) { mixed_edge1_eq(vertices, values_eq); };
    auto e3_ineq = [this](const MixedEdge1T::VertexContainer& vertices, Eigen::Ref<Eigen::VectorXd> values_ineq) {
        mixed_edge1_ineq(vertices, values_ineq);
    };
    MixedEdge1T::Ptr edge3(new MixedEdge1T(1, 2, 3, e3_precompute, e3_obj, e3_eq, e3_ineq, *v1, *v2));
    edge3->setObjectiveLsqForm(true);

    edges->addEdges({edge1o, edge2o}, {edge2e}, {edge2i}, {edge3});

    // 16x^2 + (5x + 6y)^2
    // J_x = 32*x + 10*(5x + 6y)
    // J_y = 12*(5x + 6y)
    // Eigen::MatrixXd hessian_sol1(3, 3);
    // hessian_sol1 << 0, 0, 0, 0, 32 + 50, 60, 0, 60, 12 * 6;

    // Eigen::MatrixXd hessian_sol2(3, 3);
    // hessian_sol2 << 2, 1, 2, 1, 4, 0, 2, 0, 6;

    // obj: 4*(x + 2y + 3z)^2
    // Jx = 8*(x + 2y + 3z)
    // Jy = 16*(x + 2y + 3z)
    // Jz = 24*(x + 2y + 3z)
    // Eigen::MatrixXd hessian_sol3_obj(3, 3);
    // hessian_sol3_obj << 8, 16, 24, 16, 32, 48, 24, 48, 72;

    // eq: 4*(16x^2 + (5y + 6z)^2)
    // Jx = 128x
    // Jy = 40*(5y + 6z)
    // Jz = 48*(5y + 6z)
    // Eigen::MatrixXd hessian_sol3_eq(3, 3);
    // hessian_sol3_eq << 128, 0, 0, 0, 200, 240, 0, 240, 288;

    // ineq: 4*( (x-5)^2 + (y+3)^2 + z^2 )
    // Jx = 8*(x-5)
    // Jy = 8*(y+3)
    // Jz = 8*z
    // Eigen::MatrixXd hessian_sol3_ineq(3, 3);
    // hessian_sol3_ineq << 8, 0, 0, 0, 8, 0, 0, 0, 8;

    // Eigen::MatrixXd hessian_sol_obj  = (hessian_sol1 + hessian_sol2 + hessian_sol3_obj).block(1, 1, 1, 1);
    // Eigen::MatrixXd hessian_sol_eq   = (hessian_sol1 + hessian_sol2 + hessian_sol3_eq).block(1, 1, 1, 1);
    // Eigen::MatrixXd hessian_sol_ineq = (hessian_sol1 + hessian_sol2 + hessian_sol3_ineq).block(1, 1, 1, 1);

    Eigen::MatrixXd hessian_sol_obj(1, 1);
    hessian_sol_obj << 118;

    Eigen::MatrixXd hessian_sol_eq(1, 1);
    hessian_sol_eq << 4;

    Eigen::MatrixXd hessian_sol_ineq(1, 1);
    hessian_sol_ineq << 4;

    testObjectiveAndEqualityAndInequalityHessians(hessian_sol_obj, hessian_sol_eq, hessian_sol_ineq, optim.getEqualityDimension(),
                                                  optim.getInequalityDimension(), 1, 1, 1, 1e-4);
}

TEST_F(TestHyperGraphOptimizationProblemVertexBased, hessian_lsq_modified_values)
{
    // Vertices
    vertices->addVertices({v1, v2});

    // Edges
    BaseEdge::Ptr edge2o(new Edge2T(edge2_fun, true, *v2));  // least squares edge -> error is squared

    BaseEdge::Ptr edge3o(new Edge3T(edge3_fun, false, *v1, *v2));
    BaseEdge::Ptr edge3e(new Edge3T(edge3_fun, false, *v1, *v2));
    BaseEdge::Ptr edge3i(new Edge3T(edge3_fun, false, *v1, *v2));

    edges->addEdges({edge2o, edge3o}, {edge3e}, {edge3i}, {});

    double multiplier_obj = 3;
    std::array<double, 4> multipliers_eq   = {10, 20, 20, 30};
    std::array<double, 4> multipliers_ineq = {20, 30, 40, 50};

    // 16x^2 + (5x + 6y)^2
    // J_x = 32*x + 10*(5x + 6y)
    // J_y = 12*(5x + 6y)
    // Eigen::MatrixXd hessian_sol1(3, 3);
    // hessian_sol1 << 0, 0, 0, 0, 32 + 50, 60, 0, 60, 12 * 6;

    // Eigen::MatrixXd hessian_sol2(3, 3);
    // hessian_sol2 << 2, 1, 2, 1, 4, 0, 2, 0, 6;

    // Eigen::MatrixXd hessian_sol = 10 * hessian_sol1 + 20 * hessian_sol2;

    Eigen::MatrixXd hessian_sol_obj(3, 3);
    hessian_sol_obj << 6, 3, 6, 3, 258, 180, 6, 180, 234;

    Eigen::MatrixXd hessian_sol_eq(3, 3);
    hessian_sol_eq << 20, 20, 40, 20, 40, 0, 40, 0, 60;

    Eigen::MatrixXd hessian_sol_ineq(3, 3);
    hessian_sol_ineq << 40, 30, 60, 30, 80, 0, 60, 0, 120;

    testObjectiveAndEqualityAndInequalityHessians(hessian_sol_obj, hessian_sol_eq, hessian_sol_ineq, optim.getEqualityDimension(),
                                                  optim.getInequalityDimension(), 9, 7, 7, 1e-4, multiplier_obj, multipliers_eq.data(),
                                                  multipliers_ineq.data());
}

TEST_F(TestHyperGraphOptimizationProblemVertexBased, hessian_lsq_modified_values_mixed)
{
    // Vertices
    vertices->addVertices({v1, v2});

    // Edges
    BaseEdge::Ptr edge1o(new Edge2T(edge2_fun, true, *v2));  // least squares edge -> error is squared

    BaseEdge::Ptr edge2o(new Edge3T(edge3_fun, false, *v1, *v2));
    BaseEdge::Ptr edge2e(new Edge3T(edge3_fun, false, *v1, *v2));
    BaseEdge::Ptr edge2i(new Edge3T(edge3_fun, false, *v1, *v2));

    auto e3_precompute = [this](const MixedEdge1T::VertexContainer& vertices) { mixed_edge1_precompute(vertices); };
    auto e3_obj        = [this](const MixedEdge1T::VertexContainer& vertices, Eigen::Ref<Eigen::VectorXd> values_obj) {
        mixed_edge1_obj(vertices, values_obj);
    };
    auto e3_eq = [this](const MixedEdge1T::VertexContainer& vertices, Eigen::Ref<Eigen::VectorXd> values_eq) { mixed_edge1_eq(vertices, values_eq); };
    auto e3_ineq = [this](const MixedEdge1T::VertexContainer& vertices, Eigen::Ref<Eigen::VectorXd> values_ineq) {
        mixed_edge1_ineq(vertices, values_ineq);
    };
    MixedEdge1T::Ptr edge3(new MixedEdge1T(1, 2, 3, e3_precompute, e3_obj, e3_eq, e3_ineq, *v1, *v2));
    edge3->setObjectiveLsqForm(true);

    edges->addEdges({edge1o, edge2o}, {edge2e}, {edge2i}, {edge3});

    // std::array<double, 1 + 2 + 1> multipliers_obj  = {10, 20, 20, 30};
    // std::array<double, 1 + 2 + 1> multipliers_eq   = {10, 20, 20, 40};
    // std::array<double, 1 + 2 + 1> multipliers_ineq = {10, 20, 20, 50};
    double multiplier_obj = 3;
    std::array<double, 6> multipliers_eq   = {10, 20, 20, 30, 40, 50};
    std::array<double, 7> multipliers_ineq = {1, 20, 30, 40, 50, 60, 70};

    // 16x^2 + (5x + 6y)^2
    // J_x = 32*x + 10*(5x + 6y)
    // J_y = 12*(5x + 6y)
    // Eigen::MatrixXd hessian_sol1(3, 3);
    // hessian_sol1 << 0, 0, 0, 0, 32 + 50, 60, 0, 60, 12 * 6;

    // Eigen::MatrixXd hessian_sol2(3, 3);
    // hessian_sol2 << 2, 1, 2, 1, 4, 0, 2, 0, 6;

    // obj: 4*(x + 2y + 3z)^2
    // Jx = 8*(x + 2y + 3z)
    // Jy = 16*(x + 2y + 3z)
    // Jz = 24*(x + 2y + 3z)
    // Eigen::MatrixXd hessian_sol3_obj(3, 3);
    // hessian_sol3_obj << 8, 16, 24, 16, 32, 48, 24, 48, 72;

    // eq: 4*(16x^2 + (5y + 6z)^2)
    // Jx = 128x
    // Jy = 40*(5y + 6z)
    // Jz = 48*(5y + 6z)
    // Eigen::MatrixXd hessian_sol3_eq(3, 3);
    // hessian_sol3_eq << 128, 0, 0, 0, 200, 240, 0, 240, 288;

    // ineq: 4*( (x-5)^2 + (y+3)^2 + z^2 )
    // Jx = 8*(x-5)
    // Jy = 8*(y+3)
    // Jz = 8*z
    // Eigen::MatrixXd hessian_sol3_ineq(3, 3);
    // hessian_sol3_ineq << 8, 0, 0, 0, 8, 0, 0, 0, 8;

    // Eigen::MatrixXd hessian_sol_obj  = 10 * hessian_sol1 + 20 * hessian_sol2 + 30 * hessian_sol3_obj;
    // Eigen::MatrixXd hessian_sol_eq   = 10 * hessian_sol1 + 20 * hessian_sol2 + 40 * hessian_sol3_eq;
    // Eigen::MatrixXd hessian_sol_ineq = 10 * hessian_sol1 + 20 * hessian_sol2 + 50 * hessian_sol3_ineq;

    Eigen::MatrixXd hessian_sol_obj(3, 3);
    hessian_sol_obj << 30, 51, 78, 51, 354, 324, 78, 324, 450;

    Eigen::MatrixXd hessian_sol_eq(3, 3);
    hessian_sol_eq << 20, 20, 40, 20, 40, 0, 40, 0, 60;

    Eigen::MatrixXd hessian_sol_ineq(3, 3);
    hessian_sol_ineq << 2, 20, 40, 20, 4, 0, 40, 0, 6;

    testObjectiveAndEqualityAndInequalityHessians(hessian_sol_obj, hessian_sol_eq, hessian_sol_ineq, optim.getEqualityDimension(),
                                                  optim.getInequalityDimension(), 9, 7, 7, 1e-3, multiplier_obj, multipliers_eq.data(),
                                                  multipliers_ineq.data());
}

TEST_F(TestHyperGraphOptimizationProblemVertexBased, hessian_lsq_modified_values_fixed)
{
    // Vertices
    vertices->addVertices({v1, v2});

    v1->setFixed(0, true);
    v2->setFixed(1, true);

    // Edges
    BaseEdge::Ptr edge2o(new Edge2T(edge2_fun, true, *v2));  // least squares edge -> error is squared

    BaseEdge::Ptr edge3o(new Edge3T(edge3_fun, false, *v1, *v2));
    BaseEdge::Ptr edge3e(new Edge3T(edge3_fun, false, *v1, *v2));
    BaseEdge::Ptr edge3i(new Edge3T(edge3_fun, false, *v1, *v2));

    edges->addEdges({edge2o, edge3o}, {edge3e}, {edge3i}, {});

    double multiplier_obj = 3;
    std::array<double, 4> multipliers_eq   = {10, 20, 20, 30};
    std::array<double, 4> multipliers_ineq = {20, 30, 40, 50};

    // std::array<double, 3> multipliers = {10, 20, 20};  // first edge is LSQ -> dim=1

    // 16x^2 + (5x + 6y)^2
    // J_x = 32*x + 10*(5x + 6y)
    // J_y = 12*(5x + 6y)
    // Eigen::MatrixXd hessian_sol1(3, 3);
    // hessian_sol1 << 0, 0, 0, 0, 32 + 50, 60, 0, 60, 12 * 6;

    // Eigen::MatrixXd hessian_sol2(3, 3);
    // hessian_sol2 << 2, 1, 2, 1, 4, 0, 2, 0, 6;

    // Eigen::MatrixXd hessian_sol = (10 * hessian_sol1 + 20 * hessian_sol2).block(1, 1, 1, 1);

    Eigen::MatrixXd hessian_sol_obj(1, 1);
    hessian_sol_obj << 258;

    Eigen::MatrixXd hessian_sol_eq(1, 1);
    hessian_sol_eq << 40;

    Eigen::MatrixXd hessian_sol_ineq(1, 1);
    hessian_sol_ineq << 80;

    testObjectiveAndEqualityAndInequalityHessians(hessian_sol_obj, hessian_sol_eq, hessian_sol_ineq, optim.getEqualityDimension(),
                                                  optim.getInequalityDimension(), 1, 1, 1, 1e-4, multiplier_obj, multipliers_eq.data(),
                                                  multipliers_ineq.data());
}

TEST_F(TestHyperGraphOptimizationProblemVertexBased, hessian_lsq_modified_values_fixed_mixed)
{
    // Vertices
    vertices->addVertices({v1, v2});

    v1->setFixed(0, true);
    v2->setFixed(1, true);

    // Edges
    BaseEdge::Ptr edge1o(new Edge2T(edge2_fun, true, *v2));  // least squares edge -> error is squared

    BaseEdge::Ptr edge2o(new Edge3T(edge3_fun, false, *v1, *v2));
    BaseEdge::Ptr edge2e(new Edge3T(edge3_fun, false, *v1, *v2));
    BaseEdge::Ptr edge2i(new Edge3T(edge3_fun, false, *v1, *v2));

    auto e3_precompute = [this](const MixedEdge1T::VertexContainer& vertices) { mixed_edge1_precompute(vertices); };
    auto e3_obj        = [this](const MixedEdge1T::VertexContainer& vertices, Eigen::Ref<Eigen::VectorXd> values_obj) {
        mixed_edge1_obj(vertices, values_obj);
    };
    auto e3_eq = [this](const MixedEdge1T::VertexContainer& vertices, Eigen::Ref<Eigen::VectorXd> values_eq) { mixed_edge1_eq(vertices, values_eq); };
    auto e3_ineq = [this](const MixedEdge1T::VertexContainer& vertices, Eigen::Ref<Eigen::VectorXd> values_ineq) {
        mixed_edge1_ineq(vertices, values_ineq);
    };
    MixedEdge1T::Ptr edge3(new MixedEdge1T(1, 2, 3, e3_precompute, e3_obj, e3_eq, e3_ineq, *v1, *v2));
    edge3->setObjectiveLsqForm(true);

    edges->addEdges({edge1o, edge2o}, {edge2e}, {edge2i}, {edge3});

    double multiplier_obj = 3;
    std::array<double, 4> multipliers_eq   = {10, 20, 20, 30};
    std::array<double, 4> multipliers_ineq = {20, 30, 40, 50};

    // std::array<double, 1 + 2 + 1> multipliers_obj  = {10, 20, 20, 30};
    // std::array<double, 1 + 2 + 1> multipliers_eq   = {10, 20, 20, 40};
    // std::array<double, 1 + 2 + 1> multipliers_ineq = {10, 20, 20, 50};

    // 16x^2 + (5x + 6y)^2
    // J_x = 32*x + 10*(5x + 6y)
    // J_y = 12*(5x + 6y)
    // Eigen::MatrixXd hessian_sol1(3, 3);
    // hessian_sol1 << 0, 0, 0, 0, 32 + 50, 60, 0, 60, 12 * 6;

    // Eigen::MatrixXd hessian_sol2(3, 3);
    // hessian_sol2 << 2, 1, 2, 1, 4, 0, 2, 0, 6;

    // obj: 4*(x + 2y + 3z)^2
    // Jx = 8*(x + 2y + 3z)
    // Jy = 16*(x + 2y + 3z)
    // Jz = 24*(x + 2y + 3z)
    // Eigen::MatrixXd hessian_sol3_obj(3, 3);
    // hessian_sol3_obj << 8, 16, 24, 16, 32, 48, 24, 48, 72;

    // eq: 4*(16x^2 + (5y + 6z)^2)
    // Jx = 128x
    // Jy = 40*(5y + 6z)
    // Jz = 48*(5y + 6z)
    // Eigen::MatrixXd hessian_sol3_eq(3, 3);
    // hessian_sol3_eq << 128, 0, 0, 0, 200, 240, 0, 240, 288;

    // ineq: 4*( (x-5)^2 + (y+3)^2 + z^2 )
    // Jx = 8*(x-5)
    // Jy = 8*(y+3)
    // Jz = 8*z
    // Eigen::MatrixXd hessian_sol3_ineq(3, 3);
    // hessian_sol3_ineq << 8, 0, 0, 0, 8, 0, 0, 0, 8;

    // Eigen::MatrixXd hessian_sol_obj  = (10 * hessian_sol1 + 20 * hessian_sol2 + 30 * hessian_sol3_obj).block(1, 1, 1, 1);
    // Eigen::MatrixXd hessian_sol_eq   = (10 * hessian_sol1 + 20 * hessian_sol2 + 40 * hessian_sol3_eq).block(1, 1, 1, 1);
    // Eigen::MatrixXd hessian_sol_ineq = (10 * hessian_sol1 + 20 * hessian_sol2 + 50 * hessian_sol3_ineq).block(1, 1, 1, 1);

    Eigen::MatrixXd hessian_sol_obj(1, 1);
    hessian_sol_obj << 354;

    Eigen::MatrixXd hessian_sol_eq(1, 1);
    hessian_sol_eq << 40;

    Eigen::MatrixXd hessian_sol_ineq(1, 1);
    hessian_sol_ineq << 80;

    testObjectiveAndEqualityAndInequalityHessians(hessian_sol_obj, hessian_sol_eq, hessian_sol_ineq, optim.getEqualityDimension(),
                                                  optim.getInequalityDimension(), 1, 1, 1, 1e-3, multiplier_obj, multipliers_eq.data(),
                                                  multipliers_ineq.data());
}
