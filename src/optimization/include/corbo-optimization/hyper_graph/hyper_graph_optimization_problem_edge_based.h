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

#ifndef SRC_OPTIMIZATION_INCLUDE_CORBO_OPTIMIZATION_HYPER_GRAPH_HYPER_GRAPH_OPTIMIZATION_PROBLEM_EDGE_BASED_H_
#define SRC_OPTIMIZATION_INCLUDE_CORBO_OPTIMIZATION_HYPER_GRAPH_HYPER_GRAPH_OPTIMIZATION_PROBLEM_EDGE_BASED_H_

#include <corbo-optimization/hyper_graph/hyper_graph_optimization_problem_base.h>

#include <memory>

namespace corbo {

class HyperGraphOptimizationProblemEdgeBased : public BaseHyperGraphOptimizationProblem
{
 public:
    HyperGraphOptimizationProblemEdgeBased() = default;
    HyperGraphOptimizationProblemEdgeBased(OptimizationEdgeSet::Ptr edges, VertexSetInterface::Ptr vertices)
        : BaseHyperGraphOptimizationProblem(edges, vertices)
    {
    }

    BaseHyperGraphOptimizationProblem::Ptr getInstance() const override { return std::make_shared<HyperGraphOptimizationProblemEdgeBased>(); }

    void computeGradientObjective(Eigen::Ref<Eigen::VectorXd> gradient) override;
    void computeGradientNonLsqObjective(Eigen::Ref<Eigen::VectorXd> gradient) override;

    void computeDenseJacobianLsqObjective(Eigen::Ref<Eigen::MatrixXd> jacobian, const double* multipliers = nullptr) override;
    int computeSparseJacobianLsqObjectiveNNZ() override;
    void computeSparseJacobianLsqObjectiveStructure(Eigen::Ref<Eigen::VectorXi> i_row, Eigen::Ref<Eigen::VectorXi> j_col) override;
    void computeSparseJacobianLsqObjectiveValues(Eigen::Ref<Eigen::VectorXd> values, const double* multipliers = nullptr) override;

    void computeDenseJacobianEqualities(Eigen::Ref<Eigen::MatrixXd> jacobian, const double* multipliers = nullptr) override;
    int computeSparseJacobianEqualitiesNNZ() override;
    void computeSparseJacobianEqualitiesStructure(Eigen::Ref<Eigen::VectorXi> i_row, Eigen::Ref<Eigen::VectorXi> j_col) override;
    void computeSparseJacobianEqualitiesValues(Eigen::Ref<Eigen::VectorXd> values, const double* multipliers = nullptr) override;

    void computeDenseJacobianInequalities(Eigen::Ref<Eigen::MatrixXd> jacobian, const double* multipliers = nullptr) override;
    int computeSparseJacobianInequalitiesNNZ() override;
    void computeSparseJacobianInequalitiesStructure(Eigen::Ref<Eigen::VectorXi> i_row, Eigen::Ref<Eigen::VectorXi> j_col) override;
    void computeSparseJacobianInequalitiesValues(Eigen::Ref<Eigen::VectorXd> values, const double* multipliers = nullptr) override;

    void computeDenseJacobianActiveInequalities(Eigen::Ref<Eigen::MatrixXd> jacobian, double weight = 1.0) override;
    void computeSparseJacobianActiveInequalitiesValues(Eigen::Ref<Eigen::VectorXd> values, double weight = 1.0) override;

    void computeDenseJacobians(Eigen::Ref<Eigen::VectorXd> gradient_non_lsq_obj, Eigen::Ref<Eigen::MatrixXd> jacobian_lsq_obj,
                               Eigen::Ref<Eigen::MatrixXd> jacobian_eq, Eigen::Ref<Eigen::MatrixXd> jacobian_ineq,
                               const double* multipliers_lsq_obj = nullptr, const double* multipliers_eq = nullptr,
                               const double* multipliers_ineq = nullptr, bool active_ineq = false, double active_ineq_weight = 1.0) override;

    void computeSparseJacobiansValues(Eigen::Ref<Eigen::VectorXd> values_obj, Eigen::Ref<Eigen::VectorXd> values_eq,
                                      Eigen::Ref<Eigen::VectorXd> values_ineq, const double* multipliers_obj = nullptr,
                                      const double* multipliers_eq = nullptr, const double* multipliers_ineq = nullptr, bool active_ineq = false,
                                      double active_ineq_weight = 1.0) override;

    int computeCombinedSparseJacobiansNNZ(bool objective_lsq = true, bool equality = true, bool inequality = true) override;
    void computeCombinedSparseJacobiansStructure(Eigen::Ref<Eigen::VectorXi> i_row, Eigen::Ref<Eigen::VectorXi> j_col, bool objective_lsq = true,
                                                 bool equality = true, bool inequality = true) override;
    void computeCombinedSparseJacobiansValues(Eigen::Ref<Eigen::VectorXd> values, bool objective_lsq = true, bool equality = true,
                                              bool inequality = true, const double* multipliers_obj = nullptr, const double* multipliers_eq = nullptr,
                                              const double* multipliers_ineq = nullptr) override;
    void computeCombinedSparseJacobian(Eigen::SparseMatrix<double>& jacobian, bool objective_lsq, bool equality, bool inequality,
                                       bool finite_combined_bounds, bool active_ineq = false, double weight_eq = 1.0, double weight_ineq = 1.0,
                                       double weight_bounds = 1.0, const Eigen::VectorXd* values = nullptr,
                                       const Eigen::VectorXi* col_nnz = nullptr) override;

    // useful for IPOPT (w/ hessian-approx)
    void computeGradientObjectiveAndCombinedSparseJacobiansValues(Eigen::Ref<Eigen::VectorXd> gradient, Eigen::Ref<Eigen::VectorXd> jac_values,
                                                                  bool equality = true, bool inequality = true,
                                                                  const double* multipliers_eq   = nullptr,
                                                                  const double* multipliers_ineq = nullptr) override;

    // warning: we are using H = 2* J^T*J as approximation in case of least-squares objective edges
    void computeDenseHessianObjective(Eigen::Ref<Eigen::MatrixXd> hessian, double multiplier = 1.0) override;

    int computeSparseHessianObjectiveNNZ(bool lower_part_only = false) override;
    void computeSparseHessianObjectiveStructure(Eigen::Ref<Eigen::VectorXi> i_row, Eigen::Ref<Eigen::VectorXi> j_col,
                                                bool lower_part_only = false) override;
    void computeSparseHessianObjectiveValues(Eigen::Ref<Eigen::VectorXd> values, double multiplier = 1.0, bool lower_part_only = false) override;

    void computeSparseHessianObjectiveLL(Eigen::SparseMatrix<double, Eigen::ColMajor, long long>& H, const Eigen::VectorXi* col_nnz = nullptr,
                                         bool upper_part_only = false) override;
    void computeSparseHessianObjectiveNNZperCol(Eigen::Ref<Eigen::VectorXi> col_nnz, bool upper_part_only = false) override;

    int computeSparseHessianEqualitiesNNZ(bool lower_part_only = false) override;
    void computeSparseHessianEqualitiesStructure(Eigen::Ref<Eigen::VectorXi> i_row, Eigen::Ref<Eigen::VectorXi> j_col,
                                                 bool lower_part_only = false) override;
    void computeSparseHessianEqualitiesValues(Eigen::Ref<Eigen::VectorXd> values, const double* multipliers = nullptr,
                                              bool lower_part_only = false) override;

    int computeSparseHessianInequalitiesNNZ(bool lower_part_only = false) override;
    void computeSparseHessianInequalitiesStructure(Eigen::Ref<Eigen::VectorXi> i_row, Eigen::Ref<Eigen::VectorXi> j_col,
                                                   bool lower_part_only = false) override;
    void computeSparseHessianInequalitiesValues(Eigen::Ref<Eigen::VectorXd> values, const double* multipliers = nullptr,
                                                bool lower_part_only = false) override;

    void computeSparseHessiansNNZ(int& nnz_obj, int& nnz_eq, int& nnz_ineq, bool lower_part_only = false) override;
    void computeSparseHessiansStructure(Eigen::Ref<Eigen::VectorXi> i_row_obj, Eigen::Ref<Eigen::VectorXi> j_col_obj,
                                        Eigen::Ref<Eigen::VectorXi> i_row_eq, Eigen::Ref<Eigen::VectorXi> j_col_eq,
                                        Eigen::Ref<Eigen::VectorXi> i_row_ineq, Eigen::Ref<Eigen::VectorXi> j_col_ineq,
                                        bool lower_part_only = false) override;
    void computeSparseHessiansValues(Eigen::Ref<Eigen::VectorXd> values_obj, Eigen::Ref<Eigen::VectorXd> values_eq,
                                     Eigen::Ref<Eigen::VectorXd> values_ineq, double multiplier_obj = 1.0, const double* multipliers_eq = nullptr,
                                     const double* multipliers_ineq = nullptr, bool lower_part_only = false) override;

    void computeSparseHessianLagrangian(Eigen::SparseMatrix<double, Eigen::ColMajor, long long>& H, const double* multipliers_eq,
                                        const double* multipliers_ineq, const Eigen::VectorXi* col_nnz = nullptr,
                                        bool upper_part_only = false) override;

    void computeSparseHessianLagrangianNNZperCol(Eigen::Ref<Eigen::VectorXi> col_nnz, bool upper_part_only = false) override;

    void computeSparseJacobianTwoSideBoundedLinearForm(Eigen::SparseMatrix<double, Eigen::ColMajor, long long>& A, bool include_finite_bounds,
                                                       const Eigen::VectorXi* col_nnz = nullptr) override;
    void computeSparseJacobianTwoSideBoundedLinearFormNNZPerColumn(Eigen::Ref<Eigen::VectorXi> col_nnz, bool include_finite_bounds) override;

    int computeSparseJacobianTwoSideBoundedLinearFormNNZ(bool include_finite_bounds) override;
    void computeSparseJacobianTwoSideBoundedLinearFormStructure(Eigen::Ref<Eigen::VectorXi> i_row, Eigen::Ref<Eigen::VectorXi> j_col,
                                                                bool include_finite_bounds) override;
    void computeSparseJacobianTwoSideBoundedLinearFormValues(Eigen::Ref<Eigen::VectorXd> values, bool include_finite_bounds) override;

    void computeSparseJacobianTwoSideBoundedLinearFormAndHessianLagrangian(Eigen::SparseMatrix<double, Eigen::ColMajor, long long>& H,
                                                                           const double* multipliers_eq, const double* multipliers_ineq,
                                                                           Eigen::SparseMatrix<double, Eigen::ColMajor, long long>& A,
                                                                           bool include_finite_bounds, const Eigen::VectorXi* col_nnz_H,
                                                                           const Eigen::VectorXi* col_nnz_A, bool upper_part_only_H) override;

    void computeSparseJacobianTwoSideBoundedLinearFormAndHessianObjective(Eigen::SparseMatrix<double, Eigen::ColMajor, long long>& H,
                                                                          Eigen::SparseMatrix<double, Eigen::ColMajor, long long>& A,
                                                                          bool include_finite_bounds, const Eigen::VectorXi* col_nnz_H,
                                                                          const Eigen::VectorXi* col_nnz_A, bool upper_part_only_H) override;

 protected:
    // bool _hessian_lsq_obj_approx = false; // at the moment every lsq objective edge is approximated by H = 2*J^T*J
};

FACTORY_REGISTER_HYPER_GRAPH_OPTIMIZATION_PROBLEM(HyperGraphOptimizationProblemEdgeBased)

}  // namespace corbo

#endif  // SRC_OPTIMIZATION_INCLUDE_CORBO_OPTIMIZATION_HYPER_GRAPH_HYPER_GRAPH_OPTIMIZATION_PROBLEM_EDGE_BASED_H_
