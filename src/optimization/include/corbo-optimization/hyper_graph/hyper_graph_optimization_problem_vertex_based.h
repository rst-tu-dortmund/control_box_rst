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

#ifndef SRC_OPTIMIZATION_INCLUDE_CORBO_OPTIMIZATION_HYPER_GRAPH_HYPER_GRAPH_OPTIMIZATION_PROBLEM_VERTEX_BASED_H_
#define SRC_OPTIMIZATION_INCLUDE_CORBO_OPTIMIZATION_HYPER_GRAPH_HYPER_GRAPH_OPTIMIZATION_PROBLEM_VERTEX_BASED_H_

#include <corbo-optimization/hyper_graph/hyper_graph_optimization_problem_base.h>

#include <memory>

namespace corbo {

class HyperGraphOptimizationProblemVertexBased : public BaseHyperGraphOptimizationProblem
{
 public:
    HyperGraphOptimizationProblemVertexBased() = default;
    HyperGraphOptimizationProblemVertexBased(OptimizationEdgeSet::Ptr edges, VertexSetInterface::Ptr vertices)
        : BaseHyperGraphOptimizationProblem(edges, vertices)
    {
    }

    BaseHyperGraphOptimizationProblem::Ptr getInstance() const override { return std::make_shared<HyperGraphOptimizationProblemVertexBased>(); }

    void precomputeEdgeQuantities() override;

    void computeGradientObjective(Eigen::Ref<Eigen::VectorXd> gradient) override;
    void computeGradientNonLsqObjective(Eigen::Ref<Eigen::VectorXd> gradient) override;
    void computeDenseJacobianLsqObjective(Eigen::Ref<Eigen::MatrixXd> jacobian, const double* multipliers = nullptr) override;

    int computeSparseJacobianLsqObjectiveNNZ() override;
    void computeSparseJacobianLsqObjectiveStructure(Eigen::Ref<Eigen::VectorXi> i_row, Eigen::Ref<Eigen::VectorXi> j_col) override;
    void computeSparseJacobianLsqObjectiveValues(Eigen::Ref<Eigen::VectorXd> values, const double* multipliers) override;

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

 protected:
    void precomputeConnectedMixedEdges(const VertexInterface* vertex, bool objective = true, bool equality = true, bool inequality = true);

    void computeObjectiveValuesCached(const VertexInterface* vertex, bool include_lsq_edges, bool include_nonmixed = true, bool include_mixed = true);
    void finalizeObjectiveGradient(const VertexInterface* vertex, double& gradient_coeff, bool include_lsq_edges, bool include_nonmixed = true,
                                   bool include_mixed = true, bool precompute_mixed = true);

    void computeLsqObjectiveValuesCached(const VertexInterface* vertex, bool include_nonmixed = true, bool include_mixed = true,
                                         bool precompute_mixed = true);
    void finalizeLsqObjectiveJacobian(const VertexInterface* vertex, int vtx_idx, Eigen::Ref<Eigen::MatrixXd>& jacobian,
                                      const double* multipliers = nullptr, bool include_nonmixed = true, bool include_mixed = true);
    void computeLsqObjectiveJacobianStructureForVertex(const VertexInterface* vertex, int vtx_idx, Eigen::Ref<Eigen::VectorXi> i_row,
                                                       Eigen::Ref<Eigen::VectorXi> j_col, int& nnz_idx, int row_offset = 0);
    void finalizeLsqObjectiveJacobianSparseValues(const VertexInterface* vertex, int& nnz_idx, Eigen::Ref<Eigen::VectorXd>& values,
                                                  const double* multipliers = nullptr, bool precompute_mixed = true);

    void computeEqualitiesValuesCached(const VertexInterface* vertex, bool include_nonmixed = true, bool include_mixed = true,
                                       bool precompute_mixed = true);
    void finalizeEqualitiesJacobian(const VertexInterface* vertex, int vtx_idx, Eigen::Ref<Eigen::MatrixXd>& jacobian,
                                    const double* multipliers = nullptr, bool include_nonmixed = true, bool include_mixed = true,
                                    bool precompute_mixed = true);
    void computeEqualitiesJacobianStructureForVertex(const VertexInterface* vertex, int vtx_idx, Eigen::Ref<Eigen::VectorXi> i_row,
                                                     Eigen::Ref<Eigen::VectorXi> j_col, int& nnz_idx, int row_offset = 0);
    void finalizeEqualitiesJacobianSparseValues(const VertexInterface* vertex, int& nnz_idx, Eigen::Ref<Eigen::VectorXd>& values,
                                                const double* multipliers = nullptr, bool precompute_mixed = true);

    void computeInequalitiesValuesCached(const VertexInterface* vertex, bool include_nonmixed = true, bool include_mixed = true,
                                         bool precompute_mixed = true);
    void finalizeInequalitiesJacobian(const VertexInterface* vertex, int vtx_idx, Eigen::Ref<Eigen::MatrixXd>& jacobian,
                                      const double* multipliers = nullptr, bool include_nonmixed = true, bool include_mixed = true,
                                      bool precompute_mixed = true);
    void finalizeActiveInequalitiesJacobian(const VertexInterface* vertex, int vtx_idx, Eigen::Ref<Eigen::MatrixXd>& jacobian, double weight = 1.0,
                                            bool include_nonmixed = true, bool include_mixed = true);
    void computeInequalitiesJacobianStructureForVertex(const VertexInterface* vertex, int vtx_idx, Eigen::Ref<Eigen::VectorXi> i_row,
                                                       Eigen::Ref<Eigen::VectorXi> j_col, int& nnz_idx, int row_offset = 0);
    void finalizeInequalitiesJacobianSparseValues(const VertexInterface* vertex, int& nnz_idx, Eigen::Ref<Eigen::VectorXd>& values,
                                                  const double* multipliers = nullptr, bool precompute_mixed = true);
    void finalizeActiveInequalitiesJacobianSparseValues(const VertexInterface* vertex, int& nnz_idx, Eigen::Ref<Eigen::VectorXd>& values,
                                                        double weight = 1.0);
};

FACTORY_REGISTER_HYPER_GRAPH_OPTIMIZATION_PROBLEM(HyperGraphOptimizationProblemVertexBased)

}  // namespace corbo

#endif  // SRC_OPTIMIZATION_INCLUDE_CORBO_OPTIMIZATION_HYPER_GRAPH_HYPER_GRAPH_OPTIMIZATION_PROBLEM_VERTEX_BASED_H_
