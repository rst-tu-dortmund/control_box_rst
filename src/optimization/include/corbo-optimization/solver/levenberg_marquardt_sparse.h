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

#ifndef SRC_OPTIMIZATION_INCLUDE_CORBO_OPTIMIZATION_SOLVER_LEVENBERG_MARQUARDT_SPARSE_H_
#define SRC_OPTIMIZATION_INCLUDE_CORBO_OPTIMIZATION_SOLVER_LEVENBERG_MARQUARDT_SPARSE_H_

#include <corbo-optimization/solver/nlp_solver_interface.h>

#include <Eigen/Sparse>

#ifdef CHOLMOD
#include <Eigen/CholmodSupport>
#endif

#include <memory>

namespace corbo {

/**
 * @brief Levenberg-Marquardt Solver (Sparse matrices version).
 *
 * @ingroup optimization solver
 *
 * This solver implements the Levenberg-Marquardt algorithm also implemented in [1].
 * The underlying linear system is solved using Eigens Dense Matrix algebra
 * (http://eigen.tuxfamily.org/dox/group__TutorialLinearAlgebra.html).
 *
 * The solver requires the optimization problem to be defined in Least-Squares
 * form. Bounds, inequality and equality constraints are transformed to objective functions
 * with quadratic penalties.
 * For this purpose scalar weights for all 3 approximations are defined in order to
 * configure their influence in the overall objective.
 * Theoretically, the weights should be infinity but this would cause an
 * ill-conditioned Hessian matrix and consequently slow convergence.
 * This solver implementation also provides a simple weight adapation strategy
 * which starts with user-defined weight values (setPenaltyWeights()).
 * In each invocation of the solve() method the weights are scaled by additional
 * factors (see setWeightAdapation()) up to specified maximum values.
 * Hereby, it is weight_new = weight_old * factor. The weights are resetted
 * to their original values in case solve() is called with \c new_run == true.
 *
 * [1] R. Kuemmerle, G. Grisetti, H. Strasdat, K. Konolige, W. Burgard.
 *     "g2o: A General Framework for Graph Optimization", ICRA, 2011.
 *
 * @author Christoph Rösmann (christoph.roesmann@tu-dortmund.de)
 */
class LevenbergMarquardtSparse : public NlpSolverInterface
{
 public:
    using Ptr  = std::shared_ptr<LevenbergMarquardtSparse>;
    using UPtr = std::unique_ptr<LevenbergMarquardtSparse>;

    // implements interface method
    NlpSolverInterface::Ptr getInstance() const override { return std::make_shared<LevenbergMarquardtSparse>(); }

    // implements interface method
    bool isLsqSolver() const override { return true; }

    // implements interface method
    bool initialize(OptimizationProblemInterface* problem = nullptr) override;
    // implements interface method
    SolverStatus solve(OptimizationProblemInterface& problem, bool new_structure = true, bool new_run = true, double* obj_value = nullptr) override;

    //! Define the number of solver iterations
    void setIterations(int iterations) { _iterations = iterations; }
    //! Define penalty weights (equality constraints, inequality constraints, bounds)
    void setPenaltyWeights(double weight_eq, double weight_ineq, double weight_bounds);
    //! Set parameters for weight adaptation (refer to the class description); set factors to 1 in order to disable adaptation.
    void setWeightAdapation(double factor_eq, double factor_ineq, double factor_bounds, double max_eq, double max_ineq, double max_bounds);

    void clear() override;

#ifdef MESSAGE_SUPPORT
    // implements interface method
    void toMessage(corbo::messages::NlpSolver& message) const override;
    // implements interface method
    void fromMessage(const corbo::messages::NlpSolver& message, std::stringstream* issues = nullptr) override;
#endif

 protected:
    //! Compute overall value vector including constraint approximation
    void computeValues(OptimizationProblemInterface& problem);

    //! Reset weights to their original values
    void resetWeights();
    //! Perform single weight adapation step
    void adaptWeights();

 private:
    // parameters
    int _iterations = 10;

    double _weight_init_eq     = 2;
    double _weight_init_ineq   = 2;
    double _weight_init_bounds = 2;

    double _weight_adapt_factor_eq     = 1;
    double _weight_adapt_factor_ineq   = 1;
    double _weight_adapt_factor_bounds = 1;

    double _weight_adapt_max_eq     = 500;
    double _weight_adapt_max_ineq   = 500;
    double _weight_adapt_max_bounds = 500;

    // internal states
    int _param_dim         = 0;
    int _obj_dim           = 0;
    int _eq_dim            = 0;
    int _ineq_dim          = 0;
    int _finite_bounds_dim = 0;
    int _val_dim           = 0;

    Eigen::VectorXd _values;
    Eigen::SparseMatrix<double> _jacobian;
    Eigen::SparseMatrix<double> _hessian;
    Eigen::VectorXd _delta;
    Eigen::VectorXd _rhs;

    double _weight_eq     = _weight_init_eq;
    double _weight_ineq   = _weight_init_ineq;
    double _weight_bounds = _weight_init_bounds;

/**
 * @brief Eigens sparse solver wrapper.
 * Check http://eigen.tuxfamily.org/dox/group__TopicSparseSystems.html for further information and different solvers.
 * The second template parameter specifies, whether the upper or lower triangular part should be used.
 * If CHOLMOD was found by CMake, the cholmod supernodal llt version is used rathar than Eigens simplicial ldlt.
 * Define \c FORCE_EIGEN_SOLVER if the default LDLT should sill be used even if CHOLMOD is found.
 */
#if defined(CHOLMOD) && !defined(FORCE_EIGEN_SOLVER)
    Eigen::CholmodSupernodalLLT<Eigen::SparseMatrix<double>, Eigen::Upper> _sparse_solver;
#else
    // Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>, Eigen::Upper> _sparse_solver;
    Eigen::SimplicialLLT<Eigen::SparseMatrix<double>, Eigen::Upper> _sparse_solver;
#endif
};

FACTORY_REGISTER_NLP_SOLVER(LevenbergMarquardtSparse)

}  // namespace corbo

#endif  // SRC_OPTIMIZATION_INCLUDE_CORBO_OPTIMIZATION_SOLVER_LEVENBERG_MARQUARDT_SPARSE_H_
