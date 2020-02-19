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

#ifndef SRC_OPTIMIZATION_INCLUDE_CORBO_OPTIMIZATION_SOLVER_QP_SOLVER_INTERFACE_H_
#define SRC_OPTIMIZATION_INCLUDE_CORBO_OPTIMIZATION_SOLVER_QP_SOLVER_INTERFACE_H_

#include <corbo-core/factory.h>
#include <corbo-core/types.h>
#include <corbo-optimization/optimization_problem_interface.h>
#include <corbo-optimization/types.h>

#ifdef MESSAGE_SUPPORT
#include <corbo-communication/messages/optimization/qp_solvers.pb.h>
#endif

#include <memory>

namespace corbo {

/**
 * @brief Generic interface for QP solver implementations
 *
 * @ingroup optimization solver
 *
 * This class can be used to generically define solver back-ends
 * for optimization problems with quadratic objectives and linear constraints.
 *
 * Since we are currently working with OSQP, the interface is primarly defined to
 * perfectly match the OSQP interface.
 * This means we require SparseMatrix indices vectors to be defined in long long format.
 * Furthermore, the problem definition is as follows
 * \f[
 *    \min 0.5 x^T P x + q^T x \\
 *    \text{subject to:} \\
 *    lbA <= A x <= ubB \\
 *    lb <= x <= ub
 * \f]
 *
 * The gradient of the Lagrangian is defined as follows:
 * grad L(x, lambda, mu) = grad f(x) +  A^T lambda + mu
 *
 * Some QP solvers do not allow the special treatment of bounds.
 * The particular solver implementation creates a new matrix Atile = [A; I] for that purpose.
 * However, in order to avoid this step, isSupportingSimpleBounds() can be checked a-priori.
 *
 * Note, bounds can be set to corbo_DBL_INF and -corbo_DBL_INF respectively.
 *
 * @see OptimizationProblemInterface
 *
 * @author Christoph Rösmann (christoph.roesmann@tu-dortmund.de)
 */
class QpSolverInterface
{
 public:
    using Ptr  = std::shared_ptr<QpSolverInterface>;
    using UPtr = std::unique_ptr<QpSolverInterface>;

    using SparseMatrix = Eigen::SparseMatrix<double, Eigen::ColMajor, long long>;

    //! Virtual destructor
    virtual ~QpSolverInterface() = default;

    //! Return a newly created instance of the current solver
    virtual Ptr getInstance() const = 0;

    virtual bool isSupportingSimpleBounds() = 0;

    /**
     * @brief Initialize the qp solver
     *
     */
    virtual bool initialize() { return true; }

    /**
     * @brief Solve full QP
     * @remarks We require that all arguments are non-const in order to facilitate the usage of some lowlevel c-APIs (e.g. OSQP).
     * @param[in]   P                Quadratic form matrix,
     * @param[in]   q                Linear part
     * @param[in]   A                Linear constraint matrix
     * @param[in]   lbA              Lower bounds of linear constraints
     * @param[in]   ubA              Upper bounds of linear constraints
     * @param[in]   lb               Lower bounds of parameters
     * @param[in]   ub               Upper bounds of parameters
     * @param[in]   new_structure    If true, the solver must reallocte the internal memory
     *
     * @return Solver status
     */
    virtual SolverStatus solve(SparseMatrix& P, Eigen::Ref<Eigen::VectorXd> q, SparseMatrix& A, Eigen::Ref<Eigen::VectorXd> lbA,
                               Eigen::Ref<Eigen::VectorXd> ubA, Eigen::Ref<Eigen::VectorXd> lb, Eigen::Ref<Eigen::VectorXd> ub,
                               bool new_structure = true, bool zero_x_warmstart = false) = 0;

    /**
     * @brief Solve QP without simple bounds
     * @remarks We require that all arguments are non-const in order to facilitate the usage of some lowlevel c-APIs (e.g. OSQP).
     * @param[in]   P                Quadratic form matrix,
     * @param[in]   q                Linear part
     * @param[in]   A                Linear constraint matrix
     * @param[in]   lbA              Lower bounds of linear constraints
     * @param[in]   ubA              Upper bounds of linear constraints
     * @param[in]   new_structure    If true, the solver must reallocte the internal memory
     *
     * @return Solver status
     */
    virtual SolverStatus solve(SparseMatrix& P, Eigen::Ref<Eigen::VectorXd> q, SparseMatrix& A, Eigen::Ref<Eigen::VectorXd> lbA,
                               Eigen::Ref<Eigen::VectorXd> ubA, bool new_structure = true, bool zero_x_warmstart = false, bool update_P = true,
                               bool update_q = true, bool update_A = true, bool update_bounds = true) = 0;

    virtual Eigen::Ref<Eigen::VectorXd> getPrimalSolution() = 0;
    virtual Eigen::Ref<Eigen::VectorXd> getDualSolution()   = 0;

    virtual void updatePrimalSolutionWarmStart(const Eigen::Ref<const Eigen::VectorXd>& x) = 0;
    virtual void updateDualSolutionWarmStart(const Eigen::Ref<const Eigen::VectorXd>& y)   = 0;

    virtual void enforceNewStructure(bool new_structure = true) = 0;

    //! clear internal caches
    virtual void clear() = 0;

#ifdef MESSAGE_SUPPORT
    // implement for import / export support
    virtual void toMessage(corbo::messages::QpSolver& message) const {}
    virtual void fromMessage(const corbo::messages::QpSolver& message, std::stringstream* issues = nullptr) {}
#endif
};

using QpSolverFactory = Factory<QpSolverInterface>;
#define FACTORY_REGISTER_QP_SOLVER(type) FACTORY_REGISTER_OBJECT(type, QpSolverInterface)

}  // namespace corbo

#endif  // SRC_OPTIMIZATION_INCLUDE_CORBO_OPTIMIZATION_SOLVER_QP_SOLVER_INTERFACE_H_
