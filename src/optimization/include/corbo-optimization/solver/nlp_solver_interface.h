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

#ifndef SRC_OPTIMIZATION_INCLUDE_CORBO_OPTIMIZATION_SOLVER_NLP_SOLVER_INTERFACE_H_
#define SRC_OPTIMIZATION_INCLUDE_CORBO_OPTIMIZATION_SOLVER_NLP_SOLVER_INTERFACE_H_

#include <corbo-core/factory.h>
#include <corbo-core/types.h>
#include <corbo-optimization/optimization_problem_interface.h>
#include <corbo-optimization/types.h>

#ifdef MESSAGE_SUPPORT
#include <corbo-communication/messages/optimization/nlp_solvers.pb.h>
#endif

#include <memory>

namespace corbo {

/**
 * @brief Generic interface for solver implementations
 *
 * @ingroup optimization solver
 *
 * This class can be used to generically define solver back-ends
 * for optimization problems with objectives, equality constraints,
 * inequality constraints, box constraints and their current parameter state:
 *
 * \f[
 *    \min f(x), \ f : \mathbb{R}^n \to \mathbb{R}^m \\
 *    \text{subject to:} \\
 *    ceq(x) = 0 \\
 *    c(x) < 0 \\
 *    lb < x < ub
 * \f]
 *
 * The optimization problem is usually defined in terms of
 * an OptimizationProblemInterface instance and the solver
 * needs to check whether 1st and/or 2nd order derivatives
 * in dense or sparse form are required.
 *
 * @see OptimizationProblemInterface
 *
 * @author Christoph Rösmann (christoph.roesmann@tu-dortmund.de)
 */
class NlpSolverInterface
{
 public:
    using Ptr  = std::shared_ptr<NlpSolverInterface>;
    using UPtr = std::unique_ptr<NlpSolverInterface>;

    //! Virtual destructor
    virtual ~NlpSolverInterface() {}

    //! Return a newly created instance of the current solver
    virtual Ptr getInstance() const = 0;

    //! Return true if the solver onyl supports costs in lsq form
    virtual bool isLsqSolver() const = 0;

    /**
     * @brief Initialize the solver w.r.t. a given optimization problem
     *
     * The structure of the optimization problem might be still subject to change
     * when calling this method. However, this method can be used to preallocate
     * some memory and to check general features of the provided optimization problem,
     * e.g. if Jacobian and Hessian matrices are available.
     * @param[in]   problem     Optimization problem definition with possibly incomplete structure [optional]
     * @return true if the sovler seems to be suited, false otherwise.
     */
    virtual bool initialize(OptimizationProblemInterface* problem = nullptr) { return true; }

    /**
     * @brief Solve the provided optimization problem
     * @param[in,out]   problem        Optimization problem definition,
     *                                 pre-solving: object contains the initial parameter set
     *                                 post-solving: object contains the optimized parameter set
     * @param[in]       new_structure  true if one of the dimensions or non-zero patterns has changed.
     * @param[in]       new_run        true if a new optimization run (not just an iteration) is performed,
     *                                 e.g. to reset possible weight adaptation strategies etc.
     * @param[out]      obj_value      Returns the resulting objective function value [optional]
     * @return Solver status as indicated above
     */
    virtual SolverStatus solve(OptimizationProblemInterface& problem, bool new_structure, bool new_run = true, double* obj_value = nullptr) = 0;

    //! Clear internal caches
    virtual void clear() = 0;

#ifdef MESSAGE_SUPPORT
    // implement for import / export support
    virtual void toMessage(corbo::messages::NlpSolver& message) const {}
    virtual void fromMessage(const corbo::messages::NlpSolver& message, std::stringstream* issues = nullptr) {}
#endif
};

using NlpSolverFactory = Factory<NlpSolverInterface>;
#define FACTORY_REGISTER_NLP_SOLVER(type) FACTORY_REGISTER_OBJECT(type, NlpSolverInterface)

}  // namespace corbo

#endif  // SRC_OPTIMIZATION_INCLUDE_CORBO_OPTIMIZATION_SOLVER_NLP_SOLVER_INTERFACE_H_
