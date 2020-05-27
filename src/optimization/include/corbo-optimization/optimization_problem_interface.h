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

#ifndef SRC_OPTIMIZATION_INCLUDE_CORBO_OPTIMIZATION_OPTIMIZATION_PROBLEM_INTERFACE_H_
#define SRC_OPTIMIZATION_INCLUDE_CORBO_OPTIMIZATION_OPTIMIZATION_PROBLEM_INTERFACE_H_

#include <corbo-core/macros.h>
#include <corbo-core/types.h>

#include <Eigen/Sparse>

#include <functional>
#include <memory>

namespace corbo {

/**
 * @brief Generic interface for optimization problem definitions
 *
 * @ingroup optimization
 *
 * This class can be used to generically define optimization problems
 * with objectives, equality constraints, inequality constraints,
 * box constraints and their current parameter state.
 *
 * \f[
 *    \min f(x), \ f : \mathbb{R}^n \to \mathbb{R}^m \\
 *    \text{subject to:} \\
 *    ceq(x) = 0 \\
 *    c(x) < 0 \\
 *    lb < x < ub
 * \f]
 *
 * First order and second order derivatives can be provided but
 * their actual usage depends on the SolverInterface back-end
 * which solves the optimization problem. Some solvers may require
 * the implementation of those methods.
 *
 * Solvers determine a parameter increment in each of their iterations
 * so this class provides a method applyIncrement() to modify the underlying
 * optimization parameter set. In order to convinentially store and restore
 * the current parameter set from the solver backend side,
 * appropriate backup methods must be provided by subclasses.
 *
 * @see SolverInterface
 *
 * @author Christoph Rösmann (christoph.roesmann@tu-dortmund.de)
 */
class OptimizationProblemInterface
{
 public:
    using Ptr  = std::shared_ptr<OptimizationProblemInterface>;
    using UPtr = std::unique_ptr<OptimizationProblemInterface>;

    virtual ~OptimizationProblemInterface() {}

    virtual void clear() {}

    /** @name Specify the dimension of the optimization problem  */
    //@{

    //! Get dimension of the objective (should be zero or one)
    virtual int getNonLsqObjectiveDimension() = 0;
    //! Total dimension of least-squares objective function terms
    virtual int getLsqObjectiveDimension() = 0;
    //! Get dimension of the objective (should be zero or one, includes Lsq objectives if present)
    virtual int getObjectiveDimension() = 0;

    //! Total dimension of equality constraints
    virtual int getEqualityDimension() = 0;
    //! Total dimension of general inequality constraints
    virtual int getInequalityDimension() = 0;
    //! Effictive dimension of the optimization parameter set (changeable, non-fixed part)
    virtual int getParameterDimension() = 0;

    //@}

    /** @name Generic interface to interact with the underlying parameter/variables representation  */
    //@{

    /**
     * @brief Apply increment to the current parameter set
     *
     * Usually solvers determine a new step/increment \f$ x += \Delta x \f$
     * in each iteration. This method applies the parameter update.
     * @remarks Make sure to update only non-fixed parameters according to parameterDimension().
     * @remarks It is highly recommended to override this method for particular problems in order to speed up parameter updates
     * @param[in]   increment   Parameter update [parameterDimension() x 1]
     */
    virtual void applyIncrement(const Eigen::Ref<const Eigen::VectorXd>& increment);

    /**
     * @brief Apply increment to the current parameter set (single element overload)
     *
     * @param[in]   idx         Value index according vector size [parameterDimension() x 1]
     * @param[in]   increment   Parameter update [parameterDimension() x 1]
     */
    virtual void applyIncrement(int idx, double increment);

    /**
     * @brief Return specific value of the parameter vector
     * @param[in] idx  Value index according vector size [parameterDimension() x 1]
     * @return value
     */
    virtual double getParameterValue(int idx) = 0;

    /**
     * @brief Set specific value of the parameter vector
     * @param[in] idx  Value index according vector size [parameterDimension() x 1]
     * @param[in] x    New value to be set
     */
    virtual void setParameterValue(int idx, double x) = 0;

    /**
     * @brief Return deep copy of the complete parameter vector
     * @param[out] x  Current parameter vector of size [parameterDimension() x 1] (must be pre-allocated!)
     */
    virtual void getParameterVector(Eigen::Ref<Eigen::VectorXd> x);

    /**
     * @brief Set complete parameter vector
     * @remarks The default implementation requires parameterDimension() to be valid, hence no automatic resizing
     * @param[in] x    New parameter vector of size [parameterDimension() x 1]
     */
    virtual void setParameterVector(const Eigen::Ref<const Eigen::VectorXd>& x);

    //!< Push the current parameter set to some backup stack
    virtual void backupParameters() = 0;
    //!< Restore parameter set from the last backup and keep backup if desired
    virtual void restoreBackupParameters(bool keep_backup) = 0;
    //!< Discard last backup (or all)
    virtual void discardBackupParameters(bool all = false) = 0;

    //@}

    /** @name Specify properties of the optimization problem  */
    //@{

    /**
     * @brief Check if the underlying problem is defined in the least squares form
     *
     * Least-squares problems are defined as \f$ \min f(x)^T f(x) \f$ and the
     * function values and Jacobian are computed for \f$ f(x) \f$
     * rather than for \f$ f(x)^T f(x) \f$. Specialiezed least-squares solvers
     * require the optimization problem to be defined in this particular form.
     *
     * In summary, the solver computes the squared l2-norm of the
     * objective function value (scalar) if set to true.
     *
     * @return true if all objectives are in least-squares form, false otherwise.
     */
    virtual bool isLeastSquaresProblem() const = 0;

    //@}

    /** @name Access bounds  */
    //@{

    /**
     * @brief Return specific lower bound value of a parameter
     * @param[in] idx  Value index according vector size [parameterDimension() x 1]
     * @return lower bound value
     */
    virtual double getLowerBound(int idx) = 0;

    /**
     * @brief Return specific upper bound of a parameter
     * @param[in] idx  Value index according vector size [parameterDimension() x 1]
     * @return upper bound value
     */
    virtual double getUpperBound(int idx) = 0;

    /**
     * @brief Set specific lower bound of a parameter
     * @param[in] idx  Value index according vector size [parameterDimension() x 1]
     * @param[in] lb   New lower bound to be set
     */
    virtual void setLowerBound(int idx, double lb) = 0;

    /**
     * @brief Set specific upper bound of a parameter
     * @param[in] idx  Value index according vector size [parameterDimension() x 1]
     * @param[in] ub   New upper bound to be set
     */
    virtual void setUpperBound(int idx, double ub) = 0;

    /**
     * @brief Set lower and upper bound vector
     * @param[in] lb  Lower bounds [parameterDimension() x 1]
     * @param[in] ub  Upper bounds [parameterDimension() x 1]
     */
    virtual void setBounds(const Eigen::Ref<const Eigen::VectorXd>& lb, const Eigen::Ref<const Eigen::VectorXd>& ub);

    /**
     * @brief Get lower and upper bound vector
     * @param[out] lb  Lower bounds [parameterDimension() x 1] (must be pre-allocated!)
     * @param[out] ub  Upper bounds [parameterDimension() x 1] (must be pre-allocated!)
     */
    virtual void getBounds(Eigen::Ref<Eigen::VectorXd> lb, Eigen::Ref<Eigen::VectorXd> ub);

    //@}

    /** @name Specify main equations of the optimization problem  */
    //@{

    /**
     * @brief Compute the objective function values f(x) for the current parameter set
     * @param[out]  values      The resulting value vector [objectiveDimension() x 1].
     *                          Warning: the size must be pre-allocated!
     */
    virtual void computeValuesLsqObjective(Eigen::Ref<Eigen::VectorXd> values) = 0;

    virtual double computeValueNonLsqObjective() = 0;

    virtual double computeValueObjective();

    /**
     * @brief Compute the equality constraint values ceq(x) for the current parameter set
     * @param[out]  values      The resulting value vector [equalityDimension() x 1].
     *                                        Warning: the size must be pre-allocated!
     */
    virtual void computeValuesEquality(Eigen::Ref<Eigen::VectorXd> values) = 0;

    /**
     * @brief Compute the inequality constraint values c(x) for the current parameter set
     * @param[out]  values      The resulting value vector [inequalityDimension() x 1].
     *                                        Warning: the size must be pre-allocated!
     */
    virtual void computeValuesInequality(Eigen::Ref<Eigen::VectorXd> values) = 0;

    virtual void computeValues(double& non_lsq_obj_value, Eigen::Ref<Eigen::VectorXd> lsq_obj_values, Eigen::Ref<Eigen::VectorXd> eq_values,
                               Eigen::Ref<Eigen::VectorXd> ineq_values)
    {
        non_lsq_obj_value = computeValueNonLsqObjective();
        computeValuesLsqObjective(lsq_obj_values);
        computeValuesEquality(eq_values);
        computeValuesInequality(ineq_values);
    }

    //@}

    /** @name Methods for dealing with bounds  */
    //@{

    /**
     * @brief Dimension of the set of finite bounds (combined such that each ub and lb component define a single dimension)
     *
     * This method returns the number of parameters that exhibit either finite lower or finite upper bounds.
     * In particular, it counts how often (is_finite(lb[i]) || is_finite(ub[i])) evaluates to true.
     * The maximum value is parameterDimension().
     *
     */
    virtual int finiteCombinedBoundsDimension();

    /**
     * @brief Dimension of the set of finite bounds (individual bounds ub and lb)
     *
     * This method retuns the total number of finite bounds, in particular:
     * numFinite(lb) + numFinite(ub) with numFinite counting values within interval [-CORBO_INF_DBL, CORBO_INF_DBL].
     * If all bounds are finite, the (maximum) value is 2*parameterDimension().
     *
     * @returns number of finite lower bounds + number of finite upper bounds
     */
    virtual int finiteBoundsDimension();

    //@}

    /**
     * @brief Compute the values of the active inequality constraints (elementwise max(0, c(x)))
     *
     * Return a vector [inequalityDimension() x 1] in which elements are set to max(0, c(x))
     *
     * Weight scales the values to "weight*max(0, c(x))".
     *
     * @param[out]  values      The resulting value vector [inequalityDimension() x 1].
     *                                        Warning: the size must be pre-allocated!
     * @param[in]   weight       Optionally provide a weight factor
     */
    virtual void computeValuesActiveInequality(Eigen::Ref<Eigen::VectorXd> values, double weight = 1.0);

    /**
     * @brief Compute the distance to finite bound values (combined lower and upper)
     *
     * The distance is given as (lb-x) if (x<lb) and (x-ub) if (x>ub), otherwise 0.
     *
     * @warning This dimension of the resulting vector corresponds to finite bounds only,
     *          hence the bound bust be finite and the value is non-zero if x < lb || x > ub,
     *          this method should be used according with computeFiniteCombinedBoundsJacobian().
     * @param[out]  values      The resulting value vector [finiteCombinedBoundsDimension() x 1].
     *                          Warning: the size must be pre-allocated!
     */
    virtual void computeDistanceFiniteCombinedBounds(Eigen::Ref<Eigen::VectorXd> values);

    /**
     * @brief Compute the distance between parameters and bounds
     *
     * This method might be used for sqp bound computation in the form
     *
     *      lb-x < dx < ub-x
     *
     * @warning This dimension of the resulting vector corresponds to finite bounds only,
     *          hence the bound bust be finite and the value is non-zero if x < lb || x > ub.
     * @param[out]  lb_minus_x      The resulting value vector [finiteCombinedBoundsDimension() x 1]; Warning: the size must be pre-allocated!
     * @param[out]  ub_minus_x      The resulting value vector [finiteCombinedBoundsDimension() x 1]; Warning: the size must be pre-allocated!
     *
     */
    virtual void computeLowerAndUpperBoundDiff(Eigen::Ref<Eigen::VectorXd> lb_minus_x, Eigen::Ref<Eigen::VectorXd> ub_minus_x);

    /**
     * @brief Return bound and parameter vectors only for finite boudns
     *
     * The values are returned if lb < inf or ub < inf.
     * The Dimension of each vector is finiteCombinedBoundsDimension().
     *
     * @warning This dimension of the resulting vector corresponds to finite bounds only,
     *          this method should be used according with computeFiniteCombinedBoundsJacobian().
     * @param[out]  lb_finite_bounds          Lower bounds [finiteCombinedBoundsDimension() x 1], Warning: the size must be pre-allocated!
     * @param[out]  ub_finite_bounds          Upper bounds [finiteCombinedBoundsDimension() x 1], Warning: the size must be pre-allocated!
     * @param[out]  x_finite_bounds           related parameters [finiteCombinedBoundsDimension() x 1], Warning: the size must be pre-allocated!
     */
    virtual void getParametersAndBoundsFinite(Eigen::Ref<Eigen::VectorXd> lb_finite_bounds, Eigen::Ref<Eigen::VectorXd> ub_finite_bounds,
                                              Eigen::Ref<Eigen::VectorXd> x_finite_bounds);

    virtual void computeGradientObjective(Eigen::Ref<Eigen::VectorXd> gradient);

    virtual void computeGradientNonLsqObjective(Eigen::Ref<Eigen::VectorXd> gradient);

    /**
     * @brief Compute the objective Jacobian Jf(x) for the current parameter set
     * @param[out]  jacobian    The resulting Jacobian matrix [objectiveDimension() x parameterDimension()].
     *                                          Warning: the matrix must be pre-allocated with ZEROS!
     * @param[in]   multipliers    Optionally provide a vector of multipliers to scale each cost function term [objectiveDimension() x 1]
     */
    virtual void computeDenseJacobianLsqObjective(Eigen::Ref<Eigen::MatrixXd> jacobian, const double* multipliers = nullptr);
    virtual int computeSparseJacobianLsqObjectiveNNZ();
    virtual void computeSparseJacobianLsqObjectiveStructure(Eigen::Ref<Eigen::VectorXi> i_row, Eigen::Ref<Eigen::VectorXi> j_col);
    virtual void computeSparseJacobianLsqObjectiveValues(Eigen::Ref<Eigen::VectorXd> values, const double* multipliers = nullptr);
    virtual void computeSparseJacobianLsqObjective(Eigen::SparseMatrix<double>& jacobian, const double* multipliers = nullptr);

    /**
     * @brief Compute the equality constraint Jacobian Jceq(x) for the current parameter set
     * @param[out]  jacobian      The resulting Jacobian matrix [equalityDimension() x parameterDimension()].
     *                                           Warning: the matrix must be pre-allocated with ZEROS!
     * @param[in]   multipliers    Optionally provide a vector of multipliers to scale each cost function term [equalityDimension() x 1]
     */
    virtual void computeDenseJacobianEqualities(Eigen::Ref<Eigen::MatrixXd> jacobian, const double* multipliers = nullptr);
    virtual int computeSparseJacobianEqualitiesNNZ();
    virtual void computeSparseJacobianEqualitiesStructure(Eigen::Ref<Eigen::VectorXi> i_row, Eigen::Ref<Eigen::VectorXi> j_col);
    virtual void computeSparseJacobianEqualitiesValues(Eigen::Ref<Eigen::VectorXd> values, const double* multipliers = nullptr);
    virtual void computeSparseJacobianEqualities(Eigen::SparseMatrix<double>& jacobian, const double* multipliers = nullptr);

    /**
     * @brief Compute the inequality constraint Jacobian Jc(x) for the current parameter set
     * @param[out]  jacobian      The resulting Jacobian matrix [inequalityDimension() x parameterDimension()].
     *                                           Warning: the matrix must be pre-allocated with ZEROS!
     * @param[in]   multipliers    Optionally provide a vector of multipliers to scale each cost function term [inequalityDimension() x 1]
     */
    virtual void computeDenseJacobianInequalities(Eigen::Ref<Eigen::MatrixXd> jacobian, const double* multipliers = nullptr);
    virtual int computeSparseJacobianInequalitiesNNZ();
    virtual void computeSparseJacobianInequalitiesStructure(Eigen::Ref<Eigen::VectorXi> i_row, Eigen::Ref<Eigen::VectorXi> j_col);
    virtual void computeSparseJacobianInequalitiesValues(Eigen::Ref<Eigen::VectorXd> values, const double* multipliers = nullptr);
    virtual void computeSparseJacobianInequalities(Eigen::SparseMatrix<double>& jacobian, const double* multipliers = nullptr);

    /**
     * @brief Compute the Jacobian Jc(x) with non-zeros for active constraints c(x)>= 0 and zeros for inactive ones
     * @param[out]  jacobian      The resulting Jacobian matrix [inequalityDimension() x parameterDimension()].
     *                            Warning: the matrix must be pre-allocated with ZEROS!
     * @param[in]   weight        Optionally provide a weight factor
     */
    virtual void computeDenseJacobianActiveInequalities(Eigen::Ref<Eigen::MatrixXd> jacobian, double weight = 1.0);
    virtual void computeSparseJacobianActiveInequalitiesValues(Eigen::Ref<Eigen::VectorXd> values, double weight = 1.0);
    virtual void computeSparseJacobianActiveInequalities(Eigen::SparseMatrix<double>& jacobian, double weight = 1.0);

    /**
     * @brief Compute the Jacobian for finite combined bounds
     * @details refer to computeDistanceFiniteCombinedBounds()
     * @param[out]  jacobian       The resulting Jacobian matrix [finiteCombinedBoundsDimension() x parameterDimension()].
     *                             Warning: the matrix must be pre-allocated with ZEROS!
     * @param[in]   weight         Optionally provide a weight factor
     */
    virtual void computeDenseJacobianFiniteCombinedBounds(Eigen::Ref<Eigen::MatrixXd> jacobian, double weight = 1.0);
    virtual int computeSparseJacobianFiniteCombinedBoundsNNZ();
    virtual void computeSparseJacobianFiniteCombinedBoundsStructure(Eigen::Ref<Eigen::VectorXi> i_row, Eigen::Ref<Eigen::VectorXi> j_col);
    virtual void computeSparseJacobianFiniteCombinedBoundsValues(Eigen::Ref<Eigen::VectorXd> values, double weight = 1.0);
    virtual void computeSparseJacobianFiniteCombinedBounds(Eigen::SparseMatrix<double>& jacobian, double weight = 1.0);

    /**
     * @brief Compute the Jacobian for finite combined bounds
     * @details refer to getParametersAndBoundsFinite()
     *          The resulting matrix is just the identity matrix with omitted rows for non-finite bounds.
     * @param[out]  jacobian       The resulting Jacobian matrix [finiteCombinedBoundsDimension() x parameterDimension()].
     */
    virtual void computeDenseJacobianFiniteCombinedBoundsIdentity(Eigen::Ref<Eigen::MatrixXd> jacobian);

    /**
     * @brief Compute the objective and constraint Jacobians at once
     *
     * Use this method to speed up Jacobian computation if possible. The default implementation
     * just calls the individual Jacobian computations.
     *
     * @remarks This method does not compute Jacobians for bounds due to their simplicity, compute them manually.
     *
     * @param[out]  jacobian_obj     The resulting objective Jacobian matrix [objectiveDimension() x parameterDimension()].
     *                                                 Warning: the matrix must be pre-allocated with ZEROS!
     * @param[out]  jacobian_eq      The resulting equality constraint Jacobian matrix [equalityDimension() x parameterDimension()].
     *                                                 Warning: the matrix must be pre-allocated with ZEROS!
     * @param[out]  jacobian_ineq   The resulting inequality constraint Jacobian matrix [inequalityDimension() x parameterDimension()].
     *                                                 Warning: the matrix must be pre-allocated with ZEROS!
     * @param[in]   multipliers          Optionally provide a vector of multipliers to scale each cost function term
     *                                                 [objectiveDimension() + equalityDimension() + inequalityDimension()  x 1]
     */
    virtual void computeDenseJacobians(Eigen::Ref<Eigen::VectorXd> gradient_non_lsq_obj, Eigen::Ref<Eigen::MatrixXd> jacobian_lsq_obj,
                                       Eigen::Ref<Eigen::MatrixXd> jacobian_eq, Eigen::Ref<Eigen::MatrixXd> jacobian_ineq,
                                       const double* multipliers_lsq_obj = nullptr, const double* multipliers_eq = nullptr,
                                       const double* multipliers_ineq = nullptr, bool active_ineq = false, double active_ineq_weight = 1.0);
    virtual void computeSparseJacobiansNNZ(int& nnz_lsq_obj, int& nnz_eq, int& nnz_ineq);
    virtual void computeSparseJacobiansStructure(Eigen::Ref<Eigen::VectorXi> i_row_obj, Eigen::Ref<Eigen::VectorXi> j_col_obj,
                                                 Eigen::Ref<Eigen::VectorXi> i_row_eq, Eigen::Ref<Eigen::VectorXi> j_col_eq,
                                                 Eigen::Ref<Eigen::VectorXi> i_row_ineq, Eigen::Ref<Eigen::VectorXi> j_col_ineq);
    virtual void computeSparseJacobiansValues(Eigen::Ref<Eigen::VectorXd> values_obj, Eigen::Ref<Eigen::VectorXd> values_eq,
                                              Eigen::Ref<Eigen::VectorXd> values_ineq, const double* multipliers_obj = nullptr,
                                              const double* multipliers_eq = nullptr, const double* multipliers_ineq = nullptr,
                                              bool active_ineq = false, double active_ineq_weight = 1.0);
    virtual void computeSparseJacobians(Eigen::SparseMatrix<double>& jacobian_lsq_obj, Eigen::SparseMatrix<double>& jacobian_eq,
                                        Eigen::SparseMatrix<double>& jacobian_ineq, const double* multipliers_lsq_obj = nullptr,
                                        const double* multipliers_eq = nullptr, const double* multipliers_ineq = nullptr, bool active_ineq = false,
                                        double active_ineq_weight = 1.0);

    virtual void computeCombinedSparseJacobian(Eigen::SparseMatrix<double>& jacobian, bool objective_lsq, bool equality, bool inequality,
                                               bool finite_combined_bounds, bool active_ineq = false, double weight_eq = 1.0,
                                               double weight_ineq = 1.0, double weight_bounds = 1.0, const Eigen::VectorXd* values = nullptr,
                                               const Eigen::VectorXi* col_nnz = nullptr);

    virtual int computeCombinedSparseJacobiansNNZ(bool objective_lsq = true, bool equality = true, bool inequality = true);
    virtual void computeCombinedSparseJacobiansStructure(Eigen::Ref<Eigen::VectorXi> i_row, Eigen::Ref<Eigen::VectorXi> j_col,
                                                         bool objective_lsq = true, bool equality = true, bool inequality = true);
    virtual void computeCombinedSparseJacobiansValues(Eigen::Ref<Eigen::VectorXd> values, bool objective_lsq = true, bool equality = true,
                                                      bool inequality = true, const double* multipliers_obj = nullptr,
                                                      const double* multipliers_eq = nullptr, const double* multipliers_ineq = nullptr);

    // useful for IPOPT (w/ hessian-approx)
    virtual void computeGradientObjectiveAndCombinedSparseJacobiansValues(Eigen::Ref<Eigen::VectorXd> gradient,
                                                                          Eigen::Ref<Eigen::VectorXd> jac_values, bool equality = true,
                                                                          bool inequality = true, const double* multipliers_eq = nullptr,
                                                                          const double* multipliers_ineq = nullptr);

    /**
     * @brief Compute the objective Hessian Hf(x) for the current parameter set
     * @remarks In case objectiveDimension() > 1: all individual hessians are summed up.
     * @param[in]    jacobian        Jacobian matrix which is used for the Hessian computation [objectiveDimension() x parameterDimension()]
     * @param[out]   hessian         The resulting Hessian matrix [parameterDimension() x parameterDimension()].
     *                               Warning: The matrix must be pre-allocated, the resulting hessian is ADDED to the
     *                               previous values in hessian, e.g. hessian += dense_hessian_objective.
     * @param[in]   multipliers      Optionally provide a vector of multipliers to scale each cost function term [objectiveDimension() x 1]
     * @param[in]   jacob_scaled     If multiplier is valid and selected, specify whether the \c jacobian is already scaled with \c multiplier.
     */
    virtual void computeDenseHessianObjective(const Eigen::Ref<const Eigen::MatrixXd>& jacobian, Eigen::Ref<Eigen::MatrixXd> hessian,
                                              const double* multipliers = nullptr, bool jacob_scaled = true);

    virtual void computeDenseHessianObjective(Eigen::Ref<Eigen::MatrixXd> hessian, double multiplier = 1.0);
    virtual int computeSparseHessianObjectiveNNZ(bool lower_part_only = false);
    virtual void computeSparseHessianObjectiveStructure(Eigen::Ref<Eigen::VectorXi> i_row, Eigen::Ref<Eigen::VectorXi> j_col,
                                                        bool lower_part_only = false);
    virtual void computeSparseHessianObjectiveValues(Eigen::Ref<Eigen::VectorXd> values, double multiplier = 1.0, bool lower_part_only = false);
    virtual void computeSparseHessianObjective(Eigen::SparseMatrix<double>& hessian, double multiplier = 1.0);

    virtual void computeSparseHessianObjectiveLL(Eigen::SparseMatrix<double, Eigen::ColMajor, long long>& H, const Eigen::VectorXi* col_nnz = nullptr,
                                                 bool upper_part_only = false);
    virtual void computeSparseHessianObjectiveNNZperCol(Eigen::Ref<Eigen::VectorXi> col_nnz, bool upper_part_only = false);

    virtual void computeDenseHessianEqualities(Eigen::Ref<Eigen::MatrixXd> hessian, const double* multipliers = nullptr);
    virtual int computeSparseHessianEqualitiesNNZ(bool lower_part_only = false);
    virtual void computeSparseHessianEqualitiesStructure(Eigen::Ref<Eigen::VectorXi> i_row, Eigen::Ref<Eigen::VectorXi> j_col,
                                                         bool lower_part_only = false);
    virtual void computeSparseHessianEqualitiesValues(Eigen::Ref<Eigen::VectorXd> values, const double* multipliers = nullptr,
                                                      bool lower_part_only = false);
    virtual void computeSparseHessianEqualities(Eigen::SparseMatrix<double>& hessian, const double* multipliers = nullptr);

    virtual void computeDenseHessianInequalities(Eigen::Ref<Eigen::MatrixXd> hessian, const double* multipliers = nullptr);
    virtual int computeSparseHessianInequalitiesNNZ(bool lower_part_only = false);
    virtual void computeSparseHessianInequalitiesStructure(Eigen::Ref<Eigen::VectorXi> i_row, Eigen::Ref<Eigen::VectorXi> j_col,
                                                           bool lower_part_only = false);
    virtual void computeSparseHessianInequalitiesValues(Eigen::Ref<Eigen::VectorXd> values, const double* multipliers = nullptr,
                                                        bool lower_part_only = false);
    virtual void computeSparseHessianInequalities(Eigen::SparseMatrix<double>& hessian, const double* multipliers = nullptr);

    //    virtual void computeDenseJacobianHessianObjective(Eigen::Ref<Eigen::MatrixXd> jacobian, Eigen::Ref<Eigen::MatrixXd> hessian,
    //                                                      const double* multipliers = nullptr);
    //    virtual void computeDenseJacobianHessianEqualities(Eigen::Ref<Eigen::MatrixXd> jacobian, Eigen::Ref<Eigen::MatrixXd> hessian,
    //                                                       const double* multipliers = nullptr);
    //    virtual void computeDenseJacobianHessianInequalities(Eigen::Ref<Eigen::MatrixXd> jacobian, Eigen::Ref<Eigen::MatrixXd> hessian,
    //                                                         const double* multipliers = nullptr);

    virtual void computeDenseHessians(Eigen::Ref<Eigen::MatrixXd> hessian_obj, Eigen::Ref<Eigen::MatrixXd> hessian_eq,
                                      Eigen::Ref<Eigen::MatrixXd> hessian_ineq, double multiplier_obj = 1.0, const double* multipliers_eq = nullptr,
                                      const double* multipliers_ineq = nullptr);
    virtual void computeSparseHessians(Eigen::SparseMatrix<double>& hessian_obj, Eigen::SparseMatrix<double>& hessian_eq,
                                       Eigen::SparseMatrix<double>& hessian_ineq, double multiplier_obj = 1.0, const double* multipliers_eq = nullptr,
                                       const double* multipliers_ineq = nullptr);
    virtual void computeSparseHessiansNNZ(int& nnz_obj, int& nnz_eq, int& nnz_ineq, bool lower_part_only = false);
    virtual void computeSparseHessiansStructure(Eigen::Ref<Eigen::VectorXi> i_row_obj, Eigen::Ref<Eigen::VectorXi> j_col_obj,
                                                Eigen::Ref<Eigen::VectorXi> i_row_eq, Eigen::Ref<Eigen::VectorXi> j_col_eq,
                                                Eigen::Ref<Eigen::VectorXi> i_row_ineq, Eigen::Ref<Eigen::VectorXi> j_col_ineq,
                                                bool lower_part_only = false);
    virtual void computeSparseHessiansValues(Eigen::Ref<Eigen::VectorXd> values_obj, Eigen::Ref<Eigen::VectorXd> values_eq,
                                             Eigen::Ref<Eigen::VectorXd> values_ineq, double multiplier_obj = 1.0,
                                             const double* multipliers_eq = nullptr, const double* multipliers_ineq = nullptr,
                                             bool lower_part_only = false);

    /**
     * @brief Compute the hessian of the lagrangian L(x) = f(x) + lambda1 * c(x) + lambda2 * ceq(x)
     */
    virtual void computeSparseHessianLagrangian(Eigen::SparseMatrix<double, Eigen::ColMajor, long long>& H, const double* multipliers_eq,
                                                const double* multipliers_ineq, const Eigen::VectorXi* col_nnz = nullptr,
                                                bool upper_part_only = false);
    virtual void computeSparseHessianLagrangianNNZperCol(Eigen::Ref<Eigen::VectorXi> col_nnz,
                                                         bool upper_part_only);  // TODO(roesmann): lower_part_only? does it make sense for NNZ?

    /** @name Methods specialized for QP forms  */
    //@{

    /**
     * @brief Compute the jacobian A for the linear form lbA <= A x <= lbB
     *
     * Especially some QP solvers require either
     *
     *      lbA <= A x<= ubA
     *
     *  or
     *
     *     lbA <= A x <= ubA
     *     lb <= x <= ub
     *
     * The first formulatin includes bounds and hence the bottom rows of A are the identity.
     * However, in case lb and ub are infinity, A might have obsolete rows.
     *
     * If parameter \c include_bounds is \c true: only finite bounds are considered and
     * matrix A is of size [getEqualityDimension() + getInequalityDimension() + finiteCombinedBounds() x parameterDimension()].
     *
     * If parameter \c include_bounds is \c false: the size of matrix A is [getEqualityDimension() + getInequalityDimension() x parameterDimension()].
     *
     * @remarks this version defines the sparse matrix with long long indices (since it is required by osqp for example)
     *
     * @param[out]  A                          Combined jacobian (see description above)
     * @param[in]   include_finite_bounds      Specify whether bounds should be included in A (check dimension of A as descirbed above)
     * @param[in]   col_nnz                    Estimate for number of nnz per column to speed up insertion
     */
    virtual void computeSparseJacobianTwoSideBoundedLinearForm(Eigen::SparseMatrix<double, Eigen::ColMajor, long long>& A, bool include_finite_bounds,
                                                               const Eigen::VectorXi* col_nnz = nullptr);
    virtual void computeSparseJacobianTwoSideBoundedLinearFormNNZPerColumn(Eigen::Ref<Eigen::VectorXi> col_nnz, bool include_finite_bounds);

    virtual int computeSparseJacobianTwoSideBoundedLinearFormNNZ(bool include_finite_bounds);
    virtual void computeSparseJacobianTwoSideBoundedLinearFormStructure(Eigen::Ref<Eigen::VectorXi> i_row, Eigen::Ref<Eigen::VectorXi> j_col,
                                                                        bool include_finite_bounds);
    virtual void computeSparseJacobianTwoSideBoundedLinearFormValues(Eigen::Ref<Eigen::VectorXd> values, bool include_finite_bounds);

    /**
     * @brief Compute the Jacobian and Hessian of the lagrangian
     *
     * This method combines computeSparseJacobianTwoSideBoundedLinearForm()
     * and computeSparseHessianLagrangian() in order to speed up computation.
     * Often, Hessian computations rely on first order derivatives and so we do not need to compute the Jacobian twice.
     */
    virtual void computeSparseJacobianTwoSideBoundedLinearFormAndHessianLagrangian(
        Eigen::SparseMatrix<double, Eigen::ColMajor, long long>& H, const double* multipliers_eq, const double* multipliers_ineq,
        Eigen::SparseMatrix<double, Eigen::ColMajor, long long>& A, bool include_finite_bounds, const Eigen::VectorXi* col_nnz_H = nullptr,
        const Eigen::VectorXi* col_nnz_A = nullptr, bool upper_part_only_H = false);

    /**
     * @brief Compute lower and upper bounds lbA and ubA for the linear form lbA <= A x <= ubA
     *
     * Especially some QP solvers require either
     *
     *      lbA <= A x<= ubA
     *
     *  or
     *
     *     lbA <= A x <= ubA
     *     lb <= x <= ub
     *
     * The first formulatin includes bounds and hence the bottom rows of A are the identity.
     * However, in case lb and ub are infinity, A might have obsolete rows.
     *
     * If parameter \c include_bounds is \c true: only finite bounds are considered and
     * vectors lbA and ubA are of size [getEqualityDimension() + getInequalityDimension() + finiteCombinedBounds() x 1].
     *
     * If parameter \c include_bounds is \c false: the vectors are of size [getEqualityDimension() + getInequalityDimension() x 1].
     *
     *
     * @param[out]  A                                   Combined jacobian (see description above)
     * @param[in]   include_finite_bounds    Specify whether bounds should be included in A (check dimension of A as descirbed above)
     */
    virtual void computeBoundsForTwoSideBoundedLinearForm(Eigen::Ref<Eigen::VectorXd> lbA, Eigen::Ref<Eigen::VectorXd> ubA,
                                                          bool include_finite_bounds);

    //@}

    virtual void computeSparseJacobianTwoSideBoundedLinearFormAndHessianObjective(Eigen::SparseMatrix<double, Eigen::ColMajor, long long>& H,
                                                                                  Eigen::SparseMatrix<double, Eigen::ColMajor, long long>& A,
                                                                                  bool include_finite_bounds, const Eigen::VectorXi* col_nnz_H,
                                                                                  const Eigen::VectorXi* col_nnz_A, bool upper_part_only_H);

    //! Check if a function taking the parameter value and unfixed-idx is true for all unfixed parameter values.
    virtual bool checkIfAllUnfixedParam(std::function<bool(double, int)> fun);

 protected:
    bool _warn_if_not_specialized = true;
};

}  // namespace corbo

#endif  // SRC_OPTIMIZATION_INCLUDE_CORBO_OPTIMIZATION_OPTIMIZATION_PROBLEM_INTERFACE_H_
