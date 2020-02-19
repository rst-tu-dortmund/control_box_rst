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

#ifndef SRC_OPTIMIZATION_INCLUDE_CORBO_OPTIMIZATION_SIMPLE_OPTIMIZATION_PROBLEM_H_
#define SRC_OPTIMIZATION_INCLUDE_CORBO_OPTIMIZATION_SIMPLE_OPTIMIZATION_PROBLEM_H_

#include <corbo-optimization/optimization_problem_interface.h>

#include <functional>

namespace corbo {

/**
 * @brief Simple optimization problem formulation
 *
 * @ingroup optimization
 *
 * This class defines a standard optimization problem in which
 * the optimization vector x is represented as a simple vector.
 *
 * For the mathematical description refer to the OptimizationProblemInterface.
 *
 * @remarks This class is abstract and needs to be implemented in an appropriate subclass;
 *          alternatively, refer to SimpleOptimizationProblemWithCallbacks for direct usage.
 *
 * @see OptimizationProblemInterface SimpleOptimizationProblemWithCallbacks SolverInterface
 *
 * @author Christoph Rösmann (christoph.roesmann@tu-dortmund.de)
 */
class SimpleOptimizationProblem : public OptimizationProblemInterface
{
 public:
    SimpleOptimizationProblem() { _warn_if_not_specialized = false; }
    explicit SimpleOptimizationProblem(int parameter_dim)
    {
        _warn_if_not_specialized = false;
        resizeParameterVector(parameter_dim);
    }

    /** @name Specify the dimension of the optimization problem  */
    //@{

    //! Total dimension of objective function terms
    int getNonLsqObjectiveDimension() override = 0;
    int getLsqObjectiveDimension() override    = 0;
    int getObjectiveDimension() override       = 0;
    //! Total dimension of equality constraints
    int getEqualityDimension() override = 0;
    //! Total dimension of general inequality constraints
    int getInequalityDimension() override = 0;

    //@}

    /** @name Specify main equations of the optimization problem  */
    //@{

    // implements interface method
    void computeValuesLsqObjective(Eigen::Ref<Eigen::VectorXd> values) override = 0;

    // implements interface method
    double computeValueNonLsqObjective() override = 0;

    // implements interface method
    void computeValuesEquality(Eigen::Ref<Eigen::VectorXd> values) override = 0;

    // implements interface method
    void computeValuesInequality(Eigen::Ref<Eigen::VectorXd> values) override = 0;

    //@}

    /** @name Configure parameter vector  */
    //@{

    /**
     * @brief Resize the dimension of the parameter vector
     * @remarks Values are resetted to zero and bounds to +-inf.
     * @param parameter_dim
     */
    void resizeParameterVector(int parameter_dim);

    //@}

    /** @name Access parameter vector and bounds  */
    //@{

    void setX(const Eigen::Ref<const Eigen::VectorXd>& x) { _x = x; }
    Eigen::VectorXd& getXRef() { return _x; }
    const Eigen::VectorXd& getX() const { return _x; }

    void setLowerBounds(const Eigen::Ref<const Eigen::VectorXd>& lb) { _lb = lb; }
    Eigen::VectorXd& getLowerBoundsRef() { return _lb; }
    const Eigen::VectorXd& getLowerBounds() const { return _lb; }

    void setUpperBounds(const Eigen::Ref<const Eigen::VectorXd>& ub) { _ub = ub; }
    Eigen::VectorXd& getUpperBoundsRef() { return _ub; }
    const Eigen::VectorXd& getUpperBounds() const { return _ub; }

    // implements interface method
    double getParameterValue(int idx) override { return _x[idx]; }
    // implements interface method
    void setParameterValue(int idx, double x) override { _x[idx] = x; }

    // implements interface method
    void getBounds(Eigen::Ref<Eigen::VectorXd> lb, Eigen::Ref<Eigen::VectorXd> ub) override;
    // implements interface method
    void setBounds(const Eigen::Ref<const Eigen::VectorXd>& lb, const Eigen::Ref<const Eigen::VectorXd>& ub) override;
    // implements interface method
    double getLowerBound(int idx) override { return _lb[idx]; }
    // implements interface method
    double getUpperBound(int idx) override { return _ub[idx]; }
    // implements interface method
    void setLowerBound(int idx, double lb) override { _lb[idx] = lb; }
    // implements interface method
    void setUpperBound(int idx, double ub) override { _ub[idx] = ub; }

    //! Same as getX() but less efficient (overrides interface method, but here x must not be pre-allocated)
    void getParameterVector(Eigen::Ref<Eigen::VectorXd> x) override { x = _x; }
    //! Same as setX() (overrides interface method)
    void setParameterVector(const Eigen::Ref<const Eigen::VectorXd>& x) override { setX(x); }

    //@}

    /** @name Interface implementations  */
    //@{

    // implements interface method
    int getParameterDimension() override { return _x.size(); }

    // implements interface method
    void applyIncrement(const Eigen::Ref<const Eigen::VectorXd>& increment) override { _x += increment; }
    // implements interface method
    void applyIncrement(int idx, double increment) override { _x[idx] += increment; }

    // implements interface method
    void backupParameters() override { _x_backup.push_back(_x); }
    // implements interface method
    void restoreBackupParameters(bool keep_backup) override;
    // implements interface method
    void discardBackupParameters(bool all = false) override
    {
        _x_backup.pop_back();
        if (all)
        {
            while (!_x_backup.empty()) _x_backup.pop_back();
        }
    }

    //@}

    /** @name Specify properties of the optimization problem  */
    //@{

    // implements interface method
    bool isLeastSquaresProblem() const override = 0;

    //@}

 private:
    Eigen::VectorXd _x;
    Eigen::VectorXd _lb;
    Eigen::VectorXd _ub;

    std::vector<Eigen::VectorXd> _x_backup;
};

/**
 * @brief Simple optimization problem formulation (callback based configuration)
 *
 * @ingroup optimization
 *
 * This class defines a standard optimization problem in which
 * the optimization vector x is represented as a simple vector.
 * In comparison to SimpleOptimizationProblem, the objective function and constraints
 * are integrated via callbacks rather then implementing interface methods.
 *
 * For the mathematical description refer to the OptimizationProblemInterface.
 *
 * @see SimpleOptimizationProblem OptimizationProblemInterface SolverInterface
 *
 * @author Christoph Rösmann (christoph.roesmann@tu-dortmund.de)
 */
class SimpleOptimizationProblemWithCallbacks : public SimpleOptimizationProblem
{
 public:
    //! Default constructor (do not forget to initialize the parameter vector dimension manually)
    SimpleOptimizationProblemWithCallbacks() {}
    //! Construct Optimization Problem with a given parameter vector dimension
    explicit SimpleOptimizationProblemWithCallbacks(int param_dim) : SimpleOptimizationProblem(param_dim) {}

    /** @name Set callbacks for the objective function and constraints  */
    //@{

    /**
     * @brief Set objective function callback
     *
     * The prototype is as follows:
     *
     *     void objective(const Eigen::VectorXd& x, Eigen::Ref<Eigen::VectorXd> values);
     *
     * in which \c x contains the current parameter values [parameterDimension() x 1] and \c values should be
     * filled with the corresponding objective functin values (output argument).
     * Please note, \c values is assumed to be pre-allocated with size [obj_dim x 1].
     *
     * Argument \c lsq_form indicates that the solver should take the squared l2-norm of the objective function vector:
     * e.g., \f$ \mathbf{f}^T \mathbf{f} \f$. Hereby, \f$ \mathbf{f} \f$ denotes the objective function value vector with size [obj_dim x 1].
     *
     * @param[in] obj_fun   Set callback for the objective function computation
     * @param[in] obj_dim   Dimension of the objective function vector
     * @param[in] lsq_form  If true, the solver takes the squared l2-norm of the obejctive function value.
     */
    void setObjectiveFunction(std::function<void(const Eigen::VectorXd&, Eigen::Ref<Eigen::VectorXd>)> obj_fun, int obj_dim, bool lsq_form = false);

    /**
     * @brief Set equality constraint callback
     *
     * The prototype is as follows:
     *
     *     void equality(const Eigen::VectorXd& x, Eigen::Ref<Eigen::VectorXd> values);
     *
     * in which \c x contains the current parameter values [parameterDimension() x 1] and \c values should be
     * filled with the corresponding equality constraint values (output argument).
     * Please note, \c values is assumed to be pre-allocated with size [eq_dim x 1].
     *
     * @param[in] eq_fun   Set callback for the equality constraint computation
     * @param[in] eq_dim   Dimension of the equality constraint vector
     */
    void setEqualityConstraint(std::function<void(const Eigen::VectorXd&, Eigen::Ref<Eigen::VectorXd>)> eq_fun, int eq_dim);

    /**
     * @brief Set inequality constraint callback
     *
     * The prototype is as follows:
     *
     *     void inequality(const Eigen::VectorXd& x, Eigen::Ref<Eigen::VectorXd> values);
     *
     * in which \c x contains the current parameter values [parameterDimension() x 1] and \c values should be
     * filled with the corresponding inequality constraint values (output argument).
     * Please note, \c values is assumed to be pre-allocated with size [ineq_dim x 1].
     *
     * @param[in] ineq_fun   Set callback for the inequality constraint computation
     * @param[in] ineq_dim   Dimension of the inequality constraint vector
     */
    void setInequalityConstraint(std::function<void(const Eigen::VectorXd&, Eigen::Ref<Eigen::VectorXd>)> ineq_fun, int ineq_dim);

    //@}

    // implements interface method
    int getNonLsqObjectiveDimension() override { return (!_lsq_form && _obj_dim > 0) ? 1 : 0; }
    // implements interface method
    int getLsqObjectiveDimension() override { return _lsq_form ? _obj_dim : 0; }
    // implements interface method
    int getObjectiveDimension() override { return _obj_dim > 0 ? 1 : 0; }
    // implements interface method
    int getEqualityDimension() override { return _eq_dim; }
    // implements interface method
    int getInequalityDimension() override { return _ineq_dim; }

    // implements interface method
    bool isLeastSquaresProblem() const override { return _lsq_form; }

    // implmements interface method
    double computeValueNonLsqObjective() override
    {
        if (!_obj_fun || _lsq_form) return 0;
        Eigen::VectorXd values(getNonLsqObjectiveDimension());
        _obj_fun(getX(), values);
        return values.sum();
    }

    // implements interface method
    void computeValuesLsqObjective(Eigen::Ref<Eigen::VectorXd> values) override
    {
        if (_obj_fun) _obj_fun(getX(), values);
    }

    // implements interface method
    void computeValuesEquality(Eigen::Ref<Eigen::VectorXd> values) override
    {
        if (_eq_fun) _eq_fun(getX(), values);
    }

    // implements interface method
    void computeValuesInequality(Eigen::Ref<Eigen::VectorXd> values) override
    {
        if (_ineq_fun) _ineq_fun(getX(), values);
    }

 private:
    int _obj_dim   = 0;
    int _eq_dim    = 0;
    int _ineq_dim  = 0;
    bool _lsq_form = false;

    std::function<void(const Eigen::VectorXd&, Eigen::Ref<Eigen::VectorXd>)> _obj_fun;
    std::function<void(const Eigen::VectorXd&, Eigen::Ref<Eigen::VectorXd>)> _eq_fun;
    std::function<void(const Eigen::VectorXd&, Eigen::Ref<Eigen::VectorXd>)> _ineq_fun;
};

}  // namespace corbo

#endif  // SRC_OPTIMIZATION_INCLUDE_CORBO_OPTIMIZATION_SIMPLE_OPTIMIZATION_PROBLEM_H_
