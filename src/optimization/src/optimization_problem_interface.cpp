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

#include <corbo-optimization/optimization_problem_interface.h>

#include <corbo-numerics/finite_differences.h>

namespace corbo {

double OptimizationProblemInterface::computeValueObjective()
{
    double value = computeValueNonLsqObjective();
    if (getLsqObjectiveDimension() > 0)
    {
        Eigen::VectorXd values(getLsqObjectiveDimension());
        computeValuesLsqObjective(values);
        value += values.squaredNorm();
    }
    return value;
}

void OptimizationProblemInterface::getParameterVector(Eigen::Ref<Eigen::VectorXd> x)
{
    PRINT_DEBUG_COND_ONCE(_warn_if_not_specialized, "OptimizationProblemInterface::getParameterVector(): default implementation might be slow.");
    assert(x.size() == getParameterDimension());
    for (int i = 0; i < getParameterDimension(); ++i) x[i] = getParameterValue(i);
}

void OptimizationProblemInterface::setParameterVector(const Eigen::Ref<const Eigen::VectorXd>& x)
{
    PRINT_DEBUG_COND_ONCE(_warn_if_not_specialized, "OptimizationProblemInterface::setParameterVector(): default implementation might be slow.");
    for (int i = 0; i < getParameterDimension(); ++i) setParameterValue(i, x[i]);
}

void OptimizationProblemInterface::applyIncrement(const Eigen::Ref<const Eigen::VectorXd>& increment)
{
    PRINT_DEBUG_COND_ONCE(_warn_if_not_specialized, "OptimizationProblemInterface::applyIncrement(): default implementation might be slow.");
    assert(increment.size() == getParameterDimension());
    for (int i = 0; i < increment.size(); ++i) setParameterValue(i, getParameterValue(i) + increment[i]);
}

void OptimizationProblemInterface::applyIncrement(int idx, double increment)
{
    PRINT_DEBUG_COND_ONCE(_warn_if_not_specialized, "OptimizationProblemInterface::applyIncrement(): default implementation might be slow.");
    setParameterValue(idx, getParameterValue(idx) + increment);
}

void OptimizationProblemInterface::setBounds(const Eigen::Ref<const Eigen::VectorXd>& lb, const Eigen::Ref<const Eigen::VectorXd>& ub)
{
    PRINT_DEBUG_COND_ONCE(_warn_if_not_specialized, "OptimizationProblemInterface::setBounds(): default implementation might be slow.");
    assert(lb.size() == getParameterDimension());
    assert(ub.size() == getParameterDimension());
    for (int i = 0; i < getParameterDimension(); ++i)
    {
        setLowerBound(i, lb[i]);
        setUpperBound(i, ub[i]);
    }
}

void OptimizationProblemInterface::getBounds(Eigen::Ref<Eigen::VectorXd> lb, Eigen::Ref<Eigen::VectorXd> ub)
{
    PRINT_DEBUG_COND_ONCE(_warn_if_not_specialized,
                          "OptimizationProblemInterface::finiteCombinedBoundsDimension(): default implementation might be slow.");
    assert(lb.size() == getParameterDimension());
    assert(ub.size() == getParameterDimension());
    for (int i = 0; i < getParameterDimension(); ++i)
    {
        lb[i] = getLowerBound(i);
        ub[i] = getUpperBound(i);
    }
}

int OptimizationProblemInterface::finiteCombinedBoundsDimension()
{
    PRINT_WARNING_COND_ONCE(_warn_if_not_specialized,
                            "OptimizationProblemInterface::finiteCombinedBoundsDimension(): default implementation might be slow.");
    int dim = 0;
    for (int i = 0; i < getParameterDimension(); ++i)
    {
        if (getLowerBound(i) > -CORBO_INF_DBL || getUpperBound(i) < CORBO_INF_DBL) ++dim;
    }
    return dim;
}

int OptimizationProblemInterface::finiteBoundsDimension()
{
    PRINT_DEBUG_COND_ONCE(_warn_if_not_specialized, "OptimizationProblemInterface::finiteBoundsDimension(): default implementation might be slow.");
    int dim = 0;
    for (int i = 0; i < getParameterDimension(); ++i)
    {
        if (getLowerBound(i) > -CORBO_INF_DBL) ++dim;
        if (getUpperBound(i) < CORBO_INF_DBL) ++dim;
    }
    return dim;
}

void OptimizationProblemInterface::computeValuesActiveInequality(Eigen::Ref<Eigen::VectorXd> values, double weight)
{
    PRINT_DEBUG_COND_ONCE(_warn_if_not_specialized,
                          "OptimizationProblemInterface::computeValuesActiveInequality(): default implementation might be slow.");
    computeValuesInequality(values);
    for (int i = 0; i < values.size(); ++i)
    {
        if (values[i] < 0)
            values[i] = 0;
        else
            values[i] *= weight;
    }
}

void OptimizationProblemInterface::computeDistanceFiniteCombinedBounds(Eigen::Ref<Eigen::VectorXd> values)
{
    PRINT_DEBUG_COND_ONCE(_warn_if_not_specialized,
                          "OptimizationProblemInterface::computeDistanceFiniteCombinedBounds(): default implementation might be slow.");
    double lb, ub, x;
    int idx = 0;
    for (int i = 0; i < getParameterDimension(); ++i)
    {
        lb = getLowerBound(i);
        ub = getUpperBound(i);
        if (lb > -CORBO_INF_DBL || ub < CORBO_INF_DBL)
        {
            assert(idx < finiteBoundsDimension());
            x = getParameterValue(i);
            if (x < lb)
                values[idx] = lb - x;
            else if (x > ub)
                values[idx] = x - ub;
            else
                values[idx] = 0;
            ++idx;
        }
    }
}

void OptimizationProblemInterface::computeLowerAndUpperBoundDiff(Eigen::Ref<Eigen::VectorXd> lb_minus_x, Eigen::Ref<Eigen::VectorXd> ub_minus_x)
{
    getBounds(lb_minus_x, ub_minus_x);
    Eigen::VectorXd x(getParameterDimension());
    getParameterVector(x);
    lb_minus_x -= x;
    ub_minus_x -= x;
}

void OptimizationProblemInterface::getParametersAndBoundsFinite(Eigen::Ref<Eigen::VectorXd> lb_finite_bounds,
                                                                Eigen::Ref<Eigen::VectorXd> ub_finite_bounds,
                                                                Eigen::Ref<Eigen::VectorXd> x_finite_bounds)
{
    PRINT_DEBUG_COND_ONCE(_warn_if_not_specialized,
                          "OptimizationProblemInterface::computeDistanceFiniteCombinedBounds(): default implementation might be slow.");
    double lb, ub;
    int idx = 0;
    for (int i = 0; i < getParameterDimension(); ++i)
    {
        lb = getLowerBound(i);
        ub = getUpperBound(i);
        if (lb > -CORBO_INF_DBL || ub < CORBO_INF_DBL)
        {
            assert(idx < finiteBoundsDimension());
            lb_finite_bounds(idx) = lb;
            ub_finite_bounds(idx) = ub;
            x_finite_bounds(idx)  = getParameterValue(i);
            ++idx;
        }
    }
}

void OptimizationProblemInterface::computeGradientObjective(Eigen::Ref<Eigen::VectorXd> gradient)
{
    PRINT_DEBUG_COND_ONCE(_warn_if_not_specialized,
                          "OptimizationProblemInterface::computeGradientObjective(): default implementation might be slow.");
    if (getObjectiveDimension() < 1) return;
    int dim_x = getParameterDimension();
    assert(gradient.size() == dim_x);

    // TODO(roesmann) generic interface
    auto inc  = [this](int idx, double inc) { applyIncrement(idx, inc); };
    auto eval = [this](Eigen::Ref<Eigen::VectorXd> values) { values[0] = computeValueObjective(); };
    CentralDifferences::jacobian(inc, eval, gradient.transpose());
}

void OptimizationProblemInterface::computeGradientNonLsqObjective(Eigen::Ref<Eigen::VectorXd> gradient)
{
    PRINT_DEBUG_COND_ONCE(_warn_if_not_specialized,
                          "OptimizationProblemInterface::computeGradientNonLsqObjective(): default implementation might be slow.");
    if (getNonLsqObjectiveDimension() < 1) return;
    int dim_x = getParameterDimension();
    assert(gradient.size() == dim_x);

    // TODO(roesmann) generic interface
    auto inc  = [this](int idx, double inc) { applyIncrement(idx, inc); };
    auto eval = [this](Eigen::Ref<Eigen::VectorXd> values) { values[0] = computeValueNonLsqObjective(); };
    CentralDifferences::jacobian(inc, eval, gradient.transpose());
}

void OptimizationProblemInterface::computeDenseJacobianLsqObjective(Eigen::Ref<Eigen::MatrixXd> jacobian, const double* multipliers)
{
    PRINT_DEBUG_COND_ONCE(_warn_if_not_specialized,
                          "OptimizationProblemInterface::computeDenseJacobianLsqObjective(): default implementation might be slow.");
    int dim_obj = getLsqObjectiveDimension();
    if (dim_obj < 1) return;
    int dim_x = getParameterDimension();
    assert(jacobian.rows() == dim_obj);
    assert(jacobian.cols() == dim_x);

    // TODO(roesmann) generic interface
    auto inc  = [this](int idx, double inc) { applyIncrement(idx, inc); };
    auto eval = [this](Eigen::Ref<Eigen::VectorXd> values) { computeValuesLsqObjective(values); };
    CentralDifferences::jacobian(inc, eval, jacobian);

    if (multipliers) jacobian.array().colwise() *= Eigen::Map<const Eigen::ArrayXd>(multipliers, dim_obj);
}

int OptimizationProblemInterface::computeSparseJacobianLsqObjectiveNNZ()
{
    PRINT_DEBUG_COND_ONCE(_warn_if_not_specialized,
                          "OptimizationProblemInterface::computeSparseJacobianLsqObjectiveNNZ(): default implementation might be slow.");
    // the default implementation is not aware of any sparse structure, hence we assume worst case
    return getLsqObjectiveDimension() * getParameterDimension();
}
void OptimizationProblemInterface::computeSparseJacobianLsqObjectiveStructure(Eigen::Ref<Eigen::VectorXi> i_row, Eigen::Ref<Eigen::VectorXi> j_col)
{
    PRINT_DEBUG_COND_ONCE(_warn_if_not_specialized,
                          "OptimizationProblemInterface::computeSparseJacobianLsqObjectiveStructure(): default implementation might be slow.");
    if (getLsqObjectiveDimension() < 1) return;
    // worst-case implementation
    int nz_idx = 0;
    for (int i = 0; i < getLsqObjectiveDimension(); ++i)
    {
        for (int j = 0; j < getParameterDimension(); ++j)
        {
            i_row[nz_idx] = i;
            j_col[nz_idx] = j;
            ++nz_idx;
        }
    }
}
void OptimizationProblemInterface::computeSparseJacobianLsqObjectiveValues(Eigen::Ref<Eigen::VectorXd> values, const double* multipliers)
{
    PRINT_DEBUG_COND_ONCE(_warn_if_not_specialized,
                          "OptimizationProblemInterface::computeSparseJacobianLsqObjectiveValues(): default implementation might be slow.");
    if (getLsqObjectiveDimension() < 1) return;
    // worst-case implementation
    int dim_obj = getLsqObjectiveDimension();
    int dim_x   = getParameterDimension();
    Eigen::MatrixXd dense_jacob(dim_obj, dim_x);
    computeDenseJacobianLsqObjective(dense_jacob, multipliers);

    int nz_idx = 0;
    for (int i = 0; i < getLsqObjectiveDimension(); ++i)
    {
        for (int j = 0; j < getParameterDimension(); ++j)
        {
            values[nz_idx] = dense_jacob(i, j);
            ++nz_idx;
        }
    }
}

void OptimizationProblemInterface::computeSparseJacobianLsqObjective(Eigen::SparseMatrix<double>& jacobian, const double* multipliers)
{
    PRINT_DEBUG_COND_ONCE(_warn_if_not_specialized,
                          "OptimizationProblemInterface::computeSparseJacobianLsqObjective(): default implementation might be slow.");
    int dim_obj = getLsqObjectiveDimension();
    if (dim_obj < 1) return;
    int dim_x = getParameterDimension();
    Eigen::MatrixXd dense_jacob(dim_obj, dim_x);
    computeDenseJacobianLsqObjective(dense_jacob, multipliers);
    jacobian = dense_jacob.sparseView();
}

void OptimizationProblemInterface::computeDenseJacobianEqualities(Eigen::Ref<Eigen::MatrixXd> jacobian, const double* multipliers)
{
    PRINT_DEBUG_COND_ONCE(_warn_if_not_specialized,
                          "OptimizationProblemInterface::computeDenseJacobianEqualities(): default implementation might be slow.");
    int dim_eq = getEqualityDimension();
    int dim_x  = getParameterDimension();
    assert(jacobian.rows() == dim_eq);
    assert(jacobian.cols() == dim_x);

    // TODO(roesmann) generic interface
    CentralDifferences diff;
    auto inc  = [this](int idx, double inc) { applyIncrement(idx, inc); };
    auto eval = [this](Eigen::Ref<Eigen::VectorXd> values) { computeValuesEquality(values); };
    diff.computeJacobian(inc, eval, jacobian);

    if (multipliers) jacobian.array().colwise() *= Eigen::Map<const Eigen::ArrayXd>(multipliers, dim_eq);
}

int OptimizationProblemInterface::computeSparseJacobianEqualitiesNNZ()
{
    PRINT_DEBUG_COND_ONCE(_warn_if_not_specialized,
                          "OptimizationProblemInterface::computeSparseJacobianEqualitiesNNZ(): default implementation might be slow.");
    // the default implementation is not aware of any sparse structure, hence we assume worst case
    return getEqualityDimension() * getParameterDimension();
}
void OptimizationProblemInterface::computeSparseJacobianEqualitiesStructure(Eigen::Ref<Eigen::VectorXi> i_row, Eigen::Ref<Eigen::VectorXi> j_col)
{
    PRINT_DEBUG_COND_ONCE(_warn_if_not_specialized,
                          "OptimizationProblemInterface::computeSparseJacobianEqualitiesStructure(): default implementation might be slow.");
    // worst-case implementation
    int nz_idx = 0;
    for (int i = 0; i < getEqualityDimension(); ++i)
    {
        for (int j = 0; j < getParameterDimension(); ++j)
        {
            i_row[nz_idx] = i;
            j_col[nz_idx] = j;
            ++nz_idx;
        }
    }
}
void OptimizationProblemInterface::computeSparseJacobianEqualitiesValues(Eigen::Ref<Eigen::VectorXd> values, const double* multipliers)
{
    PRINT_DEBUG_COND_ONCE(_warn_if_not_specialized,
                          "OptimizationProblemInterface::computeSparseJacobianEqualitiesValues(): default implementation might be slow.");
    // worst-case implementation
    int dim_eq = getEqualityDimension();
    int dim_x  = getParameterDimension();
    Eigen::MatrixXd dense_jacob(dim_eq, dim_x);
    computeDenseJacobianEqualities(dense_jacob, multipliers);

    int nz_idx = 0;
    for (int i = 0; i < getEqualityDimension(); ++i)
    {
        for (int j = 0; j < getParameterDimension(); ++j)
        {
            values[nz_idx] = dense_jacob(i, j);
            ++nz_idx;
        }
    }
}

void OptimizationProblemInterface::computeSparseJacobianEqualities(Eigen::SparseMatrix<double>& jacobian, const double* multipliers)
{
    PRINT_DEBUG_COND_ONCE(_warn_if_not_specialized,
                          "OptimizationProblemInterface::computeSparseJacobianEqualities(): default implementation might be slow.");
    int dim_eq = getEqualityDimension();
    int dim_x  = getParameterDimension();
    Eigen::MatrixXd dense_jacob(dim_eq, dim_x);
    computeDenseJacobianEqualities(dense_jacob, multipliers);
    jacobian = dense_jacob.sparseView();
}

void OptimizationProblemInterface::computeDenseJacobianInequalities(Eigen::Ref<Eigen::MatrixXd> jacobian, const double* multipliers)
{
    PRINT_DEBUG_COND_ONCE(_warn_if_not_specialized,
                          "OptimizationProblemInterface::computeDenseJacobianInequalities(): default implementation might be slow.");
    int dim_ineq = getInequalityDimension();
    int dim_x    = getParameterDimension();
    assert(jacobian.rows() == dim_ineq);
    assert(jacobian.cols() == dim_x);

    // TODO(roesmann) generic interface
    CentralDifferences diff;
    auto inc  = [this](int idx, double inc) { applyIncrement(idx, inc); };
    auto eval = [this](Eigen::Ref<Eigen::VectorXd> values) { computeValuesInequality(values); };
    diff.computeJacobian(inc, eval, jacobian);

    if (multipliers) jacobian.array().colwise() *= Eigen::Map<const Eigen::ArrayXd>(multipliers, dim_ineq);
}

int OptimizationProblemInterface::computeSparseJacobianInequalitiesNNZ()
{
    PRINT_DEBUG_COND_ONCE(_warn_if_not_specialized,
                          "OptimizationProblemInterface::computeSparseJacobianInequalitiesNNZ(): default implementation might be slow.");
    // the default implementation is not aware of any sparse structure, hence we assume worst case
    return getInequalityDimension() * getParameterDimension();
}
void OptimizationProblemInterface::computeSparseJacobianInequalitiesStructure(Eigen::Ref<Eigen::VectorXi> i_row, Eigen::Ref<Eigen::VectorXi> j_col)
{
    PRINT_DEBUG_COND_ONCE(_warn_if_not_specialized,
                          "OptimizationProblemInterface::computeSparseJacobianInequalitiesStructure(): default implementation might be slow.");
    // worst-case implementation
    int nz_idx = 0;
    for (int i = 0; i < getInequalityDimension(); ++i)
    {
        for (int j = 0; j < getParameterDimension(); ++j)
        {
            i_row[nz_idx] = i;
            j_col[nz_idx] = j;
            ++nz_idx;
        }
    }
}
void OptimizationProblemInterface::computeSparseJacobianInequalitiesValues(Eigen::Ref<Eigen::VectorXd> values, const double* multipliers)
{
    PRINT_DEBUG_COND_ONCE(_warn_if_not_specialized,
                          "OptimizationProblemInterface::computeSparseJacobianInequalitiesValues(): default implementation might be slow.");
    // worst-case implementation
    int dim_ineq = getInequalityDimension();
    int dim_x    = getParameterDimension();
    Eigen::MatrixXd dense_jacob(dim_ineq, dim_x);
    computeDenseJacobianInequalities(dense_jacob, multipliers);

    int nz_idx = 0;
    for (int i = 0; i < getInequalityDimension(); ++i)
    {
        for (int j = 0; j < getParameterDimension(); ++j)
        {
            values[nz_idx] = dense_jacob(i, j);
            ++nz_idx;
        }
    }
}

void OptimizationProblemInterface::computeSparseJacobianInequalities(Eigen::SparseMatrix<double>& jacobian, const double* multipliers)
{
    PRINT_DEBUG_COND_ONCE(_warn_if_not_specialized,
                          "OptimizationProblemInterface::computeSparseJacobianInequalities(): default implementation might be slow.");
    int dim_ineq = getInequalityDimension();
    int dim_x    = getParameterDimension();
    Eigen::MatrixXd dense_jacob(dim_ineq, dim_x);
    computeDenseJacobianInequalities(dense_jacob, multipliers);
    jacobian = dense_jacob.sparseView();
}

void OptimizationProblemInterface::computeDenseJacobianActiveInequalities(Eigen::Ref<Eigen::MatrixXd> jacobian, double weight)
{
    PRINT_DEBUG_COND_ONCE(_warn_if_not_specialized,
                          "OptimizationProblemInterface::computeDenseJacobianActiveInequalities(): default implementation might be slow.");
    int dim_ineq = getInequalityDimension();
    int dim_x    = getParameterDimension();
    assert(jacobian.rows() == dim_ineq);
    assert(jacobian.cols() == dim_x);

    // TODO(roesmann) generic interface
    CentralDifferences diff;
    auto inc  = [this](int idx, double inc) { applyIncrement(idx, inc); };
    auto eval = [this](Eigen::Ref<Eigen::VectorXd> values) { computeValuesActiveInequality(values); };
    // TODO(roesmann) inefficient! We can better iterate all values, check if active, and compute n x 1 subjacobians each
    diff.computeJacobian(inc, eval, jacobian);

    if (weight != 1) jacobian *= weight;
}
void OptimizationProblemInterface::computeSparseJacobianActiveInequalitiesValues(Eigen::Ref<Eigen::VectorXd> values, double weight)
{
    PRINT_DEBUG_COND_ONCE(_warn_if_not_specialized,
                          "OptimizationProblemInterface::computeSparseJacobianActiveInequalitiesValues(): default implementation might be slow.");
    // worst-case implementation
    int dim_ineq = getInequalityDimension();
    int dim_x    = getParameterDimension();
    Eigen::MatrixXd dense_jacob(dim_ineq, dim_x);
    computeDenseJacobianActiveInequalities(dense_jacob, weight);

    int nz_idx = 0;
    for (int i = 0; i < getInequalityDimension(); ++i)
    {
        for (int j = 0; j < getParameterDimension(); ++j)
        {
            values[nz_idx] = dense_jacob(i, j);
            ++nz_idx;
        }
    }
}

void OptimizationProblemInterface::computeSparseJacobianActiveInequalities(Eigen::SparseMatrix<double>& jacobian, double weight)
{
    PRINT_DEBUG_COND_ONCE(_warn_if_not_specialized,
                          "OptimizationProblemInterface::computeSparseJacobianActiveInequalities(): default implementation might be slow.");
    int dim_ineq = getInequalityDimension();
    int dim_x    = getParameterDimension();
    Eigen::MatrixXd dense_jacob(dim_ineq, dim_x);
    computeDenseJacobianActiveInequalities(dense_jacob, weight);
    jacobian = dense_jacob.sparseView();
}

void OptimizationProblemInterface::computeDenseJacobians(Eigen::Ref<Eigen::VectorXd> gradient_non_lsq_obj,
                                                         Eigen::Ref<Eigen::MatrixXd> jacobian_lsq_obj, Eigen::Ref<Eigen::MatrixXd> jacobian_eq,
                                                         Eigen::Ref<Eigen::MatrixXd> jacobian_ineq, const double* multipliers_lsq_obj,
                                                         const double* multipliers_eq, const double* multipliers_ineq, bool active_ineq,
                                                         double active_ineq_weight)
{
    PRINT_DEBUG_COND_ONCE(_warn_if_not_specialized, "OptimizationProblemInterface::computeDenseJacobians(): default implementation might be slow.");
    int non_lsq_obj_dim = getNonLsqObjectiveDimension();
    int lsq_obj_dim     = getLsqObjectiveDimension();
    int eq_dim          = getEqualityDimension();
    int ineq_dim        = getInequalityDimension();
    if (non_lsq_obj_dim)
    {
        computeGradientNonLsqObjective(gradient_non_lsq_obj);
    }
    if (lsq_obj_dim > 0)
    {
        computeDenseJacobianLsqObjective(jacobian_lsq_obj, multipliers_lsq_obj);
    }
    if (eq_dim > 0)
    {
        computeDenseJacobianEqualities(jacobian_eq, multipliers_eq);
    }
    if (ineq_dim > 0)
    {
        if (active_ineq)
            computeDenseJacobianActiveInequalities(jacobian_ineq, active_ineq_weight);
        else
            computeDenseJacobianInequalities(jacobian_ineq, multipliers_ineq);
    }
}

void OptimizationProblemInterface::computeSparseJacobians(Eigen::SparseMatrix<double>& jacobian_lsq_obj, Eigen::SparseMatrix<double>& jacobian_eq,
                                                          Eigen::SparseMatrix<double>& jacobian_ineq, const double* multipliers_lsq_obj,
                                                          const double* multipliers_eq, const double* multipliers_ineq, bool active_ineq,
                                                          double active_ineq_weight)
{
    PRINT_DEBUG_COND_ONCE(_warn_if_not_specialized, "OptimizationProblemInterface::computeSparseJacobians(): default implementation might be slow.");
    computeSparseJacobianLsqObjective(jacobian_lsq_obj, multipliers_lsq_obj);
    computeSparseJacobianEqualities(jacobian_eq, multipliers_eq);
    if (active_ineq)
        computeSparseJacobianActiveInequalities(jacobian_ineq, active_ineq_weight);
    else
        computeSparseJacobianInequalities(jacobian_ineq, multipliers_ineq);
}

void OptimizationProblemInterface::computeSparseJacobiansNNZ(int& nnz_lsq_obj, int& nnz_eq, int& nnz_ineq)
{
    nnz_lsq_obj = computeSparseJacobianLsqObjectiveNNZ();
    nnz_eq      = computeSparseJacobianEqualitiesNNZ();
    nnz_ineq    = computeSparseJacobianInequalitiesNNZ();
}

void OptimizationProblemInterface::computeSparseJacobiansStructure(Eigen::Ref<Eigen::VectorXi> i_row_obj, Eigen::Ref<Eigen::VectorXi> j_col_obj,
                                                                   Eigen::Ref<Eigen::VectorXi> i_row_eq, Eigen::Ref<Eigen::VectorXi> j_col_eq,
                                                                   Eigen::Ref<Eigen::VectorXi> i_row_ineq, Eigen::Ref<Eigen::VectorXi> j_col_ineq)
{
    computeSparseJacobianLsqObjectiveStructure(i_row_obj, j_col_obj);
    computeSparseJacobianEqualitiesStructure(i_row_eq, j_col_eq);
    computeSparseJacobianInequalitiesStructure(i_row_ineq, j_col_ineq);
}

void OptimizationProblemInterface::computeSparseJacobiansValues(Eigen::Ref<Eigen::VectorXd> values_obj, Eigen::Ref<Eigen::VectorXd> values_eq,
                                                                Eigen::Ref<Eigen::VectorXd> values_ineq, const double* multipliers_obj,
                                                                const double* multipliers_eq, const double* multipliers_ineq, bool active_ineq,
                                                                double active_ineq_weight)
{
    PRINT_DEBUG_COND_ONCE(_warn_if_not_specialized,
                          "OptimizationProblemInterface::computeSparseJacobiansValues(): default implementation might be slow.");
    computeSparseJacobianLsqObjectiveValues(values_obj, multipliers_obj);
    computeSparseJacobianEqualitiesValues(values_eq, multipliers_eq);
    if (active_ineq)
        computeSparseJacobianActiveInequalitiesValues(values_ineq, active_ineq_weight);
    else
        computeSparseJacobianInequalitiesValues(values_ineq, multipliers_ineq);
}

void OptimizationProblemInterface::computeDenseJacobianFiniteCombinedBounds(Eigen::Ref<Eigen::MatrixXd> jacobian, double weight)
{
    PRINT_DEBUG_COND_ONCE(_warn_if_not_specialized,
                          "OptimizationProblemInterface::computeDenseJacobianFiniteCombinedBounds(): default implementation might be slow.");
    jacobian.setZero();

    int row_idx = 0;
    for (int i = 0; i < getParameterDimension(); ++i)
    {
        // check if finite
        double lb = getLowerBound(i);
        double ub = getUpperBound(i);
        if (lb > -CORBO_INF_DBL || ub < CORBO_INF_DBL)
        {
            double x = getParameterValue(i);
            if (x < lb)
            {
                jacobian(row_idx, i) = -weight;
            }
            else if (x > ub)
            {
                jacobian(row_idx, i) = weight;
            }
            ++row_idx;  // increase row index in jacobian
        }
    }
}

int OptimizationProblemInterface::computeSparseJacobianFiniteCombinedBoundsNNZ()
{
    // we have a single value per row
    return finiteCombinedBoundsDimension();
}

void OptimizationProblemInterface::computeSparseJacobianFiniteCombinedBoundsStructure(Eigen::Ref<Eigen::VectorXi> i_row,
                                                                                      Eigen::Ref<Eigen::VectorXi> j_col)
{
    PRINT_DEBUG_COND_ONCE(
        _warn_if_not_specialized,
        "OptimizationProblemInterface::computeSparseJacobianFiniteCombinedBoundsStructure(): default implementation might be slow.");
    // we have a single value per row
    int jac_row_idx = 0;
    for (int i = 0; i < getParameterDimension(); ++i)
    {
        // check if finite
        if (getLowerBound(i) > -CORBO_INF_DBL || getUpperBound(i) < CORBO_INF_DBL)
        {
            i_row[jac_row_idx] = jac_row_idx;
            j_col[jac_row_idx] = i;
            ++jac_row_idx;  // increase row index in jacobian
        }
    }
}
void OptimizationProblemInterface::computeSparseJacobianFiniteCombinedBoundsValues(Eigen::Ref<Eigen::VectorXd> values, double weight)
{
    int jac_row_idx = 0;
    for (int i = 0; i < getParameterDimension(); ++i)
    {
        // check if finite
        double lb = getLowerBound(i);
        double ub = getUpperBound(i);
        if (lb > -CORBO_INF_DBL || ub < CORBO_INF_DBL)
        {
            double x = getParameterValue(i);
            if (x < lb)
            {
                values[jac_row_idx] = -weight;
            }
            else if (x > ub)
            {
                values[jac_row_idx] = weight;
            }
            else
            {
                values[jac_row_idx] = 0.0;
            }
            ++jac_row_idx;  // increase row index in jacobian
        }
    }
}

void OptimizationProblemInterface::computeSparseJacobianFiniteCombinedBounds(Eigen::SparseMatrix<double>& jacobian, double weight)
{
    PRINT_DEBUG_COND_ONCE(_warn_if_not_specialized,
                          "OptimizationProblemInterface::computeSparseJacobianFiniteCombinedBounds(): default implementation might be slow.");
    // jacobian.reserve(finiteCombinedBoundsDimension());
    jacobian.setZero();

    // we have a single value per row
    int jac_row_idx = 0;
    for (int i = 0; i < getParameterDimension(); ++i)
    {
        // check if finite
        double lb = getLowerBound(i);
        double ub = getUpperBound(i);
        if (lb > -CORBO_INF_DBL || ub < CORBO_INF_DBL)
        {
            double x = getParameterValue(i);
            if (x < lb)
            {
                jacobian.insert(jac_row_idx, i) = -weight;
            }
            else if (x > ub)
            {
                jacobian.insert(jac_row_idx, i) = weight;
            }
            // else
            //{
            //    jacobian.insert(jac_row_idx, i) = 0.0;
            //}
            ++jac_row_idx;  // increase row index in jacobian
        }
    }
    // jacobian.makeCompressed();
}

void OptimizationProblemInterface::computeDenseJacobianFiniteCombinedBoundsIdentity(Eigen::Ref<Eigen::MatrixXd> jacobian)
{
    PRINT_DEBUG_COND_ONCE(_warn_if_not_specialized,
                          "OptimizationProblemInterface::computeDenseJacobianFiniteCombinedBoundsIdentity(): default implementation might be slow.");
    jacobian.setZero();

    int row_idx = 0;
    for (int i = 0; i < getParameterDimension(); ++i)
    {
        // check if finite
        double lb = getLowerBound(i);
        double ub = getUpperBound(i);
        if (lb > -CORBO_INF_DBL || ub < CORBO_INF_DBL)
        {
            jacobian(row_idx, i) = 1;
            ++row_idx;  // increase row index in jacobian
        }
    }
}

void OptimizationProblemInterface::computeCombinedSparseJacobian(Eigen::SparseMatrix<double>& jacobian, bool objective_lsq, bool equality,
                                                                 bool inequality, bool finite_combined_bounds, bool active_ineq, double weight_eq,
                                                                 double weight_ineq, double weight_bounds, const Eigen::VectorXd* values,
                                                                 const Eigen::VectorXi* col_nnz)
{
    PRINT_DEBUG_COND_ONCE(_warn_if_not_specialized,
                          "OptimizationProblemInterface::computeCombinedSparseJacobian(): default implementation might be slow.");
    int no_cols = getParameterDimension();
    int no_rows = 0;
    // int nnz = 0;
    int dim_obj    = 0;
    int dim_eq     = 0;
    int dim_ineq   = 0;
    int dim_bounds = 0;
    if (objective_lsq)
    {
        dim_obj = getLsqObjectiveDimension();
        no_rows += dim_obj;
        // nnz += computeSparseJacobianObjectiveNNZ();
    }
    if (equality)
    {
        dim_eq = getEqualityDimension();
        no_rows += dim_eq;
        // nnz += computeSparseJacobianEqualitiesNNZ();
    }
    if (inequality)
    {
        dim_ineq = getInequalityDimension();
        no_rows += dim_ineq;
        // nnz += computeSparseJacobianInequalitiesNNZ();
    }
    if (finite_combined_bounds)
    {
        dim_bounds = finiteCombinedBoundsDimension();
        no_rows += dim_bounds;
        // nnz += computeSparseJacobianFiniteCombinedBoundsNNZ();
    }

    assert(jacobian.rows() == no_rows);
    assert(jacobian.cols() == no_cols);

    jacobian.setZero();
    // jacobian.reserve(nnz);

    // default implementation, just concatenate submatrices
    // this is inefficient due to the temporary matrices
    // TODO(roesmann) maybe we should use row-major format for vertical concatenation?
    // e.g., https://stackoverflow.com/questions/41756428/concatenate-sparse-matrix-eigen

    int row_offset = 0;

    if (objective_lsq && dim_obj > 0)
    {
        Eigen::SparseMatrix<double> temp_jacob(dim_obj, no_cols);
        computeSparseJacobianLsqObjective(temp_jacob);
        for (int k = 0; k < temp_jacob.outerSize(); ++k)
        {
            for (Eigen::SparseMatrix<double>::InnerIterator it(temp_jacob, k); it; ++it)
            {
                jacobian.insert(it.row(), it.col()) = it.value();
            }
        }
        row_offset += dim_obj;
    }

    if (equality && dim_eq > 0)
    {
        Eigen::SparseMatrix<double> temp_jacob(dim_eq, no_cols);
        computeSparseJacobianEqualities(temp_jacob);
        for (int k = 0; k < temp_jacob.outerSize(); ++k)
        {
            for (Eigen::SparseMatrix<double>::InnerIterator it(temp_jacob, k); it; ++it)
            {
                jacobian.insert(row_offset + it.row(), it.col()) = it.value() * weight_eq;
            }
        }
        row_offset += dim_eq;
    }

    if (inequality && dim_ineq > 0)
    {
        Eigen::SparseMatrix<double> temp_jacob(dim_ineq, no_cols);
        if (active_ineq)
            computeSparseJacobianActiveInequalities(temp_jacob);
        else
            computeSparseJacobianInequalities(temp_jacob);
        for (int k = 0; k < temp_jacob.outerSize(); ++k)
        {
            for (Eigen::SparseMatrix<double>::InnerIterator it(temp_jacob, k); it; ++it)
            {
                jacobian.insert(row_offset + it.row(), it.col()) = it.value() * weight_ineq;
            }
        }
        row_offset += dim_ineq;
    }

    if (finite_combined_bounds && dim_bounds > 0)
    {
        Eigen::SparseMatrix<double> temp_jacob(dim_bounds, no_cols);
        computeSparseJacobianFiniteCombinedBounds(temp_jacob);
        for (int k = 0; k < temp_jacob.outerSize(); ++k)
        {
            for (Eigen::SparseMatrix<double>::InnerIterator it(temp_jacob, k); it; ++it)
            {
                jacobian.insert(row_offset + it.row(), it.col()) = it.value() * weight_bounds;
            }
        }
    }
}

int OptimizationProblemInterface::computeCombinedSparseJacobiansNNZ(bool objective_lsq, bool equality, bool inequality)
{
    PRINT_DEBUG_COND_ONCE(_warn_if_not_specialized,
                          "OptimizationProblemInterface::computeCombinedSparseJacobiansNNZ(): default implementation might be slow.");
    int nnz = 0;
    // worst case nnz
    if (objective_lsq) nnz += getLsqObjectiveDimension() * getParameterDimension();
    if (equality) nnz += getEqualityDimension() * getParameterDimension();
    if (inequality) nnz += getInequalityDimension() * getParameterDimension();
    return nnz;
}
void OptimizationProblemInterface::computeCombinedSparseJacobiansStructure(Eigen::Ref<Eigen::VectorXi> i_row, Eigen::Ref<Eigen::VectorXi> j_col,
                                                                           bool objective_lsq, bool equality, bool inequality)
{
    PRINT_DEBUG_COND_ONCE(_warn_if_not_specialized,
                          "OptimizationProblemInterface::computeCombinedSparseJacobiansStructure(): default implementation might be slow.");
    // worst-case implementation
    int nz_idx  = 0;
    int row_idx = 0;
    if (objective_lsq)
    {
        for (int i = 0; i < getLsqObjectiveDimension(); ++i, ++row_idx)
        {
            for (int j = 0; j < getParameterDimension(); ++j)
            {
                i_row[nz_idx] = row_idx;
                j_col[nz_idx] = j;
                ++nz_idx;
            }
        }
    }

    if (equality)
    {
        for (int i = 0; i < getEqualityDimension(); ++i, ++row_idx)
        {
            for (int j = 0; j < getParameterDimension(); ++j)
            {
                i_row[nz_idx] = row_idx;
                j_col[nz_idx] = j;
                ++nz_idx;
            }
        }
    }

    if (inequality)
    {
        for (int i = 0; i < getInequalityDimension(); ++i, ++row_idx)
        {
            for (int j = 0; j < getParameterDimension(); ++j)
            {
                i_row[nz_idx] = row_idx;
                j_col[nz_idx] = j;
                ++nz_idx;
            }
        }
    }
}

void OptimizationProblemInterface::computeCombinedSparseJacobiansValues(Eigen::Ref<Eigen::VectorXd> values, bool objective_lsq, bool equality,
                                                                        bool inequality, const double* multipliers_obj, const double* multipliers_eq,
                                                                        const double* multipliers_ineq)
{
    PRINT_DEBUG_COND_ONCE(_warn_if_not_specialized,
                          "OptimizationProblemInterface::computeCombinedSparseJacobiansValues(): default implementation might be slow.");
    int dim_x  = getParameterDimension();
    int nz_idx = 0;

    // worst-case implementation

    if (objective_lsq)
    {
        int dim_obj = getLsqObjectiveDimension();
        Eigen::MatrixXd dense_jacob(dim_obj, dim_x);
        computeDenseJacobianLsqObjective(dense_jacob, multipliers_obj);

        for (int i = 0; i < dim_obj; ++i)
        {
            for (int j = 0; j < getParameterDimension(); ++j)
            {
                values[nz_idx] = dense_jacob(i, j);
                ++nz_idx;
            }
        }
    }

    if (equality)
    {
        int dim_eq = getEqualityDimension();
        Eigen::MatrixXd dense_jacob(dim_eq, dim_x);
        computeDenseJacobianEqualities(dense_jacob, multipliers_eq);

        for (int i = 0; i < dim_eq; ++i)
        {
            for (int j = 0; j < getParameterDimension(); ++j)
            {
                values[nz_idx] = dense_jacob(i, j);
                ++nz_idx;
            }
        }
    }

    if (inequality)
    {
        int dim_ineq = getInequalityDimension();
        Eigen::MatrixXd dense_jacob(dim_ineq, dim_x);
        computeDenseJacobianInequalities(dense_jacob, multipliers_ineq);

        for (int i = 0; i < dim_ineq; ++i)
        {
            for (int j = 0; j < getParameterDimension(); ++j)
            {
                values[nz_idx] = dense_jacob(i, j);
                ++nz_idx;
            }
        }
    }
}

void OptimizationProblemInterface::computeGradientObjectiveAndCombinedSparseJacobiansValues(Eigen::Ref<Eigen::VectorXd> gradient,
                                                                                            Eigen::Ref<Eigen::VectorXd> jac_values, bool equality,
                                                                                            bool inequality, const double* multipliers_eq,
                                                                                            const double* multipliers_ineq)
{
    PRINT_DEBUG_COND_ONCE(_warn_if_not_specialized,
                          "OptimizationProblemInterface::computeDenseHessianObjective(): default implementation might be slow.");
    computeGradientObjective(gradient);
    computeCombinedSparseJacobiansValues(jac_values, false, equality, inequality, nullptr, multipliers_eq, multipliers_ineq);
}

void OptimizationProblemInterface::computeSparseJacobianTwoSideBoundedLinearForm(Eigen::SparseMatrix<double, Eigen::ColMajor, long long>& A,
                                                                                 bool include_finite_bounds, const Eigen::VectorXi* col_nnz)
{
    PRINT_DEBUG_COND_ONCE(_warn_if_not_specialized,
                          "OptimizationProblemInterface::computeCombinedSparseJacobian(): default implementation might be slow.");

    int no_cols = getParameterDimension();
    // int nnz = 0;
    int dim_eq     = getEqualityDimension();
    int dim_ineq   = getInequalityDimension();
    int dim_bounds = 0;

    int no_rows = dim_eq + dim_ineq;

    if (include_finite_bounds)
    {
        dim_bounds = finiteCombinedBoundsDimension();
        no_rows += dim_bounds;
        // nnz += computeSparseJacobianFiniteCombinedBoundsNNZ();
    }

    assert(A.rows() == no_rows);
    assert(A.cols() == no_cols);

    A.setZero();
    if (col_nnz) A.reserve(*col_nnz);
    // jacobian.reserve(nnz);

    // default implementation, just concatenate submatrices
    // this is inefficient due to the temporary matrices
    // TODO(roesmann) maybe we should use row-major format for vertical concatenation?
    // e.g., https://stackoverflow.com/questions/41756428/concatenate-sparse-matrix-eigen

    int row_offset = 0;

    if (dim_eq > 0)
    {
        Eigen::SparseMatrix<double> temp_jacob(dim_eq, no_cols);
        computeSparseJacobianEqualities(temp_jacob);
        for (int k = 0; k < temp_jacob.outerSize(); ++k)
        {
            for (Eigen::SparseMatrix<double>::InnerIterator it(temp_jacob, k); it; ++it)
            {
                A.insert(row_offset + it.row(), it.col()) = it.value();
            }
        }
        row_offset += dim_eq;
    }

    if (dim_ineq > 0)
    {
        Eigen::SparseMatrix<double> temp_jacob(dim_ineq, no_cols);
        computeSparseJacobianInequalities(temp_jacob);
        for (int k = 0; k < temp_jacob.outerSize(); ++k)
        {
            for (Eigen::SparseMatrix<double>::InnerIterator it(temp_jacob, k); it; ++it)
            {
                A.insert(row_offset + it.row(), it.col()) = it.value();
            }
        }
        row_offset += dim_ineq;
    }

    if (include_finite_bounds && dim_bounds > 0)
    {
        // TODO(roesmann): inefficient
        Eigen::MatrixXd temp_jacob_dense(dim_bounds, no_cols);
        computeDenseJacobianFiniteCombinedBoundsIdentity(temp_jacob_dense);
        Eigen::SparseMatrix<double> temp_jacob = temp_jacob_dense.sparseView();
        for (int k = 0; k < temp_jacob.outerSize(); ++k)
        {
            for (Eigen::SparseMatrix<double>::InnerIterator it(temp_jacob, k); it; ++it)
            {
                A.insert(row_offset + it.row(), it.col()) = it.value();
            }
        }
    }
}

void OptimizationProblemInterface::computeSparseJacobianTwoSideBoundedLinearFormNNZPerColumn(Eigen::Ref<Eigen::VectorXi> col_nnz,
                                                                                             bool include_finite_bounds)
{
    PRINT_DEBUG_COND_ONCE(
        _warn_if_not_specialized,
        "OptimizationProblemInterface::computeSparseJacobianTwoSideBoundedLinearFormNNZPerColumn(): default implementation might be slow.");

    assert(col_nnz.size() == getParameterDimension());

    // worst case nnz
    int dim_eq     = getEqualityDimension();
    int dim_ineq   = getInequalityDimension();
    int dim_bounds = include_finite_bounds ? finiteCombinedBoundsDimension() : 0;
    col_nnz.setConstant(dim_eq + dim_ineq + dim_bounds);  // full matrix
}

int OptimizationProblemInterface::computeSparseJacobianTwoSideBoundedLinearFormNNZ(bool include_finite_bounds)
{
    PRINT_DEBUG_COND_ONCE(_warn_if_not_specialized,
                          "OptimizationProblemInterface::computeSparseJacobianTwoSideBoundedLinearFormNNZ(): default implementation might be slow.");
    int nnz = 0;
    // worst case nnz
    nnz += getEqualityDimension() * getParameterDimension();
    nnz += getInequalityDimension() * getParameterDimension();
    if (include_finite_bounds) nnz += finiteCombinedBoundsDimension();
    return nnz;
}

void OptimizationProblemInterface::computeSparseJacobianTwoSideBoundedLinearFormStructure(Eigen::Ref<Eigen::VectorXi> i_row,
                                                                                          Eigen::Ref<Eigen::VectorXi> j_col,
                                                                                          bool include_finite_bounds)
{
    PRINT_DEBUG_COND_ONCE(
        _warn_if_not_specialized,
        "OptimizationProblemInterface::computeSparseJacobianTwoSideBoundedLinearFormStructure(): default implementation might be slow.");
    // worst-case implementation
    int nz_idx  = 0;
    int row_idx = 0;

    for (int i = 0; i < getEqualityDimension(); ++i, ++row_idx)
    {
        for (int j = 0; j < getParameterDimension(); ++j)
        {
            i_row[nz_idx] = row_idx;
            j_col[nz_idx] = j;
            ++nz_idx;
        }
    }

    for (int i = 0; i < getInequalityDimension(); ++i, ++row_idx)
    {
        for (int j = 0; j < getParameterDimension(); ++j)
        {
            i_row[nz_idx] = row_idx;
            j_col[nz_idx] = j;
            ++nz_idx;
        }
    }

    if (include_finite_bounds)
    {
        for (int i = 0; i < getParameterDimension(); ++i)
        {
            // check if finite
            double lb = getLowerBound(i);
            double ub = getUpperBound(i);
            if (lb > -CORBO_INF_DBL || ub < CORBO_INF_DBL)
            {
                i_row[nz_idx] = row_idx;
                j_col[nz_idx] = i;
                ++nz_idx;
                ++row_idx;  // increase row index in jacobian
            }
        }
    }
}

void OptimizationProblemInterface::computeSparseJacobianTwoSideBoundedLinearFormValues(Eigen::Ref<Eigen::VectorXd> values, bool include_finite_bounds)
{
    PRINT_DEBUG_COND_ONCE(_warn_if_not_specialized,
                          "OptimizationProblemInterface::computeCombinedSparseJacobiansValues(): default implementation might be slow.");
    int dim_x  = getParameterDimension();
    int nz_idx = 0;

    // worst-case implementation

    int dim_eq = getEqualityDimension();
    Eigen::MatrixXd dense_jacob_eq(dim_eq, dim_x);
    computeDenseJacobianEqualities(dense_jacob_eq);

    for (int i = 0; i < dim_eq; ++i)
    {
        for (int j = 0; j < getParameterDimension(); ++j)
        {
            values[nz_idx] = dense_jacob_eq(i, j);
            ++nz_idx;
        }
    }

    int dim_ineq = getInequalityDimension();
    Eigen::MatrixXd dense_jacob_ineq(dim_ineq, dim_x);
    computeDenseJacobianInequalities(dense_jacob_ineq);

    for (int i = 0; i < dim_ineq; ++i)
    {
        for (int j = 0; j < getParameterDimension(); ++j)
        {
            values[nz_idx] = dense_jacob_ineq(i, j);
            ++nz_idx;
        }
    }

    int dim_bounds = finiteCombinedBoundsDimension();
    if (include_finite_bounds && dim_bounds > 0)
    {
        values.tail(dim_bounds).setOnes();
    }
}

void OptimizationProblemInterface::computeBoundsForTwoSideBoundedLinearForm(Eigen::Ref<Eigen::VectorXd> lbA, Eigen::Ref<Eigen::VectorXd> ubA,
                                                                            bool include_finite_bounds)
{
    int dim_eq     = getEqualityDimension();
    int dim_ineq   = getInequalityDimension();
    int dim_bounds = include_finite_bounds ? finiteCombinedBoundsDimension() : 0;

    assert(lbA.size() == dim_eq + dim_ineq + dim_bounds);
    assert(ubA.size() == dim_eq + dim_ineq + dim_bounds);

    if (dim_eq > 0)
    {
        computeValuesEquality(lbA.head(dim_eq));
        lbA.head(dim_eq) *= -1;
        ubA.head(dim_eq) = lbA.head(dim_eq);
    }

    if (dim_ineq > 0)
    {
        lbA.segment(dim_eq, dim_ineq).setConstant(-CORBO_INF_DBL);  // should be smaller than -OSQP_INFTY
        computeValuesInequality(ubA.segment(dim_eq, dim_ineq));
        ubA.segment(dim_eq, dim_ineq) *= -1;
    }

    if (dim_bounds > 0)
    {
        double lb, ub, x;
        int idx = dim_eq + dim_ineq;
        for (int i = 0; i < getParameterDimension(); ++i)
        {
            lb = getLowerBound(i);
            ub = getUpperBound(i);
            if (lb > -CORBO_INF_DBL || ub < CORBO_INF_DBL)
            {
                assert(idx < finiteBoundsDimension());
                x        = getParameterValue(i);
                lbA(idx) = lb - x;
                ubA(idx) = x - ub;
                ++idx;
            }
        }
    }
}

void OptimizationProblemInterface::computeDenseHessianObjective(const Eigen::Ref<const Eigen::MatrixXd>& jacobian,
                                                                Eigen::Ref<Eigen::MatrixXd> hessian, const double* multipliers, bool jacob_scaled)
{
    PRINT_ERROR("OptimizationProblemInterface::computeDenseHessianObjective(): NOT_YET_IMPLEMENTED");
}

void OptimizationProblemInterface::computeDenseHessianObjective(Eigen::Ref<Eigen::MatrixXd> hessian, double multiplier)
{
    PRINT_DEBUG_COND_ONCE(_warn_if_not_specialized,
                          "OptimizationProblemInterface::computeDenseHessianObjective(): default implementation might be slow.");
    if (getObjectiveDimension() == 0) return;

    int dim_x = getParameterDimension();
    assert(hessian.rows() == dim_x);
    assert(hessian.cols() == dim_x);

    // TODO(roesmann) generic interface
    CentralDifferences diff;
    auto inc  = [this](int idx, double inc) { applyIncrement(idx, inc); };
    auto eval = [this](Eigen::Ref<Eigen::VectorXd> values) { values[0] = computeValueObjective(); };
    diff.computeHessian(inc, eval, 1, hessian, nullptr);
    if (multiplier != 1.0) hessian *= multiplier;
}

int OptimizationProblemInterface::computeSparseHessianObjectiveNNZ(bool lower_part_only)
{
    PRINT_DEBUG_COND_ONCE(_warn_if_not_specialized,
                          "OptimizationProblemInterface::computeSparseHessianObjectiveNNZ(): default implementation might be slow.");
    // the default implementation is not aware of any sparse structure, hence we assume worst case
    if (lower_part_only) return (int)(0.5 * getParameterDimension() * (getParameterDimension() + 1));
    return getParameterDimension() * getParameterDimension();
}
void OptimizationProblemInterface::computeSparseHessianObjectiveStructure(Eigen::Ref<Eigen::VectorXi> i_row, Eigen::Ref<Eigen::VectorXi> j_col,
                                                                          bool lower_part_only)
{
    PRINT_DEBUG_COND_ONCE(_warn_if_not_specialized,
                          "OptimizationProblemInterface::computeSparseHessianObjectiveStructure(): default implementation might be slow.");
    if (getObjectiveDimension() == 0) return;

    // worst-case implementation
    int nz_idx = 0;
    for (int col = 0; col < getParameterDimension(); ++col)
    {
        int row_start = lower_part_only ? col : 0;
        for (int row = row_start; row < getParameterDimension(); ++row)
        {
            i_row[nz_idx] = row;
            j_col[nz_idx] = col;
            ++nz_idx;
        }
    }
}
void OptimizationProblemInterface::computeSparseHessianObjectiveValues(Eigen::Ref<Eigen::VectorXd> values, double multiplier, bool lower_part_only)
{
    PRINT_DEBUG_COND_ONCE(_warn_if_not_specialized,
                          "OptimizationProblemInterface::computeSparseHessianObjectiveValues(): default implementation might be slow.");
    if (getObjectiveDimension() == 0) return;

    // worst-case implementation
    int dim_x = getParameterDimension();
    Eigen::MatrixXd dense_hessian(dim_x, dim_x);
    computeDenseHessianObjective(dense_hessian);

    int nz_idx = 0;
    for (int col = 0; col < getParameterDimension(); ++col)
    {
        int row_start = lower_part_only ? col : 0;
        for (int row = row_start; row < getParameterDimension(); ++row)
        {
            values[nz_idx] = dense_hessian(row, col);
            ++nz_idx;
        }
    }

    if (multiplier != 1.0) values *= multiplier;
}

void OptimizationProblemInterface::computeSparseHessianObjective(Eigen::SparseMatrix<double>& hessian, double multiplier)
{
    PRINT_DEBUG_COND_ONCE(_warn_if_not_specialized,
                          "OptimizationProblemInterface::computeSparseHessianObjective(): default implementation might be slow.");
    if (getObjectiveDimension() == 0) return;

    int dim_x = getParameterDimension();
    Eigen::MatrixXd dense_hessian(dim_x, dim_x);
    computeDenseHessianObjective(dense_hessian);
    hessian = dense_hessian.sparseView();

    if (multiplier != 1.0) hessian *= multiplier;
}

void OptimizationProblemInterface::computeSparseHessianObjectiveLL(Eigen::SparseMatrix<double, Eigen::ColMajor, long long>& H,
                                                                   const Eigen::VectorXi* col_nnz, bool upper_part_only)
{
    PRINT_DEBUG_COND_ONCE(_warn_if_not_specialized,
                          "OptimizationProblemInterface::computeSparseHessiansValues(): default implementation might be slow.");

    int dim_x = getParameterDimension();

    H.setZero();

    // TODO(roesmann): can be parallelized

    if (upper_part_only)
    {
        // we cannot use .selfadjointView<>() with mixed datatypes (SparseMatrix with long long and int indices)
        Eigen::SparseMatrix<double> H_obj(dim_x, dim_x);
        computeSparseHessianObjective(H_obj, 1.0);

        Eigen::SparseMatrix<double, Eigen::ColMajor, long long> H_temp(dim_x, dim_x);
        H_temp                            = H_obj;
        H.selfadjointView<Eigen::Upper>() = H_temp.selfadjointView<Eigen::Upper>();
    }
    else
    {
        Eigen::SparseMatrix<double> H_obj(dim_x, dim_x);
        computeSparseHessianObjective(H_obj, 1.0);
        H = H_obj;  // different format (TODO(roesmann): templates...)
    }
}

void OptimizationProblemInterface::computeSparseHessianObjectiveNNZperCol(Eigen::Ref<Eigen::VectorXi> col_nnz, bool upper_part_only)
{
    PRINT_DEBUG_COND_ONCE(_warn_if_not_specialized,
                          "OptimizationProblemInterface::computeSparseHessianLagrangianNNZperCol(): default implementation might be slow.");

    assert(col_nnz.size() == getParameterDimension());

    if (upper_part_only)
    {
        for (int i = 0; i < col_nnz.size(); ++i)
        {
            col_nnz(i) = i + 1;  // worst case nnz
        }
    }
    else
        col_nnz.setConstant(getParameterDimension());  // worst case nnz
}

void OptimizationProblemInterface::computeDenseHessianEqualities(Eigen::Ref<Eigen::MatrixXd> hessian, const double* multipliers)
{
    PRINT_DEBUG_COND_ONCE(_warn_if_not_specialized,
                          "OptimizationProblemInterface::computeDenseHessianEqualities(): default implementation might be slow.");
    int dim_eq = getEqualityDimension();
    int dim_x  = getParameterDimension();

    if (dim_eq == 0) return;

    assert(hessian.rows() == dim_x);
    assert(hessian.cols() == dim_x);

    // TODO(roesmann) generic interface
    CentralDifferences diff;
    auto inc  = [this](int idx, double inc) { applyIncrement(idx, inc); };
    auto eval = [this](Eigen::Ref<Eigen::VectorXd> values) { computeValuesEquality(values); };
    diff.computeHessian(inc, eval, dim_eq, hessian, multipliers);
}

int OptimizationProblemInterface::computeSparseHessianEqualitiesNNZ(bool lower_part_only)
{
    PRINT_DEBUG_COND_ONCE(_warn_if_not_specialized,
                          "OptimizationProblemInterface::computeSparseHessianEqualitiesNNZ(): default implementation might be slow.");
    // the default implementation is not aware of any sparse structure, hence we assume worst case
    if (lower_part_only) return (int)(0.5 * getParameterDimension() * (getParameterDimension() + 1));
    return getParameterDimension() * getParameterDimension();
}
void OptimizationProblemInterface::computeSparseHessianEqualitiesStructure(Eigen::Ref<Eigen::VectorXi> i_row, Eigen::Ref<Eigen::VectorXi> j_col,
                                                                           bool lower_part_only)
{
    PRINT_DEBUG_COND_ONCE(_warn_if_not_specialized,
                          "OptimizationProblemInterface::computeSparseHessianEqualitiesStructure(): default implementation might be slow.");
    // worst-case implementation
    int nz_idx = 0;
    for (int col = 0; col < getParameterDimension(); ++col)
    {
        int row_start = lower_part_only ? col : 0;
        for (int row = row_start; row < getParameterDimension(); ++row)
        {
            i_row[nz_idx] = row;
            j_col[nz_idx] = col;
            ++nz_idx;
        }
    }
}
void OptimizationProblemInterface::computeSparseHessianEqualitiesValues(Eigen::Ref<Eigen::VectorXd> values, const double* multipliers,
                                                                        bool lower_part_only)
{
    PRINT_DEBUG_COND_ONCE(_warn_if_not_specialized,
                          "OptimizationProblemInterface::computeSparseHessianEqualitiesValues(): default implementation might be slow.");
    // worst-case implementation
    int dim_x = getParameterDimension();
    Eigen::MatrixXd dense_hessian(dim_x, dim_x);
    computeDenseHessianEqualities(dense_hessian, multipliers);

    int nz_idx = 0;
    for (int col = 0; col < getParameterDimension(); ++col)
    {
        int row_start = lower_part_only ? col : 0;
        for (int row = row_start; row < getParameterDimension(); ++row)
        {
            values[nz_idx] = dense_hessian(row, col);
            ++nz_idx;
        }
    }
}

void OptimizationProblemInterface::computeSparseHessianEqualities(Eigen::SparseMatrix<double>& hessian, const double* multipliers)
{
    PRINT_DEBUG_COND_ONCE(_warn_if_not_specialized,
                          "OptimizationProblemInterface::computeSparseHessianEqualities(): default implementation might be slow.");
    int dim_x = getParameterDimension();
    Eigen::MatrixXd dense_hessian(dim_x, dim_x);
    computeDenseHessianEqualities(dense_hessian, multipliers);
    hessian = dense_hessian.sparseView();
}

void OptimizationProblemInterface::computeDenseHessianInequalities(Eigen::Ref<Eigen::MatrixXd> hessian, const double* multipliers)
{
    PRINT_DEBUG_COND_ONCE(_warn_if_not_specialized,
                          "OptimizationProblemInterface::computeDenseHessianInequalities(): default implementation might be slow.");
    int dim_ineq = getInequalityDimension();
    int dim_x    = getParameterDimension();

    if (dim_ineq == 0) return;

    assert(hessian.rows() == dim_x);
    assert(hessian.cols() == dim_x);

    // TODO(roesmann) generic interface
    CentralDifferences diff;
    auto inc  = [this](int idx, double inc) { applyIncrement(idx, inc); };
    auto eval = [this](Eigen::Ref<Eigen::VectorXd> values) { computeValuesInequality(values); };
    diff.computeHessian(inc, eval, dim_ineq, hessian, multipliers);
}

int OptimizationProblemInterface::computeSparseHessianInequalitiesNNZ(bool lower_part_only)
{
    PRINT_DEBUG_COND_ONCE(_warn_if_not_specialized,
                          "OptimizationProblemInterface::computeSparseHessianInequalitiesNNZ(): default implementation might be slow.");
    // the default implementation is not aware of any sparse structure, hence we assume worst case
    if (lower_part_only) return (int)(0.5 * getParameterDimension() * (getParameterDimension() + 1));
    return getParameterDimension() * getParameterDimension();
}
void OptimizationProblemInterface::computeSparseHessianInequalitiesStructure(Eigen::Ref<Eigen::VectorXi> i_row, Eigen::Ref<Eigen::VectorXi> j_col,
                                                                             bool lower_part_only)
{
    PRINT_DEBUG_COND_ONCE(_warn_if_not_specialized,
                          "OptimizationProblemInterface::computeSparseHessianInequalitiesStructure(): default implementation might be slow.");
    // worst-case implementation
    int nz_idx = 0;
    for (int col = 0; col < getParameterDimension(); ++col)
    {
        int row_start = lower_part_only ? col : 0;
        for (int row = row_start; row < getParameterDimension(); ++row)
        {
            i_row[nz_idx] = row;
            j_col[nz_idx] = col;
            ++nz_idx;
        }
    }
}
void OptimizationProblemInterface::computeSparseHessianInequalitiesValues(Eigen::Ref<Eigen::VectorXd> values, const double* multipliers,
                                                                          bool lower_part_only)
{
    PRINT_DEBUG_COND_ONCE(_warn_if_not_specialized,
                          "OptimizationProblemInterface::computeSparseHessianInequalitiesValues(): default implementation might be slow.");
    // worst-case implementation
    int dim_x = getParameterDimension();
    Eigen::MatrixXd dense_hessian(dim_x, dim_x);
    computeDenseHessianInequalities(dense_hessian, multipliers);

    int nz_idx = 0;
    for (int col = 0; col < getParameterDimension(); ++col)
    {
        int row_start = lower_part_only ? col : 0;
        for (int row = row_start; row < getParameterDimension(); ++row)
        {
            values[nz_idx] = dense_hessian(row, col);
            ++nz_idx;
        }
    }
}

void OptimizationProblemInterface::computeSparseHessianInequalities(Eigen::SparseMatrix<double>& hessian, const double* multipliers)
{
    PRINT_DEBUG_COND_ONCE(_warn_if_not_specialized,
                          "OptimizationProblemInterface::computeSparseHessianInequalities(): default implementation might be slow.");
    int dim_x = getParameterDimension();
    Eigen::MatrixXd dense_hessian(dim_x, dim_x);
    computeDenseHessianInequalities(dense_hessian, multipliers);
    hessian = dense_hessian.sparseView();
}

// void OptimizationProblemInterface::computeDenseJacobianHessianObjective(Eigen::Ref<Eigen::MatrixXd> jacobian, Eigen::Ref<Eigen::MatrixXd> hessian,
//                                                                        const double* multipliers)
//{
//    computeDenseJacobianLsqObjective(jacobian, multipliers);
//    computeDenseHessianObjective(hessian, multipliers);
//}
// void OptimizationProblemInterface::computeDenseJacobianHessianEqualities(Eigen::Ref<Eigen::MatrixXd> jacobian, Eigen::Ref<Eigen::MatrixXd> hessian,
//                                                                         const double* multipliers)
//{
//    computeDenseJacobianEqualities(jacobian, multipliers);
//    computeDenseHessianEqualities(hessian, multipliers);
//}
// void OptimizationProblemInterface::computeDenseJacobianHessianInequalities(Eigen::Ref<Eigen::MatrixXd> jacobian, Eigen::Ref<Eigen::MatrixXd>
// hessian,
//                                                                           const double* multipliers)
//{
//    computeDenseJacobianInequalities(jacobian, multipliers);
//    computeDenseHessianInequalities(hessian, multipliers);
//}

void OptimizationProblemInterface::computeDenseHessians(Eigen::Ref<Eigen::MatrixXd> hessian_obj, Eigen::Ref<Eigen::MatrixXd> hessian_eq,
                                                        Eigen::Ref<Eigen::MatrixXd> hessian_ineq, double multiplier_obj, const double* multipliers_eq,
                                                        const double* multipliers_ineq)
{
    PRINT_DEBUG_COND_ONCE(_warn_if_not_specialized, "OptimizationProblemInterface::computeDenseHessians(): default implementation might be slow.");
    int obj_dim  = getObjectiveDimension();
    int eq_dim   = getEqualityDimension();
    int ineq_dim = getInequalityDimension();
    int idx      = 0;
    if (obj_dim > 0)
    {
        computeDenseHessianObjective(hessian_obj, multiplier_obj);
        idx += obj_dim;
    }
    if (eq_dim > 0)
    {
        computeDenseHessianEqualities(hessian_eq, multipliers_eq);
        idx += eq_dim;
    }
    if (ineq_dim > 0)
    {
        computeDenseHessianInequalities(hessian_ineq, multipliers_ineq);
    }
}

void OptimizationProblemInterface::computeSparseHessians(Eigen::SparseMatrix<double>& hessian_obj, Eigen::SparseMatrix<double>& hessian_eq,
                                                         Eigen::SparseMatrix<double>& hessian_ineq, double multiplier_obj,
                                                         const double* multipliers_eq, const double* multipliers_ineq)
{
    PRINT_DEBUG_COND_ONCE(_warn_if_not_specialized, "OptimizationProblemInterface::computeSparseHessians(): default implementation might be slow.");
    computeSparseHessianObjective(hessian_obj, multiplier_obj);
    computeSparseHessianEqualities(hessian_eq, multipliers_eq);
    computeSparseHessianInequalities(hessian_ineq, multipliers_ineq);
}

void OptimizationProblemInterface::computeSparseHessiansNNZ(int& nnz_obj, int& nnz_eq, int& nnz_ineq, bool lower_part_only)
{
    PRINT_DEBUG_COND_ONCE(_warn_if_not_specialized,
                          "OptimizationProblemInterface::computeSparseHessiansNNZ(): default implementation might be slow.");
    nnz_obj  = computeSparseHessianObjectiveNNZ(lower_part_only);
    nnz_eq   = computeSparseHessianEqualitiesNNZ(lower_part_only);
    nnz_ineq = computeSparseHessianInequalitiesNNZ(lower_part_only);
}

void OptimizationProblemInterface::computeSparseHessiansStructure(Eigen::Ref<Eigen::VectorXi> i_row_obj, Eigen::Ref<Eigen::VectorXi> j_col_obj,
                                                                  Eigen::Ref<Eigen::VectorXi> i_row_eq, Eigen::Ref<Eigen::VectorXi> j_col_eq,
                                                                  Eigen::Ref<Eigen::VectorXi> i_row_ineq, Eigen::Ref<Eigen::VectorXi> j_col_ineq,
                                                                  bool lower_part_only)
{
    PRINT_DEBUG_COND_ONCE(_warn_if_not_specialized,
                          "OptimizationProblemInterface::computeSparseHessiansValues(): default implementation might be slow.");
    computeSparseHessianObjectiveStructure(i_row_obj, j_col_obj, lower_part_only);
    computeSparseHessianEqualitiesStructure(i_row_eq, j_col_eq, lower_part_only);
    computeSparseHessianInequalitiesStructure(i_row_ineq, j_col_ineq, lower_part_only);
}

void OptimizationProblemInterface::computeSparseHessiansValues(Eigen::Ref<Eigen::VectorXd> values_obj, Eigen::Ref<Eigen::VectorXd> values_eq,
                                                               Eigen::Ref<Eigen::VectorXd> values_ineq, double multiplier_obj,
                                                               const double* multipliers_eq, const double* multipliers_ineq, bool lower_part_only)
{
    PRINT_DEBUG_COND_ONCE(_warn_if_not_specialized,
                          "OptimizationProblemInterface::computeSparseHessiansValues(): default implementation might be slow.");
    computeSparseHessianObjectiveValues(values_obj, multiplier_obj, lower_part_only);
    computeSparseHessianEqualitiesValues(values_eq, multipliers_eq, lower_part_only);
    computeSparseHessianInequalitiesValues(values_ineq, multipliers_ineq, lower_part_only);
}

void OptimizationProblemInterface::computeSparseHessianLagrangian(Eigen::SparseMatrix<double, Eigen::ColMajor, long long>& H,
                                                                  const double* multipliers_eq, const double* multipliers_ineq,
                                                                  const Eigen::VectorXi* col_nnz, bool upper_part_only)
{
    PRINT_DEBUG_COND_ONCE(_warn_if_not_specialized,
                          "OptimizationProblemInterface::computeSparseHessiansValues(): default implementation might be slow.");

    int dim_x    = getParameterDimension();
    int dim_eq   = getEqualityDimension();
    int dim_ineq = getInequalityDimension();

    H.setZero();

    // TODO(roesmann): can be parallelized

    Eigen::SparseMatrix<double> H_obj(dim_x, dim_x), H_eq(dim_x, dim_x), H_ineq(dim_x, dim_x);
    computeSparseHessianObjective(H_obj, 1.0);
    if (dim_eq > 0) computeSparseHessianEqualities(H_eq, multipliers_eq);
    if (dim_ineq > 0) computeSparseHessianInequalities(H_ineq, multipliers_ineq);

    if (upper_part_only)
    {
        // we cannot use .selfadjointView<>() with mixed datatypes (SparseMatrix with long long and int indices)
        Eigen::SparseMatrix<double, Eigen::ColMajor, long long> H_temp(dim_x, dim_x);

        if (dim_eq > 0 && dim_ineq > 0)
        {
            H_temp = H_obj + H_eq + H_ineq;
        }
        else if (dim_eq > 0)
        {
            H_temp = H_obj + H_eq;
        }
        else if (dim_ineq > 0)
        {
            H_temp = H_obj + H_ineq;
        }
        else
        {
            // if (upper_part_only)
            // H.selfadjointView<Eigen::Upper>() = H_obj.selfadjointView<Eigen::Upper>();
            // else
            H_temp = H_obj;
        }
        H.selfadjointView<Eigen::Upper>() = H_temp.selfadjointView<Eigen::Upper>();  // we note above!
    }
    else
    {
        if (dim_eq > 0 && dim_ineq > 0)
        {
            H = H_obj + H_eq + H_ineq;
        }
        else if (dim_eq > 0)
        {
            H = H_obj + H_eq;
        }
        else if (dim_ineq > 0)
        {
            H = H_obj + H_ineq;
        }
        else
        {
            H = H_obj;
        }
    }
}

void OptimizationProblemInterface::computeSparseHessianLagrangianNNZperCol(Eigen::Ref<Eigen::VectorXi> col_nnz, bool upper_part_only)
{
    PRINT_DEBUG_COND_ONCE(_warn_if_not_specialized,
                          "OptimizationProblemInterface::computeSparseHessianLagrangianNNZperCol(): default implementation might be slow.");

    assert(col_nnz.size() == getParameterDimension());

    if (upper_part_only)
    {
        for (int i = 0; i < col_nnz.size(); ++i)
        {
            col_nnz(i) = i + 1;  // worst case nnz
        }
    }
    else
        col_nnz.setConstant(getParameterDimension());  // worst case nnz
}

void OptimizationProblemInterface::computeSparseJacobianTwoSideBoundedLinearFormAndHessianLagrangian(
    Eigen::SparseMatrix<double, Eigen::ColMajor, long long>& H, const double* multipliers_eq, const double* multipliers_ineq,
    Eigen::SparseMatrix<double, Eigen::ColMajor, long long>& A, bool include_finite_bounds, const Eigen::VectorXi* col_nnz_H,
    const Eigen::VectorXi* col_nnz_A, bool upper_part_only_H)
{
    computeSparseJacobianTwoSideBoundedLinearForm(A, include_finite_bounds, col_nnz_A);
    computeSparseHessianLagrangian(H, multipliers_eq, multipliers_ineq, col_nnz_H, upper_part_only_H);
}

void OptimizationProblemInterface::computeSparseJacobianTwoSideBoundedLinearFormAndHessianObjective(
    Eigen::SparseMatrix<double, Eigen::ColMajor, long long>& H, Eigen::SparseMatrix<double, Eigen::ColMajor, long long>& A,
    bool include_finite_bounds, const Eigen::VectorXi* col_nnz_H, const Eigen::VectorXi* col_nnz_A, bool upper_part_only_H)
{
    computeSparseJacobianTwoSideBoundedLinearForm(A, include_finite_bounds, col_nnz_A);
    computeSparseHessianObjectiveLL(H, col_nnz_H, upper_part_only_H);
}

bool OptimizationProblemInterface::checkIfAllUnfixedParam(std::function<bool(double, int)> fun)
{
    PRINT_DEBUG_COND_ONCE(_warn_if_not_specialized, "OptimizationProblemInterface::checkIfAllUnfixedParam(): default implementation might be slow.");

    for (int i = 0; i < getParameterDimension(); ++i)
    {
        if (!fun(getParameterValue(i), i)) return false;
    }
    return true;
}

}  // namespace corbo
