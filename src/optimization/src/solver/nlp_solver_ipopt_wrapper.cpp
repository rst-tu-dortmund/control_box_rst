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

#ifdef IPOPT

#include <corbo-optimization/solver/nlp_solver_ipopt_wrapper.h>

#include <corbo-optimization/solver/nlp_solver_ipopt.h>

#include <cassert>
#include <iostream>

namespace corbo {

// constructor
IpoptWrapper::IpoptWrapper(SolverIpopt* solver) : _solver(solver) {}

// destructor
IpoptWrapper::~IpoptWrapper() {}

// returns the size of the problem
bool IpoptWrapper::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g, Index& nnz_h_lag, IndexStyleEnum& index_style)
{
    // Number of optimization variables
    n = _problem->getParameterDimension();

    // Number of equality and inequality constraints (not bounds)
    m = _problem->getEqualityDimension() + _problem->getInequalityDimension();

    // Number of nonzeros of the jacobian
    nnz_jac_g = _solver->_nnz_jac_constraints;

    if (_solver->_cache_first_order_derivatives)
    {
        _solver->_grad_f_cache.resize(n);
        _solver->_jac_constr_cache.resize(nnz_jac_g);
    }

    // Number of nonzeros of the hessian of the lagrangian, but we
    // only need the lower left corner (since it is symmetric)
    // nnz_h_lag = 0.5 * n * (n+1);
    //    assert(_solver_interface->_h_lag_i_row.size() == _solver_interface->_h_lag_j_col.size());
    nnz_h_lag = _solver->_nnz_h_lagrangian;

    // use the C style indexing (0-based)
    index_style = TNLP::C_STYLE;

    return true;
}

// returns the variable bounds
bool IpoptWrapper::get_bounds_info(Index n, Number* x_l, Number* x_u, Index m, Number* g_l, Number* g_u)
{
    // here, the n and m we gave IPOPT in get_nlp_info are passed back to us.
    // If desired, we could assert to make sure they are what we think they are.
    assert(n == _problem->getParameterDimension());
    assert(m == _problem->getEqualityDimension() + _problem->getInequalityDimension());

    // apply bounds
    Eigen::Map<Eigen::VectorXd> x_l_map(x_l, n);
    Eigen::Map<Eigen::VectorXd> x_u_map(x_u, n);
    _problem->getBounds(x_l_map, x_u_map);

    // equality and inequality constraints are in the form g_l <= g(x) <= g_u
    if (m > 0)
    {
        Eigen::Map<Eigen::VectorXd> g_l_map(g_l, m);
        Eigen::Map<Eigen::VectorXd> g_u_map(g_u, m);
        g_l_map.head(_problem->getEqualityDimension()).fill(0);
        g_u_map.head(_problem->getEqualityDimension()).fill(0);
        g_l_map.tail(_problem->getInequalityDimension()).fill(-CORBO_INF_DBL);
        g_u_map.tail(_problem->getInequalityDimension()).fill(0);
    }

    return true;
}

// returns the initial point for the problem
bool IpoptWrapper::get_starting_point(Index n, bool init_x, Number* x, bool init_z, Number* z_L, Number* z_U, Index m, bool init_lambda,
                                      Number* lambda)
{
    // Here, we assume we always have starting values for x
    assert(init_x == true);

    // initialize to the given starting point
    Eigen::Map<Eigen::VectorXd> x_map(x, n);
    _problem->getParameterVector(x_map);

    if (init_lambda)  // is true if warmstart is enabled
    {
        Eigen::Map<Eigen::VectorXd> lambda_map(lambda, m);
        lambda_map = _solver->_lambda_cache;
    }

    if (init_z)  // is true if warmstart is enabled
    {
        Eigen::Map<Eigen::VectorXd> zl_map(z_L, n);
        zl_map = _solver->_zl_cache;
        Eigen::Map<Eigen::VectorXd> zu_map(z_U, n);
        zu_map = _solver->_zu_cache;
    }

    return true;
}

// returns the value of the objective function
bool IpoptWrapper::eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
{
    assert(n == _problem->getParameterDimension());

    if (new_x)
    {
        // in the current solver implementation, we must copy the optim vector back
        // which is not really efficient. This could be changed in the future...
        Eigen::Map<const Eigen::VectorXd> x_map(x, n);
        _problem->setParameterVector(x_map);

        if (_solver->_cache_first_order_derivatives) precompute1stOrderDerivatives();
    }

    obj_value = _problem->computeValueObjective();
    return true;
}

// return the gradient of the objective function grad_{x} f(x)
bool IpoptWrapper::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
{
    assert(n == _problem->getParameterDimension());

    if (new_x)
    {
        // in the current solver implementation, we must copy the optim vector back
        // which is not really efficient. This could be changed in the future...
        Eigen::Map<const Eigen::VectorXd> x_map(x, n);
        _problem->setParameterVector(x_map);

        if (_solver->_cache_first_order_derivatives) precompute1stOrderDerivatives();
    }

    Eigen::Map<Eigen::VectorXd> grad_f_map(grad_f, n);
    if (_solver->_cache_first_order_derivatives)
        grad_f_map = _solver->_grad_f_cache;
    else
        _problem->computeGradientObjective(grad_f_map);

    return true;
}

// return the value of the constraints: g(x)
bool IpoptWrapper::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
{
    assert(n == _problem->getParameterDimension());
    assert(m == _problem->getEqualityDimension() + _problem->getInequalityDimension());

    if (new_x)
    {
        // in the current solver implementation, we must copy the optim vector back
        // which is not really efficient. This could be changed in the future...
        Eigen::Map<const Eigen::VectorXd> x_map(x, n);
        _problem->setParameterVector(x_map);

        if (_solver->_cache_first_order_derivatives) precompute1stOrderDerivatives();
    }

    Eigen::Map<Eigen::VectorXd> g_map(g, m);
    _problem->computeValuesEquality(g_map.head(_problem->getEqualityDimension()));
    _problem->computeValuesInequality(g_map.tail(_problem->getInequalityDimension()));

    return true;
}

// return the structure or values of the jacobian
bool IpoptWrapper::eval_jac_g(Index n, const Number* x, bool new_x, Index m, Index nele_jac, Index* iRow, Index* jCol, Number* values)
{
    assert(nele_jac == _solver->_nnz_jac_constraints);

    if (values == NULL)
    {
        // return the structure of the jacobian
        Eigen::Map<Eigen::VectorXi> i_row_map(iRow, nele_jac);
        Eigen::Map<Eigen::VectorXi> j_col_map(jCol, nele_jac);
        _problem->computeCombinedSparseJacobiansStructure(i_row_map, j_col_map, false, true, true);
        //   std::copy(_solver_interface->_jac_constr_i_row.begin(), _solver_interface->_jac_constr_i_row.end(), iRow);
        //   std::copy(_solver_interface->_jac_constr_j_col.begin(), _solver_interface->_jac_constr_j_col.end(), jCol);
    }
    else
    {
        if (new_x)
        {
            // in the current solver implementation, we must copy the optim vector back
            // which is not really efficient. This could be changed in the future...
            Eigen::Map<const Eigen::VectorXd> x_map(x, n);
            _problem->setParameterVector(x_map);

            if (_solver->_cache_first_order_derivatives) precompute1stOrderDerivatives();
        }
        Eigen::Map<Eigen::VectorXd> values_map(values, nele_jac);
        if (_solver->_cache_first_order_derivatives)
            values_map = _solver->_jac_constr_cache;
        else
            _problem->computeCombinedSparseJacobiansValues(values_map, false, true, true);

        // return the values of the jacobian of the constraints
        // std::copy(_solver_interface->_jac_constr_values.begin(), _solver_interface->_jac_constr_values.end(), values);
    }

    return true;
}

// return the structure or values of the hessian
bool IpoptWrapper::eval_h(Index n, const Number* x, bool new_x, Number obj_factor, Index m, const Number* lambda, bool new_lambda, Index nele_hess,
                          Index* iRow, Index* jCol, Number* values)
{
    assert(nele_hess == _solver->_nnz_h_lagrangian);
    assert(nele_hess = _solver->_nnz_hes_obj + _solver->_nnz_hes_eq + _solver->_nnz_hes_ineq);

    if (values == NULL)
    {
        // return the structure. This is a symmetric matrix, fill the lower left  triangle only.
        //        std::copy(_solver_interface->_h_lag_i_row.begin(), _solver_interface->_h_lag_i_row.end(), iRow);
        //        std::copy(_solver_interface->_h_lag_j_col.begin(), _solver_interface->_h_lag_j_col.end(), jCol);
        Eigen::Map<Eigen::VectorXi> i_row_map_obj(iRow, _solver->_nnz_hes_obj);
        Eigen::Map<Eigen::VectorXi> j_col_map_obj(jCol, _solver->_nnz_hes_obj);
        Eigen::Map<Eigen::VectorXi> i_row_map_eq(iRow + _solver->_nnz_hes_obj, _solver->_nnz_hes_eq);
        Eigen::Map<Eigen::VectorXi> j_col_map_eq(jCol + _solver->_nnz_hes_obj, _solver->_nnz_hes_eq);
        Eigen::Map<Eigen::VectorXi> i_row_map_ineq(iRow + _solver->_nnz_hes_obj + _solver->_nnz_hes_eq, _solver->_nnz_hes_ineq);
        Eigen::Map<Eigen::VectorXi> j_col_map_ineq(jCol + _solver->_nnz_hes_obj + _solver->_nnz_hes_eq, _solver->_nnz_hes_ineq);
        _problem->computeSparseHessiansStructure(i_row_map_obj, j_col_map_obj, i_row_map_eq, j_col_map_eq, i_row_map_ineq, j_col_map_ineq, true);
    }
    else
    {
        //        // return the values. This is a symmetric matrix, fill the lower left
        //        // triangle only

        if (new_x)
        {
            //            // in the current solver implementation, we must copy the optim vector back
            //            // which is not really efficient. This could be changed in the future...
            Eigen::Map<const Eigen::VectorXd> x_map(x, n);
            _problem->setParameterVector(x_map);

            if (_solver->_cache_first_order_derivatives) precompute1stOrderDerivatives();
        }

        //        _solver_interface->computeHessianLagrangianValues(obj_factor, lambda);
        Eigen::Map<Eigen::VectorXd> values_obj_map(values, _solver->_nnz_hes_obj);
        Eigen::Map<Eigen::VectorXd> values_eq_map(values + _solver->_nnz_hes_obj, _solver->_nnz_hes_eq);
        Eigen::Map<Eigen::VectorXd> values_ineq_map(values + _solver->_nnz_hes_obj + _solver->_nnz_hes_eq, _solver->_nnz_hes_ineq);
        _problem->computeSparseHessiansValues(values_obj_map, values_eq_map, values_ineq_map, obj_factor, lambda,
                                              lambda + _problem->getEqualityDimension(), true);
        //        std::copy(_solver_interface->_h_lag_values.begin(), _solver_interface->_h_lag_values.end(), values);
    }

    return true;
}

void IpoptWrapper::finalize_solution(SolverReturn status, Index n, const Number* x, const Number* z_L, const Number* z_U, Index m, const Number* g,
                                     const Number* lambda, Number obj_value, const IpoptData* ip_data, IpoptCalculatedQuantities* ip_cq)
{
    Eigen::Map<const Eigen::VectorXd> x_map(x, n);
    _problem->setParameterVector(x_map);

    _solver->_last_obj_value = obj_value;

    Eigen::Map<const Eigen::VectorXd> lambda_map(lambda, m);
    _solver->_lambda_cache = lambda_map;

    Eigen::Map<const Eigen::VectorXd> zl_map(z_L, n);
    _solver->_zl_cache = zl_map;

    Eigen::Map<const Eigen::VectorXd> zu_map(z_U, n);
    _solver->_zu_cache = zu_map;
}

void IpoptWrapper::precompute1stOrderDerivatives()
{
    _problem->computeGradientObjectiveAndCombinedSparseJacobiansValues(_solver->_grad_f_cache, _solver->_jac_constr_cache, true, true);
}

}  // namespace corbo

#endif  // IPOPT
