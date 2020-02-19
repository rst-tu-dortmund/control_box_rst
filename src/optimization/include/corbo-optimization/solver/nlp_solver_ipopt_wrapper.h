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

#ifndef SRC_OPTIMIZATION_INCLUDE_CORBO_OPTIMIZATION_SOLVER_NLP_SOLVER_IPOPT_WRAPPER_H_
#define SRC_OPTIMIZATION_INCLUDE_CORBO_OPTIMIZATION_SOLVER_NLP_SOLVER_IPOPT_WRAPPER_H_

#include <corbo-optimization/optimization_problem_interface.h>

#define HAVE_CSTDDEF
#include <IpTNLP.hpp>
#undef HAVE_CSTDDEF

// #include <cstddef>
// #define HAVE_CSTDDEF

namespace corbo {

class SolverIpopt;  // forward declaration

class IpoptWrapper : public Ipopt::TNLP
{
 private:
    using Index                     = Ipopt::Index;
    using Number                    = Ipopt::Number;
    using SolverReturn              = Ipopt::SolverReturn;
    using IpoptData                 = Ipopt::IpoptData;
    using IpoptCalculatedQuantities = Ipopt::IpoptCalculatedQuantities;

 public:
    /** default constructor */
    explicit IpoptWrapper(SolverIpopt* solver);

    /** default destructor */
    virtual ~IpoptWrapper();

    void setOptimizationProblem(OptimizationProblemInterface& problem) { _problem = &problem; }

    /**@name Overloaded from TNLP */
    //@{

    /** Method to return some info about the nlp */
    bool get_nlp_info(Index& n, Index& m, Index& nnz_jac_g, Index& nnz_h_lag, IndexStyleEnum& index_style) override;

    /** Method to return the bounds for my problem */
    bool get_bounds_info(Index n, Number* x_l, Number* x_u, Index m, Number* g_l, Number* g_u) override;

    /** Method to return the starting point for the algorithm */
    bool get_starting_point(Index n, bool init_x, Number* x, bool init_z, Number* z_L, Number* z_U, Index m, bool init_lambda,
                            Number* lambda) override;

    /** Method to return the objective value */
    bool eval_f(Index n, const Number* x, bool new_x, Number& obj_value) override;

    /** Method to return the gradient of the objective */
    bool eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f) override;

    /** Method to return the constraint residuals */
    bool eval_g(Index n, const Number* x, bool new_x, Index m, Number* g) override;

    /** Method to return:
     *   1) The structure of the jacobian (if "values" is NULL)
     *   2) The values of the jacobian (if "values" is not NULL)
     */
    bool eval_jac_g(Index n, const Number* x, bool new_x, Index m, Index nele_jac, Index* iRow, Index* jCol, Number* values) override;

    /** Method to return:
     *   1) The structure of the hessian of the lagrangian (if "values" is NULL)
     *   2) The values of the hessian of the lagrangian (if "values" is not NULL)
     */
    bool eval_h(Index n, const Number* x, bool new_x, Number obj_factor, Index m, const Number* lambda, bool new_lambda, Index nele_hess, Index* iRow,
                Index* jCol, Number* values) override;

    //@}

    /** @name Solution Methods */
    //@{

    /** This method is called when the algorithm is complete so the TNLP can store/write the solution */
    void finalize_solution(SolverReturn status, Index n, const Number* x, const Number* z_L, const Number* z_U, Index m, const Number* g,
                           const Number* lambda, Number obj_value, const IpoptData* ip_data, IpoptCalculatedQuantities* ip_cq) override;
    //@}

    void precompute1stOrderDerivatives();

 private:
    OptimizationProblemInterface* _problem = nullptr;
    SolverIpopt* _solver                   = nullptr;
};

}  // namespace corbo

#endif  // SRC_OPTIMIZATION_INCLUDE_CORBO_OPTIMIZATION_SOLVER_NLP_SOLVER_IPOPT_WRAPPER_H_
