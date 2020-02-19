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

#ifdef OSQP

#include <corbo-optimization/solver/qp_solver_osqp.h>

#include <corbo-core/console.h>

namespace corbo {

SolverOsqp::SolverOsqp()
{
    // Problem settings
    _settings = std::unique_ptr<OSQPSettings>(new OSQPSettings);

    // Define Solver settings as default
    osqp_set_default_settings(_settings.get());

    // disable verbosity by default
    _settings->verbose = 0;
}

SolverOsqp::~SolverOsqp()
{
    // Cleanup
    if (_work) osqp_cleanup(_work);
}

bool SolverOsqp::initialize()
{
    if (_initialized) return true;

    assert(CORBO_INF_DBL >= OSQP_INFTY);

    _initialized = true;
    return true;
}

SolverStatus SolverOsqp::solve(SparseMatrix& P, Eigen::Ref<Eigen::VectorXd> q, SparseMatrix& A, Eigen::Ref<Eigen::VectorXd> lbA,
                               Eigen::Ref<Eigen::VectorXd> ubA, Eigen::Ref<Eigen::VectorXd> lb, Eigen::Ref<Eigen::VectorXd> ub, bool new_structure,
                               bool zero_x_warmstart)
{
    PRINT_ERROR("This method is currently not implemented.");
    return SolverStatus::Error;
}

SolverStatus SolverOsqp::solve(SparseMatrix& P, Eigen::Ref<Eigen::VectorXd> q, SparseMatrix& A, Eigen::Ref<Eigen::VectorXd> lbA,
                               Eigen::Ref<Eigen::VectorXd> ubA, bool new_structure, bool zero_x_warmstart, bool update_P, bool update_q,
                               bool update_A, bool update_bounds)
{
    if (!_initialized)
    {
        if (!initialize()) return SolverStatus::Error;
    }

    // TODO(roesmann): we copy now the solution vectors, but just for rapid prototyping (to avoid this complete SQP cleanup mess etc.). We will change
    // this later

    // P and A must be incompressed mode
    if (!P.isCompressed()) P.makeCompressed();
    if (!A.isCompressed()) A.makeCompressed();

    csc* P_csc = csc_matrix(P.rows(), P.cols(), P.nonZeros(), P.valuePtr(), P.innerIndexPtr(), P.outerIndexPtr());
    csc* A_csc = csc_matrix(A.rows(), A.cols(), A.nonZeros(), A.valuePtr(), A.innerIndexPtr(), A.outerIndexPtr());

    // Populate data
    std::unique_ptr<OSQPData> data;
    if (new_structure || _force_new_structure || !_settings->warm_start || !_work)
    {
        // cleanup previous structure
        if (_work) osqp_cleanup(_work);

        // TODO(roesmann): whenever I move  "_data    = std::unique_ptr<OSQPData>(new OSQPData);" to initialize()
        // I get a segfault / issue in the assembler code of OSQP....
        data    = std::unique_ptr<OSQPData>(new OSQPData);
        data->n = P.rows();
        data->m = A.rows();
        data->P = P_csc;
        data->q = q.data();
        data->A = A_csc;
        data->l = lbA.data();
        data->u = ubA.data();

        // Setup workspace
        c_int exitflag = osqp_setup(&_work, data.get(), _settings.get());  // allocates memory!

        if (!_work || exitflag != 0)
        {
            PRINT_ERROR("Osqp setup failed.");
            return SolverStatus::Error;
        }

        _force_new_structure = false;
    }
    else
    {
        if (!_work)
        {
            PRINT_ERROR_NAMED("No previous workspace found. Cannot solve...");
            return SolverStatus::Error;
        }

        if (zero_x_warmstart)
        {
            if (_zero.size() != _work->data->n) _zero.setZero(_work->data->n);
            updatePrimalSolutionWarmStart(_zero);
        }

        int dim_constr = _work->data->m;

        if (update_P && update_A && dim_constr > 0)
        {
            if (osqp_update_P_A(_work, P.valuePtr(), OSQP_NULL, P.nonZeros(), A.valuePtr(), OSQP_NULL, A.nonZeros()) != 0)
            {
                PRINT_ERROR("SolverOSQP: Cannot update P and A due to dimensions mismatch. Maybe a new sparsity pattern?");
                return SolverStatus::Error;
            }
        }
        else if (update_P)  // P must be upper triangular!
        {
            if (osqp_update_P(_work, P.valuePtr(), OSQP_NULL, P.nonZeros()) != 0)
            {
                PRINT_ERROR("SolverOSQP: Cannot update P due to dimensions mismatch. Maybe a new sparsity pattern?");
                return SolverStatus::Error;
            }
        }
        else if (update_A && dim_constr > 0)
        {
            if (osqp_update_A(_work, A.valuePtr(), OSQP_NULL, A.nonZeros()) != 0)
            {
                PRINT_ERROR_NAMED("SolverOSQP: Cannot update A due to dimensions mismatch. Maybe a new sparsity pattern?");
                return SolverStatus::Error;
            }
        }

        if (update_q)
        {
            if (osqp_update_lin_cost(_work, q.data()) != 0)
            {
                PRINT_ERROR_NAMED("Cannot update q due to dimensions mismatch.");
                return SolverStatus::Error;
            }
        }

        if (update_bounds && dim_constr > 0)
        {
            if (osqp_update_bounds(_work, lbA.data(), ubA.data()) != 0)
            {
                PRINT_ERROR_NAMED("Cannot update lbA abd ubA due to dimensions mismatch.");
                return SolverStatus::Error;
            }
        }
    }

    // Solve Problem
    // c_int exit_flag =
    osqp_solve(_work);

    // WARNING the lagrange multipliers in OSQP are defiend for L = f(x) + lambda * c(x)

    //    if (exit_flag == 0)
    //    {
    //    }
    //    else
    //    {
    //    }

    if (data)
    {
        if (data->P) c_free(data->P);
        if (data->A) c_free(data->A);
    }

    // TODO(roesmann): osqp supports warm-starting and update of just the P matrix (hessian), so we should later go for that

    // TODO(roesmann): we could cache and efficiently update the local sparse identity hessian, since
    // both row and column vectors are just 0,1,2,3,... and the value vector 1,1,1,1,1
    return convertOsqpExitFlagToSolverStatus(_work->info->status_val);
}

Eigen::Ref<Eigen::VectorXd> SolverOsqp::getPrimalSolution()
{
    assert(_work);
    return Eigen::Map<Eigen::VectorXd>(_work->solution->x, _work->data->n);
}

Eigen::Ref<Eigen::VectorXd> SolverOsqp::getDualSolution()
{
    assert(_work);
    return Eigen::Map<Eigen::VectorXd>(_work->solution->y, _work->data->m);
}

void SolverOsqp::updatePrimalSolutionWarmStart(const Eigen::Ref<const Eigen::VectorXd>& x)
{
    if (!_work)
    {
        PRINT_ERROR_NAMED("No previous workspace found. Cannot update primal solution...");
        return;
    }
    if (x.size() == _work->data->n)
        osqp_warm_start_x(_work, x.data());
    else
        PRINT_ERROR_NAMED("Dimensions mismatch");
}
void SolverOsqp::updateDualSolutionWarmStart(const Eigen::Ref<const Eigen::VectorXd>& y)
{
    if (!_work)
    {
        PRINT_ERROR_NAMED("No previous workspace found. Cannot update dual solution...");
        return;
    }
    if (y.size() == _work->data->m)
    {
        if (y.size() > 0) osqp_warm_start_y(_work, y.data());
    }
    else
        PRINT_ERROR_NAMED("Dimensions mismatch");
}

SolverStatus SolverOsqp::convertOsqpExitFlagToSolverStatus(c_int status) const
{
    switch (status)
    {
        case OSQP_SOLVED:
        case OSQP_SOLVED_INACCURATE:
        {
            return SolverStatus::Converged;
        }
        case OSQP_MAX_ITER_REACHED:
        case OSQP_TIME_LIMIT_REACHED:
        case OSQP_SIGINT:
        {
            return SolverStatus::EarlyTerminated;
        }
        case OSQP_PRIMAL_INFEASIBLE:
        case OSQP_DUAL_INFEASIBLE:
        case OSQP_NON_CVX:
        {
            return SolverStatus::Infeasible;
        }
        case OSQP_UNSOLVED:
        {
            return SolverStatus::Error;
        }
        default:
        {
        }
    };
    // OSQP_DUAL_INFEASIBLE_INACCURATE
    // OSQP_PRIMAL_INFEASIBLE_INACCURATE
    return SolverStatus::Error;
}

void SolverOsqp::clear()
{
    if (_work)
    {
        osqp_cleanup(_work);
        _work = nullptr;
    }
    _force_new_structure = true;  // we could also check for _work in solve()
    _initialized         = false;
}

#ifdef MESSAGE_SUPPORT
void SolverOsqp::toMessage(corbo::messages::SolverOsqp& message) const
{
    message.set_max_iter(_settings->max_iter);
    message.set_eps_abs(_settings->eps_abs);
    message.set_eps_rel(_settings->eps_rel);
    message.set_eps_prim_inf(_settings->eps_prim_inf);
    message.set_eps_dual_inf(_settings->eps_dual_inf);
    message.set_alpha(_settings->alpha);
    message.set_verbose(_settings->verbose);
    message.set_scaled_termination(_settings->scaled_termination);
    message.set_check_termination(_settings->check_termination);
    message.set_warm_start(_settings->warm_start);
    message.set_time_limit(_settings->time_limit);

    message.set_rho(_settings->rho);
    message.set_sigma(_settings->sigma);
    message.set_scaling(_settings->scaling);
    message.set_adaptive_rho(_settings->adaptive_rho);

    message.set_delta(_settings->delta);
    message.set_polish(_settings->polish);
    message.set_polish_refine_iter(_settings->polish_refine_iter);

    switch (_settings->linsys_solver)
    {
        case QDLDL_SOLVER:
        {
            message.set_linear_solver(messages::SolverOsqp_LinearSolver_QDLDL);
            break;
        }
        case MKL_PARDISO_SOLVER:
        {
            message.set_linear_solver(messages::SolverOsqp_LinearSolver_MKL_PARADISO);
            break;
        }
    }
}

void SolverOsqp::fromMessage(const corbo::messages::SolverOsqp& message, std::stringstream* issues)
{
    _settings->max_iter           = message.max_iter();
    _settings->eps_abs            = message.eps_abs();
    _settings->eps_rel            = message.eps_rel();
    _settings->eps_prim_inf       = message.eps_prim_inf();
    _settings->eps_dual_inf       = message.eps_dual_inf();
    _settings->alpha              = message.alpha();
    _settings->verbose            = message.verbose();
    _settings->scaled_termination = message.scaled_termination();
    _settings->check_termination  = message.check_termination();
    _settings->warm_start         = message.warm_start();
    _settings->time_limit         = message.time_limit();

    _settings->rho          = message.rho();
    _settings->sigma        = message.sigma();
    _settings->scaling      = message.scaling();
    _settings->adaptive_rho = message.adaptive_rho();

    _settings->delta              = message.delta();
    _settings->polish             = message.polish();
    _settings->polish_refine_iter = message.polish_refine_iter();

    switch (message.linear_solver())
    {
        case messages::SolverOsqp_LinearSolver_QDLDL:
        {
            _settings->linsys_solver = QDLDL_SOLVER;
            break;
        }
        case messages::SolverOsqp_LinearSolver_MKL_PARADISO:
        {
            _settings->linsys_solver = MKL_PARDISO_SOLVER;
            break;
        }
    }
}
#endif

}  // namespace corbo

#endif  // OSQP
