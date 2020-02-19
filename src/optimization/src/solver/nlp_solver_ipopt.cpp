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

#include <corbo-optimization/solver/nlp_solver_ipopt.h>

#include <corbo-core/console.h>

namespace corbo {

bool SolverIpopt::initialize(OptimizationProblemInterface* /*problem*/)
{
    if (_initialized) return true;

    _ipopt_nlp = new IpoptWrapper(this);
    _ipopt_app = IpoptApplicationFactory();

    if (!_param_msg_received)  // TODO(roesmann) we need a default paramerter mechanism with protobuf!
    {
        setRelTolerance(1e-3);

        setMuStrategyAdaptive(true);
        setWarmStartInitPoint(true);
        setNlpAutoScaling(true);
        setPrintLevel(2);  // Default = 5

        /*
        #ifdef __unix__
                if (!setLinearSolver(LinearSolver::MA57)) setLinearSolver(LinearSolver::MUMPS);
        #else
                setLinearSolver(LinearSolver::MUMPS);  // the default solver on systems other than linux is MUMPS
                                                       // since setLinearSolver() returns true, even if MA57 is not available
        #endif
        */
        // _ipopt_app->Options()->SetStringValue("output_file", "ipopt.out");

        // _ipopt_app->Options()->SetStringValue("hessian_approximation", "exact"); // or "limited-memory"

        // _ipopt_app->Options()->SetStringValue("derivative_test", "second-order"); // "none", "first-order", "second-order", "only-second-order"
        // _ipopt_app->Options()->SetNumericValue("derivative_test_tol", 1e-3);

        Ipopt::ApplicationReturnStatus status;
        status = _ipopt_app->Initialize();
        if (status != Ipopt::Solve_Succeeded)
        {
            PRINT_INFO("SolverIPOPT(): Error during IPOPT initialization!");
            return false;
        }
    }

    static bool copyright_printed = false;
    if (!copyright_printed)
    {
        _ipopt_app->PrintCopyrightMessage();
        copyright_printed = true;
    }

    _initialized = true;
    return true;
}

SolverStatus SolverIpopt::solve(OptimizationProblemInterface& problem, bool new_structure, bool new_run, double* obj_value)
{
    if (!_initialized)
    {
        if (!initialize(&problem)) return SolverStatus::Error;
    }

    _ipopt_nlp->setOptimizationProblem(problem);

    if (new_structure)
    {
        _nnz_jac_constraints = problem.computeCombinedSparseJacobiansNNZ(false, true, true);

        problem.computeSparseHessiansNNZ(_nnz_hes_obj, _nnz_hes_eq, _nnz_hes_ineq, true);
        _nnz_h_lagrangian = _nnz_hes_obj + _nnz_hes_eq + _nnz_hes_ineq;

        _lambda_cache.resize(problem.getEqualityDimension() + problem.getInequalityDimension());
        _lambda_cache.setZero();

        _zl_cache.resize(problem.getParameterDimension());
        _zl_cache.setZero();

        _zu_cache.resize(problem.getParameterDimension());
        _zu_cache.setZero();

        // set max number of iterations
        _ipopt_app->Options()->SetIntegerValue("max_iter", _iterations);  // max_cpu_time // TODO(roesmann) parameter for number of iterations

        // Set bfgs method
        //         if (cfg->optim.solver.nonlin_prog.hessian.hessian_method == HessianMethod::BFGS)
        //             _ipopt_app->Options()->SetStringValue("hessian_approximation", "limited-memory");
        //         else
        //            _ipopt_app->Options()->SetStringValue("hessian_approximation", "exact");
    }

    if (_max_cpu_time > 0)
        _ipopt_app->Options()->SetNumericValue("max_cpu_time", _max_cpu_time);
    else if (_max_cpu_time == 0)
        _ipopt_app->Options()->SetNumericValue("max_cpu_time", 10e6);

    Ipopt::ApplicationReturnStatus ipopt_status;
    if (new_structure)
        ipopt_status = _ipopt_app->OptimizeTNLP(_ipopt_nlp);
    else
        ipopt_status = _ipopt_app->ReOptimizeTNLP(_ipopt_nlp);

    if (obj_value) *obj_value = _last_obj_value;

    return convertIpoptToNlpSolverStatus(ipopt_status);
}

bool SolverIpopt::setLinearSolver(LinearSolver solver_type)
{
    if (!Ipopt::IsValid(_ipopt_app)) return false;

    bool success = false;

    switch (solver_type)
    {
        case LinearSolver::MUMPS:
            success = _ipopt_app->Options()->SetStringValue("linear_solver", "mumps");
            break;
        case LinearSolver::MA27:
            success = _ipopt_app->Options()->SetStringValue("linear_solver", "ma27");
            break;
        case LinearSolver::MA57:
            success = _ipopt_app->Options()->SetStringValue("linear_solver", "ma57");
            break;
        case LinearSolver::MA77:
            success = _ipopt_app->Options()->SetStringValue("linear_solver", "ma77");
            break;
        case LinearSolver::MA86:
            success = _ipopt_app->Options()->SetStringValue("linear_solver", "ma86");
            break;
        case LinearSolver::MA97:
            success = _ipopt_app->Options()->SetStringValue("linear_solver", "ma97");
            break;
        default:
            success = false;
    }

    if (success)
    {
        _current_lin_solver = solver_type;
        return true;
    }
    return false;
}

SolverIpopt::LinearSolver SolverIpopt::getLinearSolver() const
{
    if (!Ipopt::IsValid(_ipopt_app)) return LinearSolver::NO_SOLVER;

    std::string linear_solver;
    if (!_ipopt_app->Options()->GetStringValue("linear_solver", linear_solver, "")) return LinearSolver::NO_SOLVER;

    if (!linear_solver.compare("mumps")) return LinearSolver::MUMPS;
    if (!linear_solver.compare("ma27")) return LinearSolver::MA27;
    if (!linear_solver.compare("ma57")) return LinearSolver::MA57;
    if (!linear_solver.compare("ma77")) return LinearSolver::MA77;
    if (!linear_solver.compare("ma86")) return LinearSolver::MA86;
    if (!linear_solver.compare("ma97")) return LinearSolver::MA97;

    return LinearSolver::NO_SOLVER;
}

bool SolverIpopt::setLinearSolverByName(const std::string& solver_name)
{
    if (!Ipopt::IsValid(_ipopt_app)) return false;
    return _ipopt_app->Options()->SetStringValue("linear_solver", solver_name);
}

std::string SolverIpopt::getLinearSolverByName()
{
    if (!Ipopt::IsValid(_ipopt_app)) return "";
    std::string linear_solver;
    if (!_ipopt_app->Options()->GetStringValue("linear_solver", linear_solver, "")) return "";
    return linear_solver;
}

bool SolverIpopt::setRelTolerance(double tolerance) { return _ipopt_app->Options()->SetNumericValue("tol", tolerance); }
double SolverIpopt::getRelTolerance() const
{
    double val = CORBO_INF_DBL;
    _ipopt_app->Options()->GetNumericValue("tol", val, "");
    return val;
}

bool SolverIpopt::setDualInfTolerance(double tolerance) { return _ipopt_app->Options()->SetNumericValue("dual_inf_tol", tolerance); }
double SolverIpopt::getDualInfTolerance() const
{
    double val = CORBO_INF_DBL;
    _ipopt_app->Options()->GetNumericValue("dual_inf_tol", val, "");
    return val;
}

bool SolverIpopt::setConstrViolTolerance(double tolerance) { return _ipopt_app->Options()->SetNumericValue("constr_viol_tol", tolerance); }
double SolverIpopt::getConstrViolTolerance() const
{
    double val = CORBO_INF_DBL;
    _ipopt_app->Options()->GetNumericValue("constr_viol_tol", val, "");
    return val;
}

bool SolverIpopt::setComplInfTolerance(double tolerance) { return _ipopt_app->Options()->SetNumericValue("compl_inf_tol", tolerance); }
double SolverIpopt::getComplInfTolerance() const
{
    double val = CORBO_INF_DBL;
    _ipopt_app->Options()->GetNumericValue("compl_inf_tol", val, "");
    return val;
}

bool SolverIpopt::setMuStrategyAdaptive(bool enabled)
{
    if (enabled) return _ipopt_app->Options()->SetStringValue("mu_strategy", "adaptive");
    return _ipopt_app->Options()->SetStringValue("mu_strategy", "monotone");
}
bool SolverIpopt::isMuStrategyAdaptive() const
{
    std::string opt;
    _ipopt_app->Options()->GetStringValue("mu_strategy", opt, "");
    return opt.compare("adaptive") == 0 ? true : false;
}

bool SolverIpopt::setHessianApproxExact(bool enabled)
{
    if (enabled) return _ipopt_app->Options()->SetStringValue("hessian_approximation", "exact");
    return _ipopt_app->Options()->SetStringValue("hessian_approximation", "limited-memory");
}
bool SolverIpopt::isHessianApproxExact() const
{
    std::string opt;
    _ipopt_app->Options()->GetStringValue("hessian_approximation", opt, "");
    return opt.compare("exact") == 0 ? true : false;
}

bool SolverIpopt::setWarmStartInitPoint(bool enabled)
{
    if (enabled) return _ipopt_app->Options()->SetStringValue("warm_start_init_point", "yes");
    return _ipopt_app->Options()->SetStringValue("warm_start_init_point", "no");
}
bool SolverIpopt::isWarmStartInitPoint() const
{
    std::string opt;
    _ipopt_app->Options()->GetStringValue("warm_start_init_point", opt, "");
    return opt.compare("yes") == 0 ? true : false;
}

bool SolverIpopt::setMehrotraAlgorithm(bool enabled)
{
    if (enabled) return _ipopt_app->Options()->SetStringValue("mehrotra_algorithm", "yes");
    return _ipopt_app->Options()->SetStringValue("mehrotra_algorithm", "no");
}
bool SolverIpopt::isMehrotraAlgorithm() const
{
    std::string opt;
    _ipopt_app->Options()->GetStringValue("mehrotra_algorithm", opt, "");
    return opt.compare("yes") == 0 ? true : false;
}

bool SolverIpopt::setPrintLevel(int print_level) { return _ipopt_app->Options()->SetIntegerValue("print_level", print_level); }  // level 0-5
int SolverIpopt::getPrintLevel() const
{
    int val = -1;
    _ipopt_app->Options()->GetIntegerValue("print_level", val, "");
    return val;
}

bool SolverIpopt::setNlpAutoScaling(bool enabled)
{
    if (enabled) return _ipopt_app->Options()->SetStringValue("nlp_scaling_method", "gradient-based");
    return _ipopt_app->Options()->SetStringValue("nlp_scaling_method", "none");
}
bool SolverIpopt::isNlpAutoScaling() const
{
    std::string opt;
    _ipopt_app->Options()->GetStringValue("nlp_scaling_method", opt, "");
    return opt.compare("gradient-based") == 0 ? true : false;
}

bool SolverIpopt::setCheckDerivativesForNan(bool enabled)
{
    if (enabled) return _ipopt_app->Options()->SetStringValue("check_derivatives_for_naninf", "yes");
    return _ipopt_app->Options()->SetStringValue("check_derivatives_for_naninf", "no");
}
bool SolverIpopt::isCheckDerivativesForNan() const
{
    std::string opt;
    _ipopt_app->Options()->GetStringValue("check_derivatives_for_naninf", opt, "");
    return opt.compare("yes") == 0 ? true : false;
}

bool SolverIpopt::setDerivativeTest(bool first_order, bool second_order)
{
    _ipopt_app->Options()->SetNumericValue("derivative_test_perturbation", 6e-3);
    _ipopt_app->Options()->SetNumericValue("derivative_test_tol", 1e-3);

    if (first_order && second_order)
        return _ipopt_app->Options()->SetStringValue("derivative_test", "second-order");
    else if (first_order)
        return _ipopt_app->Options()->SetStringValue("derivative_test", "first-order");
    else if (second_order)
        return _ipopt_app->Options()->SetStringValue("derivative_test", "only-second-order");
    return _ipopt_app->Options()->SetStringValue("derivative_test", "none");
}
void SolverIpopt::isDerivativeTest(bool& first_order, bool& second_order) const
{
    std::string opt;
    _ipopt_app->Options()->GetStringValue("derivative_test", opt, "");
    if (opt.compare("second-order") == 0)
    {
        first_order  = true;
        second_order = true;
    }
    else if (opt.compare("first-order") == 0)
    {
        first_order  = true;
        second_order = false;
    }
    else if (opt.compare("fonly-second-order") == 0)
    {
        first_order  = false;
        second_order = true;
    }
    else
    {
        first_order  = false;
        second_order = false;
    }
}

bool SolverIpopt::setIpoptOptionString(const std::string& param, const std::string& option)
{
    if (!Ipopt::IsValid(_ipopt_app)) return false;
    return _ipopt_app->Options()->SetStringValue(param, option);
}

bool SolverIpopt::setIpoptOptionInt(const std::string& param, int option)
{
    if (!Ipopt::IsValid(_ipopt_app)) return false;
    return _ipopt_app->Options()->SetIntegerValue(param, option);
}

bool SolverIpopt::setIpoptOptionNumeric(const std::string& param, double option)
{
    if (!Ipopt::IsValid(_ipopt_app)) return false;
    return _ipopt_app->Options()->SetNumericValue(param, option);
}

SolverStatus SolverIpopt::convertIpoptToNlpSolverStatus(Ipopt::ApplicationReturnStatus ipopt_status) const
{
    switch (ipopt_status)
    {
        case Ipopt::Solve_Succeeded:
        case Ipopt::Solved_To_Acceptable_Level:
        {
            return SolverStatus::Converged;
        }
        case Ipopt::Search_Direction_Becomes_Too_Small:
        case Ipopt::User_Requested_Stop:
        case Ipopt::Feasible_Point_Found:
        case Ipopt::Maximum_Iterations_Exceeded:
        {
            return SolverStatus::EarlyTerminated;
        }
        case Ipopt::Infeasible_Problem_Detected:
        case Ipopt::Diverging_Iterates:
        case Ipopt::Restoration_Failed:
        case Ipopt::Not_Enough_Degrees_Of_Freedom:
        case Ipopt::Invalid_Problem_Definition:
        {
            return SolverStatus::Infeasible;
        }
        default:
        {
        }
    };
    return SolverStatus::Error;
}

void SolverIpopt::clear()
{
    _last_obj_value = -1;

    _nnz_jac_constraints = 0;
    _nnz_h_lagrangian    = 0;
    _nnz_hes_obj         = 0;
    _nnz_hes_eq          = 0;
    _nnz_hes_ineq        = 0;

    // the initialize() function does not fill any problem dependent caches, so let's keep it initialized...
    // _initialized = false;
    // _ipopt_nlp = nullptr;
    // _ipopt_app = nullptr;
}

#ifdef MESSAGE_SUPPORT
void SolverIpopt::toMessage(corbo::messages::NlpSolver& message) const
{
    message.mutable_solver_ipopt()->set_iterations(getIterations());
    message.mutable_solver_ipopt()->set_max_cpu_time(getMaxCpuTime());

    switch (getLinearSolver())
    {
        case LinearSolver::MUMPS:
        {
            message.mutable_solver_ipopt()->set_linear_solver(corbo::messages::SolverIpopt::MUMPS);
            break;
        }
        case LinearSolver::MA27:
        {
            message.mutable_solver_ipopt()->set_linear_solver(corbo::messages::SolverIpopt::MA27);
            break;
        }
        case LinearSolver::MA57:
        {
            message.mutable_solver_ipopt()->set_linear_solver(corbo::messages::SolverIpopt::MA57);
            break;
        }
        case LinearSolver::MA77:
        {
            message.mutable_solver_ipopt()->set_linear_solver(corbo::messages::SolverIpopt::MA77);
            break;
        }
        case LinearSolver::MA86:
        {
            message.mutable_solver_ipopt()->set_linear_solver(corbo::messages::SolverIpopt::MA86);
            break;
        }
        case LinearSolver::MA97:
        {
            message.mutable_solver_ipopt()->set_linear_solver(corbo::messages::SolverIpopt::MA97);
            break;
        }
        default:
            PRINT_ERROR("SolverIpopt::toMessage(): selected linear solver is currently not included in message export.");
    }

    message.mutable_solver_ipopt()->set_rel_tolerance(getRelTolerance());
    message.mutable_solver_ipopt()->set_adaptive_mu_strategy(isMuStrategyAdaptive());
    message.mutable_solver_ipopt()->set_exact_hessian_approx(isHessianApproxExact());
    message.mutable_solver_ipopt()->set_mehrotra_algorithm(isMehrotraAlgorithm());
    message.mutable_solver_ipopt()->set_warm_start_init_point(isWarmStartInitPoint());
    message.mutable_solver_ipopt()->set_nlp_auto_scaling(isNlpAutoScaling());
    message.mutable_solver_ipopt()->set_print_level(getPrintLevel());
    message.mutable_solver_ipopt()->set_cache_first_order_derivatives(_cache_first_order_derivatives);
    message.mutable_solver_ipopt()->set_check_derivatives_for_naninf(isCheckDerivativesForNan());
    bool check_first_order  = false;
    bool check_second_order = false;
    isDerivativeTest(check_first_order, check_second_order);
    message.mutable_solver_ipopt()->set_derivative_test_first_order(check_first_order);
    message.mutable_solver_ipopt()->set_derivative_test_second_order(check_second_order);
}

void SolverIpopt::fromMessage(const corbo::messages::NlpSolver& message, std::stringstream* issues)
{
    if (!_initialized) initialize();  // we need valid objects // TODO(roesmann): this is not nice and not consistent in the project

    setIterations(message.solver_ipopt().iterations());
    setMaxCpuTime(message.solver_ipopt().max_cpu_time());

    switch (message.solver_ipopt().linear_solver())
    {
        case corbo::messages::SolverIpopt::MUMPS:
        {
            if (!setLinearSolver(LinearSolver::MUMPS) && issues)
                *issues << "LinearSolver::MUMPS: not supported with current IPOPT installation." << std::endl;
            break;
        }
        case corbo::messages::SolverIpopt::MA27:
        {
            if (!setLinearSolver(LinearSolver::MA27) && issues)
                *issues << "LinearSolver::MA27: not supported with current IPOPT installation." << std::endl;
            break;
        }
        case corbo::messages::SolverIpopt::MA57:
        {
            if (!setLinearSolver(LinearSolver::MA57) && issues)
                *issues << "LinearSolver::MA57: not supported with current IPOPT installation." << std::endl;
            break;
        }
        case corbo::messages::SolverIpopt::MA77:
        {
            if (!setLinearSolver(LinearSolver::MA77) && issues)
                *issues << "LinearSolver::MA77: not supported with current IPOPT installation." << std::endl;
            break;
        }
        case corbo::messages::SolverIpopt::MA86:
        {
            if (!setLinearSolver(LinearSolver::MA86) && issues)
                *issues << "LinearSolver::MA86: not supported with current IPOPT installation" << std::endl;
            break;
        }
        case corbo::messages::SolverIpopt::MA97:
        {
            if (!setLinearSolver(LinearSolver::MA97) && issues)
                *issues << "LinearSolver::MA97: not supported with current IPOPT installation." << std::endl;
            break;
        }
        default:
        {
            if (issues) *issues << "SolverIpopt::fromMessage(): selected linear solver is currently not supported." << std::endl;
        }
    }

    if (!setRelTolerance(message.solver_ipopt().rel_tolerance()) && issues) *issues << "setRelTolerance() failed." << std::endl;
    if (!setMuStrategyAdaptive(message.solver_ipopt().adaptive_mu_strategy()) && issues) *issues << "setMuStrategyAdaptive() failed." << std::endl;
    if (!setHessianApproxExact(message.solver_ipopt().exact_hessian_approx()) && issues) *issues << "setHessianApproxExact() failed." << std::endl;
    if (!setMehrotraAlgorithm(message.solver_ipopt().mehrotra_algorithm()) && issues) *issues << "setMehrotraAlgorithm() failed." << std::endl;
    if (!setWarmStartInitPoint(message.solver_ipopt().warm_start_init_point()) && issues) *issues << "setWarmStartInitPoint() failed." << std::endl;
    if (!setNlpAutoScaling(message.solver_ipopt().nlp_auto_scaling()) && issues) *issues << "setNlpAutoScaling() failed." << std::endl;
    if (!setPrintLevel(message.solver_ipopt().print_level()) && issues) *issues << "setPrintLevel() failed." << std::endl;

    if (!setCheckDerivativesForNan(message.solver_ipopt().check_derivatives_for_naninf()) && issues)
        *issues << "setCheckDerivativesForNan() failed." << std::endl;

    bool check_first_order  = message.solver_ipopt().derivative_test_first_order();
    bool check_second_order = message.solver_ipopt().derivative_test_second_order();
    if (!setDerivativeTest(check_first_order, check_second_order) && issues) *issues << "setDerivativeTest() failed." << std::endl;

    _cache_first_order_derivatives = message.solver_ipopt().cache_first_order_derivatives();

    Ipopt::ApplicationReturnStatus status;
    status = _ipopt_app->Initialize();
    if (status != Ipopt::Solve_Succeeded)
    {
        if (issues)
            *issues
                << "SolverIpopt::fromMessage(): Error during IPOPT initialization! Maybe some parameters are not supported with current installation";
        PRINT_WARNING(
            "SolverIpopt::fromMessage(): Error during IPOPT initialization! Maybe some parameters are not supported with current installation");
    }

    _param_msg_received = true;
}
#endif

}  // namespace corbo

#endif  // IPOPT
