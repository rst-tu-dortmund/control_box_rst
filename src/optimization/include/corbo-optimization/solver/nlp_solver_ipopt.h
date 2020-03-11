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

#ifndef SRC_OPTIMIZATION_INCLUDE_CORBO_OPTIMIZATION_SOLVER_NLP_SOLVER_IPOPT_H_
#define SRC_OPTIMIZATION_INCLUDE_CORBO_OPTIMIZATION_SOLVER_NLP_SOLVER_IPOPT_H_

#include <corbo-optimization/solver/nlp_solver_interface.h>

#include <memory>

#ifdef IPOPT
#include <corbo-optimization/solver/nlp_solver_ipopt_wrapper.h>
#include <IpIpoptApplication.hpp>
#endif

namespace corbo {

#ifdef IPOPT

/**
 * @brief Interface to the external interior point solver IPOPT
 *
 * @ingroup optimization solver
 *
 * @author Christoph Rösmann (christoph.roesmann@tu-dortmund.de)
 */
class SolverIpopt : public NlpSolverInterface
{
    friend class IpoptWrapper;

 public:
    using Ptr  = std::shared_ptr<SolverIpopt>;
    using UPtr = std::unique_ptr<SolverIpopt>;

    enum class LinearSolver { MUMPS, MA27, MA57, MA77, MA86, MA97, NO_SOLVER };

    // implements interface method
    NlpSolverInterface::Ptr getInstance() const override { return std::make_shared<SolverIpopt>(); }

    // implements interface method
    bool isLsqSolver() const override { return false; }

    // implements interface method
    bool initialize(OptimizationProblemInterface* problem = nullptr) override;
    // implements interface method
    SolverStatus solve(OptimizationProblemInterface& problem, bool new_structure = true, bool new_run = true, double* obj_value = nullptr) override;

    /**@name Set solver properties */
    //@{

    bool setLinearSolver(LinearSolver solver_type);
    LinearSolver getLinearSolver() const;

    // alternative method to directly pass sovler_name to ipopt
    bool setLinearSolverByName(const std::string& solver_name);
    std::string getLinearSolverByName();

    bool setRelTolerance(double tolerance);
    double getRelTolerance() const;

    bool setDualInfTolerance(double tolerance);
    double getDualInfTolerance() const;

    bool setConstrViolTolerance(double tolerance);
    double getConstrViolTolerance() const;

    bool setComplInfTolerance(double tolerance);
    double getComplInfTolerance() const;

    bool setMuStrategyAdaptive(bool enabled);
    bool isMuStrategyAdaptive() const;

    bool setHessianApproxExact(bool enabled);
    bool isHessianApproxExact() const;

    bool setWarmStartInitPoint(bool enabled);
    bool isWarmStartInitPoint() const;

    bool setMehrotraAlgorithm(bool enabled);
    bool isMehrotraAlgorithm() const;

    bool setPrintLevel(int print_level);  // level 0-5
    int getPrintLevel() const;

    bool setNlpAutoScaling(bool enabled);
    bool isNlpAutoScaling() const;

    void setIterations(int iterations) { _iterations = iterations; }
    int getIterations() const { return _iterations; }

    void setMaxCpuTime(double max_cpu_time) { _max_cpu_time = max_cpu_time; }
    double getMaxCpuTime() const { return _max_cpu_time; }

    bool setCheckDerivativesForNan(bool enabled);
    bool isCheckDerivativesForNan() const;

    bool setDerivativeTest(bool first_order, bool second_order);
    void isDerivativeTest(bool& first_order, bool& second_order) const;

    // generic setter methods
    bool setIpoptOptionString(const std::string& param, const std::string& option);
    bool setIpoptOptionInt(const std::string& param, int option);
    bool setIpoptOptionNumeric(const std::string& param, double option);

    void setCacheFirstOrderDerivatives(bool active) { _cache_first_order_derivatives = active; }

    //@}

    void clear() override;

    // void setIterations(int iterations) { _iterations = iterations; }

#ifdef MESSAGE_SUPPORT
    // implements interface method
    void toMessage(corbo::messages::NlpSolver& message) const override;
    // implements interface method
    void fromMessage(const corbo::messages::NlpSolver& message, std::stringstream* issues = nullptr) override;
#endif

 protected:
    SolverStatus convertIpoptToNlpSolverStatus(Ipopt::ApplicationReturnStatus ipopt_status) const;

 private:
    Ipopt::SmartPtr<IpoptWrapper> _ipopt_nlp;
    Ipopt::SmartPtr<Ipopt::IpoptApplication> _ipopt_app;

    int _nnz_jac_constraints = 0;
    int _nnz_h_lagrangian    = 0;
    int _nnz_hes_obj         = 0;
    int _nnz_hes_eq          = 0;
    int _nnz_hes_ineq        = 0;

    Eigen::VectorXd _lambda_cache;  //! store constraint mupltipliers between subsequent IPOPT calls
    Eigen::VectorXd _zl_cache;      //! store lower bound multipliers between subsequent IPOPT calls
    Eigen::VectorXd _zu_cache;      //! store upper bound multipliers between subsequent IPOPT calls

    Eigen::VectorXd _grad_f_cache;
    Eigen::VectorXd _jac_constr_cache;

    // std::vector<double> _lambda; // backup lagrange multiplier values between calls;

    double _max_cpu_time = -1;   // use ipopt default (0: deactivate) [but should be deactivated by default]
    int _iterations      = 100;  // ipopt max iterations

    double _last_obj_value = -1;

    bool _cache_first_order_derivatives = false;

    LinearSolver _current_lin_solver = LinearSolver::NO_SOLVER;

    bool _initialized = false;

    int _param_msg_received = false;
};

FACTORY_REGISTER_NLP_SOLVER(SolverIpopt)

#else  // IPOPT

/**
 * @brief Interface to the external interior point solver IPOPT
 *
 * @ingroup optimization solver
 *
 * @author Christoph Rösmann (christoph.roesmann@tu-dortmund.de)
 */
class SolverIpopt : public NlpSolverInterface
{
 public:
    using Ptr  = std::shared_ptr<SolverIpopt>;
    using UPtr = std::unique_ptr<SolverIpopt>;

    enum class LinearSolver { MUMPS, MA27, MA57, MA77, MA86, MA97, NO_SOLVER };

    // implements interface method
    NlpSolverInterface::Ptr getInstance() const override { return std::make_shared<SolverIpopt>(); }

    // implements interface method
    bool isLsqSolver() const override { return false; }

    // implements interface method
    bool initialize(OptimizationProblemInterface* problem = nullptr) override
    {
        printWarning();
        return false;
    }
    // implements interface method
    SolverStatus solve(OptimizationProblemInterface& problem, bool new_structure = true, bool new_run = true, double* obj_value = nullptr) override
    {
        printWarning();
        return SolverStatus::Error;
    }

    void clear() override {}

#ifdef MESSAGE_SUPPORT
    // implements interface method
    void toMessage(corbo::messages::NlpSolver& /*message*/) const override { printWarning(); }
    // implements interface method
    void fromMessage(const corbo::messages::NlpSolver& /*message*/, std::stringstream* issues) override
    {
        if (issues) *issues << "SolverIPOPT cannot be selected since it is not installed properly." << std::endl;
        printWarning();
    }
#endif

 protected:
    void printWarning() const { PRINT_WARNING("SolverIPOPT cannot be selected since it is not installed properly."); }
};

#endif  // IPOPT

}  // namespace corbo

#endif  // SRC_OPTIMIZATION_INCLUDE_CORBO_OPTIMIZATION_SOLVER_NLP_SOLVER_IPOPT_H_
