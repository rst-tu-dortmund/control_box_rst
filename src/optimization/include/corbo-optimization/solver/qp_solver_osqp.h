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

#ifndef SRC_OPTIMIZATION_INCLUDE_CORBO_OPTIMIZATION_SOLVER_QP_SOLVER_OSQP_H_
#define SRC_OPTIMIZATION_INCLUDE_CORBO_OPTIMIZATION_SOLVER_QP_SOLVER_OSQP_H_

#include <corbo-optimization/solver/qp_solver_interface.h>

#include <memory>

#ifdef OSQP
#include <osqp.h>
#endif

namespace corbo {

#ifdef OSQP

/**
 * @brief Interface to the external solver OSQP
 *
 * @ingroup optimization solver
 *
 * @author Christoph Rösmann (christoph.roesmann@tu-dortmund.de)
 */
class SolverOsqp : public QpSolverInterface
{
 public:
    using Ptr  = std::shared_ptr<SolverOsqp>;
    using UPtr = std::unique_ptr<SolverOsqp>;

    // Default constructor
    SolverOsqp();

    // Destructor
    ~SolverOsqp();

    // implements interface method
    QpSolverInterface::Ptr getInstance() const override { return std::make_shared<SolverOsqp>(); }

    // implements interface method
    bool initialize() override;

    // implements interface method
    bool isSupportingSimpleBounds() override { return false; }

    // implements interface method
    SolverStatus solve(SparseMatrix& P, Eigen::Ref<Eigen::VectorXd> q, SparseMatrix& A, Eigen::Ref<Eigen::VectorXd> lbA,
                       Eigen::Ref<Eigen::VectorXd> ubA, Eigen::Ref<Eigen::VectorXd> lb, Eigen::Ref<Eigen::VectorXd> ub, bool new_structure = true,
                       bool zero_x_warmstart = false) override;

    // implements interface method
    SolverStatus solve(SparseMatrix& P, Eigen::Ref<Eigen::VectorXd> q, SparseMatrix& A, Eigen::Ref<Eigen::VectorXd> lbA,
                       Eigen::Ref<Eigen::VectorXd> ubA, bool new_structure = true, bool zero_x_warmstart = false, bool update_P = true,
                       bool update_q = true, bool update_A = true, bool update_bounds = true) override;

    Eigen::Ref<Eigen::VectorXd> getPrimalSolution() override;
    Eigen::Ref<Eigen::VectorXd> getDualSolution() override;

    void updatePrimalSolutionWarmStart(const Eigen::Ref<const Eigen::VectorXd>& x) override;
    void updateDualSolutionWarmStart(const Eigen::Ref<const Eigen::VectorXd>& y) override;

    void enforceNewStructure(bool new_structure = true) override { _force_new_structure = new_structure; }

    /**@name Set solver properties */
    //@{
    OSQPSettings* getOsqpSettings() { return _settings.get(); }
    //@}

    void clear() override;

#ifdef MESSAGE_SUPPORT
    //! Export to message
    void toMessage(corbo::messages::SolverOsqp& message) const;
    //! Import form message
    void fromMessage(const corbo::messages::SolverOsqp& message, std::stringstream* issues = nullptr);

    // implements interface method
    void toMessage(corbo::messages::QpSolver& message) const override { toMessage(*message.mutable_solver_osqp()); }
    // implements interface method
    void fromMessage(const corbo::messages::QpSolver& message, std::stringstream* issues = nullptr) override
    {
        fromMessage(message.solver_osqp(), issues);
    }
#endif

 protected:
    SolverStatus convertOsqpExitFlagToSolverStatus(c_int status) const;

 private:
    std::unique_ptr<OSQPSettings> _settings;
    // std::unique_ptr<OSQPData> _data;
    OSQPWorkspace* _work = nullptr;

    Eigen::VectorXd _zero;

    bool _force_new_structure = false;

    bool _initialized = false;
};

FACTORY_REGISTER_QP_SOLVER(SolverOsqp)

#else  // OSQP

/**
 * @brief Interface to the external solver OSQP
 *
 * @ingroup optimization solver
 *
 * @author Christoph Rösmann (christoph.roesmann@tu-dortmund.de)
 */
class SolverOsqp : public QpSolverInterface
{
 public:
    using Ptr  = std::shared_ptr<SolverOsqp>;
    using UPtr = std::unique_ptr<SolverOsqp>;

    // implements interface method
    QpSolverInterface::Ptr getInstance() const override { return std::make_shared<SolverOsqp>(); }

    // implements interface method
    bool isSupportingSimpleBounds() override { return false; }

    // implements interface method
    SolverStatus solve(SparseMatrix& P, Eigen::Ref<Eigen::VectorXd> q, SparseMatrix& A, Eigen::Ref<Eigen::VectorXd> lbA,
                       Eigen::Ref<Eigen::VectorXd> ubA, Eigen::Ref<Eigen::VectorXd> lb, Eigen::Ref<Eigen::VectorXd> ub, bool new_structure = true,
                       bool zero_x_warmstart = false) override
    {
        printWarning();
        return SolverStatus::Error;
    }

    // implements interface method
    SolverStatus solve(SparseMatrix& P, Eigen::Ref<Eigen::VectorXd> q, SparseMatrix& A, Eigen::Ref<Eigen::VectorXd> lbA,
                       Eigen::Ref<Eigen::VectorXd> ubA, bool new_structure = true, bool zero_x_warmstart = false, bool update_P = true,
                       bool update_q = true, bool update_A = true, bool update_bounds = true) override
    {
        printWarning();
        return SolverStatus::Error;
    }

    Eigen::Ref<Eigen::VectorXd> getPrimalSolution() override { return _dummy; }
    Eigen::Ref<Eigen::VectorXd> getDualSolution() override { return _dummy; }

    void updatePrimalSolutionWarmStart(const Eigen::Ref<const Eigen::VectorXd>& x) override {}
    void updateDualSolutionWarmStart(const Eigen::Ref<const Eigen::VectorXd>& y) override {}

    void enforceNewStructure(bool new_structure = true) override { printWarning(); }

    void clear() override {}

#ifdef MESSAGE_SUPPORT
    // implements interface method
    void toMessage(corbo::messages::QpSolver& /*message*/) const override { printWarning(); }
    // implements interface method
    void fromMessage(const corbo::messages::QpSolver& /*message*/, std::stringstream* issues) override
    {
        if (issues) *issues << "SolverOsqp cannot be selected since it is not installed properly." << std::endl;
        printWarning();
    }
#endif

 protected:
    void printWarning() const { PRINT_WARNING("SolverOsqp cannot be selected since it is not installed properly."); }

 private:
    Eigen::VectorXd _dummy;
};

#endif  // OSQP

}  // namespace corbo

#endif  // SRC_OPTIMIZATION_INCLUDE_CORBO_OPTIMIZATION_SOLVER_QP_SOLVER_OSQP_H_
