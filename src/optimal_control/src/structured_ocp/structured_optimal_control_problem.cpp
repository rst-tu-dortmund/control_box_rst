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

#include <corbo-optimal-control/structured_ocp/structured_optimal_control_problem.h>

#include <corbo-optimal-control/structured_ocp/discretization_grids/discretization_grid_interface.h>

#include <corbo-optimization/hyper_graph/generic_edge.h>

#include <corbo-communication/utilities.h>
#include <corbo-numerics/algebraic_riccati_continuous.h>
#include <corbo-numerics/algebraic_riccati_discrete.h>
#include <corbo-numerics/integrator_interface.h>
#include <corbo-numerics/matrix_utilities.h>

#include <memory>

namespace corbo {

StructuredOptimalControlProblem::StructuredOptimalControlProblem() {}

StructuredOptimalControlProblem::StructuredOptimalControlProblem(DiscretizationGridInterface::Ptr grid, SystemDynamicsInterface::Ptr dynamics,
                                                                 BaseHyperGraphOptimizationProblem::Ptr optim_prob, NlpSolverInterface::Ptr solver)
    : _grid(grid), _optim_prob(optim_prob), _dynamics(dynamics), _solver(solver)
{
    _optim_prob->setGraph(_edges, _grid);
}

bool StructuredOptimalControlProblem::providesFutureControls() const { return true; }

bool StructuredOptimalControlProblem::providesFutureStates() const { return _grid ? _grid->providesStateTrajectory() : false; }

bool StructuredOptimalControlProblem::initialize()
{
    if (!_optim_prob)
    {
        PRINT_ERROR("StructuredOptimalControlProblem::initialize(): no hyper-graph optimization problem strategy specified.");
        return false;
    }
    if (!_solver || !_solver->initialize(_optim_prob.get()))
    {
        PRINT_ERROR("StructuredOptimalControlProblem::initialize(): no solver specified or solver initialization failed.");
        return false;
    }

    if (_u_prev.size() == 0)
    {
        // Default setting for previous control and dt if nothing specified
        _u_prev.setZero(_dynamics->getInputDimension());
        _u_prev_dt = _grid->getInitialDt();
    }

    return true;
}

bool StructuredOptimalControlProblem::compute(const StateVector& x, ReferenceTrajectoryInterface& xref, ReferenceTrajectoryInterface& uref,
                                              ReferenceTrajectoryInterface* sref, const Time& t, bool new_run, SignalTargetInterface* signal_target,
                                              ReferenceTrajectoryInterface* xinit, ReferenceTrajectoryInterface* uinit, const std::string& ns)
{
    if (!_grid)
    {
        PRINT_ERROR("StructuredOptimalControlProblem::compute(): no discretization grid specified.");
        return false;
    }

    if (!_dynamics)
    {
        PRINT_ERROR("StructuredOptimalControlProblem::compute(): no system dynamics model specified.");
        return false;
    }
    if (!_optim_prob)
    {
        PRINT_ERROR("StructuredOptimalControlProblem::compute(): no hyper-graph optimization strategy specified.");
        return false;
    }
    if (!_solver)
    {
        PRINT_ERROR("StructuredOptimalControlProblem::compute(): no solver specified.");
        return false;
    }
    if (!_functions.stage_cost && !_functions.final_stage_cost)
    {
        PRINT_WARNING("StructuredOptimalControlProblem::compute(): no cost function specified.");
    }

    bool success = false;

    // reset ts caches
    _ts_x_cache.reset();
    _ts_u_cache.reset();
    _ts_dt_cache = 0;

    if (_statistics) _statistics->clear();

    Time t1 = Time::now();

    GridUpdateResult grid_udpate_result =
        _grid->update(x, xref, uref, _functions, *_edges, _dynamics, new_run, t, sref, &_u_prev, _u_prev_dt, xinit, uinit);

    if (grid_udpate_result.vertices_updated)
    {
        _optim_prob->precomputeVertexQuantities();
    }
    if (grid_udpate_result.updated())
    {
        _optim_prob->precomputeEdgeQuantities();
    }

    assert(_optim_prob->getGraph().checkGraphConsistency());

    Time t2 = Time::now();

    SolverStatus status = _solver->solve(*_optim_prob, grid_udpate_result.updated(), new_run, &_objective_value);
    if (status == SolverStatus::Converged || status == SolverStatus::EarlyTerminated)
        success = true;
    else if (_increase_n_if_infeas)
    {
        PRINT_WARNING("infeasible solution found. Increasing n for next ocp iteration.");
        _grid->setN(_grid->getN() + 1);
    }

    Time t3 = Time::now();

    if (_statistics)
    {
        _statistics->preparation_time = t2 - t1;
        _statistics->solving_time     = t3 - t2;
    }

    // PRINT_INFO("CPU time of only the solving phase: " << (t2 - t1).toSec() * 1000.0 << " ms.");

    return success;
}

bool StructuredOptimalControlProblem::getFirstControlInput(ControlVector& u0) const
{
    if (!_grid) return false;
    if (_grid->getFirstControlInput(u0))
        return true;
    else
        return false;
}

void StructuredOptimalControlProblem::setBounds(const Eigen::VectorXd& x_lb, const Eigen::VectorXd& x_ub, const Eigen::VectorXd& u_lb,
                                                const Eigen::VectorXd& u_ub)
{
    _functions.x_lb = x_lb;
    _functions.x_ub = x_ub;
    _functions.u_lb = u_lb;
    _functions.u_ub = u_ub;
    // TODO(roesmann): ocp modified?! we just changed bounds and not dimensions
}

void StructuredOptimalControlProblem::setStateBounds(const Eigen::VectorXd& x_lb, const Eigen::VectorXd& x_ub)
{
    _functions.x_lb = x_lb;
    _functions.x_ub = x_ub;
    // TODO(roesmann): ocp modified?! we just changed bounds and not dimensions
}

void StructuredOptimalControlProblem::setControlBounds(const Eigen::VectorXd& u_lb, const Eigen::VectorXd& u_ub)
{
    _functions.u_lb = u_lb;
    _functions.u_ub = u_ub;
    // TODO(roesmann): ocp modified?! we just changed bounds and not dimensions
}

void StructuredOptimalControlProblem::getTimeSeries(TimeSeries::Ptr x_sequence, TimeSeries::Ptr u_sequence, double t_max)
{
    if (!_grid)
    {
        PRINT_ERROR_NAMED("No grid loaded.");
        return;
    }
    _grid->getStateAndControlTimeSeries(x_sequence, u_sequence);
}

void StructuredOptimalControlProblem::reset()
{
    if (_grid) _grid->clear();
    if (_optim_prob) _optim_prob->clear();
    if (_dynamics) _dynamics->reset();
    if (_solver) _solver->clear();
    if (_statistics) _statistics->clear();

    _u_prev.setZero();  // TODO(roesmann): should we set this to zero here=

    _ocp_modified    = true;
    _objective_value = -1;
}

#ifdef MESSAGE_SUPPORT
void StructuredOptimalControlProblem::toMessage(corbo::messages::OptimalControlProblem& message) const
{
    messages::StructuredOptimalControlProblem* msg = message.mutable_structured_optimal_control_problem();

    if (_dynamics) _dynamics->toMessage(*msg->mutable_system_model());
    if (_optim_prob) _optim_prob->toMessage(*msg->mutable_hyper_graph_strategy());
    if (_solver) _solver->toMessage(*msg->mutable_nlp_solver());
    if (_grid) _grid->toMessage(*msg->mutable_discretization_grid());
    if (_functions.stage_cost) _functions.stage_cost->toMessage(*msg->mutable_stage_cost());
    if (_functions.final_stage_cost) _functions.final_stage_cost->toMessage(*msg->mutable_final_stage_cost());
    if (_functions.stage_equalities) _functions.stage_equalities->toMessage(*msg->mutable_stage_equalities());
    if (_functions.stage_inequalities) _functions.stage_inequalities->toMessage(*msg->mutable_stage_inequalities());
    if (_functions.final_stage_constraints) _functions.final_stage_constraints->toMessage(*msg->mutable_final_stage_constraints());
    if (_functions.stage_preprocessor) _functions.stage_preprocessor->toMessage(*msg->mutable_stage_preprocessors());

    msg->set_increase_n_if_infeas(_increase_n_if_infeas);

    // bounds
    msg->mutable_x_min()->Resize(_functions.x_lb.size(), 0);
    Eigen::Map<Eigen::VectorXd>(msg->mutable_x_min()->mutable_data(), _functions.x_lb.size()) = _functions.x_lb;
    msg->mutable_x_max()->Resize(_functions.x_ub.size(), 0);
    Eigen::Map<Eigen::VectorXd>(msg->mutable_x_max()->mutable_data(), _functions.x_ub.size()) = _functions.x_ub;
    msg->mutable_u_min()->Resize(_functions.u_lb.size(), 0);
    Eigen::Map<Eigen::VectorXd>(msg->mutable_u_min()->mutable_data(), _functions.u_lb.size()) = _functions.u_lb;
    msg->mutable_u_max()->Resize(_functions.u_ub.size(), 0);
    Eigen::Map<Eigen::VectorXd>(msg->mutable_u_max()->mutable_data(), _functions.u_ub.size()) = _functions.u_ub;
}

void StructuredOptimalControlProblem::fromMessage(const corbo::messages::OptimalControlProblem& message, std::stringstream* issues)
{
    _dynamics.reset();
    _solver.reset();
    _grid.reset();
    _optim_prob.reset();

    reset();
    //    _graph.clear();

    const messages::StructuredOptimalControlProblem& msg = message.structured_optimal_control_problem();

    // system model
    if (!msg.has_system_model())
    {
        if (issues) *issues << "StructuredOptimalControlProblem: no system model specified.\n";
        return;
    }
    // construct object
    std::string type;
    util::get_oneof_field_type_expand_isolated(msg.system_model(), "system_dynamics", type, false, 1);
    SystemDynamicsInterface::Ptr dynamics = SystemDynamicsFactory::instance().create(type);
    // import parameters
    if (dynamics)
    {
        dynamics->fromMessage(msg.system_model(), issues);
        setSystemDynamics(dynamics);
    }
    else
    {
        if (issues) *issues << "StructuredOptimalControlProblem: unknown system model specified.\n";
        return;
    }
    int dim_x = dynamics->getStateDimension();
    int dim_u = dynamics->getInputDimension();

    // grid
    if (!msg.has_discretization_grid())
    {
        if (issues) *issues << "StructuredOptimalControlProblem: no discretization grid specified.\n";
        return;
    }
    // construct object
    util::get_oneof_field_type(msg.discretization_grid(), "discretization_grid", type, false);
    DiscretizationGridInterface::Ptr grid = DiscretizationGridFactory::instance().create(type);
    // import parameters
    if (grid)
    {
        grid->fromMessage(msg.discretization_grid(), issues);
        setDiscretizationGrid(grid);
    }
    else
    {
        if (issues) *issues << "StructuredOptimalControlProblem: unknown discretization grid specified.\n";
        return;
    }

    // stage cost
    if (msg.has_stage_cost() && !msg.stage_cost().has_no_function())
    {
        // construct object
        util::get_oneof_field_type_expand_isolated(msg.stage_cost(), "stage_cost", type, false, 1);
        StageCost::Ptr stage_cost = StageCostFactory::instance().create(type);
        // import parameters
        if (stage_cost)
        {
            if (!stage_cost->fromMessage(msg.stage_cost(), issues)) return;
            if (!stage_cost->checkParameters(_dynamics->getStateDimension(), _dynamics->getInputDimension(), issues)) return;
            setStageCost(stage_cost);
        }
        else
        {
            if (issues) *issues << "StructuredOptimalControlProblem: unknown stage_cost specified.\n";
            return;
        }
    }
    else
        setStageCost({});

    // final stage cost
    if (msg.has_final_stage_cost() && !msg.final_stage_cost().has_no_function())
    {
        // construct object
        util::get_oneof_field_type_expand_isolated(msg.final_stage_cost(), "final_stage_cost", type, false, 1);
        FinalStageCost::Ptr final_stage_cost = FinalStageCostFactory::instance().create(type);
        // import parameters
        if (final_stage_cost)
        {
            if (!final_stage_cost->fromMessage(msg.final_stage_cost(), issues)) return;
            if (!final_stage_cost->checkParameters(_dynamics->getStateDimension(), _dynamics->getInputDimension(), issues)) return;
            setFinalStageCost(final_stage_cost);
        }
        else
        {
            if (issues) *issues << "StructuredOptimalControlProblem: unknown final_stage_cost specified.\n";
            return;
        }
    }
    else
        setFinalStageCost({});

    // stage equalities
    if (msg.has_stage_equalities() && !msg.stage_equalities().has_no_function())
    {
        // construct object
        util::get_oneof_field_type_expand_isolated(msg.stage_equalities(), "stage_equalities", type, false, 1);
        StageEqualityConstraint::Ptr stage_equalities = StageEqualitiesFactory::instance().create(type);
        // import parameters
        if (stage_equalities)
        {
            if (!stage_equalities->fromMessage(msg.stage_equalities(), issues)) return;
            if (!stage_equalities->checkParameters(_dynamics->getStateDimension(), _dynamics->getInputDimension(), issues)) return;

            setStageEqualityConstraint(stage_equalities);
        }
        else
        {
            if (issues) *issues << "StructuredOptimalControlProblem: unknown stage_equalities specified.\n";
            return;
        }
    }
    else
        setStageEqualityConstraint({});

    // stage inequalities
    if (msg.has_stage_inequalities() && !msg.stage_inequalities().has_no_function())
    {
        // construct object
        util::get_oneof_field_type_expand_isolated(msg.stage_inequalities(), "stage_inequalities", type, false, 1);
        StageInequalityConstraint::Ptr stage_inequalities = StageInequalitiesFactory::instance().create(type);
        // import parameters
        if (stage_inequalities)
        {
            if (!stage_inequalities->fromMessage(msg.stage_inequalities(), issues)) return;
            if (!stage_inequalities->checkParameters(_dynamics->getStateDimension(), _dynamics->getInputDimension(), issues)) return;

            setStageInequalityConstraint(stage_inequalities);
        }
        else
        {
            if (issues) *issues << "StructuredOptimalControlProblem: unknown stage_inequalities specified.\n";
            return;
        }
    }
    else
        setStageInequalityConstraint({});

    // final stage constraints
    if (msg.has_final_stage_constraints() && !msg.final_stage_constraints().has_no_function())
    {
        // construct object
        util::get_oneof_field_type_expand_isolated(msg.final_stage_constraints(), "final_stage_constraints", type, false, 1);
        FinalStageConstraint::Ptr final_stage_constraint = FinalStageConstraintFactory::instance().create(type);
        // import parameters
        if (final_stage_constraint)
        {
            if (!final_stage_constraint->fromMessage(msg.final_stage_constraints(), issues)) return;
            if (!final_stage_constraint->checkParameters(_dynamics->getStateDimension(), _dynamics->getInputDimension(), _functions.final_stage_cost,
                                                         issues))
                return;
            setFinalStageConstraint(final_stage_constraint);
        }
        else
        {
            if (issues) *issues << "StructuredOptimalControlProblem: unknown final_stage_constraints specified.\n";
            return;
        }
    }
    else
        setFinalStageConstraint({});

    // stage preprocessor
    if (msg.has_stage_preprocessors() && !msg.stage_preprocessors().has_no_function())
    {
        // construct object
        util::get_oneof_field_type_expand_isolated(msg.stage_preprocessors(), "stage_preprocessors", type, false, 1);

        StagePreprocessor::Ptr stage_preprocessor = StagePreprocessorFactory::instance().create(type);
        // import parameters
        if (stage_preprocessor)
        {
            if (!stage_preprocessor->fromMessage(msg.stage_preprocessors(), issues)) return;
            setStagePreprocessor(stage_preprocessor);
        }
        else
        {
            if (issues) *issues << "StructuredOptimalControlProblem: unknown stage_processors specified.\n";
            return;
        }
    }
    else
        setStagePreprocessor({});

    // bounds
    assert(dynamics);  // should be loaded before, otherwise return ...
    Eigen::VectorXd x_lb, x_ub, u_lb, u_ub;
    if (msg.x_min_size() > 0 && msg.x_min_size() != dim_x && issues)
    {
        *issues << "StructuredOptimalControlProblem: x_min size does not match dimension of the system dynamics object" << dim_x << ".\n";
    }
    else
    {
        x_lb = Eigen::Map<const Eigen::VectorXd>(msg.x_min().data(), msg.x_min_size());
    }

    if (msg.x_max_size() > 0 && msg.x_max_size() != dim_x && issues)
    {
        *issues << "StructuredOptimalControlProblem: x_max size does not match dimension of the system dynamics object " << dim_x << ".\n";
    }
    else
    {
        x_ub = Eigen::Map<const Eigen::VectorXd>(msg.x_max().data(), msg.x_max_size());
    }

    if (msg.u_min_size() > 0 && msg.u_min_size() != dim_u && issues)
    {
        *issues << "StructuredOptimalControlProblem: u_min size does not match dimension of the system dynamics object" << dim_u << ".\n";
    }
    else
    {
        u_lb = Eigen::Map<const Eigen::VectorXd>(msg.u_min().data(), msg.u_min_size());
    }

    if (msg.u_max_size() > 0 && msg.u_max_size() != dim_u && issues)
    {
        *issues << "StructuredOptimalControlProblem: u_max size does not match dimension of the system dynamics object " << dim_u << ".\n";
    }
    else
    {
        u_ub = Eigen::Map<const Eigen::VectorXd>(msg.u_max().data(), msg.u_max_size());
    }
    setBounds(x_lb, x_ub, u_lb, u_ub);

    // hyper-graph optimization problem
    if (!msg.has_hyper_graph_strategy())
    {
        if (issues) *issues << "StructuredOptimalControlProblem: no hyper_graph_strategy specified.\n";
        return;
    }
    // construct object
    util::get_oneof_field_type(msg.hyper_graph_strategy(), "hyper_graph_optimization_problem", type, false);
    BaseHyperGraphOptimizationProblem::Ptr optim_prob = HyperGraphOptimizationProblemFactory::instance().create(type);
    // import parameters
    if (optim_prob)
    {
        optim_prob->fromMessage(msg.hyper_graph_strategy(), issues);
        setHyperGraphOptimizationProblem(optim_prob);
    }
    else
    {
        if (issues) *issues << "StructuredOptimalControlProblem: unknown hyper_graph_strategy specified.\n";
        return;
    }

    // solver
    if (!msg.has_nlp_solver())
    {
        if (issues) *issues << "StructuredOptimalControlProblem: no solver specified.\n";
        return;
    }
    // construct object
    util::get_oneof_field_type(msg.nlp_solver(), "solver", type, false);
    NlpSolverInterface::Ptr solver = NlpSolverFactory::instance().create(type);
    // import parameters
    if (solver)
    {
        solver->fromMessage(msg.nlp_solver(), issues);
        setSolver(solver);
    }
    else
    {
        if (issues) *issues << "StructuredOptimalControlProblem: unknown solver specified.\n";
        return;
    }

    _increase_n_if_infeas = msg.increase_n_if_infeas();
}
#endif

}  // namespace corbo
