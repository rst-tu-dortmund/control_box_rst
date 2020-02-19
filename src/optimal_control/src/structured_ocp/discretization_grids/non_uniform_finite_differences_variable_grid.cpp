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

#include <corbo-optimal-control/structured_ocp/discretization_grids/non_uniform_finite_differences_variable_grid.h>

#include <corbo-optimal-control/structured_ocp/edges/finite_differences_collocation_edges.h>
#include <corbo-optimal-control/structured_ocp/edges/misc_edges.h>

#include <corbo-communication/utilities.h>
#include <corbo-core/console.h>

#include <algorithm>
#include <cmath>
#include <memory>

namespace corbo {

void NonUniformFiniteDifferencesVariableGrid::setDtBounds(double dt_lb, double dt_ub)
{
    _dt_lb = dt_lb;
    _dt_ub = dt_ub;
}

void NonUniformFiniteDifferencesVariableGrid::setGridAdaptTimeBasedSingleStep(int n_max, double dt_hyst_ratio)
{
    _grid_adapt    = GridAdaptStrategy::TimeBasedSingleStep;
    _n_max         = n_max;
    _dt_hyst_ratio = dt_hyst_ratio;
}

void NonUniformFiniteDifferencesVariableGrid::setGridAdaptRedundantControls(int n_max, int num_backup_nodes, double epsilon)
{
    _grid_adapt             = GridAdaptStrategy::RedundantControls;
    _n_max                  = n_max;
    _redundant_ctrl_backup  = num_backup_nodes;
    _redundant_ctrl_epsilon = epsilon;
}

void NonUniformFiniteDifferencesVariableGrid::createEdges(NlpFunctions& nlp_fun, OptimizationEdgeSet& edges, SystemDynamicsInterface::Ptr dynamics)
{
    assert(isValid());

    // clear edges first
    // TODO(roesmann): we could implement a more efficient strategy without recreating the whole edgeset everytime
    edges.clear();

    int n = getN();

    std::vector<BaseEdge::Ptr> cost_terms, eq_terms, ineq_terms;
    for (int k = 0; k < n - 1; ++k)
    {
        VectorVertex& x_next  = (k < n - 2) ? _x_seq[k + 1] : _xf;
        VectorVertex& u_prev  = (k > 0) ? _u_seq[k - 1] : _u_prev;
        ScalarVertex& dt_prev = (k > 0) ? _dt_seq[k - 1] : _u_prev_dt;

        cost_terms.clear();
        eq_terms.clear();
        ineq_terms.clear();
        nlp_fun.getNonIntegralStageFunctionEdges(k, _x_seq[k], _u_seq[k], _dt_seq[k], u_prev, dt_prev, cost_terms, eq_terms, ineq_terms);
        for (BaseEdge::Ptr& edge : cost_terms) edges.addObjectiveEdge(edge);
        for (BaseEdge::Ptr& edge : eq_terms) edges.addEqualityEdge(edge);
        for (BaseEdge::Ptr& edge : ineq_terms) edges.addInequalityEdge(edge);

        if (nlp_fun.stage_cost && nlp_fun.stage_cost->hasIntegralTerms(k))
        {
            if (_cost_integration == CostIntegrationRule::TrapezoidalRule)
            {
                TrapezoidalIntegralCostEdge::Ptr edge =
                    std::make_shared<TrapezoidalIntegralCostEdge>(_x_seq[k], _u_seq[k], x_next, _dt_seq[k], nlp_fun.stage_cost, k);
                edges.addObjectiveEdge(edge);
            }
            else if (_cost_integration == CostIntegrationRule::LeftSum)
            {
                LeftSumCostEdge::Ptr edge = std::make_shared<LeftSumCostEdge>(_x_seq[k], _u_seq[k], _dt_seq[k], nlp_fun.stage_cost, k);
                edges.addObjectiveEdge(edge);
            }
            else
                PRINT_ERROR_NAMED("Cost integration rule not implemented");
        }

        if (nlp_fun.stage_equalities && nlp_fun.stage_equalities->hasIntegralTerms(k))
        {
            if (_cost_integration == CostIntegrationRule::TrapezoidalRule)
            {
                TrapezoidalIntegralEqualityDynamicsEdge::Ptr edge = std::make_shared<TrapezoidalIntegralEqualityDynamicsEdge>(
                    dynamics, _x_seq[k], _u_seq[k], x_next, _dt_seq[k], nlp_fun.stage_equalities, k);
                edge->setFiniteDifferencesCollocationMethod(_fd_eval);
                edges.addEqualityEdge(edge);
            }
            else if (_cost_integration == CostIntegrationRule::LeftSum)
            {
                LeftSumEqualityEdge::Ptr edge = std::make_shared<LeftSumEqualityEdge>(_x_seq[k], _u_seq[k], _dt_seq[k], nlp_fun.stage_equalities, k);
                edges.addEqualityEdge(edge);

                // system dynamics edge
                FDCollocationEdge::Ptr sys_edge = std::make_shared<FDCollocationEdge>(dynamics, _x_seq[k], _u_seq[k], x_next, _dt_seq[k]);
                sys_edge->setFiniteDifferencesCollocationMethod(_fd_eval);
                edges.addEqualityEdge(sys_edge);
            }
            else
                PRINT_ERROR_NAMED("Cost integration rule not implemented");
        }
        else
        {
            // just the system dynamics edge
            FDCollocationEdge::Ptr edge = std::make_shared<FDCollocationEdge>(dynamics, _x_seq[k], _u_seq[k], x_next, _dt_seq[k]);
            edge->setFiniteDifferencesCollocationMethod(_fd_eval);
            edges.addEqualityEdge(edge);
        }

        if (nlp_fun.stage_inequalities && nlp_fun.stage_inequalities->hasIntegralTerms(k))
        {
            if (_cost_integration == CostIntegrationRule::TrapezoidalRule)
            {
                TrapezoidalIntegralInequalityEdge::Ptr edge =
                    std::make_shared<TrapezoidalIntegralInequalityEdge>(_x_seq[k], _u_seq[k], x_next, _dt_seq[k], nlp_fun.stage_inequalities, k);
                edges.addInequalityEdge(edge);
            }
            else if (_cost_integration == CostIntegrationRule::LeftSum)
            {
                LeftSumInequalityEdge::Ptr edge =
                    std::make_shared<LeftSumInequalityEdge>(_x_seq[k], _u_seq[k], _dt_seq[k], nlp_fun.stage_inequalities, k);
                edges.addInequalityEdge(edge);
            }
            else
                PRINT_ERROR_NAMED("Cost integration rule not implemented");
        }

        if (_dt_eq_constraint && k > 0)
        {
            TwoScalarEqualEdge::Ptr edge = std::make_shared<TwoScalarEqualEdge>(_dt_seq[k - 1], _dt_seq[k]);
            edges.addEqualityEdge(edge);
        }
    }

    // check if we have a separate unfixed final state
    if (!_xf.isFixed())
    {
        // set final state cost
        BaseEdge::Ptr cost_edge = nlp_fun.getFinalStateCostEdge(n - 1, _xf);
        if (cost_edge) edges.addObjectiveEdge(cost_edge);

        // set final state constraint
        BaseEdge::Ptr constr_edge = nlp_fun.getFinalStateConstraintEdge(n - 1, _xf);
        if (constr_edge)
        {
            if (nlp_fun.final_stage_constraints->isEqualityConstraint())
                edges.addEqualityEdge(constr_edge);
            else
                edges.addInequalityEdge(constr_edge);
        }
    }
}

bool NonUniformFiniteDifferencesVariableGrid::adaptGrid(bool new_run, NlpFunctions& nlp_fun)
{
    // do not adapt grid in a new run
    if (new_run && !_adapt_first_iter) return false;

    bool changed = false;
    switch (_grid_adapt)
    {
        case GridAdaptStrategy::NoGridAdapt:
        {
            break;
        }
        case GridAdaptStrategy::TimeBasedSingleStep:
        {
            changed = adaptGridTimeBasedSingleStep(nlp_fun);
            break;
        }
        case GridAdaptStrategy::RedundantControls:
        {
            changed = adaptGridRedundantControls(nlp_fun);
            break;
        }
        default:
        {
            PRINT_ERROR_NAMED("selected grid adaptation strategy not implemented.");
        }
    }
    return changed;
}

bool NonUniformFiniteDifferencesVariableGrid::adaptGridTimeBasedSingleStep(NlpFunctions& nlp_fun)
{
    assert(isValid());
    // TODO(roesmann): should we implement a linear resampling according to the uniform grid?

    // iterate sequence once and check time differences
    bool changed = false;

    int n = getN();

    for (int i = 0; i < (int)_dt_seq.size(); ++i)
    {
        double dt = _dt_seq[i].value();

        if (dt > _dt_ref * (1.0 + _dt_hyst_ratio) && n < _n_max)
        {
            double new_dt = 0.5 * dt;

            // insert tupel
            Eigen::VectorXd new_x = 0.5 * (_x_seq[i].values() + _x_seq[i + 1].values());
            _x_seq.emplace(_x_seq.begin() + i + 1, new_x, nlp_fun.x_lb, nlp_fun.x_ub);
            _u_seq.emplace(_u_seq.begin() + i + 1, _u_seq[i].values(), nlp_fun.u_lb, nlp_fun.u_ub);
            _dt_seq.emplace(_dt_seq.begin() + i + 1, new_dt, _dt_lb, _dt_ub);

            changed = true;

            break;  // TODO(roesmann): better strategy, instead of just searching the first dt?

            ++i;  // skip newly inserted node
        }
        else if (dt < _dt_ref * (1.0 - _dt_hyst_ratio) && n > _n_min)
        {
            // if (i < (int)_dt_seq.size()-2)
            {
                _dt_seq[i + 1].value() += dt;
                _x_seq.erase(_x_seq.begin() + i);
                _u_seq.erase(_u_seq.begin() + i);
                _dt_seq.erase(_dt_seq.begin() + i);

                // if i== 0: fix new start state for optimization
                if (i == 0 && !_x_seq.empty()) _x_seq.front().setFixed(true);

                changed = true;

                break;  // TODO(roesmann): better strategy, instead of just searching the first dt?
            }
        }
    }
    if (changed) setModified(true);
    _n_adapt = getN();
    return changed;
}

bool NonUniformFiniteDifferencesVariableGrid::adaptGridRedundantControls(NlpFunctions& nlp_fun)
{
    assert(isValid());
    // find approx(u_k, u_{k+1})

    if (getN() < 3) return false;

    int num_interv = (int)_u_seq.size();

    std::vector<std::size_t> non_unique_indices;
    for (std::size_t idx = 0; idx < num_interv - 1; ++idx)  // never delete the last control
    {
        // if ( _u_seq[idx-1].values().isApprox(_u_seq[idx].values(), 1e-1 ) )
        // also delete if the time diff of the interval is sufficiently small
        if (_dt_seq[idx].value() < 1e-6)
        {
            non_unique_indices.emplace_back(idx);
            // non_unique_indices.emplace_back(idx+1); // since we have zero dt, at least two samples are obsolete
            // ++idx;
            continue;
        }
        // PRINT_INFO("control threshold: " << _controls_similar_threshold);
        // if ( (_ctrl_seq[idx+1].controls()-_ctrl_seq[idx].controls()).isZero(_controls_similar_threshold) )
        //     non_unique_indices.emplace_back(idx);
        if (((_u_seq[idx + 1].values() - _u_seq[idx].values()).cwiseAbs().array() <= _redundant_ctrl_epsilon).all())
            non_unique_indices.emplace_back(idx);
    }

    //   PRINT_INFO("number non-unique: " << non_unique_indices.size());

    bool changed = false;

    int backup_diff = (int)non_unique_indices.size() - _redundant_ctrl_backup;

    if (backup_diff < 0)
    {
        changed     = true;
        backup_diff = std::abs(backup_diff);
        for (int i = 0; i < backup_diff && getN() < _n_max; ++i)
        {
            // add new sample
            int dt_max_idx = 0;
            if (getN() > 2)
            {
                // TODO(roesmann): maybe rewrite the following lines, this is modified from an old more completex data structure...
                // insert inbetween largest gap (largest DT)
                auto end_it = _dt_seq.end();
                std::advance(end_it, -1);  // we do not want to a new point to the goal (therefore -2)
                auto max_elem = std::max_element(_dt_seq.begin(), end_it,
                                                 [](const ScalarVertex& dtv1, const ScalarVertex& dtv2) { return dtv1.value() < dtv2.value(); });
                if (max_elem == end_it)
                {
                    PRINT_INFO_NAMED("Invalid time max element in resizeTrajectoryRedundantControls(). break...");
                    break;
                }
                dt_max_idx = std::distance(_dt_seq.begin(), max_elem);
            }
            assert(dt_max_idx + 1 < num_interv);

            double new_dt               = 0.5 * _dt_seq[dt_max_idx].value();
            _dt_seq[dt_max_idx].value() = new_dt;

            // insert tupel
            _x_seq.emplace(_x_seq.begin() + dt_max_idx + 1, 0.5 * (_x_seq[dt_max_idx].values() + _x_seq[dt_max_idx + 1].values()), nlp_fun.x_lb,
                           nlp_fun.x_ub);
            _u_seq.emplace(_u_seq.begin() + dt_max_idx + 1, _u_seq[dt_max_idx].values(), nlp_fun.u_lb, nlp_fun.u_ub);
            _dt_seq.emplace(_dt_seq.begin() + dt_max_idx + 1, new_dt, _dt_lb, _dt_ub);
            // PRINT_INFO("inserted");
        }
    }
    else if (backup_diff > 0)
    {
        changed     = true;
        auto idx_it = non_unique_indices.rbegin();  // erase starting from the last one (reverse iterator)
        for (int i = 0; i < backup_diff && getN() > _n_min; ++i)
        {
            int k = (int)*idx_it;
            if (k >= getN() - 2) --k;

            assert(k + 1 < num_interv);

            _dt_seq[k].value() += _dt_seq[k + 1].value();

            _x_seq.erase(_x_seq.begin() + k + 1);
            _u_seq.erase(_u_seq.begin() + k + 1);
            _dt_seq.erase(_dt_seq.begin() + k + 1);

            ++idx_it;
            // PRINT_INFO("removed");
        }
    }
    if (changed) setModified(true);
    _n_adapt = getN();
    return changed;
}

#ifdef MESSAGE_SUPPORT
void NonUniformFiniteDifferencesVariableGrid::fromMessage(const messages::NonUniformFiniteDifferencesVariableGrid& message, std::stringstream* issues)
{
    if (message.n() < 2 && issues) *issues << "NonUniformFiniteDifferencesVariableGrid: Number of states must be greater than or equal 2.\n";
    if (message.dt_ref() <= 0 && issues) *issues << "NonUniformFiniteDifferencesVariableGrid: Dt must be greater than 0.0.\n";

    setNRef(message.n());
    setDtRef(message.dt_ref());

    switch (message.cost_integration_rule())
    {
        case messages::NonUniformFiniteDifferencesVariableGrid::CostIntegrationRule::
            NonUniformFiniteDifferencesVariableGrid_CostIntegrationRule_LeftSum:
        {
            setCostIntegrationRule(CostIntegrationRule::LeftSum);
            break;
        }
        case messages::NonUniformFiniteDifferencesVariableGrid::CostIntegrationRule::
            NonUniformFiniteDifferencesVariableGrid_CostIntegrationRule_TrapezoidalRule:
        {
            setCostIntegrationRule(CostIntegrationRule::TrapezoidalRule);
            break;
        }
        default:
        {
            PRINT_ERROR_NAMED("Selected cost integration rule not implemented");
        }
    };

    // fd collocation method
    // construct object
    std::string type;
    util::get_oneof_field_type(message.fd_collocation(), "fd_collocation", type, false);
    FiniteDifferencesCollocationInterface::Ptr fd_eval = create_from_factory<FiniteDifferencesCollocationInterface>(type);
    // import parameters
    if (fd_eval)
    {
        fd_eval->fromMessage(message.fd_collocation(), issues);
        setFiniteDifferencesCollocationMethod(fd_eval);
    }
    else
    {
        if (issues) *issues << "NonUniformFiniteDifferencesVariableGrid: unknown finite differences collocation method specified.\n";
        return;
    }

    // xf fixed states
    // if (grid_msg.xf_fixed_size() != _p && issues) *issues << "FullDiscretizationGrid: xf_fixed size does not match state dimension " << _p <<
    // ".\n";
    // TODO(roesmann): we need a warning in the gui if xf_fixed has the wrong dimension.
    //                 maybe we could add a "verifyParameters()" method to all interfaces, so predictive control asks the ocp, the ocp the grid and so
    //                 on
    _xf_fixed = Eigen::Map<const Eigen::Matrix<bool, -1, 1>>(message.xf_fixed().data(), message.xf_fixed_size());

    // dt bounds
    setDtBounds(message.dt_min(), message.dt_max());

    _warm_start = message.warm_start();

    // auto resize
    if (message.has_grid_adapt_strategy())
    {
        if (message.grid_adapt_strategy().has_no_grid_adapt())
        {
            disableGridAdaptation();
        }
        else if (message.grid_adapt_strategy().has_time_based_single_step())
        {
            setGridAdaptTimeBasedSingleStep(message.grid_adapt_strategy().time_based_single_step().n_max(),
                                            message.grid_adapt_strategy().time_based_single_step().dt_hyst_ratio());
        }
        else if (message.grid_adapt_strategy().has_redundant_controls())
        {
            setGridAdaptRedundantControls(message.grid_adapt_strategy().redundant_controls().n_max(),
                                          message.grid_adapt_strategy().redundant_controls().num_backup(),
                                          message.grid_adapt_strategy().redundant_controls().epsilon());
        }
    }
    _adapt_first_iter = message.grid_adapt_strategy().adapt_first_iter();

    // others
    _n_min            = message.n_min();
    _dt_eq_constraint = message.dt_eq_constraint();
}

void NonUniformFiniteDifferencesVariableGrid::toMessage(messages::NonUniformFiniteDifferencesVariableGrid& message) const
{
    message.set_n(getNRef());
    message.set_dt_ref(getDtRef());

    switch (_cost_integration)
    {
        case CostIntegrationRule::LeftSum:
        {
            message.set_cost_integration_rule(messages::NonUniformFiniteDifferencesVariableGrid::CostIntegrationRule::
                                                  NonUniformFiniteDifferencesVariableGrid_CostIntegrationRule_LeftSum);
            break;
        }
        case CostIntegrationRule::TrapezoidalRule:
        {
            message.set_cost_integration_rule(messages::NonUniformFiniteDifferencesVariableGrid::CostIntegrationRule::
                                                  NonUniformFiniteDifferencesVariableGrid_CostIntegrationRule_TrapezoidalRule);
            break;
        }
        default:
        {
            PRINT_ERROR_NAMED("Selected cost integration rule not implemented");
        }
    };

    // fd collocation method
    if (_fd_eval) _fd_eval->fromMessage(*message.mutable_fd_collocation());

    // xf fixed states
    if (_xf_fixed.size() > 0)
    {
        message.mutable_xf_fixed()->Resize(_xf_fixed.size(), false);
        Eigen::Map<Eigen::Matrix<bool, -1, 1>>(message.mutable_xf_fixed()->mutable_data(), _xf_fixed.size()) = _xf_fixed;
    }

    // dt bounds
    message.set_dt_min(_dt_lb);
    message.set_dt_max(_dt_ub);

    message.set_warm_start(_warm_start);

    // auto resize
    switch (_grid_adapt)
    {
        case GridAdaptStrategy::NoGridAdapt:
        {
            message.mutable_grid_adapt_strategy()->mutable_no_grid_adapt();
            break;
        }
        case GridAdaptStrategy::TimeBasedSingleStep:
        {
            message.mutable_grid_adapt_strategy()->mutable_time_based_single_step()->set_n_max(_n_max);
            message.mutable_grid_adapt_strategy()->mutable_time_based_single_step()->set_dt_hyst_ratio(_dt_hyst_ratio);
            break;
        }
        case GridAdaptStrategy::RedundantControls:
        {
            message.mutable_grid_adapt_strategy()->mutable_redundant_controls()->set_n_max(_n_max);
            message.mutable_grid_adapt_strategy()->mutable_redundant_controls()->set_num_backup(_redundant_ctrl_backup);
            message.mutable_grid_adapt_strategy()->mutable_redundant_controls()->set_epsilon(_redundant_ctrl_epsilon);
            break;
        }
        default:
        {
            PRINT_ERROR_NAMED("exporting of the selected auto resize strategy not implemented.");
        }
    }
    message.mutable_grid_adapt_strategy()->set_adapt_first_iter(_adapt_first_iter);

    // others
    message.set_n_min(_n_min);
    message.set_dt_eq_constraint(_dt_eq_constraint);
}
#endif

}  // namespace corbo
