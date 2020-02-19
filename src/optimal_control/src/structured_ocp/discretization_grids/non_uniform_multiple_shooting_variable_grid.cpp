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

#include <corbo-optimal-control/structured_ocp/discretization_grids/non_uniform_multiple_shooting_variable_grid.h>

#include <corbo-optimal-control/structured_ocp/edges/misc_edges.h>
#include <corbo-optimal-control/structured_ocp/edges/multiple_shooting_edges.h>

#include <corbo-communication/utilities.h>
#include <corbo-core/console.h>

#include <algorithm>
#include <cmath>
#include <memory>

namespace corbo {

void NonUniformMultipleShootingVariableGrid::setDtBounds(double dt_lb, double dt_ub)
{
    _dt_lb = dt_lb;
    _dt_ub = dt_ub;
}

void NonUniformMultipleShootingVariableGrid::setGridAdaptTimeBasedSingleStep(int n_max, double dt_hyst_ratio)
{
    _grid_adapt    = GridAdaptStrategy::TimeBasedSingleStep;
    _n_max         = n_max;
    _dt_hyst_ratio = dt_hyst_ratio;
}

void NonUniformMultipleShootingVariableGrid::setGridAdaptRedundantControls(int n_max, int num_backup_nodes, double epsilon)
{
    _grid_adapt             = GridAdaptStrategy::RedundantControls;
    _n_max                  = n_max;
    _redundant_ctrl_backup  = num_backup_nodes;
    _redundant_ctrl_epsilon = epsilon;
}

void NonUniformMultipleShootingVariableGrid::createEdges(NlpFunctions& nlp_fun, OptimizationEdgeSet& edges, SystemDynamicsInterface::Ptr dynamics)
{
    assert(isValid());
    assert(_integrator);

    // clear edges first
    // TODO(roesmann): we could implement a more efficient strategy without recreating the whole edgeset everytime
    edges.clear();

    int n = getN();

    std::vector<BaseEdge::Ptr> cost_terms, eq_terms, ineq_terms;
    int k = 0;
    for (int i = 0; i < (int)_intervals.size(); ++i)
    {
        VectorVertex& s                   = _intervals[i].s;
        VectorVertex& s_next              = (i < (int)_intervals.size() - 1) ? _intervals[i + 1].s : _xf;
        std::vector<VectorVertex>& u_seq  = _intervals[i].u_seq;
        VectorVertex& u_prev              = (i == 0) ? _u_prev : _intervals[i - 1].u_seq.back();  // TODO(roesmann) add also to other grids...
        std::vector<ScalarVertex>& dt_seq = _intervals[i].dt_seq;
        ScalarVertex& dt_prev             = (i == 0) ? _u_prev_dt : _intervals[i - 1].dt_seq.back();
        assert(u_seq.size() == dt_seq.size());

        if (u_seq.size() == 1)
        {
            // full discretization
            cost_terms.clear();
            eq_terms.clear();
            ineq_terms.clear();
            nlp_fun.getNonIntegralStageFunctionEdges(k, s, u_seq.front(), dt_seq.front(), u_prev, dt_prev, cost_terms, eq_terms, ineq_terms);
            for (BaseEdge::Ptr& edge : cost_terms) edges.addObjectiveEdge(edge);
            for (BaseEdge::Ptr& edge : eq_terms) edges.addEqualityEdge(edge);
            for (BaseEdge::Ptr& edge : ineq_terms) edges.addInequalityEdge(edge);

            if (nlp_fun.hasIntegralTerms(k))
            {
                // system dynamics and integral cost/constraints
                MultipleShootingEdgeSingleControl::Ptr edge = std::make_shared<MultipleShootingEdgeSingleControl>(
                    dynamics, nlp_fun.stage_cost, nlp_fun.stage_equalities, nlp_fun.stage_inequalities, k, s, u_seq.front(), dt_seq.front(), s_next);
                edge->setIntegrator(_integrator);
                edges.addMixedEdge(edge);
            }
            else
            {
                // system dynamics
                MSVariableDynamicsOnlyEdge::Ptr edge =
                    std::make_shared<MSVariableDynamicsOnlyEdge>(dynamics, s, u_seq.front(), s_next, dt_seq.front());
                edge->setIntegrator(_integrator);
                edges.addEqualityEdge(edge);
            }
            //            if (_dt_eq_constraint && i > 0)
            //            {
            //                TwoScalarEqualEdge::Ptr edge = std::make_shared<TwoScalarEqualEdge>(_intervals[i - 1].dt_seq.back(), dt_seq.front());
            //                edges.addEqualityEdge(edge);
            //            }
        }
        else
        {
            // TODO(roesmann): maybe we are forgetting to evaluate non-integral stage functions at non-shooting grid points?? -> Check (
            // it seems so according to a quick simulation...: quadratic form cost (non-integral)

            int num_controls = (int)u_seq.size();

            if (nlp_fun.hasIntegralTerms(k) || _intermediate_x_constraints)
            {
                MultipleShootingEdge::Ptr edge =
                    std::make_shared<MultipleShootingEdge>(dynamics, nlp_fun.stage_cost, nlp_fun.stage_equalities, nlp_fun.stage_inequalities, k,
                                                           num_controls + 1, num_controls,  // TODO(roesmann): here could use a different nc
                                                           (int)dt_seq.size(), _intermediate_x_constraints);
                edge->setIntegrator(_integrator);

                edge->setVertex(0, s);
                edge->setVertex(1, s_next);
                int vert_idx = 2;
                for (int i = 0; i < (int)u_seq.size(); ++i)
                {
                    edge->setVertex(vert_idx++, u_seq[i]);
                }
                for (int i = 0; i < (int)dt_seq.size(); ++i)
                {
                    edge->setVertex(vert_idx++, dt_seq[i]);
                }
                edge->finalize();
                edges.addMixedEdge(edge);
            }
            else
            {
                // dynamics only
                MSDynamicsOnlyMultiControlsMultiDtsEdge::Ptr edge =
                    std::make_shared<MSDynamicsOnlyMultiControlsMultiDtsEdge>(dynamics, s, s_next, (int)u_seq.size());
                edge->setIntegrator(_integrator);

                for (int i = 0; i < (int)u_seq.size(); ++i)
                {
                    edge->setVertex(2 + 2 * i, u_seq[i]);
                    edge->setVertex(2 + 2 * i + 1, dt_seq[i]);
                }
                edge->finalize();
                edges.addEqualityEdge(edge);
            }

            // add non-integral terms at shooting nodes
            cost_terms.clear();
            eq_terms.clear();
            ineq_terms.clear();
            nlp_fun.getNonIntegralStageFunctionEdges(k, s, u_seq.front(), dt_seq.front(), u_prev, dt_prev, cost_terms, eq_terms, ineq_terms);
            for (BaseEdge::Ptr& edge : cost_terms) edges.addObjectiveEdge(edge);
            for (BaseEdge::Ptr& edge : eq_terms) edges.addEqualityEdge(edge);
            for (BaseEdge::Ptr& edge : ineq_terms) edges.addInequalityEdge(edge);

            std::vector<BaseEdge::Ptr> function_terms;

            // also add terms that depend only on controls and dts to the interval variables (skip first ones since they are considered above)
            for (int j = 1; j < (int)u_seq.size(); ++j)
            {
                // VectorVertex& ukk = i < (int)u_seq.size() ? u_seq[i] : *interval.controls.back();    // ith control of the shooting interval -> for
                // move blocking
                if (nlp_fun.stage_cost)
                {
                    cost_terms.clear();
                    getStageFunctionEdgesIntermediateControlDtMultipleShooting(k + j, u_seq[j], dt_seq[j], u_prev, dt_prev, *nlp_fun.stage_cost,
                                                                               cost_terms);
                    for (BaseEdge::Ptr& edge : cost_terms) edges.addObjectiveEdge(edge);
                }
                if (nlp_fun.stage_equalities)
                {
                    eq_terms.clear();
                    getStageFunctionEdgesIntermediateControlDtMultipleShooting(k + j, u_seq[j], dt_seq[j], u_prev, dt_prev, *nlp_fun.stage_equalities,
                                                                               eq_terms);
                    for (BaseEdge::Ptr& edge : eq_terms) edges.addEqualityEdge(edge);
                }
                if (nlp_fun.stage_inequalities)
                {
                    ineq_terms.clear();
                    getStageFunctionEdgesIntermediateControlDtMultipleShooting(k + j, u_seq[j], dt_seq[j], u_prev, dt_prev,
                                                                               *nlp_fun.stage_inequalities, ineq_terms);
                    for (BaseEdge::Ptr& edge : ineq_terms) edges.addInequalityEdge(edge);
                }

                //                if (_dt_eq_constraint && j > 0)
                //                {
                //                    TwoScalarEqualEdge::Ptr edge = std::make_shared<TwoScalarEqualEdge>(dt_seq[j - 1], dt_seq[j]);
                //                    edges.addEqualityEdge(edge);
                //                }
            }
        }
        k += (int)u_seq.size();
    }

    if (_dt_eq_constraint)
    {
        for (int i = 0; i < (int)_intervals.size(); ++i)
        {
            for (int j = 0; j < (int)_intervals[i].dt_seq.size(); ++j)
            {
                if (i == 0 && j == 0) continue;
                ScalarVertex& dt1 = (j == 0) ? _intervals[i - 1].dt_seq.back() : _intervals[i].dt_seq[j - 1];
                ScalarVertex& dt2 = _intervals[i].dt_seq[j];

                TwoScalarEqualEdge::Ptr edge = std::make_shared<TwoScalarEqualEdge>(dt1, dt2);
                edges.addEqualityEdge(edge);
            }
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

    // add control deviation edges for last control
    cost_terms.clear();
    eq_terms.clear();
    ineq_terms.clear();
    nlp_fun.getFinalControlDeviationEdges(n, _u_ref, _intervals.back().u_seq.back(), _intervals.back().dt_seq.back(), cost_terms, eq_terms,
                                          ineq_terms);
    for (BaseEdge::Ptr& edge : cost_terms) edges.addObjectiveEdge(edge);
    for (BaseEdge::Ptr& edge : eq_terms) edges.addEqualityEdge(edge);
    for (BaseEdge::Ptr& edge : ineq_terms) edges.addInequalityEdge(edge);
}

void NonUniformMultipleShootingVariableGrid::getStageFunctionEdgesIntermediateControlDtMultipleShooting(int k, VectorVertex& uk, ScalarVertex& dt,
                                                                                                        VectorVertex& u_prev, ScalarVertex& dt_prev,
                                                                                                        const StageFunction& stage_fun,
                                                                                                        std::vector<BaseEdge::Ptr>& edges)
{
    int dim = stage_fun.getNonIntegralControlTermDimension(k);
    if (dim > 0)
    {
        using Edge = UnaryVectorVertexEdge<StageFunction, &StageFunction::computeNonIntegralControlTerm>;
        Edge::Ptr edge =
            std::make_shared<Edge>(dim, k, stage_fun, uk, stage_fun.isLinearNonIntegralControlTerm(k), stage_fun.isLsqFormNonIntegralControlTerm(k));
        edges.push_back(edge);
    }

    dim = stage_fun.getNonIntegralDtTermDimension(k);
    if (dim > 0)
    {
        using Edge = UnaryScalarVertexEdge<StageFunction, &StageFunction::computeNonIntegralDtTerm>;
        Edge::Ptr edge =
            std::make_shared<Edge>(dim, k, stage_fun, dt, stage_fun.isLinearNonIntegralDtTerm(k), stage_fun.isLsqFormNonIntegralDtTerm(k));
        edges.push_back(edge);
    }

    dim = stage_fun.getNonIntegralControlDeviationTermDimension(k);
    if (dim > 0)
    {
        using Edge     = TernaryVectorScalarVertexEdge<StageFunction, &StageFunction::computeNonIntegralControlDeviationTerm>;
        Edge::Ptr edge = std::make_shared<Edge>(dim, k, stage_fun, uk, u_prev, dt_prev, false, false);
        edges.push_back(edge);
    }
}
bool NonUniformMultipleShootingVariableGrid::adaptGrid(bool new_run, NlpFunctions& nlp_fun)
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

bool NonUniformMultipleShootingVariableGrid::adaptGridTimeBasedSingleStep(NlpFunctions& nlp_fun)
{
    assert(isValid());
    // TODO(roesmann): should we implement a linear resampling according to the uniform grid?

    // iterate sequence once and check time differences
    bool changed = false;

    if (_full_discretization)
    {
        int n          = getN();
        int num_interv = (int)_intervals.size();  // -> full discretization

        for (int i = 0; i < num_interv; ++i)
        {
            double dt = _intervals[i].dt_seq.front().value();

            if (dt > _dt_ref * (1.0 + _dt_hyst_ratio) && n < _n_max)
            {
                double new_dt = 0.5 * dt;

                // insert tupel
                // insert tupel
                auto new_interv_it         = _intervals.insert(_intervals.begin() + i + 1, ShootingInterval());
                ShootingInterval& interval = *new_interv_it;

                interval.s.set(0.5 * (_intervals[i].s.values() + _intervals[i + 2].s.values()), nlp_fun.x_lb,
                               nlp_fun.x_ub);  // we need to consider i+2 since we inserted at i+1
                interval.u_seq.emplace_back(_intervals[i].u_seq.front().values(), nlp_fun.u_lb, nlp_fun.u_ub);
                interval.dt_seq.emplace_back(new_dt, _dt_lb, _dt_ub);

                changed = true;

                break;  // TODO(roesmann): better strategy, instead of just searching the first dt?

                ++i;  // skip newly inserted node
            }
            else if (dt < _dt_ref * (1.0 - _dt_hyst_ratio) && n > _n_min)
            {
                // if (i < (int)_dt_seq.size()-2)
                {
                    _intervals[i + 1].dt_seq.front().value() += dt;
                    _intervals.erase(_intervals.begin() + i);

                    // if i== 0: fix new start state for optimization
                    if (i == 0 && !_intervals.empty()) _intervals.front().s.setFixed(true);

                    changed = true;

                    break;  // TODO(roesmann): better strategy, instead of just searching the first dt?
                }
            }
        }
    }
    else
    {
        PRINT_WARNING_ONCE("method for shooting grids with more than 1 control per interval not yet implemented.");
        // TODO(roesmann): probably we must adhere to _num_interv_ref while resampling, so recalculate how many controls are in every shooting
        // interval
    }
    if (changed) setModified(true);
    _n_adapt = getN();
    return changed;
}

bool NonUniformMultipleShootingVariableGrid::adaptGridRedundantControls(NlpFunctions& nlp_fun)
{
    assert(isValid());
    // find approx(u_k, u_{k+1})

    if (getN() < 3) return false;

    bool changed = false;

    if (_full_discretization)
    {
        int num_interv = (int)_intervals.size();

        std::vector<std::size_t> non_unique_indices;
        for (std::size_t idx = 0; idx < num_interv - 1; ++idx)  // never delete the last control
        {
            // if ( _u_seq[idx-1].values().isApprox(_u_seq[idx].values(), 1e-1 ) )
            // also delete if the time diff of the interval is sufficiently small
            if (_intervals[idx].dt_seq.front().value() < 1e-6)
            {
                non_unique_indices.emplace_back(idx);
                // non_unique_indices.emplace_back(idx+1); // since we have zero dt, at least two samples are obsolete
                // ++idx;
                continue;
            }
            // PRINT_INFO("control threshold: " << _controls_similar_threshold);
            // if ( (_ctrl_seq[idx+1].controls()-_ctrl_seq[idx].controls()).isZero(_controls_similar_threshold) )
            //     non_unique_indices.emplace_back(idx);
            if (((_intervals[idx + 1].u_seq.front().values() - _intervals[idx].u_seq.front().values()).cwiseAbs().array() <= _redundant_ctrl_epsilon)
                    .all())
                non_unique_indices.emplace_back(idx);
        }

        //   PRINT_INFO("number non-unique: " << non_unique_indices.size());

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
                    auto end_it = _intervals.end();
                    std::advance(end_it, -1);  // we do not want to a new point to the goal (therefore -2)
                    auto max_elem =
                        std::max_element(_intervals.begin(), end_it, [](const ShootingInterval& interv1, const ShootingInterval& interv2) {
                            return interv1.dt_seq.front().value() < interv2.dt_seq.front().value();
                        });
                    if (max_elem == end_it)
                    {
                        PRINT_INFO_NAMED("Invalid time max element. break...");
                        break;
                    }
                    dt_max_idx = std::distance(_intervals.begin(), max_elem);
                }
                if (dt_max_idx + 1 >= num_interv)
                {
                    PRINT_ERROR("Number of backup parameters too high. Not yet implemented.");  // TODO(roesmann): we cannot interpolate, we must
                                                                                                // extrapolate...
                }
                else
                {
                    double new_dt                                 = 0.5 * _intervals[dt_max_idx].dt_seq.front().value();
                    _intervals[dt_max_idx].dt_seq.front().value() = new_dt;

                    // insert tupel
                    Eigen::VectorXd new_s = 0.5 * (_intervals[dt_max_idx].s.values() + _intervals[dt_max_idx + 1].s.values());

                    auto new_interv_it         = _intervals.insert(_intervals.begin() + dt_max_idx + 1, ShootingInterval());
                    ShootingInterval& interval = *new_interv_it;

                    interval.s.set(new_s, nlp_fun.x_lb, nlp_fun.x_ub);
                    interval.u_seq.emplace_back(_intervals[dt_max_idx].u_seq.front().values(), nlp_fun.u_lb, nlp_fun.u_ub);
                    interval.dt_seq.emplace_back(new_dt, _dt_lb, _dt_ub);
                    // PRINT_INFO("inserted");
                }
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

                _intervals[k].dt_seq.front().value() += _intervals[k + 1].dt_seq.front().value();

                _intervals.erase(_intervals.begin() + k + 1);

                ++idx_it;
                // PRINT_INFO("removed");
            }
        }
    }
    else
    {
        PRINT_WARNING_ONCE("method for shooting grids with more than 1 control per interval not yet implemented.");
        // TODO(roesmann): probably we must adhere to _num_interv_ref while resampling, so recalculate how many controls are in every shooting
        // interval
    }

    if (changed) setModified(true);
    _n_adapt = getN();
    return changed;
}

#ifdef MESSAGE_SUPPORT
void NonUniformMultipleShootingVariableGrid::fromMessage(const messages::NonUniformMultipleShootingVariableGrid& message, std::stringstream* issues)
{
    if (message.n() < 2 && issues) *issues << "MultipleShootingGrid: Number of states must be greater than or equal 2.\n";
    if (message.dt_ref() <= 0 && issues) *issues << "MultipleShootingGrid: Dt must be greater than 0.0.\n";

    setNRef(message.n());
    setDtRef(message.dt_ref());
    setNumControlsPerShootingInterval(message.num_u_per_interval(), message.intermediate_x_constraints());
    setWarmStart(message.warm_start());

    // fd collocation method
    // construct object
    std::string type;
    util::get_oneof_field_type(message.integrator(), "explicit_integrator", type, false);
    NumericalIntegratorExplicitInterface::Ptr integrator = create_from_factory<NumericalIntegratorExplicitInterface>(type);
    // import parameters
    if (integrator)
    {
        integrator->fromMessage(message.integrator(), issues);
        setNumericalIntegrator(integrator);
    }
    else
    {
        if (issues) *issues << "MultipleShootingGrid: unknown finite differences collocation method specified.\n";
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

void NonUniformMultipleShootingVariableGrid::toMessage(messages::NonUniformMultipleShootingVariableGrid& message) const
{
    message.set_n(getNRef());
    message.set_dt_ref(getDtRef());
    message.set_num_u_per_interval(_num_u_per_interv_ref);
    message.set_intermediate_x_constraints(_intermediate_x_constraints);
    message.set_warm_start(_warm_start);

    // fd collocation method
    if (_integrator) _integrator->fromMessage(*message.mutable_integrator());

    // xf fixed states
    if (_xf_fixed.size() > 0)
    {
        message.mutable_xf_fixed()->Resize(_xf_fixed.size(), false);
        Eigen::Map<Eigen::Matrix<bool, -1, 1>>(message.mutable_xf_fixed()->mutable_data(), _xf_fixed.size()) = _xf_fixed;
    }

    // dt bounds
    message.set_dt_min(_dt_lb);
    message.set_dt_max(_dt_ub);

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
