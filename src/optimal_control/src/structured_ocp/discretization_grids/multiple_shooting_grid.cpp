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

#include <corbo-optimal-control/structured_ocp/discretization_grids/multiple_shooting_grid.h>

#include <corbo-optimal-control/structured_ocp/edges/multiple_shooting_edges.h>

#include <corbo-communication/utilities.h>
#include <corbo-core/console.h>

#include <algorithm>
#include <cmath>
#include <memory>

namespace corbo {

void MultipleShootingGrid::createEdges(NlpFunctions& nlp_fun, OptimizationEdgeSet& edges, SystemDynamicsInterface::Ptr dynamics)
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
        VectorVertex& s                  = _intervals[i].s;
        VectorVertex& s_next             = (i < (int)_intervals.size() - 1) ? _intervals[i + 1].s : _xf;
        std::vector<VectorVertex>& u_seq = _intervals[i].u_seq;
        VectorVertex& u_prev             = (i == 0) ? _u_prev : _intervals[i - 1].u_seq.back();  // TODO(roesmann): add also to other grids...
        ScalarVertex& dt_prev            = (i == 0) ? _u_prev_dt : _dt;

        if (u_seq.size() == 1)
        {
            // full discretization
            cost_terms.clear();
            eq_terms.clear();
            ineq_terms.clear();
            nlp_fun.getNonIntegralStageFunctionEdges(k, s, u_seq.front(), _dt, u_prev, dt_prev, cost_terms, eq_terms, ineq_terms);
            for (BaseEdge::Ptr& edge : cost_terms) edges.addObjectiveEdge(edge);
            for (BaseEdge::Ptr& edge : eq_terms) edges.addEqualityEdge(edge);
            for (BaseEdge::Ptr& edge : ineq_terms) edges.addInequalityEdge(edge);

            if (nlp_fun.hasIntegralTerms(k))
            {
                // system dynamics and integral cost/constraints
                MultipleShootingEdgeSingleControl::Ptr edge = std::make_shared<MultipleShootingEdgeSingleControl>(
                    dynamics, nlp_fun.stage_cost, nlp_fun.stage_equalities, nlp_fun.stage_inequalities, k, s, u_seq.front(), _dt, s_next);
                edge->setIntegrator(_integrator);
                edges.addMixedEdge(edge);
            }
            else
            {
                // system dynamics
                MSVariableDynamicsOnlyEdge::Ptr edge = std::make_shared<MSVariableDynamicsOnlyEdge>(
                    dynamics, s, u_seq.front(), s_next, _dt);  // TODO(roesmann): this is version that allows a variable dt
                edge->setIntegrator(_integrator);
                edges.addEqualityEdge(edge);
            }
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
                                                           1, _intermediate_x_constraints);
                edge->setIntegrator(_integrator);

                edge->setVertex(0, s);
                edge->setVertex(1, s_next);
                int vert_idx = 2;
                for (int i = 0; i < num_controls; ++i)
                {
                    edge->setVertex(vert_idx++, u_seq[i]);
                }
                edge->setVertex(vert_idx++, _dt);
                edge->finalize();
                edges.addMixedEdge(edge);
            }
            else
            {
                // dynamics only
                MSDynamicsOnlyMultiControlsEdge::Ptr edge = std::make_shared<MSDynamicsOnlyMultiControlsEdge>(dynamics, num_controls);
                edge->setIntegrator(_integrator);

                edge->setVertex(0, s);
                edge->setVertex(1, s_next);
                edge->setVertex(2, _dt);
                for (int i = 0; i < num_controls; ++i)
                {
                    edge->setVertex(3 + i, u_seq[i]);
                }
                edge->finalize();
                edges.addEqualityEdge(edge);
            }

            // add non-integral terms at shooting nodes
            cost_terms.clear();
            eq_terms.clear();
            ineq_terms.clear();
            nlp_fun.getNonIntegralStageFunctionEdges(k, s, u_seq.front(), _dt, u_prev, dt_prev, cost_terms, eq_terms, ineq_terms);
            for (BaseEdge::Ptr& edge : cost_terms) edges.addObjectiveEdge(edge);
            for (BaseEdge::Ptr& edge : eq_terms) edges.addEqualityEdge(edge);
            for (BaseEdge::Ptr& edge : ineq_terms) edges.addInequalityEdge(edge);

            // also add terms that depend only on controls and dts to the interval variables (skip first ones since they are considered above)
            for (int i = 1; i < (int)u_seq.size(); ++i)
            {
                // VectorVertex& ukk = i < (int)u_seq.size() ? u_seq[i] : *interval.controls.back();    // ith control of the shooting interval -> for
                // move blocking

                if (nlp_fun.stage_cost)
                {
                    cost_terms.clear();
                    getStageFunctionEdgesIntermediateControlDtMultipleShooting(k + i, u_seq[i], _dt, u_prev, dt_prev, *nlp_fun.stage_cost,
                                                                               cost_terms);
                    for (BaseEdge::Ptr& edge : cost_terms) edges.addObjectiveEdge(edge);
                }
                if (nlp_fun.stage_equalities)
                {
                    eq_terms.clear();
                    getStageFunctionEdgesIntermediateControlDtMultipleShooting(k + i, u_seq[i], _dt, u_prev, dt_prev, *nlp_fun.stage_equalities,
                                                                               eq_terms);
                    for (BaseEdge::Ptr& edge : eq_terms) edges.addEqualityEdge(edge);
                }
                if (nlp_fun.stage_inequalities)
                {
                    ineq_terms.clear();
                    getStageFunctionEdgesIntermediateControlDtMultipleShooting(k + i, u_seq[i], _dt, u_prev, dt_prev, *nlp_fun.stage_inequalities,
                                                                               ineq_terms);
                    for (BaseEdge::Ptr& edge : ineq_terms) edges.addInequalityEdge(edge);
                }
            }
        }
        k += (int)u_seq.size();
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
    nlp_fun.getFinalControlDeviationEdges(n, _u_ref, _intervals.back().u_seq.back(), _dt, cost_terms, eq_terms, ineq_terms);
    for (BaseEdge::Ptr& edge : cost_terms) edges.addObjectiveEdge(edge);
    for (BaseEdge::Ptr& edge : eq_terms) edges.addEqualityEdge(edge);
    for (BaseEdge::Ptr& edge : ineq_terms) edges.addInequalityEdge(edge);
}

void MultipleShootingGrid::getStageFunctionEdgesIntermediateControlDtMultipleShooting(int k, VectorVertex& uk, ScalarVertex& dt, VectorVertex& u_prev,
                                                                                      ScalarVertex& dt_prev, const StageFunction& stage_fun,
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

#ifdef MESSAGE_SUPPORT
void MultipleShootingGrid::fromMessage(const messages::MultipleShootingGrid& message, std::stringstream* issues)
{
    if (message.n() < 2 && issues) *issues << "MultipleShootingGrid: Number of states must be greater than or equal 2.\n";
    if (message.dt() <= 0 && issues) *issues << "MultipleShootingGrid: Dt must be greater than 0.0.\n";

    setNRef(message.n());
    setDtRef(message.dt());
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
        if (issues) *issues << "MultipleShootingGrid: unknown integration method specified.\n";
        return;
    }
}

void MultipleShootingGrid::toMessage(messages::MultipleShootingGrid& message) const
{
    message.set_n(getNRef());
    message.set_dt(getDtRef());
    message.set_num_u_per_interval(_num_u_per_interv_ref);
    message.set_intermediate_x_constraints(_intermediate_x_constraints);
    message.set_warm_start(_warm_start);

    // fd collocation method
    if (_integrator) _integrator->fromMessage(*message.mutable_integrator());
}
#endif

}  // namespace corbo
