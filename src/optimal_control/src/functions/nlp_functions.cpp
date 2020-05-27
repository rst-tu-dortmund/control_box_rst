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

#include <corbo-optimal-control/functions/nlp_functions.h>

#include <corbo-optimal-control/structured_ocp/discretization_grids/discretization_grid_interface.h>

namespace corbo {

bool NlpFunctions::update(int n, double t, ReferenceTrajectoryInterface& xref, ReferenceTrajectoryInterface& uref, ReferenceTrajectoryInterface* sref,
                          bool single_dt, const Eigen::VectorXd& x0, const std::vector<double>& dts, const DiscretizationGridInterface* grid)
{
    bool dimension_modified = false;

    if (stage_preprocessor) dimension_modified |= stage_preprocessor->update(n, t, xref, uref, sref, single_dt, x0, dts, grid);
    if (stage_cost) dimension_modified |= stage_cost->update(n, t, xref, uref, sref, single_dt, x0, stage_preprocessor, dts, grid);
    if (final_stage_cost) dimension_modified |= final_stage_cost->update(n, t, xref, uref, sref, single_dt, x0, stage_preprocessor, dts, grid);
    if (stage_equalities) dimension_modified |= stage_equalities->update(n, t, xref, uref, sref, single_dt, x0, stage_preprocessor, dts, grid);
    if (stage_inequalities) dimension_modified |= stage_inequalities->update(n, t, xref, uref, sref, single_dt, x0, stage_preprocessor, dts, grid);
    if (final_stage_constraints)
        dimension_modified |= final_stage_constraints->update(n, t, xref, uref, sref, single_dt, x0, final_stage_cost, stage_preprocessor, dts, grid);

    return dimension_modified;
}

void NlpFunctions::checkAndInitializeBoundDimensions(int x_dim, int u_dim)
{
    if (x_lb.size() == 0)
        x_lb.setConstant(x_dim, -CORBO_INF_DBL);
    else if (x_lb.size() != x_dim)
        PRINT_ERROR_NAMED("Error in lower state bounds: dimensions mismatch");

    if (x_ub.size() == 0)
        x_ub.setConstant(x_dim, CORBO_INF_DBL);
    else if (x_ub.size() != x_dim)
        PRINT_ERROR_NAMED("Error in upper state bounds: dimensions mismatch");

    if (u_lb.size() == 0)
        u_lb.setConstant(u_dim, -CORBO_INF_DBL);
    else if (u_lb.size() != u_dim)
        PRINT_ERROR_NAMED("Error in lower control input bounds: dimensions mismatch");

    if (u_ub.size() == 0)
        u_ub.setConstant(u_dim, CORBO_INF_DBL);
    else if (u_ub.size() != u_dim)
        PRINT_ERROR_NAMED("Error in upper control input bounds: dimensions mismatch");
}

void NlpFunctions::getNonIntegralStageFunctionEdges(int k, VectorVertex& xk, VectorVertex& uk, ScalarVertex& dt, VectorVertex& u_prev,
                                                    ScalarVertex& u_prev_dt, const StageFunction& stage_fun, std::vector<BaseEdge::Ptr>& edges)
{
    int dim = stage_fun.getNonIntegralStateTermDimension(k);
    if (dim > 0)
    {
        using Edge = UnaryVectorVertexEdge<StageFunction, &StageFunction::computeNonIntegralStateTerm>;
        Edge::Ptr edge =
            std::make_shared<Edge>(dim, k, stage_fun, xk, stage_fun.isLinearNonIntegralStateTerm(k), stage_fun.isLsqFormNonIntegralStateTerm(k));
        edges.push_back(edge);
    }

    dim = stage_fun.getNonIntegralControlTermDimension(k);
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

    dim = stage_fun.getNonIntegralDtTermDimension(k);
    if (dim > 0)
    {
        using Edge = UnaryScalarVertexEdge<StageFunction, &StageFunction::computeNonIntegralDtTerm>;
        Edge::Ptr edge =
            std::make_shared<Edge>(dim, k, stage_fun, dt, stage_fun.isLinearNonIntegralDtTerm(k), stage_fun.isLsqFormNonIntegralDtTerm(k));
        edges.push_back(edge);
    }

    dim = stage_fun.getNonIntegralStateControlTermDimension(k);
    if (dim > 0)
    {
        using Edge     = BinaryVectorVertexEdge<StageFunction, &StageFunction::computeNonIntegralStateControlTerm>;
        Edge::Ptr edge = std::make_shared<Edge>(dim, k, stage_fun, xk, uk, false, false);
        edges.push_back(edge);
    }

    dim = stage_fun.getNonIntegralControlDeviationTermDimension(k);
    if (dim > 0)
    {
        using Edge     = TernaryVectorScalarVertexEdge<StageFunction, &StageFunction::computeNonIntegralControlDeviationTerm>;
        Edge::Ptr edge = std::make_shared<Edge>(dim, k, stage_fun, uk, u_prev, u_prev_dt, false, false);
        edges.push_back(edge);
    }

    dim = stage_fun.getNonIntegralStateControlDtTermDimension(k);
    if (dim > 0)
    {
        using Edge     = TernaryVectorScalarVertexEdge<StageFunction, &StageFunction::computeNonIntegralStateControlDtTerm>;
        Edge::Ptr edge = std::make_shared<Edge>(dim, k, stage_fun, xk, uk, dt, false, false);
        edges.push_back(edge);
    }
}

void NlpFunctions::getNonIntegralStageFunctionEdges(int k, VectorVertex& xk, VectorVertex& uk, ScalarVertex& dt, VectorVertex& u_prev,
                                                    ScalarVertex& u_prev_dt, std::vector<BaseEdge::Ptr>& cost_edges,
                                                    std::vector<BaseEdge::Ptr>& eq_edges, std::vector<BaseEdge::Ptr>& ineq_edges)
{
    if (stage_cost)
    {
        getNonIntegralStageFunctionEdges(k, xk, uk, dt, u_prev, u_prev_dt, *stage_cost, cost_edges);
    }
    if (stage_equalities)
    {
        getNonIntegralStageFunctionEdges(k, xk, uk, dt, u_prev, u_prev_dt, *stage_equalities, eq_edges);
    }
    if (stage_inequalities)
    {
        getNonIntegralStageFunctionEdges(k, xk, uk, dt, u_prev, u_prev_dt, *stage_inequalities, ineq_edges);
    }
}

void NlpFunctions::getFinalControlDeviationEdges(int n, VectorVertex& u_ref, VectorVertex& u_prev, ScalarVertex& u_prev_dt,
                                                 std::vector<BaseEdge::Ptr>& cost_edges, std::vector<BaseEdge::Ptr>& eq_edges,
                                                 std::vector<BaseEdge::Ptr>& ineq_edges)
{
    if (stage_cost)
    {
        int dim = stage_cost->getNonIntegralControlDeviationTermDimension(n);
        if (dim > 0)
        {
            using Edge     = TernaryVectorScalarVertexEdge<StageFunction, &StageFunction::computeNonIntegralControlDeviationTerm>;
            Edge::Ptr edge = std::make_shared<Edge>(dim, n, *stage_cost, u_ref, u_prev, u_prev_dt, false, false);
            cost_edges.push_back(edge);
        }
    }
    if (stage_equalities)
    {
        int dim = stage_equalities->getNonIntegralControlDeviationTermDimension(n);
        if (dim > 0)
        {
            using Edge     = TernaryVectorScalarVertexEdge<StageFunction, &StageFunction::computeNonIntegralControlDeviationTerm>;
            Edge::Ptr edge = std::make_shared<Edge>(dim, n, *stage_equalities, u_ref, u_prev, u_prev_dt, false, false);
            eq_edges.push_back(edge);
        }
    }
    if (stage_inequalities)
    {
        int dim = stage_inequalities->getNonIntegralControlDeviationTermDimension(n);
        if (dim > 0)
        {
            using Edge     = TernaryVectorScalarVertexEdge<StageFunction, &StageFunction::computeNonIntegralControlDeviationTerm>;
            Edge::Ptr edge = std::make_shared<Edge>(dim, n, *stage_inequalities, u_ref, u_prev, u_prev_dt, false, false);
            ineq_edges.push_back(edge);
        }
    }
}

BaseEdge::Ptr NlpFunctions::getFinalStateCostEdge(int k, VectorVertex& xf)
{
    if (final_stage_cost)
    {
        int dim = final_stage_cost->getNonIntegralStateTermDimension(k);
        if (dim > 0)
        {
            using Edge     = UnaryVectorVertexEdge<StageFunction, &StageFunction::computeNonIntegralStateTerm>;
            Edge::Ptr edge = std::make_shared<Edge>(dim, k, *final_stage_cost, xf, final_stage_cost->isLinearNonIntegralStateTerm(k),
                                                    final_stage_cost->isLsqFormNonIntegralStateTerm(k));
            return edge;
        }
    }
    return {};
}

BaseEdge::Ptr NlpFunctions::getFinalStateConstraintEdge(int k, VectorVertex& xf)
{
    if (final_stage_constraints)
    {
        int dim = final_stage_constraints->getNonIntegralStateTermDimension(k);
        if (dim > 0)
        {
            using Edge     = UnaryVectorVertexEdge<StageFunction, &StageFunction::computeNonIntegralStateTerm>;
            Edge::Ptr edge = std::make_shared<Edge>(dim, k, *final_stage_constraints, xf, final_stage_constraints->isLinearNonIntegralStateTerm(k),
                                                    final_stage_constraints->isLsqFormNonIntegralStateTerm(k));
            return edge;
        }
    }
    return {};
}

}  // namespace corbo
