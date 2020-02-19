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

#ifndef SRC_OPTIMAL_CONTROL_INCLUDE_CORBO_OPTIMAL_CONTROL_STRUCTURED_OCP_DISCRETIZATION_GRIDS_NON_UNIFORM_FINITE_DIFFERENCES_VARIABLE_GRID_H_
#define SRC_OPTIMAL_CONTROL_INCLUDE_CORBO_OPTIMAL_CONTROL_STRUCTURED_OCP_DISCRETIZATION_GRIDS_NON_UNIFORM_FINITE_DIFFERENCES_VARIABLE_GRID_H_

#include <corbo-optimal-control/structured_ocp/discretization_grids/non_uniform_full_discretization_grid_base.h>

#include <memory>

namespace corbo {

class NonUniformFiniteDifferencesVariableGrid : public NonUniformFullDiscretizationGridBase
{
 public:
    using Ptr  = std::shared_ptr<NonUniformFiniteDifferencesVariableGrid>;
    using UPtr = std::unique_ptr<NonUniformFiniteDifferencesVariableGrid>;

    enum class GridAdaptStrategy { NoGridAdapt, TimeBasedSingleStep, RedundantControls };

    NonUniformFiniteDifferencesVariableGrid()          = default;
    virtual ~NonUniformFiniteDifferencesVariableGrid() = default;

    //! Return a newly created shared instance of the implemented class
    DiscretizationGridInterface::Ptr getInstance() const override { return std::make_shared<NonUniformFiniteDifferencesVariableGrid>(); }

    //! Get access to the associated factory
    static Factory<DiscretizationGridInterface>& getFactory() { return Factory<DiscretizationGridInterface>::instance(); }

    void setDtBounds(double dt_lb, double dt_ub);
    void setNmin(int n_min) { _n_min = n_min; }
    void disableGridAdaptation() { _grid_adapt = GridAdaptStrategy::NoGridAdapt; }
    void setGridAdaptTimeBasedSingleStep(int n_max, double dt_hyst_ratio = 0.1);
    void setGridAdaptRedundantControls(int n_max, int num_backup_nodes = 1, double epsilon = 1e-3);

    void setDtEqConstraint(bool active) { _dt_eq_constraint = active; }

#ifdef MESSAGE_SUPPORT
    void fromMessage(const messages::NonUniformFiniteDifferencesVariableGrid& message, std::stringstream* issues);
    void toMessage(messages::NonUniformFiniteDifferencesVariableGrid& message) const;

    void fromMessage(const messages::DiscretizationGrid& message, std::stringstream* issues) override
    {
        fromMessage(message.non_uniform_fd_variable_grid(), issues);
    }
    void toMessage(messages::DiscretizationGrid& message) const override { toMessage(*message.mutable_non_uniform_fd_variable_grid()); }
#endif

 protected:
    void createEdges(NlpFunctions& nlp_fun, OptimizationEdgeSet& edges, SystemDynamicsInterface::Ptr dynamics) override;

    bool isDtFixedIntended() const override { return false; }

    bool adaptGrid(bool new_run, NlpFunctions& nlp_fun) override;
    bool adaptGridTimeBasedSingleStep(NlpFunctions& nlp_fun);
    bool adaptGridRedundantControls(NlpFunctions& nlp_fun);

    bool isMovingHorizonWarmStartActive() const override { return false; }
    bool isGridAdaptActive() const override { return true; }

    // auto resize stuff
    GridAdaptStrategy _grid_adapt  = GridAdaptStrategy::NoGridAdapt;
    int _n_max                     = 1000;
    double _dt_hyst_ratio          = 0.1;
    int _n_min                     = 2;
    double _redundant_ctrl_epsilon = 1e-2;
    int _redundant_ctrl_backup     = 1;

    bool _adapt_first_iter = false;

    bool _dt_eq_constraint = false;
};

FACTORY_REGISTER_DISCRETIZATION_GRID(NonUniformFiniteDifferencesVariableGrid)

}  // namespace corbo

#endif  // SRC_OPTIMAL_CONTROL_INCLUDE_CORBO_OPTIMAL_CONTROL_STRUCTURED_OCP_DISCRETIZATION_GRIDS_NON_UNIFORM_FINITE_DIFFERENCES_VARIABLE_GRID_H_
