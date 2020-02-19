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

#ifndef SRC_OPTIMAL_CONTROL_INCLUDE_CORBO_OPTIMAL_CONTROL_STRUCTURED_OCP_DISCRETIZATION_GRIDS_FINITE_DIFFERENCES_VARIABLE_GRID_H_
#define SRC_OPTIMAL_CONTROL_INCLUDE_CORBO_OPTIMAL_CONTROL_STRUCTURED_OCP_DISCRETIZATION_GRIDS_FINITE_DIFFERENCES_VARIABLE_GRID_H_

#include <corbo-optimal-control/structured_ocp/discretization_grids/finite_differences_grid.h>

#include <memory>

namespace corbo {

class FiniteDifferencesVariableGrid : public FiniteDifferencesGrid
{
 public:
    using Ptr  = std::shared_ptr<FiniteDifferencesVariableGrid>;
    using UPtr = std::unique_ptr<FiniteDifferencesVariableGrid>;

    enum class GridAdaptStrategy { NoGridAdapt, TimeBasedSingleStep, TimeBasedAggressiveEstimate, SimpleShrinkingHorizon };

    FiniteDifferencesVariableGrid()          = default;
    virtual ~FiniteDifferencesVariableGrid() = default;

    //! Return a newly created shared instance of the implemented class
    DiscretizationGridInterface::Ptr getInstance() const override { return std::make_shared<FiniteDifferencesVariableGrid>(); }

    //! Get access to the associated factory
    static Factory<DiscretizationGridInterface>& getFactory() { return Factory<DiscretizationGridInterface>::instance(); }

    void setDtBounds(double dt_lb, double dt_ub);
    void setNmin(int n_min) { _n_min = n_min; }
    void disableGridAdaptation() { _grid_adapt = GridAdaptStrategy::NoGridAdapt; }
    void setGridAdaptTimeBasedSingleStep(int n_max, double dt_hyst_ratio = 0.1, bool adapt_first_iter = false);
    void setGridAdaptTimeBasedAggressiveEstimate(int n_max, double dt_hyst_ratio = 0.1, bool adapt_first_iter = false);
    void setGridAdaptSimpleShrinkingHorizon(bool adapt_first_iter = false);

#ifdef MESSAGE_SUPPORT
    void fromMessage(const messages::FiniteDifferencesVariableGrid& message, std::stringstream* issues);
    void toMessage(messages::FiniteDifferencesVariableGrid& message) const;

    void fromMessage(const messages::DiscretizationGrid& message, std::stringstream* issues) override
    {
        fromMessage(message.finite_differences_variable_grid(), issues);
    }
    void toMessage(messages::DiscretizationGrid& message) const override { toMessage(*message.mutable_finite_differences_variable_grid()); }
#endif

 protected:
    bool isDtFixedIntended() const override { return false; }

    bool adaptGrid(bool new_run, NlpFunctions& nlp_fun) override;
    bool adaptGridTimeBasedSingleStep(NlpFunctions& nlp_fun);
    bool adaptGridTimeBasedAggressiveEstimate(NlpFunctions& nlp_fun);
    bool adaptGridSimpleShrinkingHorizon(NlpFunctions& nlp_fun);

    bool isMovingHorizonWarmStartActive() const override { return false; }
    bool isGridAdaptActive() const override { return true; }

    // auto resize stuff
    GridAdaptStrategy _grid_adapt = GridAdaptStrategy::NoGridAdapt;
    int _n_max                    = 1000;
    double _dt_hyst_ratio         = 0.1;
    int _n_min                    = 2;
    bool _adapt_first_iter        = false;
};

FACTORY_REGISTER_DISCRETIZATION_GRID(FiniteDifferencesVariableGrid)

}  // namespace corbo

#endif  // SRC_OPTIMAL_CONTROL_INCLUDE_CORBO_OPTIMAL_CONTROL_STRUCTURED_OCP_DISCRETIZATION_GRIDS_FINITE_DIFFERENCES_VARIABLE_GRID_H_
