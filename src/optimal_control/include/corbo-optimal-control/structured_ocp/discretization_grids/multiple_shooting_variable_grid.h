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

#ifndef SRC_OPTIMAL_CONTROL_INCLUDE_CORBO_OPTIMAL_CONTROL_STRUCTURED_OCP_DISCRETIZATION_GRIDS_MULTIPLE_SHOOTING_VARIABLE_GRID_H_
#define SRC_OPTIMAL_CONTROL_INCLUDE_CORBO_OPTIMAL_CONTROL_STRUCTURED_OCP_DISCRETIZATION_GRIDS_MULTIPLE_SHOOTING_VARIABLE_GRID_H_

#include <corbo-optimal-control/structured_ocp/discretization_grids/multiple_shooting_grid.h>

#include <memory>

namespace corbo {

class MultipleShootingVariableGrid : public MultipleShootingGrid
{
 public:
    using Ptr  = std::shared_ptr<MultipleShootingVariableGrid>;
    using UPtr = std::unique_ptr<MultipleShootingVariableGrid>;

    enum class GridAdaptStrategy { NoGridAdapt, TimeBasedSingleStep, TimeBasedAggressiveEstimate, SimpleShrinkingHorizon };

    MultipleShootingVariableGrid()          = default;
    virtual ~MultipleShootingVariableGrid() = default;

    //! Return a newly created shared instance of the implemented class
    DiscretizationGridInterface::Ptr getInstance() const override { return std::make_shared<MultipleShootingVariableGrid>(); }

    //! Get access to the associated factory
    static Factory<DiscretizationGridInterface>& getFactory() { return Factory<DiscretizationGridInterface>::instance(); }

    void setDtBounds(double dt_lb, double dt_ub);
    void setNmin(int n_min) { _n_min = n_min; }
    void disableGridAdaptation() { _grid_adapt = GridAdaptStrategy::NoGridAdapt; }
    void setGridAdaptTimeBasedSingleStep(int n_max, double dt_hyst_ratio = 0.1);
    void setGridAdaptTimeBasedAggressiveEstimate(int n_max, double dt_hyst_ratio = 0.1);
    void setGridAdaptSimpleShrinkingHorizon();

#ifdef MESSAGE_SUPPORT
    void fromMessage(const messages::MultipleShootingVariableGrid& message, std::stringstream* issues);
    void toMessage(messages::MultipleShootingVariableGrid& message) const;

    void fromMessage(const messages::DiscretizationGrid& message, std::stringstream* issues) override
    {
        fromMessage(message.multiple_shooting_variable_grid(), issues);
    }
    void toMessage(messages::DiscretizationGrid& message) const override { toMessage(*message.mutable_multiple_shooting_variable_grid()); }
#endif

 protected:
    bool isMovingHorizonWarmStartActive() const override { return false; }
    bool isGridAdaptActive() const override { return true; }
    bool isDtFixedIntended() const override { return false; }

    bool adaptGrid(bool new_run, NlpFunctions& nlp_fun) override;
    bool adaptGridTimeBasedSingleStep(NlpFunctions& nlp_fun);
    bool adaptGridTimeBasedAggressiveEstimate(NlpFunctions& nlp_fun);
    bool adaptGridSimpleShrinkingHorizon(NlpFunctions& nlp_fun);

    // auto resize stuff
    GridAdaptStrategy _grid_adapt = GridAdaptStrategy::NoGridAdapt;
    int _n_max                    = 1000;
    double _dt_hyst_ratio         = 0.1;
    int _n_min                    = 2;
    bool _adapt_first_iter        = false;
};

FACTORY_REGISTER_DISCRETIZATION_GRID(MultipleShootingVariableGrid)

}  // namespace corbo

#endif  // SRC_OPTIMAL_CONTROL_INCLUDE_CORBO_OPTIMAL_CONTROL_STRUCTURED_OCP_DISCRETIZATION_GRIDS_MULTIPLE_SHOOTING_VARIABLE_GRID_H_
