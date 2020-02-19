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

#ifndef SRC_OPTIMAL_CONTROL_INCLUDE_CORBO_OPTIMAL_CONTROL_STRUCTURED_OCP_DISCRETIZATION_GRIDS_MULTIPLE_SHOOTING_GRID_H_
#define SRC_OPTIMAL_CONTROL_INCLUDE_CORBO_OPTIMAL_CONTROL_STRUCTURED_OCP_DISCRETIZATION_GRIDS_MULTIPLE_SHOOTING_GRID_H_

#include <corbo-optimal-control/structured_ocp/discretization_grids/shooting_grid_base.h>

#include <memory>

namespace corbo {

class MultipleShootingGrid : public ShootingGridBase
{
 public:
    using Ptr  = std::shared_ptr<MultipleShootingGrid>;
    using UPtr = std::unique_ptr<MultipleShootingGrid>;

    MultipleShootingGrid()          = default;
    virtual ~MultipleShootingGrid() = default;

    //! Return a newly created shared instance of the implemented class
    DiscretizationGridInterface::Ptr getInstance() const override { return std::make_shared<MultipleShootingGrid>(); }

    //! Get access to the associated factory
    static Factory<DiscretizationGridInterface>& getFactory() { return Factory<DiscretizationGridInterface>::instance(); }

#ifdef MESSAGE_SUPPORT
    void fromMessage(const messages::MultipleShootingGrid& message, std::stringstream* issues);
    void toMessage(messages::MultipleShootingGrid& message) const;

    void fromMessage(const messages::DiscretizationGrid& message, std::stringstream* issues) override
    {
        fromMessage(message.multiple_shooting_grid(), issues);
    }
    void toMessage(messages::DiscretizationGrid& message) const override { toMessage(*message.mutable_multiple_shooting_grid()); }
#endif

 protected:
    void createEdges(NlpFunctions& nlp_fun, OptimizationEdgeSet& edges, SystemDynamicsInterface::Ptr dynamics) override;
    void getStageFunctionEdgesIntermediateControlDtMultipleShooting(int k, VectorVertex& uk, ScalarVertex& dt, VectorVertex& u_prev,
                                                                    ScalarVertex& dt_prev, const StageFunction& stage_fun,
                                                                    std::vector<BaseEdge::Ptr>& edges);

    bool isDtFixedIntended() const override { return true; }
};

FACTORY_REGISTER_DISCRETIZATION_GRID(MultipleShootingGrid)

}  // namespace corbo

#endif  // SRC_OPTIMAL_CONTROL_INCLUDE_CORBO_OPTIMAL_CONTROL_STRUCTURED_OCP_DISCRETIZATION_GRIDS_MULTIPLE_SHOOTING_GRID_H_
