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
 *  Authors: Maximilian Krämer, Christoph Rösmann
 *********************************************************************/

#ifndef SRC_OPTIMAL_CONTROL_INCLUDE_CORBO_OPTIMAL_CONTROL_STRUCTURED_OCP_DISCRETIZATION_GRIDS_FINITE_DIFFERENCES_GRID_MOVE_BLOCKING_H_
#define SRC_OPTIMAL_CONTROL_INCLUDE_CORBO_OPTIMAL_CONTROL_STRUCTURED_OCP_DISCRETIZATION_GRIDS_FINITE_DIFFERENCES_GRID_MOVE_BLOCKING_H_

#include <corbo-optimal-control/structured_ocp/discretization_grids/full_discretization_grid_move_blocking_base.h>

namespace corbo {

class FiniteDifferencesGridMoveBlocking : public FullDiscretizationGridMoveBlockingBase
{
 public:
    using Ptr  = std::shared_ptr<FiniteDifferencesGridMoveBlocking>;
    using UPtr = std::unique_ptr<FiniteDifferencesGridMoveBlocking>;

    FiniteDifferencesGridMoveBlocking()          = default;
    virtual ~FiniteDifferencesGridMoveBlocking() = default;

    //! Return a newly created shared instance of the implemented class
    DiscretizationGridInterface::Ptr getInstance() const override { return std::make_shared<FiniteDifferencesGridMoveBlocking>(); }

    //! Get access to the associated factory
    static Factory<DiscretizationGridInterface>& getFactory() { return Factory<DiscretizationGridInterface>::instance(); }

#ifdef MESSAGE_SUPPORT
    void fromMessage(const messages::FiniteDifferencesGridMoveBlocking& message, std::stringstream* issues);
    void toMessage(messages::FiniteDifferencesGridMoveBlocking& message) const;

    void fromMessage(const messages::DiscretizationGrid& message, std::stringstream* issues) override
    {
        fromMessage(message.finite_differences_grid_move_blocking(), issues);
    }
    void toMessage(messages::DiscretizationGrid& message) const override { toMessage(*message.mutable_finite_differences_grid_move_blocking()); }
#endif

 protected:
    void createEdges(NlpFunctions& nlp_fun, OptimizationEdgeSet& edges, SystemDynamicsInterface::Ptr dynamics) override;
};

FACTORY_REGISTER_DISCRETIZATION_GRID(FiniteDifferencesGridMoveBlocking)

}  // namespace corbo

#endif  // SRC_OPTIMAL_CONTROL_INCLUDE_CORBO_OPTIMAL_CONTROL_STRUCTURED_OCP_DISCRETIZATION_GRIDS_FINITE_DIFFERENCES_GRID_MOVE_BLOCKING_H_
