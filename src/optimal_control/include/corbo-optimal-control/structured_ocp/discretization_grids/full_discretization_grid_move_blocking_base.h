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

#ifndef SRC_OPTIMAL_CONTROL_INCLUDE_CORBO_OPTIMAL_CONTROL_STRUCTURED_OCP_DISCRETIZATION_GRIDS_FULL_DISCRETIZATION_GRID_MOVE_BLOCKING_BASE_H_
#define SRC_OPTIMAL_CONTROL_INCLUDE_CORBO_OPTIMAL_CONTROL_STRUCTURED_OCP_DISCRETIZATION_GRIDS_FULL_DISCRETIZATION_GRID_MOVE_BLOCKING_BASE_H_

#include <corbo-optimal-control/structured_ocp/discretization_grids/full_discretization_grid_base.h>

namespace corbo {

class FullDiscretizationGridMoveBlockingBase : public FullDiscretizationGridBase
{
 public:
    using Ptr  = std::shared_ptr<FullDiscretizationGridMoveBlockingBase>;
    using UPtr = std::unique_ptr<FullDiscretizationGridMoveBlockingBase>;

    FullDiscretizationGridMoveBlockingBase()          = default;
    virtual ~FullDiscretizationGridMoveBlockingBase() = default;

    //! Return a newly created shared instance of the implemented class
    DiscretizationGridInterface::Ptr getInstance() const override = 0;

    void getStateAndControlTimeSeries(TimeSeries::Ptr x_sequence, TimeSeries::Ptr u_sequence, double t_max = CORBO_INF_DBL) const override;

    bool isValid() const override;

    void getBlockingMatrix(Eigen::Ref<Eigen::VectorXi> B) { B = _blocking_vector; }
    void setBlockingMatrix(const Eigen::Ref<const Eigen::VectorXi>& B) { _blocking_vector = B; }

    void setWarmStartShiftU(bool shift_u) { _warm_start_shift_u = shift_u; }

#ifdef MESSAGE_SUPPORT
    void fromMessage(const messages::DiscretizationGrid& message, std::stringstream* issues) override {}
    void toMessage(messages::DiscretizationGrid& message) const override {}
#endif

 protected:
    Eigen::VectorXi _blocking_vector;
    bool _warm_start_shift_u = false;

    void initializeSequences(const Eigen::VectorXd& x0, const Eigen::VectorXd& xf, ReferenceTrajectoryInterface& uref,
                             NlpFunctions& nlp_fun) override;
    void initializeSequences(const Eigen::VectorXd& x0, const Eigen::VectorXd& xf, ReferenceTrajectoryInterface& xref,
                             ReferenceTrajectoryInterface& uref, NlpFunctions& nlp_fun) override;

    void warmStartShifting(const Eigen::VectorXd& x0) override;

    void resampleTrajectory(int n_new) override;

    void createEdges(NlpFunctions& nlp_fun, OptimizationEdgeSet& edges, SystemDynamicsInterface::Ptr dynamics) override = 0;

    void computeActiveVertices() override;
};

}  // namespace corbo

#endif  // SRC_OPTIMAL_CONTROL_INCLUDE_CORBO_OPTIMAL_CONTROL_STRUCTURED_OCP_DISCRETIZATION_GRIDS_FULL_DISCRETIZATION_GRID_MOVE_BLOCKING_BASE_H_
