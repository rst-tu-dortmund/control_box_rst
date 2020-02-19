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

#ifndef SRC_OPTIMAL_CONTROL_INCLUDE_CORBO_OPTIMAL_CONTROL_STRUCTURED_OCP_DISCRETIZATION_GRIDS_SHOOTING_GRID_BASE_H_
#define SRC_OPTIMAL_CONTROL_INCLUDE_CORBO_OPTIMAL_CONTROL_STRUCTURED_OCP_DISCRETIZATION_GRIDS_SHOOTING_GRID_BASE_H_

#include <corbo-optimal-control/structured_ocp/discretization_grids/discretization_grid_interface.h>

#include <corbo-optimization/hyper_graph/scalar_vertex.h>
#include <corbo-optimization/hyper_graph/vector_vertex.h>

#include <corbo-numerics/integrator_interface.h>

#include <memory>

namespace corbo {

class ShootingGridBase : public DiscretizationGridInterface
{
 public:
    using Ptr  = std::shared_ptr<ShootingGridBase>;
    using UPtr = std::unique_ptr<ShootingGridBase>;

    struct ShootingInterval
    {
        VectorVertex s;
        std::vector<VectorVertex> u_seq;
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    };

    ShootingGridBase()          = default;
    virtual ~ShootingGridBase() = default;

    //! Return a newly created shared instance of the implemented class
    DiscretizationGridInterface::Ptr getInstance() const override = 0;

    GridUpdateResult update(const Eigen::VectorXd& x0, ReferenceTrajectoryInterface& xref, ReferenceTrajectoryInterface& uref, NlpFunctions& nlp_fun,
                            OptimizationEdgeSet& edges, SystemDynamicsInterface::Ptr dynamics, bool new_run, const Time& t,
                            ReferenceTrajectoryInterface* sref = nullptr, const Eigen::VectorXd* prev_u = nullptr, double prev_u_dt = 0,
                            ReferenceTrajectoryInterface* xinit = nullptr, ReferenceTrajectoryInterface* uinit = nullptr) override;

    double getFirstDt() const override { return getDt(); }
    double getFinalTime() const override { return double(getN() - 1) * getDt(); }

    bool hasConstantControls() const override { return true; }
    bool hasSingleDt() const override { return true; }
    bool isTimeVariableGrid() const override { return !isDtFixedIntended(); }
    bool isUniformGrid() const override { return true; }
    bool providesStateTrajectory() const override { return true; }  // TODO(roesmann): clarify

    bool getFirstControlInput(Eigen::VectorXd& u0) override;

    void getStateAndControlTimeSeries(TimeSeries::Ptr x_sequence, TimeSeries::Ptr u_sequence, double t_max = CORBO_INF_DBL) const override;

    void clear() override;

    bool isEmpty() const override { return _intervals.empty(); }
    bool isValid() const { return !_intervals.empty() && !_intervals.front().u_seq.empty(); }

    void setN(int n, bool try_resample = true) override
    {
        clear();  // resampling not yet implemented
        setNRef(n);
    }
    void setInitialDt(double dt) override { setDtRef(dt); }
    double getInitialDt() const override { return getDtRef(); }
    int getInitialN() const override { return getNRef(); }

    int getNRef() const { return _n_ref; }
    int getN() const override;
    void setNumControlsPerShootingInterval(int num_u_per_interv) { _num_u_per_interv_ref = num_u_per_interv; }
    void setNumControlsPerShootingInterval(int num_u_per_interv, bool intermediate_x_constraints)
    {
        setNumControlsPerShootingInterval(num_u_per_interv);
        _intermediate_x_constraints = intermediate_x_constraints;
    }
    void setNRef(int n);
    void setDtRef(double dt) { _dt_ref = dt; }
    double getDtRef() const { return _dt_ref; }
    double getDt() const { return _dt.value(); }
    void setWarmStart(bool active) { _warm_start = active; }

    void setConsiderIntermediateStateConstraints(bool active) { _intermediate_x_constraints = active; }

    void setNumericalIntegrator(NumericalIntegratorExplicitInterface::Ptr integrator) { _integrator = integrator; }

    void setXfFixed(const Eigen::Matrix<bool, -1, 1>& xf_fixed)
    {
        _xf_fixed = xf_fixed;
        setModified(true);
    }

    // VertexSet methods
    std::vector<VertexInterface*>& getActiveVertices() override { return _active_vertices; }
    void getVertices(std::vector<VertexInterface*>& vertices) override;

#ifdef MESSAGE_SUPPORT
    void fromMessage(const messages::DiscretizationGrid& message, std::stringstream* issues) override {}
    void toMessage(messages::DiscretizationGrid& message) const override {}
#endif

 protected:
    void initializeSequences(const Eigen::VectorXd& x0, const Eigen::VectorXd& xf, ReferenceTrajectoryInterface& uref, NlpFunctions& nlp_fun);
    void initializeSequences(const Eigen::VectorXd& x0, const Eigen::VectorXd& xf, ReferenceTrajectoryInterface& xref,
                             ReferenceTrajectoryInterface& uref, NlpFunctions& nlp_fun);
    void warmStartShifting(const Eigen::VectorXd& x0);

    virtual bool adaptGrid(bool new_run, NlpFunctions& nlp_fun) { return false; }  // A subclass might add a grid adaptation, returns true if adapted

    int findNearestShootingInterval(const Eigen::VectorXd& x0);

    void updateBounds(const NlpFunctions& nlp_fun);

    void resampleTrajectory(int n_new, NlpFunctions& nlp_fun);

    bool checkAndInitializeXfFixedFlags(int dim_x);

    virtual void createEdges(NlpFunctions& nlp_fun, OptimizationEdgeSet& edges, SystemDynamicsInterface::Ptr dynamics) = 0;

    virtual bool isDtFixedIntended() const { return true; }
    virtual bool isMovingHorizonWarmStartActive() const { return _warm_start; }
    virtual bool isGridAdaptActive() const { return false; }

    void computeActiveVertices() override;

 protected:
    virtual bool isXfShootingNode() const { return true; }

    NumericalIntegratorExplicitInterface::Ptr _integrator;

    std::vector<ShootingInterval> _intervals;
    ScalarVertex _dt;
    PartiallyFixedVectorVertex _xf;
    std::vector<VertexInterface*> _active_vertices;

    int _num_u_per_interv_ref = 1;  // 1 -> full discretization
    int _n_ref                = 11;
    int _n_adapt              = 0;  // if adaption is on and warmstart of, we might use this n instead of n_ref (only if n_adapt > 0)
    double _dt_ref            = 0.1;

    bool _warm_start = false;
    bool _first_run  = true;

    // might be required if the last dt should be fixed or if dt is not fixed
    Eigen::Matrix<bool, -1, 1> _xf_fixed;
    double _dt_lb = 0;
    double _dt_ub = CORBO_INF_DBL;

    bool _full_discretization        = true;
    bool _intermediate_x_constraints = false;
};

}  // namespace corbo

#endif  // SRC_OPTIMAL_CONTROL_INCLUDE_CORBO_OPTIMAL_CONTROL_STRUCTURED_OCP_DISCRETIZATION_GRIDS_SHOOTING_GRID_BASE_H_
