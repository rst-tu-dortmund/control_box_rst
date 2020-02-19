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

#ifndef SRC_OPTIMAL_CONTROL_INCLUDE_CORBO_OPTIMAL_CONTROL_STRUCTURED_OCP_DISCRETIZATION_GRIDS_NON_UNIFORM_FULL_DISCRETIZATION_GRID_BASE_H_
#define SRC_OPTIMAL_CONTROL_INCLUDE_CORBO_OPTIMAL_CONTROL_STRUCTURED_OCP_DISCRETIZATION_GRIDS_NON_UNIFORM_FULL_DISCRETIZATION_GRID_BASE_H_

#include <corbo-optimal-control/structured_ocp/discretization_grids/discretization_grid_interface.h>

#include <corbo-optimization/hyper_graph/scalar_vertex.h>
#include <corbo-optimization/hyper_graph/vector_vertex.h>

#include <corbo-numerics/finite_differences_collocation.h>

#include <memory>

namespace corbo {

class NonUniformFullDiscretizationGridBase : public DiscretizationGridInterface
{
 public:
    using Ptr  = std::shared_ptr<NonUniformFullDiscretizationGridBase>;
    using UPtr = std::unique_ptr<NonUniformFullDiscretizationGridBase>;

    enum class CostIntegrationRule { LeftSum, TrapezoidalRule };

    NonUniformFullDiscretizationGridBase()          = default;
    virtual ~NonUniformFullDiscretizationGridBase() = default;

    //! Return a newly created shared instance of the implemented class
    DiscretizationGridInterface::Ptr getInstance() const override = 0;

    GridUpdateResult update(const Eigen::VectorXd& x0, ReferenceTrajectoryInterface& xref, ReferenceTrajectoryInterface& uref, NlpFunctions& nlp_fun,
                            OptimizationEdgeSet& edges, SystemDynamicsInterface::Ptr dynamics, bool new_run, const Time& t,
                            ReferenceTrajectoryInterface* sref = nullptr, const Eigen::VectorXd* prev_u = nullptr, double prev_u_dt = 0,
                            ReferenceTrajectoryInterface* xinit = nullptr, ReferenceTrajectoryInterface* uinit = nullptr) override;

    double getFirstDt() const override { return !_dt_seq.empty() ? _dt_seq.front().value() : 0; }
    double getFinalTime() const override;

    bool hasConstantControls() const override { return true; }
    bool hasSingleDt() const override { return false; }
    bool isTimeVariableGrid() const override { return !isDtFixedIntended(); }
    bool isUniformGrid() const override { return false; }
    bool providesStateTrajectory() const override { return true; }

    bool getFirstControlInput(Eigen::VectorXd& u0) override;

    void getStateAndControlTimeSeries(TimeSeries::Ptr x_sequence, TimeSeries::Ptr u_sequence, double t_max = CORBO_INF_DBL) const override;

    void clear() override;

    bool isEmpty() const override { return _u_seq.empty(); }
    bool isValid() const { return (_x_seq.size() == _u_seq.size()) && _u_seq.size() == _dt_seq.size(); }

    void setN(int n, bool try_resample = true) override
    {
        clear();  // resampling not yet implemented
        setNRef(n);
    }
    void setInitialDt(double dt) override { setDtRef(dt); }
    double getInitialDt() const override { return getDtRef(); }
    int getInitialN() const override { return getNRef(); }

    int getNRef() const { return _n_ref; }
    int getN() const override { return _x_seq.size() + 1; }
    void setNRef(int n);
    void setDtRef(double dt) { _dt_ref = dt; }
    double getDtRef() const { return _dt_ref; }
    void setWarmStart(bool active) { _warm_start = active; }

    void getDts(std::vector<double>& dts) const;

    void setFiniteDifferencesCollocationMethod(FiniteDifferencesCollocationInterface::Ptr fd_eval) { _fd_eval = fd_eval; }
    void setCostIntegrationRule(CostIntegrationRule integration) { _cost_integration = integration; }

    void setXfFixed(const Eigen::Matrix<bool, -1, 1>& xf_fixed)
    {
        _xf_fixed = xf_fixed;
        setModified(true);
    }

    // VertexSet methods
    std::vector<VertexInterface*>& getActiveVertices() override { return _active_vertices; }
    void getVertices(std::vector<VertexInterface*>& vertices) override;  // TODO(roesmann) add dt, just in case

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

    int findNearestState(const Eigen::VectorXd& x0);

    void updateBounds(const NlpFunctions& nlp_fun);

    bool checkAndInitializeXfFixedFlags(int dim_x);

    virtual void createEdges(NlpFunctions& nlp_fun, OptimizationEdgeSet& edges, SystemDynamicsInterface::Ptr dynamics) = 0;

    virtual bool isDtFixedIntended() const { return false; }
    virtual bool isMovingHorizonWarmStartActive() const { return _warm_start; }
    virtual bool isGridAdaptActive() const { return false; }

    void computeActiveVertices() override;

 protected:
    FiniteDifferencesCollocationInterface::Ptr _fd_eval = std::make_shared<CrankNicolsonDiffCollocation>();

    std::vector<VectorVertex> _x_seq;
    std::vector<VectorVertex> _u_seq;
    std::vector<ScalarVertex> _dt_seq;
    PartiallyFixedVectorVertex _xf;
    std::vector<VertexInterface*> _active_vertices;

    int _n_ref       = 11;
    int _n_adapt     = 0;  // if adaption is on and warmstart of, we might use this n instead of n_ref (only if n_adapt > 0)
    double _dt_ref   = 0.1;
    bool _warm_start = false;
    bool _first_run  = true;

    // might be required if the last dt should be fixed or if dt is not fixed
    Eigen::Matrix<bool, -1, 1> _xf_fixed;
    double _dt_lb = 0;
    double _dt_ub = CORBO_INF_DBL;

    CostIntegrationRule _cost_integration;
};

}  // namespace corbo

#endif  // SRC_OPTIMAL_CONTROL_INCLUDE_CORBO_OPTIMAL_CONTROL_STRUCTURED_OCP_DISCRETIZATION_GRIDS_NON_UNIFORM_FULL_DISCRETIZATION_GRID_BASE_H_
