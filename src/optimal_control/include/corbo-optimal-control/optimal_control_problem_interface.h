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

#ifndef SRC_OPTIMAL_CONTROL_INCLUDE_CORBO_OPTIMAL_CONTROL_OPTIMAL_CONTROL_PROBLEM_INTERFACE_H_
#define SRC_OPTIMAL_CONTROL_INCLUDE_CORBO_OPTIMAL_CONTROL_OPTIMAL_CONTROL_PROBLEM_INTERFACE_H_

#include <corbo-core/factory.h>
#include <corbo-core/reference_trajectory.h>
#include <corbo-core/signal_target_interface.h>
#include <corbo-core/time.h>
#include <corbo-core/types.h>

#include <corbo-optimal-control/statistics.h>

#ifdef MESSAGE_SUPPORT
#include <corbo-communication/messages/optimal_control/optimal_control_problem.pb.h>
#endif

#include <memory>

namespace corbo {

class OptimalControlProblemInterface
{
 public:
    using Ptr           = std::shared_ptr<OptimalControlProblemInterface>;
    using UPtr          = std::unique_ptr<OptimalControlProblemInterface>;
    using StateVector   = Eigen::VectorXd;
    using ControlVector = Eigen::VectorXd;

    virtual ~OptimalControlProblemInterface() {}

    virtual Ptr getInstance() const = 0;

    virtual void reset() = 0;

    virtual int getControlInputDimension() const = 0;
    virtual int getStateDimension() const        = 0;
    virtual int getN() const                     = 0;

    virtual bool initialize() { return true; }

    virtual bool compute(const StateVector& x, ReferenceTrajectoryInterface& xref, ReferenceTrajectoryInterface& uref,
                         ReferenceTrajectoryInterface* sref, const Time& t, bool new_run = true, SignalTargetInterface* signal_target = nullptr,
                         ReferenceTrajectoryInterface* xinit = nullptr, ReferenceTrajectoryInterface* uinit = nullptr,
                         const std::string& ns = "") = 0;

    virtual bool getFirstControlInput(ControlVector& u0) const = 0;
    virtual double getFirstDt() const                          = 0;

    virtual bool isConstantControlAction() const = 0;

    virtual double getCurrentObjectiveValue() { return -1; }

    // TODO(roesmann): someday we should change this to return a time-series, e.g. a std::pair<TS,TS>, because pre-allocating a shared pointer is ugly
    virtual void getTimeSeries(TimeSeries::Ptr x_sequence, TimeSeries::Ptr u_sequence, double t_max = CORBO_INF_DBL) = 0;
    // virtual std::pair<TimeSeries::Ptr, TimeSeries::Ptr> getTimeSeries(bool states = true, bool controls = true, double t_max = CORBO_INF_DBL) = 0;

    virtual bool providesFutureControls() const = 0;
    virtual bool providesFutureStates() const   = 0;

    virtual void setPreviousControlInput(const Eigen::Ref<const ControlVector>& u_prev, double dt) {}
    virtual void setPreviousControlInputDt(double dt) {}

    virtual OptimalControlProblemStatistics::Ptr getStatistics() const { return {}; }

#ifdef MESSAGE_SUPPORT
    virtual void toMessage(corbo::messages::OptimalControlProblem& message) const {}
    virtual void fromMessage(const corbo::messages::OptimalControlProblem& message, std::stringstream* issues = nullptr) {}
#endif
};

using OptimalControlProgramFactory = Factory<OptimalControlProblemInterface>;
#define FACTORY_REGISTER_OCP(type) FACTORY_REGISTER_OBJECT(type, OptimalControlProblemInterface)

}  // namespace corbo

#endif  // SRC_OPTIMAL_CONTROL_INCLUDE_CORBO_OPTIMAL_CONTROL_OPTIMAL_CONTROL_PROBLEM_INTERFACE_H_
