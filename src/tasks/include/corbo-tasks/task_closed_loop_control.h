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

#ifndef SRC_TASKS_INCLUDE_CORBO_TASKS_TASK_CLOSED_LOOP_CONTROL_H_
#define SRC_TASKS_INCLUDE_CORBO_TASKS_TASK_CLOSED_LOOP_CONTROL_H_

#include <corbo-tasks/task_interface.h>

#include <corbo-systems/filter_interface.h>
#include <corbo-systems/one_step_predictor.h>

#include <memory>
#include <string>

namespace corbo {

/**
 * @brief Perform closed-loop control task
 *
 * @ingroup tasks
 *
 * This class performs a closed-loop control task given
 * a state and control input reference trajectory.
 *
 * @see TaskInterface Environment
 *
 * @author Christoph Rösmann (christoph.roesmann@tu-dortmund.de)
 */
class ClosedLoopControlTask : public TaskInterface
{
 public:
    using Ptr = std::shared_ptr<ClosedLoopControlTask>;

    //! Default constructor
    ClosedLoopControlTask();

    // implements interface method
    TaskInterface::Ptr getInstance() const override { return std::make_shared<ClosedLoopControlTask>(); }

    // implements interface method
    void performTask(Environment& environment, SignalTargetInterface* signal_target = nullptr, std::string* msg = nullptr,
                     const std::string& ns = "") override;

    //! Set state reference trajectory
    void setStateReference(ReferenceTrajectoryInterface::Ptr xreference) { _xreference = xreference; }
    //! Set control input reference trajectory
    void setControlReference(ReferenceTrajectoryInterface::Ptr ureference) { _ureference = ureference; }

    // implements interface method
    bool verify(const Environment& environment, std::string* msg = nullptr) const override;

    // implements interface method
    void getAvailableSignals(const Environment& environment, SignalTargetInterface& signal_target, const std::string& ns = "") const override;

    // implements interface method
    void reset() override;

#ifdef MESSAGE_SUPPORT
    void toMessage(messages::ClosedLoopControlTask& message) const;
    void fromMessage(const messages::ClosedLoopControlTask& message, std::stringstream* issues = nullptr);

    // implements interface method
    void toMessage(messages::Task& message) const override { toMessage(*message.mutable_closed_loop_control_task()); }
    // implements interface method
    void fromMessage(const messages::Task& message, std::stringstream* issues = nullptr) override
    {
        fromMessage(message.closed_loop_control_task(), issues);
    }

#endif

 private:
    double _sim_time           = 10.0;
    double _dt                 = 0.1;
    bool _realtime_sync        = false;
    bool _use_wall_time        = false;
    bool _compensate_dead_time = false;
    bool _compensate_cpu_time  = false;
    double _min_dt             = 0;
    double _max_dt             = CORBO_INF_DBL;

    ReferenceTrajectoryInterface::Ptr _xreference;
    ReferenceTrajectoryInterface::Ptr _ureference;

    double _computation_delay = 0;
    OneStepPredictor _compensator;
    TimeValueBuffer _time_value_buffer;
    FilterInterface::Ptr _computation_delay_filter;
};

FACTORY_REGISTER_TASK(ClosedLoopControlTask)

}  // namespace corbo

#endif  // SRC_TASKS_INCLUDE_CORBO_TASKS_TASK_CLOSED_LOOP_CONTROL_H_
